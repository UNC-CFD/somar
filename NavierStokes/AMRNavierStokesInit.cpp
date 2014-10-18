#include "AMRNavierStokes.H"
#include "MappedAMRPoissonOpFactory.H"
#include "StratUtilsF_F.H"
#include "StratUtils.H"
#include "MappedFineInterp.H"
#include "computeMappedNorm.H"
#include "computeMappedSum.H"
#include "ProblemContext.H"
#include "MappedLevelBackwardEuler.H"
#include "MappedLevelCrankNicolson.H"
#include "MappedLevelTGA.H"
#include "SetValLevel.H"
#include <fstream>
#include <iomanip>

#include "AMRCCProjector.H"


// -----------------------------------------------------------------------------
// Initialize grid
// -----------------------------------------------------------------------------
void AMRNavierStokes::initialGrid (const Vector<Box>& a_new_grids)
{
    if (s_verbosity >= 5) {
        pout () << "AMRNavierStokes::initialGrid " << m_level << endl;
    }

    // Copy this level's new grids
    m_level_grids = a_new_grids;

    // From Wikipedia: Morton ordering maps multidimensional data to one
    // dimension while preserving locality of the data points.
    mortonOrdering(m_level_grids);

    // Is this level empty?
    m_is_empty = !(m_level_grids.size() > 0);

    // Set either this or the coarser level as the finest
    this->finestLevel(!m_is_empty);
    if (m_coarser_level_ptr != NULL) {
        this->crseNSPtr()->finestLevel(m_is_empty);
    }

    // Now balance the load
    const DisjointBoxLayout grids = this->loadBalance(m_level_grids);

    if (s_verbosity >= 4) {
        pout () << "New grids on level " << m_level << ": " << endl;
        for (LayoutIterator lit = grids.layoutIterator(); lit.ok(); ++lit) {
            pout() << grids[lit] << endl;
        }
    }

    // Now that we have grids on this level, regrid the levGeo.
    m_levGeoPtr->reset();
    m_levGeoPtr->regrid(grids);

    // Next, define all of the data holders.
    // Do velocity first.
    IntVect ghostVect(D_DECL(1,1,1));

    if (m_vel_new_ptr != NULL) {
        delete m_vel_new_ptr;
        m_vel_new_ptr = NULL;
    }
    m_vel_new_ptr = new LevelData<FArrayBox>(grids, CH_SPACEDIM, ghostVect);

    if (m_vel_old_ptr != NULL) {
        delete m_vel_old_ptr;
        m_vel_old_ptr = NULL;
    }
    m_vel_old_ptr = new LevelData<FArrayBox>(grids, CH_SPACEDIM, ghostVect);

    // Next, do lambda.
    if (m_lambda_new_ptr != NULL) {
        delete m_lambda_new_ptr;
        m_lambda_new_ptr = NULL;
    }
    m_lambda_new_ptr = new LevelData<FArrayBox>(grids, 1, ghostVect);

    if (m_lambda_old_ptr != NULL) {
        delete m_lambda_old_ptr;
        m_lambda_old_ptr = NULL;
    }
    m_lambda_old_ptr = new LevelData<FArrayBox>(grids, 1, ghostVect);

    // Now do scalars. (We need to allocate flux registers for the scalars too.)
    if (m_scal_new.size() != s_num_scal_comps) {
        //  initialize pointers to null
        Vector<LevelData<FArrayBox>* > tempVect(s_num_scal_comps, NULL);
        m_scal_new = tempVect;
        m_scal_old = tempVect;
    }

    if (m_scal_fluxreg_ptrs.size() < s_num_scal_comps) {
        m_scal_fluxreg_ptrs.resize(s_num_scal_comps,NULL);
    }

    for (int comp = 0; comp < s_num_scal_comps; ++comp) {
        if (m_scal_new[comp] != NULL) {
            delete m_scal_new[comp];
            m_scal_new[comp] = NULL;
        }
        m_scal_new[comp] = new LevelData<FArrayBox>(grids, 1, ghostVect);

        if (m_scal_old[comp] != NULL) {
            delete m_scal_old[comp];
            m_scal_old[comp] = NULL;
        }
        m_scal_old[comp] = new LevelData<FArrayBox>(grids, 1, ghostVect);

        if (m_scal_fluxreg_ptrs[comp] != NULL) {
            delete m_scal_fluxreg_ptrs[comp];
            m_scal_fluxreg_ptrs[comp] = NULL;
        }
        m_scal_fluxreg_ptrs[comp] = new MappedLevelFluxRegister;
    } // end loop over scalar components

    // Finally, pressure.
    m_macPressure.define(grids, 1, IntVect::Unit);
    m_ccPressure.define(grids, 1, IntVect::Unit);
    m_syncPressure.define(grids, 1, IntVect::Unit);
    m_eLambda.define(grids, 1, IntVect::Unit);
    m_gradELambda.define(grids, 1, IntVect::Unit);

    // Fill data holders with bogus data
    this->setAllBogus();
}


// -----------------------------------------------------------------------------
// Performs operations required after the grid has been defined but before
// data initialization. This will also be called after readCheckpointLevel
// during a restart procedure with argument a_restart set to true.
// -----------------------------------------------------------------------------
void AMRNavierStokes::postInitialGrid (const bool a_restart)
{
    if (s_verbosity >= 5) {
        pout () << "AMRNavierStokes::postInitialGrid " << m_level << endl;
    }
}


// -----------------------------------------------------------------------------
// Initialize data
// -----------------------------------------------------------------------------
void AMRNavierStokes::initialData ()
{
    if (s_verbosity >= 5) {
        pout () << "AMRNavierStokes::initialData " << m_level << endl;
    }

    // Do nothing if this level is empty
    if (m_is_empty) return;

    // Gather geometric info
    const Box domBox = m_problem_domain.domainBox();
    const IntVect iv0 = domBox.smallEnd();
    const RealVect dx = m_levGeoPtr->getDx();
    const DisjointBoxLayout& grids = newVel().getBoxes();
    DataIterator dit(grids);

    // Initialize Lambda
    for (dit.reset(); dit.ok(); ++dit) {
        newLambda()[dit].setVal(1.0);
    }

    // Fill velocity ICs
    for (int dir = 0; dir < SpaceDim; ++dir) {
        m_physBCPtr->setVelIC(*m_vel_new_ptr, dir, *m_levGeoPtr);
    }

    // Fill scalar ICs
    for (int comp = 0; comp < s_num_scal_comps; ++comp) {
        m_physBCPtr->setScalarIC(*m_scal_new[comp], comp, *m_levGeoPtr);
    }

    // Send velocity vector to the mapped basis
    m_levGeoPtr->sendToMappedBasis(*m_vel_new_ptr, true);

    // Do exchanges
    m_vel_new_ptr->exchange(m_oneGhostExCopier);
    m_lambda_new_ptr->exchange(m_oneGhostExCopier);
    for (int comp = 0; comp < s_num_scal_comps; ++comp) {
        m_scal_new[comp]->exchange(m_oneGhostExCopier);
    }

    // Remove background scalar from newScal.
    CH_assert(m_time == 0.0);
    for (int comp = 0; comp < s_num_scal_comps; ++comp) {
        m_physBCPtr->subtractBackgroundScalar(newScal(comp), comp, m_time, *m_levGeoPtr);
    }

    // Copy new data to old holders
    for (dit.reset(); dit.ok(); ++dit) {
        (*m_vel_old_ptr)[dit].copy((*m_vel_new_ptr)[dit]);
        (*m_lambda_old_ptr)[dit].copy((*m_lambda_new_ptr)[dit]);
        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            (*m_scal_old[comp])[dit].copy((*m_scal_new[comp])[dit]);
        }
    }

    // Initialize pressure to zero. A better initial guess of the CC pressure
    // will be generated in postInitialize().
    setValLevel(m_macPressure, 0.0);
    setValLevel(m_ccPressure, 0.0);
    setValLevel(m_syncPressure, 0.0);
    setValLevel(m_eLambda, 0.0);
    setValLevel(m_gradELambda, 0.0);

    m_ccPressureState = CCPressureState::ZERO;
    m_syncPressureState = SyncPressureState::ZERO;
    m_eLambdaState = ELambdaState::ZERO;
    m_gradELambdaState = GradELambdaState::ZERO;

    // Compute the KdV stuff (wave speed and vertical structure).
    if (m_level == 0 && !m_c0iPtr) {
        m_c0iPtr = new LevelData<FArrayBox>(grids, SpaceDim);
        m_phi0Ptr = new LevelData<FArrayBox>(grids, 1);
        this->initializeInternalWaveSpeed();
    }
}


// -----------------------------------------------------------------------------
// Solve for the internal wave speed due to stratification.
// WARNING: This must be called by level 0 first! Then work your way up the
// levels.
// -----------------------------------------------------------------------------
void AMRNavierStokes::initializeInternalWaveSpeed ()
{
    CH_TIME("AMRNavierStokes::initializeInternalWaveSpeed");

    static int maxLevelSolved = -1;
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();

    if (m_level == 0) {
        // Solve the eigenvalue problem...
        CH_TIME("Eigenvalue_problem");
        if (s_verbosity > 3) {
            pout() << "Solving the eigenproblem on level 0..." << endl;
        }

        DataIterator dit = grids.dataIterator();

        // Fill the FC background buoyancy field.
        LevelData<FluxBox> bbar(grids, 1);
        for (dit.reset(); dit.ok(); ++dit) {
            FluxBox& bbarFB = bbar[dit];

            D_TERM(m_physBCPtr->setBackgroundScalar(bbarFB[0], 0, *m_levGeoPtr, dit(), m_time);,
                   m_physBCPtr->setBackgroundScalar(bbarFB[1], 0, *m_levGeoPtr, dit(), m_time);,
                   m_physBCPtr->setBackgroundScalar(bbarFB[2], 0, *m_levGeoPtr, dit(), m_time);)
        }

        // Compute the BV freq^2.
        LevelData<FArrayBox> Nsq(grids, 1);
        computeBVFreq(Nsq, bbar, *m_levGeoPtr, 2.0);

        // Solve d^2[phi]/dz^2 + (Nsq/c^2)phi = 0 w/ homog Diri BCs.
        LevelData<FArrayBox> c0(grids, 1);
        solveVertEigenProblem(c0, *m_phi0Ptr, Nsq, *m_levGeoPtr);

        // Loop over grids and project the phase speed, c0,
        // onto the mapped coordinate lines.
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& c0iFAB = (*m_c0iPtr)[dit];
            const FArrayBox& c0FAB = c0[dit];
            const Box& valid = grids[dit];

            // Fill dx^i/dXi^j
            FArrayBox dXidxFAB(valid, SpaceDim*SpaceDim);
            m_levGeoPtr->fill_dXidx(dXidxFAB);

            FORT_PROJECTPHASESPEED(
                CHF_FRA(c0iFAB),
                CHF_CONST_FRA1(c0FAB,0),
                CHF_CONST_FRA(dXidxFAB),
                CHF_BOX(valid));
        }

        maxLevelSolved = 0;

    } else {
        MayDay::Error("Right now, I'm not allowing this to be solved on levels > 0. "
                      "Use the fill function instead.");
    }
}


// -----------------------------------------------------------------------------
// postInitialize
// -----------------------------------------------------------------------------
void AMRNavierStokes::postInitialize ()
{
    if (s_verbosity >= 5) {
        pout () << "AMRNavierStokes::postInitialize " << m_level << endl;
    }
    CH_TIME("AMRNavierStokes::postInitialize");

    // Do nothing if this level is empty
    if (m_is_empty) return;


    // Perform initial projection and set up the initial pressure.
    // This must be done over the entire hierarchy.
    // We do this from level 0 since postInitialize() is
    // called from fine->coarse.
    if (m_level == 0) {
        // Calculate the number of levels to set up.
        AMRNavierStokes* levelNSPtr = finestNSPtr();
        const int numLevels = levelNSPtr->m_level + 1;

        // Initialize data structures which haven't yet been initialized --
        // solvers and operators and such.
        levelNSPtr = this;
        for (int lev = 0; lev < numLevels; ++lev) {
            levelNSPtr->levelSetup(levelNSPtr->newVel().getBoxes());
            levelNSPtr = levelNSPtr->fineNSPtr();
        }

        // Prepare for initial projection...
        // Gather the velocities from each level
        Vector<LevelData<FArrayBox>*> amrVel = gatherNewVelStartingWith(0);
        CH_assert(amrVel.size() == numLevels);

        // Set physical BCs on velocity
        bool isViscous = (s_nu > 0.0);
        VelBCHolder velBC(m_physBCPtr->uStarFuncBC(isViscous));

        levelNSPtr = this;
        for (int lev = 0; lev < numLevels; ++lev) {
            LevelData<FArrayBox>& levelVel = *amrVel[lev];
            const RealVect& levelDx = levelNSPtr->m_levGeoPtr->getDx();
            const LevelData<FluxBox>& levelJgup = levelNSPtr->m_levGeoPtr->getFCJgup();

            velBC.setGhosts(levelVel,
                            NULL,
                            levelDx,
                            &levelJgup,
                            false, // is homogeneous?
                            m_time);

            levelNSPtr = levelNSPtr->fineNSPtr();
        }

        // Do initial velocity projection
        if (s_isIncompressible && s_initial_projection_iters > 0) {
            pout() << "Init projection on levels 0 to "
                   << numLevels-1 << ": " << std::flush;

            Vector<LevelData<FArrayBox>*> eSync = gatherSyncPressure();

            AMRCCProjector projObj;
            projObj.define(eSync, *m_physBCPtr, *m_levGeoPtr, NULL);

            for (int iter = 0; iter < s_initial_projection_iters; ++iter) {
                projObj.project(amrVel,
                                *m_levGeoPtr, //amrLevGeos,
                                0,           // lmin
                                numLevels-1, // lmax
                                0.0,         // newTime
                                1.0,         // dt
                                false,       // velIsFlux,
                                true,        // zero-out pressure
                                false);      // force homogeneous
            }

            // m_syncPressure on this and all higher levels contains init
            // projection data that should not be used to influence dynamics.
            setSyncPressureStates(SyncPressureState::INIT);

            // Reset BCs on projected velocity
            levelNSPtr = this;
            for (int lev = 0; lev < numLevels; ++lev) {
                LevelData<FArrayBox>& levelVel = *amrVel[lev];
                const RealVect& levelDx = levelNSPtr->m_levGeoPtr->getDx();
                const LevelData<FluxBox>& levelJgup = levelNSPtr->m_levGeoPtr->getFCJgup();

                velBC.setGhosts(levelVel,
                                NULL,
                                levelDx,
                                &levelJgup,
                                false, // is homogeneous?
                                m_time);

                levelNSPtr = levelNSPtr->fineNSPtr();
            }
        }
        // Now, the velocities are projected in the composite sense
        // and BCs have been set.

        // No need to initialize the VD correction here.
        // lambda should already be set to one and grad set to zero.
        // Just specify that this is valid data.
        setELambdaStates(ELambdaState::VALID);
        setGradELambdaStates(GradELambdaState::VALID);

        // Initialize the pressures.
        if (s_initial_pressure_iters > 0) {
            this->initializeGlobalPressure();
        }

        // If desired, dump out grids and processor mappings
        if (s_write_grids) {
#ifdef CH_MPI
            if (procID() == 0)
#endif
            {
                ofstream os("grids.out", ios::out);
                if (os.fail()) {
                    pout() << "cannot open grid output file " << os << endl;
                    MayDay::Error();
                }

                levelNSPtr = this;
                os << "NumLevels = " << numLevels << endl;
                for (int lev = 0; lev < numLevels; ++lev) {
                    const DisjointBoxLayout& thisLevelGrids = levelNSPtr->m_vel_new_ptr->getBoxes();
                    os << "Number of Grids on Level " << lev << ": "
                       << thisLevelGrids.size() << "\n"
                       << thisLevelGrids;
                    levelNSPtr = levelNSPtr->fineNSPtr();
                }
                os.close();
            } // end if procID == 0

        } // end if writing grids to file
    } // end if level 0

    // Write terminal output...
    if (s_write_stdout) {
        // Write header
        if (procID() == 0) {
            static bool headerIsWritten = false;
            if (!headerIsWritten) {
                std::cout << color::hiwhite << std::left << setw(6) << "iter"
                          << std::left << setw(8) << "level"
                          << std::left << setw(18) << "time"
                          << std::left << setiosflags(ios::fixed) << setprecision(8) << setw(18) << "dt"
                          << color::higreen
                          D_TERM(
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << "max|u|",
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << "max|v|",
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << "max|w|")
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << "max|b|"
                          << color::hiblue
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << "Sum[m]"
                          D_TERM(
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << "Sum[u]",
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << "Sum[v]",
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << "Sum[w]")
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << "Sum[E]"
                          << std::endl;
                const int totalWidth = 6 + 8 + 18 + 18 + SpaceDim*18 + 18 + 18 + SpaceDim*18 + 18;
                std::cout << color::white << std::string(totalWidth, '-') << color::none << std::endl;
                headerIsWritten = true;
            }
        }

        // Write level output
        {
            // First, calculate max norms
            const LevelGeometry* fineLevGeoPtr = m_levGeoPtr->getFinerPtr();
            const DisjointBoxLayout* fineGridsPtr = NULL;
            if (fineLevGeoPtr != NULL) {
                fineGridsPtr = &(fineLevGeoPtr->getBoxes());
            }

            const Real tol = 1e8;

            // Compute max|velocity|
            D_TERM(
            const Real unorm = computeUnmappedNorm(*m_vel_new_ptr, fineGridsPtr, *m_levGeoPtr, 0, 0);,
            const Real vnorm = computeUnmappedNorm(*m_vel_new_ptr, fineGridsPtr, *m_levGeoPtr, 0, 1);,
            const Real wnorm = computeUnmappedNorm(*m_vel_new_ptr, fineGridsPtr, *m_levGeoPtr, 0, 2);)

            // Compute max|buoyancy|
            Real bnorm = 0.0;
            if (s_num_scal_comps > 0) bnorm = computeUnmappedNorm(*m_scal_new[0], fineGridsPtr, *m_levGeoPtr, 0, 0);

            // Then, write results to stdout
            if (procID() == 0) {
                std::cout << color::white << "      "
                          << std::left << setw(8) << m_level
                          << std::left << setiosflags(ios::fixed) << setprecision(8) << setw(18) << m_time
                          << std::left << setiosflags(ios::fixed) << setprecision(8) << setw(18) << ((m_dt < tol)? m_dt: -123)
                          << color::green
                          D_TERM(
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << unorm,
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << vnorm,
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << wnorm)
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << bnorm
                          << color::none
                          << std::endl;
            }
        }

        // Write composite output
        static int iteration = -1;
        if (m_level == 0) {
            // Gather composite data
            Vector<LevelData<FArrayBox>*> amrVel(0), amrB(0);
            AMRNavierStokes* levelNSPtr = this;
            while(levelNSPtr != NULL) {
                amrVel.push_back(levelNSPtr->m_vel_new_ptr);
                if (s_num_scal_comps > 0) amrB.push_back(levelNSPtr->m_scal_new[0]);

                levelNSPtr = levelNSPtr->fineNSPtr();
            }

            const Real tol = 1e8;

            // Compute max|velocity|
            D_TERM(
            const Real unorm = computeUnmappedNorm(amrVel, *m_levGeoPtr, 0, 0);,
            const Real vnorm = computeUnmappedNorm(amrVel, *m_levGeoPtr, 0, 1);,
            const Real wnorm = computeUnmappedNorm(amrVel, *m_levGeoPtr, 0, 2);)

            // Compute max|buoyancy|
            Real bnorm = 0.0;
            if (s_num_scal_comps > 0) bnorm = computeMappedNorm(amrB, *m_levGeoPtr, 0, 0);

            // Compute total mass
            Real mass = 0.0;
            if (s_num_scal_comps > 0) mass = computeMappedSum(amrB, *m_levGeoPtr, 0);

            // Compute total energy
            Real globalEnergy = this->totalEnergy();

            // Compute total momenta
            RealVect totalMom(D_DECL(computeMappedSum(amrVel, *m_levGeoPtr, 0),
                                     computeMappedSum(amrVel, *m_levGeoPtr, 1),
                                     computeMappedSum(amrVel, *m_levGeoPtr, 2)));

            // Write results to stdout
            if (procID() == 0) {
                std::cout << color::hiwhite << std::left << setw(6) << ++iteration
                          << "---     "
                          << std::left << setiosflags(ios::fixed) << setprecision(8) << setw(18) << m_time
                          << std::left << setiosflags(ios::fixed) << setprecision(8) << setw(18) << ((m_dt < tol)? m_dt: -123)
                          << color::higreen
                          D_TERM(
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << unorm,
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << vnorm,
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << wnorm)
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << bnorm
                          << color::hiblue
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << mass
                          D_TERM(
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << totalMom[0],
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << totalMom[1],
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << totalMom[2])
                          << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << globalEnergy
                          << color::none
                          << "\n"
                          << std::endl;
            }

            // Save the total energy for diagnostics at a later time step
            s_totalEnergy = globalEnergy;

        } // end if level 0
    } // end if write to terminal

    pout() << endl; // For some reason, this prevents a segfault in release mode.
}



// -----------------------------------------------------------------------------
// Set up data structures associated with this level. This function is called
// from the coarsest level that needs setting up to the finest.
// -----------------------------------------------------------------------------
void AMRNavierStokes::levelSetup (const DisjointBoxLayout& a_grids)
{
    if (s_verbosity >= 5) {
        pout() << "AMRNavierStokes::levelSetup on level " << m_level << endl;
        pout() << "grids = " << a_grids << endl;
    }

    CH_assert(!m_is_empty);

    const RealVect dx = m_levGeoPtr->getDx();
    const DisjointBoxLayout* crseGridsPtr = NULL;
    AMRNavierStokes* crse_amrns_ptr = this->crseNSPtr();
    AMRNavierStokes* fine_amrns_ptr = this->fineNSPtr();
    LevelGeometry* crseLevGeoPtr = NULL;
    LevelGeometry* fineLevGeoPtr = NULL;
    IntVect nRefCrse(D_DECL(-1,-1,-1));


    // Set levGeo hierarchy ------------------------
    // Get the coarser levGeoPtr
    if (crse_amrns_ptr != NULL) {
        crseLevGeoPtr = crse_amrns_ptr->m_levGeoPtr;
    }
    m_levGeoPtr->setCoarserPtr(crseLevGeoPtr);

    // Get the finer levGeoPtr
    // TODO: Maybe we want this to be NULL and then repaired if the finer level exists.
    if (fine_amrns_ptr != NULL) {
        fineLevGeoPtr = fine_amrns_ptr->m_levGeoPtr;
    }
    m_levGeoPtr->setFinerPtr(fineLevGeoPtr);

    m_levGeoPtr->regrid(a_grids);
    // end set levGeo hierarchy --------------------


    // Set up data structures that depend on a coarser level.
    if (m_coarser_level_ptr != NULL) {
        // Grab coarse level data
        nRefCrse = crse_amrns_ptr->m_ref_ratio;
        // crseProjPtr = &(crse_amrns_ptr->m_projection);
        crseGridsPtr = &(crse_amrns_ptr->newVel().getBoxes());
        CH_assert(m_levGeoPtr->getCoarserPtr()->getBoxes() == *crseGridsPtr);

        const LevelGeometry* crseLevGeoPtr = m_levGeoPtr->getCoarserPtr();

        // Define averaging tools
        m_coarse_average.define(a_grids, *crseGridsPtr, SpaceDim, nRefCrse, IntVect::Zero);
        m_coarse_average_scal.define(a_grids, *crseGridsPtr, 1, nRefCrse, IntVect::Zero);

        // Define flux registers and initialize to zero
        crse_amrns_ptr->m_vel_flux_reg.undefine();
        crse_amrns_ptr->m_vel_flux_reg.define(a_grids,
                                               *crseGridsPtr,
                                               m_problem_domain,
                                               nRefCrse,
                                               SpaceDim);
        crse_amrns_ptr->m_vel_flux_reg.setToZero();

        if (s_advective_lambda_reflux) {
            crse_amrns_ptr->m_lambda_flux_reg.undefine();
            crse_amrns_ptr->m_lambda_flux_reg.define(a_grids,
                                                     *crseGridsPtr,
                                                     m_problem_domain,
                                                     nRefCrse,
                                                     1);
            crse_amrns_ptr->m_lambda_flux_reg.setToZero();
        }

        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            crse_amrns_ptr->m_scal_fluxreg_ptrs[comp]->undefine();
            crse_amrns_ptr->m_scal_fluxreg_ptrs[comp]->define(a_grids,
                                                              *crseGridsPtr,
                                                              m_problem_domain,
                                                              nRefCrse,
                                                              1);
            crse_amrns_ptr->m_scal_fluxreg_ptrs[comp]->setToZero();
        }

        // Define CF interpolation operator
        m_velCFInterp.define(a_grids, crseGridsPtr,
                             dx, nRefCrse, SpaceDim,
                             m_problem_domain);
    } // end if coarser level exists


    // Define viscous/diffusive operators...
    // (This needs knowledge of a coarser level, but not finer levels.)

    // If one of the operators are created, then we can use its diagonals
    // for all other operators that are needed. This will reduce memory
    // requirements.
    RefCountedPtr<LevelData<FArrayBox> > lapDiagsPtr(NULL);

    // Check if any of the scalars are diffusive
    bool isDiffusive = false;
    for (int comp = 0; comp < s_num_scal_comps; ++comp) {
        if (s_scal_coeffs[comp] > 0) isDiffusive = true;
    }

    // If so, define the Laplacian op for the scalars
    // and grab a pointer to its diagonals.
    if (isDiffusive) {
        // Define this level's op as if there is no finer level just in case
        // this level doesn't get tagged. We do this through the levGeo ptrs.
        const LevelGeometry* fineLevGeoSave = m_levGeoPtr->getFinerPtr();
        m_levGeoPtr->setFinerPtr(NULL);

        BCMethodHolder bcHolder = m_physBCPtr->diffusiveSourceFuncBC();
        m_scalarsAMRPoissonOp.define(bcHolder, *m_levGeoPtr, 0.0, 1.0,
                                     RefCountedPtr<LevelData<FArrayBox> >(NULL));   // lapDiags
        TODO(); // lapDiagsPtr = m_scalarsAMRPoissonOp.m_lapDiag;

        m_levGeoPtr->setFinerPtr(fineLevGeoSave);

        if (m_coarser_level_ptr != NULL) {
            // Define the coarser level's op
            BCMethodHolder bcHolder = crse_amrns_ptr->m_physBCPtr->diffusiveSourceFuncBC();
            crse_amrns_ptr->m_scalarsAMRPoissonOp.define(bcHolder, *crseLevGeoPtr, 0.0, 1.0,
                                                         RefCountedPtr<LevelData<FArrayBox> >(NULL));   // lapDiagsPtr
        }
    }

    // Define viscous solver here as well
    if (s_nu > 0.0) {
        // Define this level's op as if there is no finer level just in case
        // this level doesn't get tagged. We do this through the levGeo ptrs.
        const LevelGeometry* fineLevGeoSave = m_levGeoPtr->getFinerPtr();
        m_levGeoPtr->setFinerPtr(NULL);

        m_velocityAMRPoissonOp.define(*m_levGeoPtr,
                                      0.0,
                                      1.0,
                                      RefCountedPtr<LevelData<FArrayBox> >(NULL));   // lapDiagsPtr

#       ifdef USE_STRESSMETRIC
            m_stressMetric.define(m_levGeoPtr->getGeoSourcePtr());
            m_velocityAMRPoissonOp.setJgup(&m_stressMetric);
#       endif

        m_levGeoPtr->setFinerPtr(fineLevGeoSave);

        // Now, (re)define the coarser level's ops
        if (m_coarser_level_ptr != NULL) {
            // Define the coarser level's op
            crse_amrns_ptr->m_velocityAMRPoissonOp.define(*crseLevGeoPtr, 0.0, 1.0,
                                                          RefCountedPtr<LevelData<FArrayBox> >(NULL));   // lapDiagsPtr

#           ifdef USE_STRESSMETRIC
                crse_amrns_ptr->m_velocityAMRPoissonOp.setJgup(&m_stressMetric);
#           endif
        }
    }

    // Define viscous and diffusive solvers.
    if ((s_nu > 0.0) || isDiffusive) {
        this->defineViscousMGSolver(a_grids, crseGridsPtr, nRefCrse);
    }

    // The flux registers will be defined by the finer level itself,
    // since it represents data from the fine grids. However, if we
    // have no finer level, we won't need flux registers.
    if (m_finer_level_ptr == NULL) {
        m_vel_flux_reg.undefine();
        m_lambda_flux_reg.undefine();
        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            m_scal_fluxreg_ptrs[comp]->undefine();
        }
    }

    CH_assert(m_level == 0 || nRefCrse == m_levGeoPtr->getCrseRefRatio());

    // Define the velocity projector.
    // Remember: this function is always called from coarse to fine levels.
    if (s_isIncompressible) {
        const LevelData<FArrayBox>* crsePresBC = NULL;
        if (m_level > 0) crsePresBC = &(crseNSPtr()->m_ccPressure);

        m_macProjector.define(&m_macPressure, crsePresBC, *m_physBCPtr, *m_levGeoPtr);
        m_ccProjector.define(&m_ccPressure, crsePresBC, *m_physBCPtr, *m_levGeoPtr);

        setValLevel(m_macPressure, 0.0);
        setValLevel(m_ccPressure, 0.0);
        setValLevel(m_syncPressure, 0.0);
        setValLevel(m_eLambda, 0.0);
        setValLevel(m_gradELambda, 0.0);

        m_ccPressureState = CCPressureState::ZERO;
        m_syncPressureState = SyncPressureState::ZERO;
        m_eLambdaState = ELambdaState::ZERO;
        m_gradELambdaState = GradELambdaState::ZERO;

    } else {
        // Just for good measure.
        setValLevel(m_macPressure, 0.0);
        setValLevel(m_ccPressure, 0.0);
        setValLevel(m_syncPressure, 0.0);
        setValLevel(m_eLambda, 0.0);
        setValLevel(m_gradELambda, 0.0);

        m_ccPressureState = CCPressureState::VALID;
        m_syncPressureState = SyncPressureState::SYNC;
        m_eLambdaState = ELambdaState::VALID;
        m_gradELambdaState = GradELambdaState::ZERO;
    }

    // Compute the KdV stuff (wave speed and vertical structure).
    if (m_level == 0 && !m_c0iPtr) {
        m_c0iPtr = new LevelData<FArrayBox>(a_grids, SpaceDim);
        m_phi0Ptr = new LevelData<FArrayBox>(a_grids, 1);
        this->initializeInternalWaveSpeed();
    }

    TODO(); // Try to trim the edges, etc.
    // Define exchange copiers
    m_oneGhostExCopier.define(a_grids, a_grids, m_problem_domain, IntVect::Unit, true);
    // m_oneGhostExCopier.exchangeDefine(a_grids, IntVect::Unit);
    // m_oneGhostExCopier.trimEdges(a_grids, IntVect::Unit);

    // Define the exchange copiers used by the hyperbolic tracing scheme
    m_tracingGhosts = IntVect(D_DECL(ADVECT_GROW, ADVECT_GROW, ADVECT_GROW));
    // m_tracingExCopier.exchangeDefine(a_grids, m_tracingGhosts);
    // m_tracingExCopier.trimEdges(a_grids, m_tracingGhosts);
    m_tracingExCopier.define(a_grids, a_grids, m_problem_domain, m_tracingGhosts, true);
    m_tracingExCornerCopier.define(a_grids, a_grids, m_problem_domain, m_tracingGhosts, true);

    // Define the advection utilities
    m_advectUtilVel.define(m_levGeoPtr,
                           s_normalPredOrderVel,
                           s_useFourthOrderSlopesVel,
                           s_useLimitingVel,
                           s_useHighOrderLimiterVel,
                           s_useUpwindingVel);

    m_advectUtilLambda.define(m_levGeoPtr,
                              s_normalPredOrderVel,
                              s_useFourthOrderSlopesVel,
                              s_useLimitingVel,
                              s_useHighOrderLimiterVel,
                              s_useUpwindingVel);

    m_advectUtilScal.define(m_levGeoPtr,
                            s_normalPredOrderScal,
                            s_useFourthOrderSlopesScal,
                            s_useLimitingScal,
                            s_useHighOrderLimiterScal,
                            s_useUpwindingScal);
}


// -----------------------------------------------------------------------------
// Define the viscous and diffusive solvers.
// -----------------------------------------------------------------------------
void AMRNavierStokes::defineViscousMGSolver (const DisjointBoxLayout& a_grids,
                                             const DisjointBoxLayout* a_crseGridsPtr,
                                             const IntVect&           a_refCrse)
{
    if (s_verbosity >= 6) {
        pout() << "AMRNavierStokes::defineViscousMGSolver " << m_level << endl;
    }

    // Sanity check
    CH_assert(m_levGeoPtr != NULL);
    CH_assert(m_physBCPtr != NULL);

    // This will make the code easier to read.
    typedef MappedAMRLevelOpFactory< LevelData<FArrayBox> > LOFact;
    typedef MappedAMRMultiGrid<LevelData<FArrayBox> > AMRMGSolver;

    // Define some data needed by the ops and factories.
    ProblemDomain baseDomain(m_problem_domain);
    Real alpha = 1.0;
    int numSolverLevels = (a_crseGridsPtr == NULL) ? 1 : 2;
    Vector<DisjointBoxLayout> allGrids(numSolverLevels);
    Vector<IntVect> refRatios(1, a_refCrse);

    if (a_crseGridsPtr != NULL) {
        // coarser level exists:  define solver on two levels
        ProblemDomain origDomain = baseDomain;
        coarsen(baseDomain, origDomain, a_refCrse);

        allGrids[0] = *a_crseGridsPtr;
        allGrids[1] = a_grids;
        refRatios.push_back(IntVect::Unit);
    } else {
        // no coarser level:  define solver on only one level
        allGrids[0] = a_grids;
    }

    // This can easily become confusing, so I'm gonna define all ops, solvers,
    // and factories from the bottom up. Begin with the most basic solver...

    // Define the bottom solver for the viscous and diffusive level solvers.
    m_viscousBottomSolver.m_eps = s_viscous_bottom_eps;
    m_viscousBottomSolver.m_reps = s_viscous_bottom_reps;
    m_viscousBottomSolver.m_imax = s_viscous_bottom_imax;
    m_viscousBottomSolver.m_numRestarts = s_viscous_bottom_numRestarts;
    m_viscousBottomSolver.m_hang = s_viscous_bottom_hang;
    m_viscousBottomSolver.m_small = s_viscous_bottom_small;
    m_viscousBottomSolver.m_verbosity = s_viscous_bottom_verbosity;
    m_viscousBottomSolver.m_normType = s_viscous_bottom_normType;

    // Working our way up in complexity...

    // Allocate and define the viscous MappedAMRPoissonOpFactory(s) needed by LevelTGA
    Tuple<RefCountedPtr<LOFact>,SpaceDim> velTGAOpFactoryPtrs;
    if (s_viscSolverScheme != ProblemContext::HeatSolverScheme::EXPLICIT) {
        for (int idir = 0; idir < SpaceDim; ++idir) {

            velTGAOpFactoryPtrs[idir] = RefCountedPtr<LOFact>((LOFact*)(new MappedAMRPoissonOpFactory));
            CH_assert(!velTGAOpFactoryPtrs[idir].isNull());

            const FillJgupInterface* fillJgupPtr = NULL;
#           ifdef USE_STRESSMETRIC
                m_stressMetric.define(m_levGeoPtr->getGeoSourcePtr());
                fillJgupPtr = &m_stressMetric;
#           endif

            // NOTE: Changed significantly! Does not perform 2-level setup.
            BCMethodHolder bcHolder = m_physBCPtr->viscousSolveFuncBC(idir);
            MappedAMRPoissonOpFactory& factRef = (MappedAMRPoissonOpFactory&)(*velTGAOpFactoryPtrs[idir]);
            factRef.define(m_levGeoPtr,
                           alpha,
                           s_nu,
                           bcHolder,
                           s_viscous_AMRMG_maxDepth,
                           s_viscous_AMRMG_num_smooth_precond,
                           s_viscous_AMRMG_precondMode,
                           s_viscous_AMRMG_relaxMode,
                           false,  // horizontal factory?
                           fillJgupPtr);
        }
    }

    // Allocate the diffusive MappedAMRPoissonOpFactory(s) needed by LevelTGA
    Vector< RefCountedPtr<LOFact> > scalTGAOpFactoryPtrs;
    if (s_num_scal_comps > 0 && s_diffSolverScheme != ProblemContext::HeatSolverScheme::EXPLICIT) {
        scalTGAOpFactoryPtrs.resize(s_num_scal_comps);
        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            scalTGAOpFactoryPtrs[comp] = RefCountedPtr<LOFact>((LOFact*)(new MappedAMRPoissonOpFactory));
            CH_assert(!scalTGAOpFactoryPtrs[comp].isNull());

            // NOTE: Changed significantly! Does not perform 2-level setup.
            BCMethodHolder bcHolder = m_physBCPtr->diffusiveSolveFuncBC();
            MappedAMRPoissonOpFactory& factRef = (MappedAMRPoissonOpFactory&)(*scalTGAOpFactoryPtrs[comp]);
            factRef.define(m_levGeoPtr,
                           alpha,
                           s_scal_coeffs[comp],
                           bcHolder,
                           s_viscous_AMRMG_maxDepth,
                           s_viscous_AMRMG_num_smooth_precond,
                           s_viscous_AMRMG_precondMode,
                           s_viscous_AMRMG_relaxMode);
        }
    }

    // Working our way up in complexity...

    // Define the viscous AMRMG and BaseLevelHeatSolver solvers
    if (s_viscSolverScheme != ProblemContext::HeatSolverScheme::EXPLICIT) {
        for (int idir = 0; idir < SpaceDim; ++idir) {
            // First, do the AMRMG solvers needed by BaseLevelHeatSolver
            m_viscAMRMGPtrs[idir] = RefCountedPtr<AMRMGSolver>(new AMRMGSolver());
            CH_assert(!m_viscAMRMGPtrs[idir].isNull());

            m_viscAMRMGPtrs[idir]->define(baseDomain, // on either this level or coarser level
                                          *velTGAOpFactoryPtrs[idir],
                                          &m_viscousBottomSolver,
                                          numSolverLevels);
            m_viscAMRMGPtrs[idir]->m_verbosity = s_viscous_AMRMG_verbosity;
            m_viscAMRMGPtrs[idir]->m_imin = s_viscous_AMRMG_imin;
            m_viscAMRMGPtrs[idir]->setSolverParameters(s_viscous_AMRMG_num_smooth_down,
                                                       s_viscous_AMRMG_num_smooth_up,
                                                       s_viscous_AMRMG_num_smooth_bottom,
                                                       s_viscous_AMRMG_numMG,
                                                       s_viscous_AMRMG_imax,
                                                       s_viscous_AMRMG_eps,
                                                       s_viscous_AMRMG_hang,
                                                       s_viscous_AMRMG_normThresh);

            // Finally, do the BaseLevelHeatSolver object.
            switch (s_viscSolverScheme) {
                case ProblemContext::HeatSolverScheme::BACKWARD_EULER:
                {
                    typedef MappedLevelBackwardEuler HeatSolverType;
                    HeatSolverType* newSolver = new HeatSolverType(allGrids,
                                                                   refRatios,
                                                                   baseDomain,
                                                                   velTGAOpFactoryPtrs[idir],  // Not used
                                                                   m_viscAMRMGPtrs[idir]);
                    m_viscSolverPtrs[idir] = RefCountedPtr<HeatSolverType>(newSolver);
                    break;
                }
                case ProblemContext::HeatSolverScheme::CRANK_NICOLSON:
                {
                    typedef MappedLevelCrankNicolson HeatSolverType;
                    HeatSolverType* newSolver = new HeatSolverType(allGrids,
                                                                   refRatios,
                                                                   baseDomain,
                                                                   velTGAOpFactoryPtrs[idir],  // Not used
                                                                   m_viscAMRMGPtrs[idir]);
                    m_viscSolverPtrs[idir] = RefCountedPtr<HeatSolverType>(newSolver);
                    break;
                }
                case ProblemContext::HeatSolverScheme::TGA:
                {
                    typedef MappedLevelTGA HeatSolverType;
                    HeatSolverType* newSolver = new HeatSolverType(allGrids,
                                                                   refRatios,
                                                                   baseDomain,
                                                                   velTGAOpFactoryPtrs[idir],  // Not used
                                                                   m_viscAMRMGPtrs[idir]);
                    m_viscSolverPtrs[idir] = RefCountedPtr<HeatSolverType>(newSolver);
                    break;
                }
                case ProblemContext::HeatSolverScheme::EXPLICIT:
                {
                    // No solver to set up!
                    break;
                }
                default:
                    ostringstream msg;
                    msg << "Invalid viscous solver type. s_viscSolverScheme = "
                        << s_viscSolverScheme;
                    MayDay::Error(msg.str().c_str());
            }
            CH_assert(!m_viscSolverPtrs[idir].isNull());
        } // end loop over vel dirs
    }

    // Define the diffusive AMRMG and BaseLevelHeatSolver solvers
    if (s_num_scal_comps > 0 && s_diffSolverScheme != ProblemContext::HeatSolverScheme::EXPLICIT) {
        m_diffAMRMGPtrs.resize(s_num_scal_comps);
        m_diffSolverPtrs.resize(s_num_scal_comps);
        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            // First, do the AMRMG solvers needed by BaseLevelHeatSolver
            m_diffAMRMGPtrs[comp] = RefCountedPtr<AMRMGSolver>(new AMRMGSolver());
            CH_assert(!m_diffAMRMGPtrs[comp].isNull());

            m_diffAMRMGPtrs[comp]->define(baseDomain, // on either this level or coarser level
                                          *scalTGAOpFactoryPtrs[comp],
                                          &m_viscousBottomSolver,
                                          numSolverLevels);
            m_diffAMRMGPtrs[comp]->m_verbosity = s_viscous_AMRMG_verbosity;
            m_diffAMRMGPtrs[comp]->m_imin = s_viscous_AMRMG_imin;
            m_diffAMRMGPtrs[comp]->setSolverParameters(s_viscous_AMRMG_num_smooth_down,
                                                       s_viscous_AMRMG_num_smooth_up,
                                                       s_viscous_AMRMG_num_smooth_bottom,
                                                       s_viscous_AMRMG_numMG,
                                                       s_viscous_AMRMG_imax,
                                                       s_viscous_AMRMG_eps,
                                                       s_viscous_AMRMG_hang,
                                                       s_viscous_AMRMG_normThresh);

            // Finally, do the BaseLevelHeatSolver object.
            switch (s_diffSolverScheme) {
                case ProblemContext::HeatSolverScheme::BACKWARD_EULER:
                {
                    typedef MappedLevelBackwardEuler HeatSolverType;
                    HeatSolverType* newSolver = new HeatSolverType(allGrids,
                                                                   refRatios,
                                                                   baseDomain,
                                                                   scalTGAOpFactoryPtrs[comp],  // Not used
                                                                   m_diffAMRMGPtrs[comp]);
                    m_diffSolverPtrs[comp] = RefCountedPtr<HeatSolverType>(newSolver);
                    break;
                }
                case ProblemContext::HeatSolverScheme::CRANK_NICOLSON:
                {
                    typedef MappedLevelCrankNicolson HeatSolverType;
                    HeatSolverType* newSolver = new HeatSolverType(allGrids,
                                                                   refRatios,
                                                                   baseDomain,
                                                                   scalTGAOpFactoryPtrs[comp],  // Not used
                                                                   m_diffAMRMGPtrs[comp]);
                    m_diffSolverPtrs[comp] = RefCountedPtr<HeatSolverType>(newSolver);
                    break;
                }
                case ProblemContext::HeatSolverScheme::TGA:
                {
                    typedef MappedLevelTGA HeatSolverType;
                    HeatSolverType* newSolver = new HeatSolverType(allGrids,
                                                                   refRatios,
                                                                   baseDomain,
                                                                   scalTGAOpFactoryPtrs[comp],  // Not used
                                                                   m_diffAMRMGPtrs[comp]);
                    m_diffSolverPtrs[comp] = RefCountedPtr<HeatSolverType>(newSolver);
                    break;
                }
                case ProblemContext::HeatSolverScheme::EXPLICIT:
                {
                    // No solver to set up!
                    break;
                }
                default:
                    ostringstream msg;
                    msg << "Invalid diffusive solver type. s_diffSolverScheme = "
                        << s_diffSolverScheme;
                    MayDay::Error(msg.str().c_str());
            }
            CH_assert(!m_diffSolverPtrs[comp].isNull());
        } // end loop over scalars
    } // end if more than one scalar
}


// -----------------------------------------------------------------------------
// This function manages the pressure initialization after
// initialization and regridding
// -----------------------------------------------------------------------------
void AMRNavierStokes::initializeGlobalPressure ()
{
    if (s_verbosity >= 5) {
        pout() << "AMRNavierStokes::initializeGlobalPressure " << m_level << endl;
    }

    // Find out what the finest level is
    AMRNavierStokes* thisNSPtr = this;
    while (!thisNSPtr->finestLevel()) {
        thisNSPtr = thisNSPtr->fineNSPtr();
    }
    CH_assert(thisNSPtr->finestLevel());
    int finest_level = thisNSPtr->m_level;

    // Find the coarse level dt.
    Real crseDtLimit = 1e8;
    if (m_level > 0) {
        const AMRNavierStokes* crseNSPtr = this->crseNSPtr();
        const IntVect& crseRefRatio = crseNSPtr->m_ref_ratio;

        D_TERM(int maxRatio = crseRefRatio[0];,
               if (crseRefRatio[1] > maxRatio) maxRatio = crseRefRatio[1];,
               if (crseRefRatio[2] > maxRatio) maxRatio = crseRefRatio[2];)
        CH_assert(maxRatio >= 2);

        crseDtLimit = crseNSPtr->m_dt / Real(maxRatio);

        if (s_verbosity >= 4) {
            pout() << setiosflags(ios::fixed)
                   << "\tlevel " << m_level-1 << ": "
                   << "\tm_time = " << crseNSPtr->m_time
                   << "\tm_dt = " << crseNSPtr->m_dt
                   << "\tcrseDtLimit = " << crseDtLimit
                   << "\n";
        }
    }

    // Find the smallest dt over all of the levels.
    // also save dt for each level so it can be reset later
    const Real cur_time = m_time;
    Vector<Real> dtSave(finest_level+1, 1.0e8);
    Real dtInit = 10000000.0;
    thisNSPtr = this;
    for (int lev = m_level; lev <= finest_level; ++lev) {
        dtSave[lev] = thisNSPtr->m_dt;

        // TODO: Maybe we should interpolate the pressure from the coarser
        // levels to give computeDt a better idea of the current state.
        const Real dtLevel = thisNSPtr->computeDt();
        if (dtLevel < dtInit) dtInit = dtLevel;

        if (s_verbosity >= 4) {
            pout() << setiosflags(ios::fixed)
                   << "\tlevel " << lev << ": "
                   << "\tm_time = " << thisNSPtr->m_time
                   << "\tdtSave = " << dtSave[lev]
                   << "\tdtLevel = " << dtLevel
                   << "\n";
        }

        thisNSPtr = thisNSPtr->fineNSPtr();
    }

    // Limit dtInit's value based on coarser level dt.
    if (dtInit > crseDtLimit && crseDtLimit > 0.0) {
        dtInit = crseDtLimit;
    }

    // Finally, set dtInit equal to half of the finest level's dt.
    dtInit *= 0.5;
    if (s_verbosity >= 4) {
        pout() << "\tdtInit = " << dtInit << endl;
    }

    // Loop through levels and compute estimate of level pressure (Pi)
    for (int iter = 0; iter < s_initial_pressure_iters; ++iter) {
        if (s_verbosity >= 2) {
            pout() << "\nInitializing level pressures on levels "
                   << m_level << " to " << finest_level
                   << " -- pass #" << iter + 1 << endl;
        }

        int lbase = m_level;
        AMRNavierStokes* thisNSPtr = this;

        // Initialize each level's pressure field (Pi)
        for (int lev = lbase; lev <= finest_level; ++lev) {
            // Initialize the _new_ pressure fields to zero
            setValLevel(thisNSPtr->m_ccPressure, 0.0);
            thisNSPtr->m_ccPressureState = CCPressureState::ZERO;

            // Dump data to pout.*
            if (s_verbosity >= 2) {
                pout() << setiosflags(ios::fixed)
                       << "Initializing level pressure on level "
                       << thisNSPtr->m_level << " of " << finest_level
                       << " pass #" << iter + 1
                       << setiosflags(ios::scientific) << setprecision(8) << endl;
            }
            if (s_verbosity >= 4) {
                pout() << setiosflags(ios::fixed) << setprecision(8)
                       << "\ta_currentTime = m_time = " << cur_time
                       << "\n\ta_dtInit = m_dt = " << dtInit
                       << setiosflags(ios::scientific) << setprecision(8) << endl;
            }

            // Do the time-stepping procedure
            CH_assert(thisNSPtr->m_level == lev);
            CH_assert(   (thisNSPtr->m_time - TIME_EPS < cur_time)
                      && (cur_time < thisNSPtr->m_time + TIME_EPS)   );
            thisNSPtr->dt(dtInit);
            thisNSPtr->swapOldAndNewStates();

            if (s_updateScheme == ProblemContext::UpdateScheme::FiniteVolume) {
                if (s_gravityMethod == ProblemContext::GravityMethod::IMPLICIT) {
                    thisNSPtr->PPMIGTimeStep(cur_time,
                                             m_dt,
                                             false,     // updatePassiveScalars
                                             true);     // doLevelProj
                } else {
                    thisNSPtr->PPMTimeStep(cur_time,
                                           m_dt,
                                           false,     // updatePassiveScalars
                                           true);     // doLevelProj
                }

            } else if (s_updateScheme == ProblemContext::UpdateScheme::RK3) {
                thisNSPtr->RK3TimeStep(cur_time,
                                       m_dt,
                                       false,     // updatePassiveScalars
                                       true);     // doLevelProj
            } else {
                MayDay::Error("Unrecognized update scheme");
            }

            // Dump data to pout.*
            if (s_verbosity >= 2) {
                pout() << setiosflags(ios::fixed)
                       << "Finished initializing level pressure on level "
                       << thisNSPtr->m_level << "\n"
                       << setiosflags(ios::scientific) << setprecision(8) << endl;
            }

            // Move on to the next level
            thisNSPtr = thisNSPtr->fineNSPtr();
        }

        // Now reset times and dt
        // NOTE: This must be done after ALL levels are advanced
        // otherwise time interpolations may be taken in error.
        thisNSPtr = this;
        for (int lev = lbase; lev <= finest_level; ++lev) {
            thisNSPtr->resetStates(cur_time);
            thisNSPtr->dt(dtSave[lev]);
            thisNSPtr = thisNSPtr->fineNSPtr();
        }
    } // end loop over init passes
}
