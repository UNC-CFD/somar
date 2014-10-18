/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Copyright (C) 2014 Edward Santilli & Alberto Scotti
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *  USA
 *
 *  For up-to-date contact information, please visit the repository homepage,
 *  https://github.com/somarhub.
 ******************************************************************************/
#include "AMRNavierStokes.H"
#include "AMRNSF_F.H"
#include "computeMappedNorm.H"
#include "computeMappedSum.H"
#include "SetValLevel.H"
#include "MappedAMRPoissonOpFactory.H"
#include "Constants.H"
#include "AMRLESMeta.H"
#include "Printing.H"
#include "AMRCCProjector.H"
#include <iomanip>


// -----------------------------------------------------------------------------
// Things to do after a basic timestep - synchronization and diagnostics
// -----------------------------------------------------------------------------
void AMRNavierStokes::postTimeStep()
{
    if (s_verbosity >= 5) {
        pout () << "AMRNavierStokes::postTimeStep " << m_level << endl;
    }
    CH_TIME("AMRNavierStokes::postTimeStep");

    // If we are at the finest level, just print diagnostic info and scram.
    if (finestLevel()) {
        if (m_level == 0) {
            // TODO: This may be better suited in tagCells on level 0 so
            // that we know not to tag the turbulent region.
            this->syncWithSGS();
            this->syncSingleGridDiagnostics();
        }
        this->syncTermDiagnostics();

        return;
    }

    // Multi-grid and not finest level: we have work to do...
    CH_assert(!finestLevel());

    // 1.) Do refluxing and avgDown for conservation.

    // 1.0) Gather some needed objects
    AMRNavierStokes* fineAMRNavierStokesPtr = fineNSPtr();
    bool isViscous = (s_nu > 0.0);
    VelBCHolder velBC(m_physBCPtr->uStarFuncBC(isViscous));
    const RealVect& dx = m_levGeoPtr->getDx();
    const bool considerCellVols = true;

    if (!finestLevel()) {
        // 1.1) Momentum refluxing
        {
            LevelData<FArrayBox>& newVel = *m_vel_new_ptr;
            if (s_advective_momentum_reflux || s_diffusive_momentum_reflux) {
                if (s_implicit_momentum_reflux && s_nu > 0.0) {
                    // [Defer implicit refluxing until we're doing sync projection]
                } else {
                    const DisjointBoxLayout& grids = newVel.disjointBoxLayout();
                    DataIterator dit = newVel.dataIterator();

                    // Instead of sending newVel to a Cartesian basis,
                    // refluxing, then sending it back to a mapped basis, let's
                    // just create a field that can hold the refluxing gain, be
                    // converted to a mapped basis, then added to newVel.
                    LevelData<FArrayBox> velRefluxGain(grids, SpaceDim);
                    setValLevel(velRefluxGain, 0.0);

                    m_vel_flux_reg.reflux(velRefluxGain, *m_levGeoPtr);
                    m_levGeoPtr->sendToMappedBasis(velRefluxGain);
                    for (dit.reset(); dit.ok(); ++dit) {
                        newVel[dit].plus(velRefluxGain[dit], 0, 0, SpaceDim);
                    }
                    velRefluxGain.clear();

                    // Now, send the average of the new fine data down to the
                    // next coarser level.
                    MappedCoarseAverage& velAvgDown = fineAMRNavierStokesPtr->m_coarse_average;
                    velAvgDown.averageToCoarse(newVel,
                                               *fineAMRNavierStokesPtr->m_vel_new_ptr,
                                               fineAMRNavierStokesPtr->m_levGeoPtr,
                                               considerCellVols);
                }
            } else {
                // We are not refluxing. Just to the averaging.
                MappedCoarseAverage& velAvgDown = fineAMRNavierStokesPtr->m_coarse_average;
                velAvgDown.averageToCoarse(newVel,
                                           *fineAMRNavierStokesPtr->m_vel_new_ptr,
                                           fineAMRNavierStokesPtr->m_levGeoPtr,
                                           considerCellVols);
            }
        }

        // 1.2) Lambda refluxing
        if (s_advective_lambda_reflux) {
            LevelData<FArrayBox>& newLambda = *m_lambda_new_ptr;

            m_lambda_flux_reg.reflux(newLambda, *m_levGeoPtr);

            MappedCoarseAverage& lambdaAvgDown = fineAMRNavierStokesPtr->m_coarse_average_scal;
            lambdaAvgDown.averageToCoarse(newLambda,
                                          *fineAMRNavierStokesPtr->m_lambda_new_ptr,
                                          fineAMRNavierStokesPtr->m_levGeoPtr,
                                          considerCellVols);
        } else {
            // We are not refluxing. Just to the averaging.
            LevelData<FArrayBox>& newLambda = *m_lambda_new_ptr;

            MappedCoarseAverage& lambdaAvgDown = fineAMRNavierStokesPtr->m_coarse_average_scal;
            lambdaAvgDown.averageToCoarse(newLambda,
                                          *fineAMRNavierStokesPtr->m_lambda_new_ptr,
                                          fineAMRNavierStokesPtr->m_levGeoPtr,
                                          considerCellVols);

        }

        // 1.3) Scalar refluxing
        if (s_advective_scalar_reflux || s_diffusive_scalar_reflux) {
            if (s_implicit_scalar_reflux) {
                // [Defer implicit refluxing until we're doing multilevel stuff]
            } else {
                for (int comp = 0; comp < s_num_scal_comps; ++comp) {
                    LevelData<FArrayBox>& newScalars = *m_scal_new[comp];

                    m_scal_fluxreg_ptrs[comp]->reflux(newScalars, *m_levGeoPtr);

                    MappedCoarseAverage& scalAvgDown = fineAMRNavierStokesPtr->m_coarse_average_scal;
                    scalAvgDown.averageToCoarse(newScalars,
                                                *fineAMRNavierStokesPtr->m_scal_new[comp],
                                                fineAMRNavierStokesPtr->m_levGeoPtr,
                                                considerCellVols);
                }
            }
        } else {
            // We are not refluxing. Just to the averaging.
            for (int comp = 0; comp < s_num_scal_comps; ++comp) {
                LevelData<FArrayBox>& newScalars = *m_scal_new[comp];

                MappedCoarseAverage& scalAvgDown = fineAMRNavierStokesPtr->m_coarse_average_scal;
                scalAvgDown.averageToCoarse(newScalars,
                                            *fineAMRNavierStokesPtr->m_scal_new[comp],
                                            fineAMRNavierStokesPtr->m_levGeoPtr,
                                            considerCellVols);
            }
        }
    } // end if !finestLevel


    // 2.) Do multi-level sync operations if this is the coarsest level
    //     or if coarser level is not at the same time as this level.
    const Real crseTime = (m_level > 0)? m_coarser_level_ptr->time(): -1.0e300;   // Was 0.0
    if (m_level == 0 || (abs(crseTime-m_time) > TIME_EPS)) {

        // 2.0) Find lbase. This will not be altered.
        int lBase = m_level;
        AMRNavierStokes* lBaseAMRNavierStokes = this;
        if (m_level > 0) {
            --lBase;
            lBaseAMRNavierStokes = crseNSPtr();
            CH_assert(lBaseAMRNavierStokes != NULL);
            CH_assert(lBaseAMRNavierStokes->m_level == lBase);
        }

        // 2.1) Find the finest level.
        int finest_level = m_level;
        AMRNavierStokes* thisNSPtr = this;
        while (!thisNSPtr->finestLevel()) {
            thisNSPtr = thisNSPtr->fineNSPtr();
            CH_assert(thisNSPtr != NULL);
        }
        CH_assert(thisNSPtr->finestLevel());
        finest_level = thisNSPtr->m_level;

        // 2.2) Do implicit refluxing for scalars, if necessary.
        if (!finestLevel()
            && s_implicit_scalar_reflux
            && (s_advective_scalar_reflux || s_diffusive_scalar_reflux))
        {
            this->doImplicitScalarReflux();
        }

        // 2.3) Construct composite fields for velocity
        Vector<LevelData<FArrayBox>*> compVel(finest_level+1, NULL);

        if (m_level > 0) {
            // Get coarser grids
            AMRNavierStokes* crseAMRNavierStokesPtr = crseNSPtr();
            DisjointBoxLayout crseGrids = crseAMRNavierStokesPtr->newVel().getBoxes();

            // Get coarser velocity
            LevelData<FArrayBox>* crseVelPtr = new LevelData<FArrayBox>(crseGrids, SpaceDim);
            crseAMRNavierStokesPtr->velocity(*crseVelPtr, m_time);
            compVel[m_level-1] = crseVelPtr;

            if (s_verbosity > 5) {
                pout() << "Collected vel at level " << crseAMRNavierStokesPtr->m_level
                       << " as base for sync projection" << endl;
            }
        }

        thisNSPtr = this;
        for (int lev = m_level; lev <= finest_level; ++lev) {
            compVel[lev] = thisNSPtr->m_vel_new_ptr;
            thisNSPtr = thisNSPtr->fineNSPtr();
        }

        // 2.4) Do implicit momentum refluxing, if necessary
        if (!finestLevel()
            && s_implicit_momentum_reflux
            && (s_advective_momentum_reflux || s_diffusive_momentum_reflux)
            && s_nu > 0.0)
        {
            this->doImplicitMomentumReflux(compVel);

            // Re-apply physical BCs after refluxing
            thisNSPtr = this;
            for (int lev = m_level; lev < finest_level; ++lev) {
                LevelData<FArrayBox>& levelVel = thisNSPtr->newVel();
                velBC.setGhosts(levelVel,
                                NULL,
                                thisNSPtr->m_levGeoPtr->getDx(),
                                &(thisNSPtr->m_levGeoPtr->getFCJgup()),
                                false,  // not homogeneous
                                m_time);

                thisNSPtr = thisNSPtr->fineNSPtr();
            }
        }

        // 2.5) Sync projection
        int projIters = (s_isIncompressible? s_sync_projection_iters: 0);
        {
            // Collect sync pressure pointers and set an initial guess of zero.
            Vector<LevelData<FArrayBox>*> eSync = lBaseAMRNavierStokes->gatherSyncPressure();
            setValLevels(eSync, m_level, eSync.size()-1, 0.0);

            // m_syncPressure on this and all higher levels will contain sync
            // projection data. Make sure lBase contains valid sync data too.
            setSyncPressureStates(SyncPressureState::SYNC);
            if (lBaseAMRNavierStokes->m_syncPressureState == SyncPressureState::INIT) {
                if (s_verbosity > 5) {
                    pout() << "Sync pressure at lbase = " << lBase << " is in state"
                           << " INIT. Setting to zero." << endl;
                }
                setValLevel(*eSync[lBase], 0.0);

            } else if (lBaseAMRNavierStokes->m_syncPressureState != SyncPressureState::SYNC) {
                ostringstream msg;
                msg << "Sync pressure at lBase = " << lBase << " is in state "
                    << lBaseAMRNavierStokes->m_syncPressureState;
                MayDay::Warning(msg.str().c_str());
                pout() << msg << endl;
                setValLevel(*eSync[lBase], 0.0);
            }

            // Sync projection is done from this level, up.
            if (s_isIncompressible && s_sync_projection_iters > 0) {
                pout() << "Sync projection on levels " << m_level << " to "
                       << eSync.size()-1 << ": " << std::flush;

                AMRCCProjector projObj;
                projObj.define(eSync, *m_physBCPtr, *m_levGeoPtr, NULL);

                for (int iter = 0; iter < s_sync_projection_iters; ++iter) {
                    projObj.project(compVel,
                                    *m_levGeoPtr,
                                    m_level,        // lmin
                                    eSync.size()-1, // lmax
                                    m_time,         // newTime
                                    m_dt,           // dt
                                    false,          // velIsFlux,
                                    false,          // zero-out pressure
                                    false);         // force homogeneous
                }
            }
        } // end sync projection

        // 2.6) Compute volume discrepancy correction
        this->computeVDCorrection(false);

        // 2.7) Re-apply physical BCs
        thisNSPtr = this;
        for (int lev = m_level; lev <= finest_level; ++lev) {       // Changed from < finest_level
            LevelData<FArrayBox>& levelVel = thisNSPtr->newVel();
            velBC.setGhosts(levelVel,
                            NULL,
                            thisNSPtr->m_levGeoPtr->getDx(),
                            &(thisNSPtr->m_levGeoPtr->getFCJgup()),
                            false,  // not homogeneous
                            m_time);

            thisNSPtr = thisNSPtr->fineNSPtr();
        }

        // 2.8) Do sync diagnostics
        if (!finestLevel() && m_level == 0) {
            this->syncBottomDiagnostics(compVel);
        } else if (finestLevel() && m_level == 0) {
            this->syncSingleGridDiagnostics();
        }

        // 2.9) Clean up storage
        if (m_level > 0) {
            if (compVel[m_level-1] != NULL) {
                delete compVel[m_level-1];
                compVel[m_level-1] = NULL;
            }
        }
    } // End multi-level, sync operations

    // 3.) Sync with LES
    // TODO: This may be better suited in tagCells on level 0 so
    // that we know not to tag the turbulent region.
    this->syncWithSGS();

    // 4.) Write terminal output
    this->syncTermDiagnostics();
}


// -----------------------------------------------------------------------------
// Do implicit scalar refluxing
// -----------------------------------------------------------------------------
void AMRNavierStokes::doImplicitScalarReflux ()
{
    if (s_verbosity >= 6) {
        pout() << "AMRNavierStokes::doImplicitScalarReflux" << endl;
    }
    CH_TIME("AMRNavierStokes::doImplicitScalarReflux");

    CH_assert(s_implicit_scalar_reflux);
    CH_assert(s_advective_scalar_reflux || s_diffusive_scalar_reflux);
    // CH_assert(s_nu > 0.0);

    const bool considerCellVols = false;

    // Find the coarsest level and collect its data
    AMRNavierStokes* thisNSPtr = this;
    while (thisNSPtr->m_coarser_level_ptr != NULL) {
        thisNSPtr = thisNSPtr->crseNSPtr();
    }
    const ProblemDomain& coarsestDomain = thisNSPtr->m_problem_domain;

    // Find the finest level
    int finest_level = m_level;
    thisNSPtr = this;
    while (!thisNSPtr->finestLevel()) {
        thisNSPtr = thisNSPtr->fineNSPtr();
    }
    CH_assert(thisNSPtr->finestLevel());
    finest_level = thisNSPtr->m_level;

    // Create containers for diffusive solve
    Vector<LevelData<FArrayBox>*> scalRefluxCorr(finest_level+1, NULL);
    Vector<LevelData<FArrayBox>*> scalRefluxRHS(finest_level+1, NULL);

    // Climb down levels and collect data. Collect everything down to m_level-1
    // if it exists. It will be used for setting BCs.
    int startLev;
    while (thisNSPtr != NULL && thisNSPtr->m_level >= m_level-1) {
        startLev = thisNSPtr->m_level;
        const DisjointBoxLayout& levelGrids = thisNSPtr->newVel().getBoxes();

        // Recall that AMRMultiGrid can only do one component.
        // RHS has no ghost cells
        // Soln will have one layer of ghost cells
        scalRefluxRHS[startLev] = new LevelData<FArrayBox>(levelGrids, 1);
        scalRefluxCorr[startLev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);

        // AmrJ[startLev] = (LevelData<FArrayBox>*)&(thisNSPtr->m_levGeoPtr->getCCJ());

        // Set RHS and corr to zero
        // We only need to do this for m_level-1 and finest_level -- the other
        // levels will be initialized in the next block of code.
        if (startLev == m_level-1 || startLev == finest_level) {
            LevelData<FArrayBox>& levelRefluxRHS = *(scalRefluxRHS[startLev]);
            LevelData<FArrayBox>& levelRefluxCorr = *(scalRefluxCorr[startLev]);

            DataIterator levelDit = scalRefluxRHS[startLev]->dataIterator();
            for (levelDit.reset(); levelDit.ok(); ++levelDit) {
                levelRefluxRHS[levelDit()].setVal(0.0);
                levelRefluxCorr[levelDit()].setVal(0.0);
            }
        }

        thisNSPtr = thisNSPtr->crseNSPtr();
    }

    // Now, do each scalar component
    const Interval solverComps(0,0);
    for (int scalComp = 0; scalComp < s_num_scal_comps; ++scalComp) {
        // Loop over levels and establish RHS, corr
        thisNSPtr = this;
        for (int lev = m_level; lev < finest_level; ++lev) {
            CH_assert(thisNSPtr->m_level == lev);

            LevelData<FArrayBox>& levelRefluxCorr = *(scalRefluxCorr[lev]);
            LevelData<FArrayBox>& levelRefluxRHS = *(scalRefluxRHS[lev]);
            MappedLevelFluxRegister& levelFR = *(thisNSPtr->m_scal_fluxreg_ptrs[scalComp]);

            // Initialize RHS to 0
            DataIterator levelDit = levelRefluxRHS.dataIterator();
            for (levelDit.reset(); levelDit.ok(); ++levelDit) {
                levelRefluxRHS[levelDit()].setVal(0.0);
            }

            // Do refluxing
            thisNSPtr->m_scal_fluxreg_ptrs[scalComp]->reflux(levelRefluxRHS, *(thisNSPtr->m_levGeoPtr));

            // Initial guess for correction is RHS
            levelRefluxRHS.copyTo(solverComps, levelRefluxCorr, solverComps);

            thisNSPtr = thisNSPtr->fineNSPtr();
        } // end loop over levels to compute RHS

        // Define levelOp and diffusive solver
        // We need to define new AMRMultiGrid for each component
        // because they may have different coefficients

        Real nuComp = s_scal_coeffs[scalComp];
        // Only do all of this if actually diffusive
        if (nuComp > 0) {

            // Define the diffusive solver
            const int numLevels = finest_level+1;
            const Real alpha = 1.0;
            const Real beta = -nuComp * m_dt;
            BCMethodHolder diffusiveBCs = m_physBCPtr->scalarRefluxSolveBC(scalComp);

            MappedAMRPoissonOpFactory diffusiveOpFactory;
            diffusiveOpFactory.define(m_levGeoPtr,
                                      alpha,
                                      beta,
                                      diffusiveBCs,
                                      s_viscous_AMRMG_maxDepth,
                                      s_viscous_AMRMG_num_smooth_precond,
                                      s_viscous_AMRMG_precondMode,
                                      s_viscous_AMRMG_relaxMode);

            BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
            bottomSolver.m_eps = s_viscous_bottom_eps;
            bottomSolver.m_reps = s_viscous_bottom_reps;
            bottomSolver.m_imax = s_viscous_bottom_imax;
            bottomSolver.m_numRestarts = s_viscous_bottom_numRestarts;
            bottomSolver.m_hang = s_viscous_bottom_hang;
            bottomSolver.m_small = s_viscous_bottom_small;
            bottomSolver.m_verbosity = s_viscous_bottom_verbosity;
            bottomSolver.m_normType = s_viscous_bottom_normType;

            typedef MappedAMRLevelOpFactory<LevelData<FArrayBox> > LevOpFact;
            MappedAMRMultiGrid<LevelData<FArrayBox> > diffusionSolver;
            LevOpFact& castFact = (LevOpFact&)diffusiveOpFactory;

            diffusionSolver.define(coarsestDomain,
                                   castFact,
                                   &bottomSolver,
                                   numLevels);

            diffusionSolver.m_verbosity = s_viscous_AMRMG_verbosity;
            diffusionSolver.m_imin = s_viscous_AMRMG_imin;
            diffusionSolver.setSolverParameters(s_viscous_AMRMG_num_smooth_down,
                                                s_viscous_AMRMG_num_smooth_up,
                                                s_viscous_AMRMG_num_smooth_bottom,
                                                s_viscous_AMRMG_numMG,
                                                s_viscous_AMRMG_imax,
                                                s_viscous_AMRMG_eps,
                                                s_viscous_AMRMG_hang,
                                                s_viscous_AMRMG_normThresh);

            // Quick check -- want sum(RHS) to be 0
#           ifndef _NDEBUG
                const Real sumRHS = computeMappedSum(scalRefluxRHS, *m_levGeoPtr, 0, m_level);
                pout() << "Sum(RHS) for implicit reflux for "
                       << s_scal_names[scalComp] << " = "
                       << setiosflags(ios::scientific) << setprecision(8)
                       << sumRHS << endl;
#           endif

            // Set time in all operators that will need to set BCs during solve.
            Vector<MappedMGLevelOp<LevelData<FArrayBox> > * > allOps = diffusionSolver.getAllOperators();
            for (int idx = 0; idx < allOps.size(); ++idx) {
                MappedAMRPoissonOp* thisOp = dynamic_cast<MappedAMRPoissonOp*>(allOps[idx]);
                thisOp->setTime(m_time);
            }

            // Solve!
            diffusionSolver.solve(scalRefluxCorr,
                                  scalRefluxRHS,
                                  finest_level,
                                  m_level,
                                  false, // don't initialize to zero
                                  true); // force homogeneous BCs

            // Increment NS pointer to finest level.
            thisNSPtr = this;
            while (!thisNSPtr->finestLevel()) {
                thisNSPtr = thisNSPtr->fineNSPtr();
            }
            CH_assert(thisNSPtr->m_level == finest_level);

            // Now increment scalars with reflux correction.
            // Go from finest->coarsest so that we can also avgDown.
            for (int lev = finest_level; lev >= m_level; --lev) {
                LevelData<FArrayBox>& levelScal = thisNSPtr->newScal(scalComp);
                LevelData<FArrayBox>& levelCorr = *(scalRefluxCorr[lev]);

                DataIterator levelDit = levelCorr.dataIterator();
                for (levelDit.reset(); levelDit.ok(); ++levelDit) {
                    levelScal[levelDit()] += levelCorr[levelDit()];
                }

                // If a finer level exists, do avgDown as well
                if (lev < finest_level) {
                    AMRNavierStokes& fineNS = *(thisNSPtr->fineNSPtr());
                    CH_assert(fineNS.m_level = lev+1);
                    LevelData<FArrayBox>& fineScal = fineNS.newScal(scalComp);

                    MappedCoarseAverage& scalAvgDown = fineNS.m_coarse_average_scal;
                    scalAvgDown.averageToCoarse(levelScal,
                                                fineScal,
                                                fineNS.m_levGeoPtr,
                                                considerCellVols);
                }

                thisNSPtr = thisNSPtr->crseNSPtr();
            } // end loop over levels

        } else { // End if diffusive, begin if not diffusive

            // If this component is not diffusive, just add RHS to scalar

            // Increment NS pointer to finest level
            thisNSPtr = this;
            while (!thisNSPtr->finestLevel()) {
                thisNSPtr = thisNSPtr->fineNSPtr();
            }
            CH_assert(thisNSPtr->m_level == finest_level);

            // Increment scalars with reflux correction
            // Go from finest->coarsest so that we can also avgDown
            for (int lev = finest_level; lev >= m_level; --lev) {
                LevelData<FArrayBox>& levelScal = thisNSPtr->newScal(scalComp);
                LevelData<FArrayBox>& levelCorr = *(scalRefluxRHS[lev]);

                DataIterator levelDit = levelCorr.dataIterator();
                for (levelDit.reset(); levelDit.ok(); ++levelDit) {
                    levelScal[levelDit()] += levelCorr[levelDit()];
                }

                // If a finer level exists, do avgDown as well
                if (lev < finest_level) {
                    AMRNavierStokes& fineNS = *(thisNSPtr->fineNSPtr());
                    CH_assert (fineNS.m_level == lev+1);
                    LevelData<FArrayBox>& fineScal = fineNS.newScal(scalComp);

                    MappedCoarseAverage& scalAvgDown = fineNS.m_coarse_average_scal;
                    scalAvgDown.averageToCoarse(levelScal,
                                                fineScal,
                                                fineNS.m_levGeoPtr,
                                                considerCellVols);
                }

                thisNSPtr = thisNSPtr->crseNSPtr();
            } // end loop over levels
        } // end if not diffusive
    } // end loop over scalar components

    // Clean up temporary scalar storage
    for (int lev = 0; lev <= finest_level; ++lev) {
        if (scalRefluxRHS[lev] != NULL) {
            delete scalRefluxRHS[lev];
            scalRefluxRHS[lev] = NULL;
        }

        if (scalRefluxCorr[lev] != NULL) {
            delete scalRefluxCorr[lev];
            scalRefluxCorr[lev] = NULL;
        }
    }
}



// -----------------------------------------------------------------------------
// Do implicit momentum refluxing (used in sync ops)
// -----------------------------------------------------------------------------
void AMRNavierStokes::doImplicitMomentumReflux (const Vector<LevelData<FArrayBox>*>& a_compVel)
{
    if (s_verbosity >= 6) {
        pout() << "AMRNavierStokes::doImplicitMomentumReflux" << endl;
    }
    CH_TIME("AMRNavierStokes::doImplicitMomentumReflux");

    CH_assert(s_implicit_momentum_reflux);
    CH_assert(s_advective_momentum_reflux || s_diffusive_momentum_reflux);

    // Find the coarsest level and collect its data
    AMRNavierStokes* thisNSPtr = this;
    while (thisNSPtr->m_coarser_level_ptr != NULL) {
        thisNSPtr = thisNSPtr->crseNSPtr();
    }
    const ProblemDomain& coarsestDomain = thisNSPtr->m_problem_domain;

    // Find the finest level
    int finest_level = m_level;
    thisNSPtr = this;
    while (!thisNSPtr->finestLevel()) {
        thisNSPtr = thisNSPtr->fineNSPtr();
    }
    CH_assert(thisNSPtr->finestLevel());
    finest_level = thisNSPtr->m_level;
    AMRNavierStokes* finestNSPtr = thisNSPtr;

    // Loop over levels and compute RHS
    Vector<LevelData<FArrayBox>*> refluxRHS(finest_level+1,NULL);
    Vector<LevelData<FArrayBox>*> refluxCorr(finest_level+1,NULL);
    Vector<LevelData<FArrayBox>*> tempRefluxData(finest_level+1,NULL);

    thisNSPtr = this;
    int startLev = m_level;
    // if crser level exists, define it as well for BC's
    if (startLev > 0) {
      startLev = startLev-1;
      thisNSPtr = thisNSPtr->crseNSPtr();
    }

    for (int lev = startLev; lev <= finest_level; ++lev) {
        CH_assert(thisNSPtr->m_level == lev);
        const DisjointBoxLayout& levelGrids = thisNSPtr->newVel().getBoxes();

        // Recall that AMRMultiGrid can only do one component
        // RHS has no ghost cells
        // Soln will have one layer of ghost cells
        refluxRHS[lev] = new LevelData<FArrayBox>(levelGrids, 1);
        tempRefluxData[lev] = new LevelData<FArrayBox>(levelGrids, SpaceDim);
        refluxCorr[lev] = new LevelData<FArrayBox>(levelGrids, 1, IntVect::Unit);

        // Initialize rhs to 0
        DataIterator levelDit = tempRefluxData[lev]->dataIterator();
        LevelData<FArrayBox>& levelRefluxData = *(tempRefluxData[lev]);
        for (levelDit.reset(); levelDit.ok(); ++levelDit) {
            levelRefluxData[levelDit()].setVal(0.0);
        }

        // While we're here, do refluxing.
        // (Recall that startLev may be coarser than m_level for BC's,
        // however, don't want to do anything to that level)
        if  ((lev >= m_level) && (lev < finest_level)) {
            thisNSPtr->m_vel_flux_reg.reflux(levelRefluxData, *(thisNSPtr->m_levGeoPtr));
        }

        thisNSPtr = thisNSPtr->fineNSPtr();
    } // end loop over levels

    // coarse BC is 0
    if (m_level > 0) {
        LevelData<FArrayBox>& crseBCData = *refluxCorr[m_level-1];
        DataIterator crseDit = crseBCData.dataIterator();
        for (crseDit.reset(); crseDit.ok(); ++crseDit) {
            crseBCData[crseDit()].setVal(0.0);
        }
    }

    const Interval solverComps(0,0);
    for (int dir = 0; dir < SpaceDim; ++dir) {
        const Interval velComps(dir,dir);
        for (int lev = m_level; lev <= finest_level; ++lev) {
            // Copy rhs to single-component holder
            tempRefluxData[lev]->copyTo(velComps, *(refluxRHS[lev]), solverComps);

            // Initial guess for correction is RHS.
            LevelData<FArrayBox>& levelCorr = *(refluxCorr[lev]);
            DataIterator levelDit = levelCorr.dataIterator();
            for (levelDit.reset(); levelDit.ok(); ++levelDit) {
                levelCorr[levelDit()].setVal(0.0);
            }
            refluxRHS[lev]->copyTo(solverComps, *(refluxCorr[lev]), solverComps);
        } // end set initial guess

        // Define the viscous solver
        const int numLevels = finest_level+1;
        const Real alpha = 1.0;
        const Real beta = -s_nu * m_dt;
        BCMethodHolder viscousBCs = m_physBCPtr->viscousRefluxBC(dir);

        // TODO: This is an expensive solver construction. We should instead
        // create one solver and repurpose it for each solve. The only parameter
        // that varies is the BC.
        TODO();

        MappedAMRPoissonOpFactory viscousOpFactory;
        viscousOpFactory.define(&(*m_levGeoPtr),
                                alpha,
                                beta,
                                viscousBCs,
                                s_viscous_AMRMG_maxDepth,
                                s_viscous_AMRMG_num_smooth_precond,
                                s_viscous_AMRMG_precondMode,
                                s_viscous_AMRMG_relaxMode);

        BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
        bottomSolver.m_eps = s_viscous_bottom_eps;
        bottomSolver.m_reps = s_viscous_bottom_reps;
        bottomSolver.m_imax = s_viscous_bottom_imax;
        bottomSolver.m_numRestarts = s_viscous_bottom_numRestarts;
        bottomSolver.m_hang = s_viscous_bottom_hang;
        bottomSolver.m_small = s_viscous_bottom_small;
        bottomSolver.m_verbosity = s_viscous_bottom_verbosity;
        bottomSolver.m_normType = s_viscous_bottom_normType;

        MappedAMRMultiGrid<LevelData<FArrayBox> > viscousSolver;
        MappedAMRLevelOpFactory<LevelData<FArrayBox> >& viscCastFact
            = (MappedAMRLevelOpFactory<LevelData<FArrayBox> >&) viscousOpFactory;

        viscousSolver.define(coarsestDomain,
                             viscCastFact,
                             &bottomSolver,
                             numLevels);

        viscousSolver.m_verbosity = s_viscous_AMRMG_verbosity;
        viscousSolver.m_imin = s_viscous_AMRMG_imin;
        viscousSolver.setSolverParameters(s_viscous_AMRMG_num_smooth_down,
                                          s_viscous_AMRMG_num_smooth_up,
                                          s_viscous_AMRMG_num_smooth_bottom,
                                          s_viscous_AMRMG_numMG,
                                          s_viscous_AMRMG_imax,
                                          s_viscous_AMRMG_eps,
                                          s_viscous_AMRMG_hang,
                                          s_viscous_AMRMG_normThresh);

        // Set time in all operators that will need to set BCs during solve.
        Vector<MappedMGLevelOp<LevelData<FArrayBox> > * > allOps = viscousSolver.getAllOperators();
        for (int idx = 0; idx < allOps.size(); ++idx) {
            MappedAMRPoissonOp* thisOp = dynamic_cast<MappedAMRPoissonOp*>(allOps[idx]);
            thisOp->setTime(m_time);
        }

        pout() << "Implicit reflux for velocity comp " << dir << endl;

        // Solve!
        viscousSolver.solve(refluxCorr,
                            refluxRHS,
                            finest_level,
                            m_level,
                            false, // don't initialize to zero
                            true); // force homogeneous BCs

        // Copy correction to a holder with SpaceDim comps.
        thisNSPtr = this;
        for (int lev = m_level; lev <= finest_level; ++lev) {
            // Sanity check
            CH_assert(thisNSPtr != NULL);
            CH_assert(thisNSPtr->m_level == lev);

            // Create references for convenience
            LevelData<FArrayBox>& levelRefluxData = (*tempRefluxData[lev]);
            const LevelData<FArrayBox>& levelCorr = (*refluxCorr[lev]);
            DataIterator levelDit = levelCorr.dataIterator();

            // Perform the copy
            for (levelDit.reset(); levelDit.ok(); ++levelDit) {
                const DataIndex& di = levelDit();
                levelRefluxData[di].copy(levelCorr[di], 0, dir, 1);
            }

            // Move on to the next finer level
            thisNSPtr = thisNSPtr->fineNSPtr();
        } // end loop over levels
    } // end loop over directions

    // We are done with refluxRHS and refluxCorr.
    // The solution has been copied to tempRefluxData.
    for (int lev = startLev; lev <= finest_level; ++lev) {
        if (refluxRHS[lev] != NULL) {
            delete refluxRHS[lev];
            refluxRHS[lev] = NULL;
        }

        if (refluxCorr[lev] != NULL) {
            delete refluxCorr[lev];
            refluxCorr[lev] = NULL;
        }
    }

    // Increment velocity with reflux correction
    thisNSPtr = this;
    for (int lev = m_level; lev <= finest_level; ++lev) {
        // Sanity check
        CH_assert(thisNSPtr != NULL);
        CH_assert(thisNSPtr->m_level == lev);

        // Create references for convenience
        const LevelGeometry& levelLevGeo = *(thisNSPtr->m_levGeoPtr);
        LevelData<FArrayBox>& levelRefluxData = *(tempRefluxData[lev]);
        LevelData<FArrayBox>& levelVel = *(a_compVel[lev]);
        DataIterator levelDit = levelVel.dataIterator();

        // Send the correction to a mapped basis
        levelLevGeo.sendToMappedBasis(levelRefluxData, true); // TODO: Why are we using ghosts???

        // Add correction to velocity
        for (levelDit.reset(); levelDit.ok(); ++levelDit) {
            const DataIndex& di = levelDit();
            levelVel[di].plus(levelRefluxData[di], 0, 0, SpaceDim);
        }

        // Move on to the next finer level
        thisNSPtr = thisNSPtr->fineNSPtr();
    }

    // We should not use tempRefluxData anymore since we clobbered the data.
    for (int lev = startLev; lev <= finest_level; ++lev) {
        if (tempRefluxData[lev] != NULL) {
            delete tempRefluxData[lev];
            tempRefluxData[lev] = NULL;
        }
    }

    // Averaging down the levels will be done in doSyncOperations.
}


// -----------------------------------------------------------------------------
// Compute the volume discrepancy correction.
// Set a_init to true during initialization VD solves and false for sync solves.
// -----------------------------------------------------------------------------
void AMRNavierStokes::computeVDCorrection (const bool a_init)
{
    // Should we really do this?
    if (s_etaLambda == 0.0) return;

    // Find the finest level.
    AMRNavierStokes* thisNSPtr = this;
    while (!thisNSPtr->finestLevel()) {
        thisNSPtr = thisNSPtr->fineNSPtr();
    }
    CH_assert(thisNSPtr->finestLevel());
    const int finest_level = thisNSPtr->m_level;


    // Set up RHS. (We use lambda as a temp holder.)
    Vector<LevelData<FArrayBox>*> compRHS(finest_level+1, NULL);
    const Real lambdaMult = s_etaLambda / m_dt;

    thisNSPtr = this;
    for (int lev = m_level; lev <= finest_level; ++lev) {
        compRHS[lev] = thisNSPtr->m_lambda_new_ptr;

        CH_assert(m_dt > 0.0);
        CH_assert(s_etaLambda > 0.0);

        DataIterator levelDit = compRHS[lev]->dataIterator();
        for (levelDit.reset(); levelDit.ok(); ++levelDit) {
            FArrayBox& rhsFAB = (*compRHS[lev])[levelDit];

            rhsFAB -= 1.0;
            rhsFAB *= lambdaMult;
        }

        thisNSPtr = thisNSPtr->fineNSPtr();
    }

    // For debugging -- should be zero
    if (s_verbosity >= 3) {
        Real vol = 0.0;
        Real sumRHS = computeMappedSum(vol, compRHS, *m_levGeoPtr, 0, m_level);
        pout() << "\nVD solve..."
               << setiosflags(ios::scientific) << setprecision(8)
               << "\nSum(RHS) = " << sumRHS
               << endl;
    }

    // Collect eLambda pointers. These are the solution holders.
    Vector<LevelData<FArrayBox>*> compELambda = gatherELambda();
    CH_assert(compELambda.size() == finest_level + 1);
    CH_assert(m_level == 0 || m_eLambdaState == ELambdaState::VALID);

    // If this is an initialization VD solve, the CFBCs are zero.
    // NOTE: Chombo just sets the coarse level eLambda to zero directly.
    if (m_level > 0) {
        setValLevel(*compELambda[m_level-1], 0.0);
    }

    // Set initial guess.
    setValLevels(compELambda, m_level, finest_level, 0.0);

    // Define the VD solver
    AMRPressureSolver solver;
    {
        const ProblemContext* ctx = ProblemContext::getInstance();

        solver.setAMRMGParameters(ctx->VD_AMRMG_imin,
                                  ctx->VD_AMRMG_imax,
                                  ctx->VD_AMRMG_eps,
                                  ctx->VD_AMRMG_maxDepth,
                                  ctx->VD_AMRMG_num_smooth_precond,
                                  ctx->VD_AMRMG_num_smooth_down,
                                  ctx->VD_AMRMG_num_smooth_up,
                                  ctx->VD_AMRMG_num_smooth_bottom,
                                  ctx->VD_AMRMG_precondMode,
                                  ctx->VD_AMRMG_relaxMode,
                                  ctx->VD_AMRMG_numMG,
                                  ctx->VD_AMRMG_hang,
                                  ctx->VD_AMRMG_normThresh,
                                  ctx->VD_AMRMG_verbosity);

        solver.setBottomParameters(ctx->VD_bottom_imax,
                                   ctx->VD_bottom_numRestarts,
                                   ctx->VD_bottom_eps,
                                   ctx->VD_bottom_reps,
                                   ctx->VD_bottom_hang,
                                   ctx->VD_bottom_small,
                                   ctx->VD_bottom_normType,
                                   ctx->VD_bottom_verbosity);

        Vector<const LevelGeometry*> compLevGeos
                = ((const LevelGeometry*)m_levGeoPtr)->getAMRLevGeos();
        const Box& lminDomBox = compLevGeos[0]->getDomain().domainBox();
        const int numLevels = compELambda.size();
        BCMethodHolder bcHolder = m_physBCPtr->FreestreamCorrFuncBC();

        solver.define(bcHolder, *m_levGeoPtr, lminDomBox, numLevels, NULL);
    } // end define VD solver

    // Solve!
    if (s_verbosity >= 1) {
        if (a_init) {
            pout() << "init VD solve on levels " << m_level << " to " << finest_level << ": " << flush;
        } else {
            pout() << "sync VD solve on levels " << m_level << " to " << finest_level << ": " << flush;
        }
    }
    solver.solve(compELambda, compRHS, m_level, finest_level, false, false);

    // Apply all BC's here
    {
        BCMethodHolder bcHolder = m_physBCPtr->gradELambdaFuncBC();

        thisNSPtr = this;
        for (int lev = m_level; lev <= finest_level; ++lev) {
            CH_assert(thisNSPtr->m_level == lev);

            LevelData<FArrayBox>& levelELambda = (*compELambda[lev]);
            const DisjointBoxLayout& levelGrids = levelELambda.getBoxes();
            const ProblemDomain& levelDomain = levelGrids.physDomain();
            const RealVect& levelDx = thisNSPtr->m_levGeoPtr->getDx();

            DataIterator ditLev = levelELambda.dataIterator(); // TODO: May not be needed.
            for (ditLev.reset(); ditLev.ok(); ++ditLev) {
                bcHolder.setGhosts(levelELambda[ditLev],// stateFAB
                                   NULL,                // &extrapFAB
                                   levelGrids[ditLev],  // valid box
                                   levelDomain,         // ProblemDomain
                                   levelDx,             // dx
                                   ditLev(),            // DataIndex
                                   NULL,                // &JgupFAB
                                   false,               // isHomogeneous
                                   m_time);             // time
            }

            Copier exc(levelGrids, levelGrids, levelDomain, IntVect::Unit, true);
            levelELambda.exchange(exc);

            CornerCopier excc(levelGrids, levelGrids, levelDomain, IntVect::Unit, true);
            levelELambda.exchange(excc);

            thisNSPtr = thisNSPtr->fineNSPtr();
        }
    }

    // Fill gradELambda from the finest level, down.
    thisNSPtr = this;
    while (thisNSPtr->m_level < finest_level) {
        thisNSPtr = thisNSPtr->fineNSPtr();
        CH_assert(thisNSPtr != NULL);
    }

    for (signed int lev = finest_level; lev >= m_level; --lev) {
        CH_assert(thisNSPtr->m_level == lev);

        // Compute gradient.
        BCMethodHolder fluxBC = thisNSPtr->m_physBCPtr->gradELambdaFuncBC();
        const LevelData<FArrayBox>* crseELambdaPtr = ((lev > 0)? compELambda[lev-1]: NULL);

        // Recall that composite MAC gradient is the same as the level
        // gradient, since finer level is not considered to be part
        // of this level (take care of covered regions with avgDown)
        Gradient::levelGradientMAC(thisNSPtr->m_gradELambda,
                                   *compELambda[lev],
                                   crseELambdaPtr,
                                   *(thisNSPtr->m_levGeoPtr),
                                   m_time,
                                   &fluxBC);

        // Average gradient down from finer level.
        if (lev < finest_level) {
            const AMRNavierStokes* fineLevelPtr = thisNSPtr->fineNSPtr();
            const LevelGeometry* fineLevGeoPtr = fineLevelPtr->m_levGeoPtr;

            const DisjointBoxLayout& fineGrids = fineLevGeoPtr->getBoxes();
            const IntVect& nRefFine = fineLevGeoPtr->getCrseRefRatio();
            const int nComp = 1;
            const LevelData<FluxBox>& fineEdgeGrad = fineLevelPtr->m_gradELambda;

            MappedCoarseAverageFace avgDownObj(fineGrids, nComp, nRefFine);
            avgDownObj.averageToCoarse(thisNSPtr->m_gradELambda, fineEdgeGrad);
        }

        // Exchange.
        const DisjointBoxLayout& levelGrids = thisNSPtr->m_gradELambda.getBoxes();
        const ProblemDomain& levelDomain = levelGrids.physDomain();

        Copier exc(levelGrids, levelGrids, levelDomain, IntVect::Unit, true);
        thisNSPtr->m_gradELambda.exchange(exc);

        CornerCopier excc(levelGrids, levelGrids, levelDomain, IntVect::Unit, true);
        thisNSPtr->m_gradELambda.exchange(excc);

        // Move down one level
        thisNSPtr = thisNSPtr->crseNSPtr();
    }

    // Restore lambda
    for (int lev = m_level; lev <= finest_level; ++lev) {
        DataIterator levelDit = compRHS[lev]->dataIterator();
        for (levelDit.reset(); levelDit.ok(); ++levelDit) {
            FArrayBox& lambdaFAB = (*compRHS[lev])[levelDit];

            lambdaFAB /= lambdaMult;
            lambdaFAB += 1.0;
        }
    }

    // Set the states
    setELambdaStates(ELambdaState::VALID);
    setGradELambdaStates(GradELambdaState::VALID);
}


// -----------------------------------------------------------------------------
// This can be used to monitor max(divergence), sum(KE), etc.
// -----------------------------------------------------------------------------
void AMRNavierStokes::syncBottomDiagnostics(const Vector<LevelData<FArrayBox>*>& a_compVel)
{
    // Put anything you like here.
}


// -----------------------------------------------------------------------------
// This is a diagnostic -- compute max(divergence)
// -----------------------------------------------------------------------------
void AMRNavierStokes::syncSingleGridDiagnostics()
{
    if (s_verbosity >= 6) {
        pout() << "AMRNavierStokes::syncSingleGridDiagnostics" << endl;
    }
    CH_TIME("AMRNavierStokes::syncSingleGridDiagnostics");

    CH_assert(m_level == 0);

    // Collect some needed data
    const RealVect& dx = m_levGeoPtr->getDx();
    const DisjointBoxLayout& grids = newVel().getBoxes();
    DataIterator dit = grids.dataIterator();
    const bool isViscous = (s_nu > 0.0);
    LevelData<FArrayBox>& vel = newVel();

    // Exchange ghosts
    vel.exchange(m_oneGhostExCopier);

    // Set physical BC's
    VelBCHolder velBC(m_physBCPtr->uStarFuncBC(isViscous));
    velBC.setGhosts(vel,
                    NULL,
                    dx,
                    &(m_levGeoPtr->getFCJgup()),
                    false,  // not homogeneous
                    m_time);

    // Compute divergence
    LevelData<FArrayBox> thisDiv(grids, 1);
    Divergence::compDivergenceCC(thisDiv,
                                 newVel(),
                                 NULL,      // crseVelPtr
                                 NULL,      // fineVelPtr
                                 true,
                                 *m_levGeoPtr);

    // Compute max|Divergence|
    const Real maxDiv = computeMappedNorm(thisDiv, NULL, *m_levGeoPtr,
                                          0);  // norm type

    pout() << setiosflags(ios::scientific) << setprecision(15);
    pout() << "Time = " << setw(15) << m_time
           << setw(30) << " Max Div(u) = "
           << setw(23)  << maxDiv << endl;
    pout() << setiosflags(ios::scientific) << setprecision(8) << std::flush;
}


// -----------------------------------------------------------------------------
// Write diagnostic info to stdout
// -----------------------------------------------------------------------------
void AMRNavierStokes::syncTermDiagnostics()
{
    if (s_verbosity >= 6) {
        pout() << "AMRNavierStokes::syncTermDiagnostics" << endl;
    }
    CH_TIME("AMRNavierStokes::syncTermDiagnostics");

    if (!s_write_stdout) return;

    // Write level output
    {
        // Get reference to finer grids
        const DisjointBoxLayout* fineGridsPtr = NULL;
        const LevelGeometry* fineLevGeoPtr = m_levGeoPtr->getFinerPtr();
        if (fineLevGeoPtr != NULL) {
            fineGridsPtr = &(fineLevGeoPtr->getBoxes());
        }

        // Compute max|velocity|
        D_TERM(
        const Real unorm = computeUnmappedNorm(*m_vel_new_ptr, fineGridsPtr, *m_levGeoPtr, 0, 0);,
        const Real vnorm = computeUnmappedNorm(*m_vel_new_ptr, fineGridsPtr, *m_levGeoPtr, 0, 1);,
        const Real wnorm = computeUnmappedNorm(*m_vel_new_ptr, fineGridsPtr, *m_levGeoPtr, 0, 2);)

        // Compute max|buoyancy|
        Real bnorm = 0.0;
        if (s_num_scal_comps > 0) {
            bnorm = computeUnmappedNorm(*m_scal_new[0], fineGridsPtr, *m_levGeoPtr, 0, 0);
        }

        // Print results to terminal
        if (procID() == 0) {
            const Real tol = 1e8;
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
    if (m_level == 0) {
        // Collect AMR data
        Vector<LevelData<FArrayBox>*> amrVel(0), amrB(0);
        AMRNavierStokes* levelNSPtr = this;
        while(levelNSPtr != NULL) {
            amrVel.push_back(levelNSPtr->m_vel_new_ptr);
            if (s_num_scal_comps > 0) amrB.push_back(levelNSPtr->m_scal_new[0]);

            if (levelNSPtr->m_finer_level_ptr != NULL) {
                levelNSPtr = dynamic_cast<AMRNavierStokes*>(levelNSPtr->m_finer_level_ptr);
            } else {
                levelNSPtr = NULL;
            }
        }

        // Compute max|velocity|
        D_TERM(
        const Real unorm = computeUnmappedNorm(amrVel, *m_levGeoPtr, 0, 0);,
        const Real vnorm = computeUnmappedNorm(amrVel, *m_levGeoPtr, 0, 1);,
        const Real wnorm = computeUnmappedNorm(amrVel, *m_levGeoPtr, 0, 2);)

        // Compute max|buoyancy|
        Real bnorm = 0.0;
        if (s_num_scal_comps > 0) {
            bnorm = computeUnmappedNorm(amrB, *m_levGeoPtr, 0, 0);
        }

        // Compute total mass
        Real mass = 0.0;
        if (s_num_scal_comps > 0) {
            mass = computeMappedSum(amrB, *m_levGeoPtr);
        }

        // Compute total momenta
        const RealVect totalMom(D_DECL(
            computeMappedSum(amrVel, *m_levGeoPtr, 0),
            computeMappedSum(amrVel, *m_levGeoPtr, 1),
            computeMappedSum(amrVel, *m_levGeoPtr, 2)
        ));

        // Compute total energy
        const Real globalEnergy = this->totalEnergy();

        if (procID() == 0) {
            const Real tol = 1e8;

            std::cout << color::hiwhite << std::left << setw(6) << ++s_step_number
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
                      << ((s_totalEnergy < globalEnergy)? color::red: color::hiblue)
                      << std::left << setiosflags(ios::scientific) << setprecision(8) << setw(18) << globalEnergy
                      << color::none
                      << "\n"
                      << std::endl;
        }
        s_totalEnergy = globalEnergy;
    }
}


// -----------------------------------------------------------------------------
// Provides syncing with a subgrid scale model.
// This function does nothing by default. Feel free to add whatever code you
// like, but future versions of SOMAR will use this function to update the
// stress tensor via LES.
// -----------------------------------------------------------------------------
void AMRNavierStokes::syncWithSGS ()
{
    return;
}
