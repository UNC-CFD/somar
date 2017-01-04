/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2014 University of North Carolina at Chapel Hill
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
#include "SetValLevel.H"
#include "CellToEdge.H"
#include "EdgeToCell.H"
#include "ConvertFABF_F.H" // For 4th-order CellToEdge
#include "Gradient.H"
#include "timeInterp.H"
#include "ProblemContext.H"
#include "Debug.H"



// -----------------------------------------------------------------------------
// Erases the data in oldPtr and swaps the pointers.
// -----------------------------------------------------------------------------
void AMRNavierStokes::swapPointers (LevelData<FArrayBox>*& a_newPtr,
                                    LevelData<FArrayBox>*& a_oldPtr)
{
    if (s_set_bogus_values) {
        setValLevel(*a_oldPtr, s_bogus_value);
    }

    LevelData<FArrayBox>* tmpPtr = a_oldPtr;
    a_oldPtr = a_newPtr;
    a_newPtr = tmpPtr;
}


// -----------------------------------------------------------------------------
// Computes the source terms for a Runge-Kutta update.
// This assumes only one scalar (not counting lambda).
// -----------------------------------------------------------------------------
void AMRNavierStokes::RK3TimeStep (const Real a_oldTime,
                                   const Real a_dt,
                                   const bool a_updatePassiveScalars,
                                   const bool a_doLevelProj)
{
    CH_TIME("AMRNavierStokes::RK3TimeStep");

    // Sanity checks
    CH_assert(m_levGeoPtr->getDomain() == m_problem_domain);
    CH_assert(Abs(a_oldTime - (m_time - a_dt)) < TIME_EPS);
    CH_assert(s_num_scal_comps == 1);

    const DisjointBoxLayout& grids = newVel().getBoxes();
    CH_assert(grids.compatible(m_levGeoPtr->getBoxes()));

    LevelData<FArrayBox>& u = newVel();     // The current velocity state
    LevelData<FArrayBox>& b = newScal(0);   // The current buoyancy state
    LevelData<FArrayBox>* newSuPtr = NULL;  // The velocity eq's advective source + external forces
    LevelData<FArrayBox>* newSbPtr = NULL;  // The buoyancy eq's advective source + external forces
    LevelData<FArrayBox>* oldSuPtr = NULL;  // Ditto, but from the previous RK stage.
    LevelData<FArrayBox>* oldSbPtr = NULL;  // Ditto, but from the previous RK stage.
    Real stateTime = a_oldTime;             // The time the current state (u,b) is at.
    int RKStage;                            // The RK stage.

    // Initialize state
    this->fillVelocity(u, stateTime);
    u.exchange();
    {
        LevelData<FArrayBox> tmpb;
        fillScalars(tmpb, stateTime, 0);

        DataIterator dit = grids.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            b[dit].copy(tmpb[dit]);
        }

        b.exchange();
    }
    checkForValidNAN(u);
    checkForValidNAN(b);


    // Initialize sources
    newSuPtr = new LevelData<FArrayBox>(grids, SpaceDim, IntVect::Zero);
    newSbPtr = new LevelData<FArrayBox>(grids, 1       , IntVect::Zero);
    oldSuPtr = new LevelData<FArrayBox>(grids, SpaceDim, IntVect::Zero);
    oldSbPtr = new LevelData<FArrayBox>(grids, 1       , IntVect::Zero);
    if (s_set_bogus_values) {
        setValLevel(*newSuPtr, s_bogus_value);
        setValLevel(*newSbPtr, s_bogus_value);
        setValLevel(*oldSuPtr, s_bogus_value);
        setValLevel(*oldSbPtr, s_bogus_value);
    }

    // Stage 1...
    RKStage = 1;
    this->computeMOLSources(*newSuPtr, *newSbPtr, u, b, stateTime, RKStage);                    // Compute (newSu,newSb) given (u,b).
    this->updateState(u, b, stateTime, *newSuPtr, *newSbPtr, *oldSuPtr, *oldSbPtr, RKStage);    // Get new (u,b) and update stateTime.
    this->swapPointers(newSuPtr, oldSuPtr);                                                     // Clear old data and swap pointers.
    this->swapPointers(newSbPtr, oldSbPtr);

    // Stage 2...
    RKStage = 2;
    this->computeMOLSources(*newSuPtr, *newSbPtr, u, b, stateTime, RKStage);                    // Compute (newSu,newSb) given (u,b).
    this->updateState(u, b, stateTime, *newSuPtr, *newSbPtr, *oldSuPtr, *oldSbPtr, RKStage);    // Get new (u,b) and update stateTime.
    this->swapPointers(newSuPtr, oldSuPtr);                                                     // Clear old data and swap pointers.
    this->swapPointers(newSbPtr, oldSbPtr);

    // Stage 3...
    RKStage = 3;
    this->computeMOLSources(*newSuPtr, *newSbPtr, u, b, stateTime, RKStage);                    // Compute (newSu,newSb) given (u,b).
    this->updateState(u, b, stateTime, *newSuPtr, *newSbPtr, *oldSuPtr, *oldSbPtr, RKStage);    // Get new (u,b) and update stateTime.
    delete oldSbPtr; oldSbPtr = NULL;
    delete oldSuPtr; oldSuPtr = NULL;
    delete newSbPtr; newSbPtr = NULL;
    delete newSuPtr; newSuPtr = NULL;

    CH_assert(Abs(stateTime - m_time) < TIME_EPS);

    // CC pressure field now contains valid data.
    m_ccPressureState = CCPressureState::VALID;
}


// -----------------------------------------------------------------------------
// Computes
// Su = advective source + external force. Returned in the Cartesian basis.
// Sb = advective source + external force.
// -----------------------------------------------------------------------------
void AMRNavierStokes::computeMOLSources (LevelData<FArrayBox>& a_Su,
                                         LevelData<FArrayBox>& a_Sb,
                                         LevelData<FArrayBox>& a_u,
                                         LevelData<FArrayBox>& a_b,
                                         const Real            a_stateTime,
                                         const int             a_stage)
{
    CH_TIME("AMRNavierStokes::computeMOLSources");

    // Sanity checks
    CH_assert(a_stage == 1 || a_stage == 2 || a_stage == 3);
    CH_assert(m_levGeoPtr->getBoxes() == a_u.getBoxes());
    CH_assert(m_levGeoPtr->getDomain() == m_problem_domain);
    CH_assert(m_time - m_dt - TIME_EPS < a_stateTime);
    CH_assert(a_stateTime < m_time + TIME_EPS);
    CH_assert(s_num_scal_comps == 1);

    // Set up some basic values
    const RealVect& dx = m_levGeoPtr->getDx();
    const DisjointBoxLayout& grids = a_u.getBoxes();
    DataIterator dit = grids.dataIterator();
    const Box domainBox = m_problem_domain.domainBox();

    const bool isDiffusive = (s_scal_coeffs[0] > 0.0);
    const bool isViscous = (s_nu > 0.0);

    // Set RK3 stage coefficients.
    Real h;     // The effective stage dt.
    switch (a_stage) {
    case 1:
        h = m_dt * 8. / 15.;
        break;
    case 2:
        h = m_dt * 2. / 15.;
        break;
    case 3:
        h = m_dt * 1. / 3.;
        break;
    default:
        MayDay::Error("Bad RK3 stage");
    }

    // Set Su and Sb to bogus values
    if (s_set_bogus_values) {
        setValLevel(a_Su, s_bogus_value);
        setValLevel(a_Sb, s_bogus_value);
    }

    // Set all ghosts on (u,b)
    {
        VelBCHolder velBC(m_physBCPtr->tracingVelFuncBC());
        this->setGhostsVelocity(a_u, velBC, a_stateTime);
        a_u.exchange();

        BCMethodHolder scalBC = m_physBCPtr->scalarTraceFuncBC(0);
        this->setGhostsScalar(a_b, scalBC, a_stateTime, 0);
        a_b.exchange();
    }


    // Advecting velocity ------------------------------------------------------
    LevelData<FluxBox> uAD(grids, 1, IntVect::Unit);
    {
        if (s_set_bogus_values) {
            setValLevel(uAD, s_bogus_value);
        }

        // Send u to FC
        CellToEdge(a_u, uAD);

        // Set BCs
        EdgeVelBCHolder edgeVelBC(m_physBCPtr->advectingVelFuncBC(isViscous));
        edgeVelBC.setGhosts(uAD,
                            NULL,
                            dx,
                            &(m_levGeoPtr->getFCJgup()),
                            false, // inhomogeneous
                            a_stateTime);

        // MAC project
        if (s_isIncompressible) {
            Vector<LevelData<FluxBox>*> amrVel(1, &uAD);

            pout() << "Level " << m_level << " MAC proj: " << flush;
            m_macProjector.levelProject(amrVel,
                                        m_levGeoPtr,
                                        a_stateTime,
                                        h,
                                        false,       // a_advVel is not a flux
                                        true,        // isLevelSolve
                                        false);      // forceHomogSolve
        }

        // Set BCs
        uAD.exchange();
        edgeVelBC.setGhosts(uAD,
                            NULL,
                            dx,
                            &(m_levGeoPtr->getFCJgup()),
                            false,  // inhomogeneous
                            a_stateTime);

        // Scale as a flux
        m_levGeoPtr->multByJ(uAD);

        checkForValidNAN(uAD);
    }


    // Buoyancy ----------------------------------------------------------------
    // Compute advective source term...
    {
        const bool useFourthOrder = false;

        // Compute -flux.
        LevelData<FluxBox> bFlux(grids, 1);
        if (s_set_bogus_values) {
            setValLevel(bFlux, s_bogus_value);
        }

        if (useFourthOrder) {
            // Use 4th-order interpolation from CC to FC.
            // TODO: This should be its own function.

            // Copy the CC data and fill 2 exchange/CFBC ghosts
            // and 1 physical boundary ghost.
            LevelData<FArrayBox> bGrow(grids, 1, 2*IntVect::Unit);
            if (s_set_bogus_values) {
                setValLevel(bGrow, s_bogus_value);
            }
            a_b.copyTo(bGrow);

            BCMethodHolder scalBC = m_physBCPtr->scalarTraceFuncBC(0);
            this->setGhostsScalar(bGrow, scalBC, a_stateTime, 0);

            Copier exCopier;
            exCopier.exchangeDefine(grids, bGrow.ghostVect());
            exCopier.trimEdges(grids, bGrow.ghostVect());
            bGrow.exchange(exCopier);

            // Interpolate to FC.
            for (dit.reset(); dit.ok(); ++dit) {
                const FArrayBox& ccFAB = bGrow[dit];

                for (int FCdir = 0; FCdir < SpaceDim; ++FCdir) {
                    FArrayBox& fcFAB = bFlux[dit][FCdir];

                    // We will use 1 ghost at physical boundaries
                    // and 2 ghosts elsewhere.
                    Box centerDomBox = m_problem_domain.domainBox();
                    Box loDomBox = bdryLo(centerDomBox, FCdir, 1);
                    Box hiDomBox = bdryHi(centerDomBox, FCdir, 1);
                    centerDomBox.surroundingNodes(FCdir);
                    if (!m_problem_domain.isPeriodic(FCdir)) {
                        centerDomBox.grow(FCdir, -1);
                    }

                    Box entireBox = surroundingNodes(grids[dit], FCdir);
                    Box centerBox = entireBox & centerDomBox;
                    Box loBox;
                    Box hiBox;
                    int hasLo = 0;
                    int hasHi = 0;
                    if (!m_problem_domain.isPeriodic(FCdir)) {
                        loBox = entireBox & loDomBox;
                        hiBox = entireBox & hiDomBox;
                        hasLo = (loBox.isEmpty()? 0: 1);
                        hasHi = (hiBox.isEmpty()? 0: 1);
                    }

                    CH_assert(fcFAB.box().contains(centerBox));
                    CH_assert(hasLo == 0);
                    CH_assert(hasHi == 0);
                    CH_assert(ccFAB.box().contains(enclosedCells(centerBox).grow(FCdir,2)));

                    // Interpolate!
                    FORT_CELLTOEDGE4TH (
                        CHF_FRA(fcFAB),
                        CHF_CONST_FRA(ccFAB),
                        CHF_CONST_INT(FCdir),
                        CHF_BOX(loBox),
                        CHF_CONST_INT(hasLo),
                        CHF_BOX(hiBox),
                        CHF_CONST_INT(hasHi),
                        CHF_BOX(centerBox));
                }
            }

        } else {
            // 2nd-order interpolation (average) from CC to FC.
            CellToEdge(a_b, bFlux);
        }

        for (dit.reset(); dit.ok(); ++dit) {
            bFlux[dit].mult(uAD[dit], grids[dit], 0, 0, 1);
            bFlux[dit].negate();
        }

        // Compute -Div[flux]
        Divergence::levelDivergenceMAC(a_Sb, bFlux, *m_levGeoPtr);

        // Factor in the background scalar, if needed.
        if (m_physBCPtr->useBackgroundScalar()) {
            if (s_set_bogus_values) {
                setValLevel(bFlux, s_bogus_value);
            }

            // Compute flux
            for (dit.reset(); dit.ok(); ++dit) {
                for (int dir = 0; dir < SpaceDim; ++dir) {
                    m_physBCPtr->setBackgroundScalar(bFlux[dit][dir],
                                                     0,
                                                     *m_levGeoPtr,
                                                     dit(),
                                                     a_stateTime);
                }
                bFlux[dit].mult(uAD[dit], grids[dit], 0, 0, 1);
            }

            // Compute divergence
            LevelData<FArrayBox> bkgdSrc(grids, 1);
            Divergence::levelDivergenceMAC(bkgdSrc, bFlux, *m_levGeoPtr);

            // Subtract from source terms
            for (dit.reset(); dit.ok(); ++dit) {
                a_Sb[dit].plus(bkgdSrc[dit], -1.0);
            }
        }

        checkForValidNAN(a_Sb);
    }

    // Add diffusive term...
    if (isDiffusive && s_diffSolverScheme == ProblemContext::HeatSolverScheme::EXPLICIT) {
        // Set up crse level BC.
        // Remember, coarse level may be at a more advanced time than this level.
        LevelData<FArrayBox>* crseDataPtr = NULL;
        if (m_level > 0) {
            // Get coarse-level scalars for BC's if necessary
            LevelData<FArrayBox>* newCrseScalPtr = &(crseNSPtr()->newScal(0));
            LevelData<FArrayBox>* oldCrseScalPtr = &(crseNSPtr()->oldScal(0));
            Real newCrseTime = crseNSPtr()->time();
            Real oldCrseTime = newCrseTime - crseNSPtr()->dt();

            const DisjointBoxLayout& crseGrids = crseNSPtr()->newVel().getBoxes();
            CH_assert(crseGrids == m_levGeoPtr->getCoarserPtr()->getBoxes());
            crseDataPtr = new LevelData<FArrayBox>(crseGrids, 1);

            timeInterp(*crseDataPtr, a_stateTime,
                       *oldCrseScalPtr, oldCrseTime,
                       *newCrseScalPtr, newCrseTime,
                       Interval(0,0));
        }

        // Compute the diffusive source term, D[kappa G[scalar]] and add to total.
        // (The op takes care of the exchanges.)
        LevelData<FArrayBox> diffusiveSrc(grids, 1);
        this->computeDiffusiveSrc(diffusiveSrc, a_b, crseDataPtr, 0, a_stateTime);
        for (dit.reset(); dit.ok(); ++dit) {
            a_Sb[dit].plus(diffusiveSrc[dit], 1.0);
        }

        // Free memory
        delete crseDataPtr;
        crseDataPtr = NULL;

    } else if (isDiffusive && s_diffSolverScheme != ProblemContext::HeatSolverScheme::EXPLICIT) {
        LevelData<FArrayBox> diffusiveSrc(grids, 1);
        LevelData<FluxBox> diffFlux(grids, 1, IntVect::Zero); // The IntVect::Zero is important here!

        MappedLevelFluxRegister* fineScalFluxRegPtr = NULL;
        MappedLevelFluxRegister* crseScalFluxRegPtr = NULL;
        if (s_diffusive_scalar_reflux) {
            if (!finestLevel()) {
                fineScalFluxRegPtr = m_scal_fluxreg_ptrs[0];
                CH_assert(fineScalFluxRegPtr->isDefined());
            }
            if (m_level > 0) {
                crseScalFluxRegPtr = (crseNSPtr()->m_scal_fluxreg_ptrs[0]);
                CH_assert(crseScalFluxRegPtr->isDefined());
            }
        }

        LevelData<FArrayBox>* oldCrseScalPtr = NULL;
        LevelData<FArrayBox>* newCrseScalPtr = NULL;
        Real oldCrseTime = -1e8;
        Real newCrseTime = 1e8;
        if (m_level > 0) {
            newCrseScalPtr = &(crseNSPtr()->newScal(0));
            oldCrseScalPtr = &(crseNSPtr()->oldScal(0));
            newCrseTime = crseNSPtr()->time();
            oldCrseTime = newCrseTime - crseNSPtr()->dt();
        }

        int numberMGlevels = (m_level == 0) ? 0 : 1;

        m_diffSolverPtrs[0]->computeDiffusion(diffusiveSrc,
                                              a_b,          // old b
                                              a_Sb,         // source terms
                                              diffFlux,
                                              fineScalFluxRegPtr,
                                              crseScalFluxRegPtr,
                                              oldCrseScalPtr,
                                              newCrseScalPtr,
                                              a_stateTime,  // old time
                                              oldCrseTime,
                                              newCrseTime,
                                              h,            // stage dt
                                              numberMGlevels,
                                              false,        // zero-out new b?
                                              true);        // already kappa weighted?

        for (dit.reset(); dit.ok(); ++dit) {
            a_Sb[dit].plus(diffusiveSrc[dit], 1.0);
        }
    }

    // Velocity ----------------------------------------------------------------
    // Compute advective source term...
    {
        switch (s_nonlinearDifferencingForm) {
        case ProblemContext::NonlinearDifferencingForm::NONE:
            {
                setValLevel(a_Su, 0.0);
            }
            break;
        case ProblemContext::NonlinearDifferencingForm::CONSERVATIVE:
            {
                LevelData<FluxBox> momentumFlux(grids, SpaceDim);
                if (s_set_bogus_values) {
                    setValLevel(momentumFlux, s_bogus_value);
                }

                // Send u to FC
                CellToEdge(a_u, momentumFlux);

                // Project
                LevelData<FluxBox> gradPhi(grids, SpaceDim);
                gradMACPressure(gradPhi, h);
                m_levGeoPtr->divByJ(gradPhi);
                for (dit.reset(); dit.ok(); ++dit) {
                    D_TERM(momentumFlux[dit][0].plus(gradPhi[dit][0], -1.0);,
                           momentumFlux[dit][1].plus(gradPhi[dit][1], -1.0);,
                           momentumFlux[dit][2].plus(gradPhi[dit][2], -1.0);)
                }
                EdgeVelBCHolder edgeVelBC(m_physBCPtr->advectingVelFuncBC(isViscous));

                LevelData<FluxBox> uHalfFC(grids, 1);
                if (s_set_bogus_values) {
                    setValLevel(uHalfFC, s_bogus_value);
                }
                m_levGeoPtr->sendToCartesianBasis(a_u, true);
                CellToEdge(a_u, uHalfFC);
                m_levGeoPtr->sendToMappedBasis(a_u, true);

                // Loop over grids and multiply the Cartesian based uHalfFC
                // by the mapping based advecting vel (stored in momentumFlux).
                for (dit.reset(); dit.ok(); ++dit) {
                    for (int dir = 0; dir < SpaceDim; ++dir) {
                        for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                            momentumFlux[dit][dir].mult(uHalfFC[dit][dir], 0, velComp, 1);
                        }
                    }
                }

                // Compute advective forcing term = -Div[momentumFlux]
                Divergence::levelDivergenceMAC(a_Su, momentumFlux, *m_levGeoPtr, a_stateTime, NULL);
                for (dit.reset(); dit.ok(); ++dit) {
                    momentumFlux[dit] *= -1.0;
                }

                // Use these fluxes to update the momentum flux registers.
                if (s_advective_momentum_reflux) {
                    if (a_stage == 0) {
                        Real bdt = m_dt / 4.0;
                        updateVelFluxRegister(momentumFlux, bdt);
                    } else if (a_stage == 2) {
                        Real bdt = m_dt * 3.0 / 4.0;
                        updateVelFluxRegister(momentumFlux, bdt);
                    }
                }
            }
            break;
        case ProblemContext::NonlinearDifferencingForm::ADVECTIVE:
            MayDay::Error("RK3 needs to be updated");
            // {
            //     // Initialize adv_term
            //     if (s_set_bogus_values) {
            //         setValLevel(adv_term, s_bogus_value);
            //     }

            //     // Compute Av[adv_vel / J]. Use a_newVel as a temp holder.
            //     LevelData<FArrayBox>& half_vel = a_newVel;
            //     EdgeToCell(a_advVel, half_vel);
            //     m_levGeoPtr->divByJ(half_vel);

            //     // Compute the negative of the advective term, half_vel.Grad[pred_vel]
            //     Gradient::levelCCDotGradFC(adv_term, half_vel, pred_vel, *m_levGeoPtr);

            //     // Set adv_vel to -u.Grad[u] to mimic other source terms.
            //     for (dit.reset(); dit.ok(); ++dit) {
            //         adv_term[dit] *= -1.0;
            //     }

            //     if (s_advective_momentum_reflux) {
            //         // Loop over grids and multiply the Cartesian based pred_vel
            //         // by the mapping based a_advVel. This will give us the
            //         // fluxes adv_vel * pred_vel.
            //         for (dit.reset(); dit.ok(); ++dit) {
            //             for (int dir = 0; dir < SpaceDim; ++dir) {
            //                 for (int velComp = 0; velComp < SpaceDim; ++velComp) {
            //                     pred_vel[dit][dir].mult(a_advVel[dit][dir], 0, velComp, 1);
            //                 }
            //             }
            //         }

            //         // Now, pred_vel contains the momentum fluxes
            //         //                   comp 0                comp 1                comp 2
            //         // FC dir 0  (J*Uad[0]*pred_vel[0], J*Uad[0]*pred_vel[1], J*Uad[0]*pred_vel[2])
            //         // FC dir 1  (J*Uad[1]*pred_vel[0], J*Uad[1]*pred_vel[1], J*Uad[1]*pred_vel[2])
            //         // FC dir 2  (J*Uad[2]*pred_vel[0], J*Uad[2]*pred_vel[1], J*Uad[2]*pred_vel[2])

            //         // Use these fluxes to update the momentum flux registers.
            //         updateVelFluxRegister(pred_vel, a_FRscale);
            //     }
            // }
            break;
        default:
            MayDay::Error("Unknown nonlinear differencing form");
        };
    }


    // Add external forces...

    // If there is gravity, we need to add its contribution to the forcing.
    if (s_gravityMethod != ProblemContext::GravityMethod::NONE) {
        // Compute the gravitational source term.
        LevelData<FArrayBox> gravSource(grids, SpaceDim);
        this->fillGravSource(gravSource, a_stateTime,
                             false);  // add background?

        // Combine the advective and gravitational source terms.
        for (dit.reset(); dit.ok(); ++dit) {
            a_Su[dit].plus(gravSource[dit], 1.0);
        }
    }

    // Compute the tidal forcing term
    if (s_tidalU0.sum() * s_tidalOmega != 0.0) {
        // Compute the gravitational source term.
        LevelData<FArrayBox> tidalSource(grids, SpaceDim);
        this->fillTidalSource(tidalSource, a_stateTime, h);

        // Combine the advective and gravitational source terms.
        for (dit.reset(); dit.ok(); ++dit) {
            a_Su[dit].plus(tidalSource[dit], 1.0);
        }
    }

    // If we are using a sponge layer, add its contribution.
    if (m_physBCPtr->useSpongeLayer()) {
        m_levGeoPtr->sendToCartesianBasis(a_u, true);

        LevelData<FArrayBox> spongeLayerSource(grids, SpaceDim);
        m_physBCPtr->fillSpongeLayerSrcTerm(spongeLayerSource,
                                            a_u,
                                            a_stateTime,
                                            h,
                                            *m_levGeoPtr);

        for (dit.reset(); dit.ok(); ++dit) {
            a_Su[dit].plus(spongeLayerSource[dit], 1.0);
        }

        m_levGeoPtr->sendToMappedBasis(a_u, true);
    }

    // Add viscous term...
    if (isViscous && s_viscSolverScheme == ProblemContext::HeatSolverScheme::EXPLICIT) {
        // Add explicit viscous source

        m_levGeoPtr->sendToCartesianBasis(a_u, true);

        LevelData<FArrayBox> viscSource(grids, SpaceDim);
        computeViscousSrc(viscSource, a_u, a_stateTime);
        for (dit.reset(); dit.ok(); ++dit) {
            a_Su[dit].plus(viscSource[dit], 1.0);
        }

        m_levGeoPtr->sendToMappedBasis(a_u, true);

    } else if (isViscous && s_viscSolverScheme != ProblemContext::HeatSolverScheme::EXPLICIT) {
        // Add implicitly derived source term

        LevelData<FArrayBox> viscSource(grids, 1);
        LevelData<FluxBox> viscFlux(grids, 1, IntVect::Zero); // The IntVect::Zero is important here!

        const int numberMGlevels = (m_level == 0) ? 0 : 1;

        MappedLevelFluxRegister* fineFluxRegisterPtr = NULL;
        MappedLevelFluxRegister* crseFluxRegisterPtr = NULL;
        if (s_diffusive_momentum_reflux) {
            if (!finestLevel()) {
                fineFluxRegisterPtr = &m_vel_flux_reg;
                CH_assert(fineFluxRegisterPtr->isDefined());
            }
            if (m_level > 0) {
                crseFluxRegisterPtr = &(crseNSPtr()->m_vel_flux_reg);
                CH_assert(crseFluxRegisterPtr->isDefined());
            }
        }

        Real crseOldTime = -1e8;
        Real crseNewTime = 1e8;
        if (numberMGlevels > 0) {
            crseNewTime = crseNSPtr()->m_time;
            Real crseDt = crseNSPtr()->m_dt;
            crseOldTime = crseNewTime - crseDt;
        }

        LevelData<FArrayBox>* tmpCrseOldVelPtr = NULL;
        LevelData<FArrayBox>* tmpCrseNewVelPtr = NULL;
        if (m_level > 0) {
            const AMRNavierStokes* crseLevelPtr = this->crseNSPtr();
            const LevelData<FArrayBox>& crseOldVel = *(crseLevelPtr->m_vel_old_ptr);
            const LevelData<FArrayBox>& crseNewVel = *(crseLevelPtr->m_vel_new_ptr);
            const DisjointBoxLayout& crseGrids = crseNewVel.getBoxes();

            tmpCrseOldVelPtr = new LevelData<FArrayBox>(crseGrids, SpaceDim);
            tmpCrseNewVelPtr = new LevelData<FArrayBox>(crseGrids, SpaceDim);

            if (s_set_bogus_values) {
                setValLevel(*tmpCrseOldVelPtr, s_bogus_value);
                setValLevel(*tmpCrseNewVelPtr, s_bogus_value);
            }

            DataIterator crseDit = crseNewVel.dataIterator();
            for (crseDit.reset(); crseDit.ok(); ++crseDit) {
                (*tmpCrseOldVelPtr)[crseDit].copy(crseOldVel[crseDit]);
                (*tmpCrseNewVelPtr)[crseDit].copy(crseNewVel[crseDit]);
            }

            m_levGeoPtr->sendToCartesianBasis(*tmpCrseOldVelPtr);
            m_levGeoPtr->sendToCartesianBasis(*tmpCrseNewVelPtr);
        } // end convert coarse-level velocities

        m_levGeoPtr->sendToCartesianBasis(a_u, true);


        for (int comp = 0; comp < SpaceDim; comp++) {
            const Interval intvl(comp, comp);

            LevelData<FArrayBox> compOldVelocity;
            aliasLevelData(compOldVelocity, &a_u, intvl);

            LevelData<FArrayBox> compSrc;
            aliasLevelData(compSrc, &a_Su, intvl);

            if (numberMGlevels == 0) {
                // Bottom level - no coarser.
                m_viscSolverPtrs[comp]->computeDiffusion(viscSource,
                                                         compOldVelocity,
                                                         compSrc,
                                                         viscFlux,
                                                         fineFluxRegisterPtr,
                                                         NULL,         // crse flux reg
                                                         NULL,         // crse old phi ptr
                                                         NULL,         // crse new phi ptr
                                                         a_stateTime,  // old time
                                                         crseOldTime,
                                                         crseNewTime,
                                                         h,            // stage dt
                                                         numberMGlevels,
                                                         false,        // zero-out new vel?
                                                         true,         // already kappa weighted?
                                                         comp);        // flux reg start comp
            } else {
                // Coarser level exists.
                LevelData<FArrayBox> compCrseOldVel;
                aliasLevelData(compCrseOldVel,
                               tmpCrseOldVelPtr,
                               intvl);

                LevelData<FArrayBox> compCrseNewVel;
                aliasLevelData(compCrseNewVel,
                               tmpCrseNewVelPtr,
                               intvl);

                m_viscSolverPtrs[comp]->computeDiffusion(viscSource,
                                                         compOldVelocity,
                                                         compSrc,
                                                         viscFlux,
                                                         fineFluxRegisterPtr,
                                                         crseFluxRegisterPtr,
                                                         &compCrseOldVel,
                                                         &compCrseNewVel,
                                                         a_stateTime,  // old time
                                                         crseOldTime,
                                                         crseNewTime,
                                                         h,            // stage dt
                                                         numberMGlevels,
                                                         false,        // zero-out new vel?
                                                         true,         // already kappa weighted?
                                                         comp);        // flux reg start comp
            }

            // Add to viscous source
            for (dit.reset(); dit.ok(); ++dit) {
                a_Su[dit].plus(viscSource[dit], 0, comp, 1);
            }
        }

        // Restore basis
        m_levGeoPtr->sendToMappedBasis(a_u, true);
    } // end if viscous

    // Source term is needed in mapped basis.
    m_levGeoPtr->sendToMappedBasis(a_Su, true);

    checkForValidNAN(a_Su);
}


// -----------------------------------------------------------------------------
// Performs the TGA solve.
// [I - h*L/2]q = [I + h*L/2]q + h*(beta*newS + zeta*oldS)
// a_time is the time a_newSu is at.
// a_stage can be 1, 2, or 3.
// -----------------------------------------------------------------------------
void AMRNavierStokes::updateState (LevelData<FArrayBox>&       a_u,
                                   LevelData<FArrayBox>&       a_b,
                                   Real&                       a_stateTime,
                                   const LevelData<FArrayBox>& a_newSu,
                                   const LevelData<FArrayBox>& a_newSb,
                                   const LevelData<FArrayBox>& a_oldSu,
                                   const LevelData<FArrayBox>& a_oldSb,
                                   const int                   a_stage)
{
    // Sanity checks
    CH_assert(a_stage == 1 || a_stage == 2 || a_stage == 3);
    CH_assert(m_levGeoPtr->getBoxes() == a_u.getBoxes());
    CH_assert(m_levGeoPtr->getDomain() == m_problem_domain);
    CH_assert(m_time - m_dt - TIME_EPS < a_stateTime);
    CH_assert(a_stateTime < m_time + TIME_EPS);
    CH_assert(s_num_scal_comps == 1);

    // Set up some basic values
    const RealVect& dx = m_levGeoPtr->getDx();
    const DisjointBoxLayout& grids = a_u.getBoxes();
    DataIterator dit = grids.dataIterator();
    const Box domainBox = m_problem_domain.domainBox();

    const bool isViscous = (s_nu > 0.0);

    // Set RK3 stage coefficients.
    Real h;     // The effective stage dt.
    Real beta;  // The newS weight.
    Real zeta;  // The oldS weight.
    switch (a_stage) {
    case 1:
        h = m_dt * 8. / 15.;
        beta = 1.;
        zeta = 0.;
        break;
    case 2:
        h = m_dt * 2. / 15.;
        beta = 25. / 8.;
        zeta = -17. / 8.;
        break;
    case 3:
        h = m_dt * 1. / 3.;
        beta = 9. / 4.;
        zeta = -5. / 4.;
        break;
    default:
        MayDay::Error("Bad RK3 stage");
    }

    // Project
    bool useImplicitGrav = (s_gravityMethod == ProblemContext::GravityMethod::IMPLICIT);
    useImplicitGrav &= m_physBCPtr->useBackgroundScalar();

    if (useImplicitGrav) {
        // b update.
        LevelData<FArrayBox> newb(grids, 1);
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& newStateFAB = newb[dit];
            const FArrayBox& oldStateFAB = a_b[dit];
            const FArrayBox& newSFAB = a_newSb[dit];
            const FArrayBox& oldSFAB = a_oldSb[dit];

            newStateFAB.copy(oldStateFAB);
            newStateFAB.plus(newSFAB, h*beta);
            if (zeta != 0.0) newStateFAB.plus(oldSFAB, h*zeta);
        }

        // u update.
        LevelData<FArrayBox> newu(grids, SpaceDim, IntVect::Unit);
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& newStateFAB = newu[dit];
            const FArrayBox& oldStateFAB = a_u[dit];
            const FArrayBox& newSFAB = a_newSu[dit];
            const FArrayBox& oldSFAB = a_oldSu[dit];

            newStateFAB.copy(oldStateFAB);
            newStateFAB.plus(newSFAB, h*beta);
            if (zeta != 0.0) newStateFAB.plus(oldSFAB, h*zeta);
        }

        // Advecting velocity
        LevelData<FluxBox> uAD(grids, 1);
        {
            if (s_set_bogus_values) {
                setValLevel(uAD, s_bogus_value);
            }

            // Send u to FC and scale by J
            CellToEdge(a_u, uAD);

            // Set BCs
            EdgeVelBCHolder edgeVelBC(m_physBCPtr->advectingVelFuncBC(isViscous));
            edgeVelBC.setGhosts(uAD,
                                NULL,
                                dx,
                                &(m_levGeoPtr->getFCJgup()),
                                false, // inhomogeneous
                                a_stateTime);

            // MAC project
            if (s_isIncompressible) {
                LevelData<FluxBox> gradPhi(grids, 1);
                gradMACPressure(gradPhi, h);

                m_levGeoPtr->divByJ(gradPhi);
                for (dit.reset(); dit.ok(); ++dit) {
                    for (int dir = 0; dir < SpaceDim; ++dir) {
                        uAD[dit][dir].plus(gradPhi[dit][dir], -1.0);
                    }
                }
            }

            // Set BCs
            edgeVelBC.setGhosts(uAD,
                                NULL,
                                dx,
                                &(m_levGeoPtr->getFCJgup()),
                                false,  // inhomogeneous
                                a_stateTime);

            // Scale as a flux
            m_levGeoPtr->multByJ(uAD);

            checkForValidNAN(uAD);
        }

        // Project!
        doCCIGProjection(newu, newb, a_u, a_b, uAD, a_stateTime, h, true);

        // Update states
        for (dit.reset(); dit.ok(); ++dit) {
            a_u[dit].copy(newu[dit]);
            a_b[dit].copy(newb[dit]);
        }

    } else {
        // b update.
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& stateFAB = a_b[dit];
            const FArrayBox& newSFAB = a_newSb[dit];
            const FArrayBox& oldSFAB = a_oldSb[dit];

            stateFAB.plus(newSFAB, h*beta);
            if (zeta != 0.0) stateFAB.plus(oldSFAB, h*zeta);
        }

        // u update.
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& stateFAB = a_u[dit];
            const FArrayBox& newSFAB = a_newSu[dit];
            const FArrayBox& oldSFAB = a_oldSu[dit];

            stateFAB.plus(newSFAB, h*beta);
            if (zeta != 0.0) stateFAB.plus(oldSFAB, h*zeta);
        }

        // Project!
        doCCProjection(a_u, a_stateTime, h, true);
    }

    // Update stateTime
    a_stateTime += h;

    checkForValidNAN(a_b);
    checkForValidNAN(a_u);
}

