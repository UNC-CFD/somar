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
#include "timeInterp.H"
#include "Gradient.H"
#include "EdgeToCell.H"
#include <iomanip>
#include "Constants.H"
#include "ProblemContext.H"

// #define nanCheck(x) checkForValidNAN(x)
#define nanCheck(x)

// All of these should be set to true for accuracy.
// These switches should only be used to debug the algorithm.
static const bool s_useLaggedPressure = true;
static const bool s_useScalAdvectiveSource = true;
static const bool s_useVelAdvectiveSource = true;


// -----------------------------------------------------------------------------
// This does the actual computation to update the state variables.
// I've included this function since the initializeGlobalPressure() and advance()
// functions are, for the most part, the same piece of code.
// -----------------------------------------------------------------------------
void AMRNavierStokes::PPMTimeStep (const Real a_oldTime,
                                   const Real a_dt,
                                   const bool a_updatePassiveScalars,
                                   const bool a_doLevelProj)
{
    CH_TIME("AMRNavierStokes::PPMTimeStep");
    pout() << setiosflags(ios::scientific) << setprecision(8) << flush;

    // Set up some basic values
    const RealVect& dx = m_levGeoPtr->getDx();
    const DisjointBoxLayout& grids = newVel().getBoxes();
    DataIterator dit = grids.dataIterator();
    const Box domainBox = m_problem_domain.domainBox();
    const bool isViscous = (s_nu > 0.0);

    // Initialize all flux registers
    if (!finestLevel()) {
        m_vel_flux_reg.setToZero();
        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            m_scal_fluxreg_ptrs[comp]->setToZero();
        }
        m_lambda_flux_reg.setToZero();
    }

    // Sanity checks
    CH_assert(m_levGeoPtr->getBoxes() == grids);
    CH_assert(m_levGeoPtr->getDomain() == m_problem_domain);
    CH_assert(Abs(a_oldTime - (m_time - a_dt)) < TIME_EPS);
    CH_assert(s_num_scal_comps <= 1);
    CH_assert(s_num_scal_comps > 0 || s_gravityMethod == ProblemContext::GravityMethod::NONE);

    // Reference new holders
    LevelData<FArrayBox>& new_vel = newVel();
    LevelData<FArrayBox>& new_lambda = newLambda();
    LevelData<FArrayBox> new_b;
    if (s_num_scal_comps > 0) {
        aliasLevelData(new_b, &(newScal(0)), Interval(0,0));
    }


    // Compute advecting velocities
    LevelData<FArrayBox> old_vel(grids, SpaceDim, m_tracingGhosts);
    fillVelocity(old_vel, a_oldTime);
    old_vel.exchange(m_tracingExCopier);
    nanCheck(old_vel);

    LevelData<FluxBox> adv_vel(grids, 1, IntVect::Unit); // Changed from m_tracingGhosts to 1
    computeAdvectingVelocities(adv_vel, old_vel, a_oldTime, a_dt);
    nanCheck(adv_vel);

    if (a_updatePassiveScalars) {
        // Lambda update
        LevelData<FArrayBox> old_lambda;
        fillLambda(old_lambda, a_oldTime);
        old_lambda.exchange(m_tracingExCopier);
        nanCheck(old_lambda);

        LevelData<FArrayBox> dLdt(grids, 1);
        getNewLambda(dLdt, new_lambda, old_lambda, old_vel, adv_vel, a_oldTime, a_dt, a_dt);
        nanCheck(new_lambda);
    }

    if (s_num_scal_comps > 0) {
        // Scalar update
        LevelData<FArrayBox> old_b;
        fillScalars(old_b, a_oldTime, 0);
        old_b.exchange(m_tracingExCopier);
        nanCheck(old_b);

        LevelData<FArrayBox> dSdt(grids, 1);
        getNewScalar(dSdt, new_b, old_b, old_vel, adv_vel, a_oldTime, a_dt, a_dt, 0);
        nanCheck(new_b);
    }

    for (int comp = 1; comp < s_num_scal_comps; ++comp) {
        // Scalar update
        LevelData<FArrayBox> old_scal;
        fillScalars(old_scal, a_oldTime, comp);
        old_scal.exchange(m_tracingExCopier);

        LevelData<FArrayBox> dSdt(grids, 1);
        getNewScalar(dSdt, newScal(comp), old_scal, old_vel, adv_vel, a_oldTime, a_dt, a_dt, comp);
    }

    {   // Update CC velocities
        LevelData<FArrayBox> dUdt(grids, SpaceDim);
        getNewVelocity(dUdt, new_vel, old_vel, adv_vel, a_oldTime, a_dt, a_dt);
        nanCheck(new_vel);
    }

    // Project CC velocities
    doCCProjection(new_vel, a_oldTime + a_dt, a_dt, a_doLevelProj);
}


// -----------------------------------------------------------------------------
// Compute the half-time, VD corrected, projected, FC advecting velocity.
// All BCs and exchanges will be done on a_advVel, but not inputs.
// -----------------------------------------------------------------------------
void AMRNavierStokes::computeAdvectingVelocities (LevelData<FluxBox>&   a_advVel,
                                                  LevelData<FArrayBox>& a_oldVel,
                                                  const Real            a_oldTime,
                                                  const Real            a_dt)
{
    CH_TIME("AMRNavierStokes::computeAdvectingVelocities");

    const RealVect& dx = m_levGeoPtr->getDx();
    const DisjointBoxLayout& grids = newVel().getBoxes();
    DataIterator dit = grids.dataIterator();
    const bool isViscous = (s_nu > 0.0);
    const Real halfTime = a_oldTime + 0.5 * a_dt;

    // We need an advecting velocity to do the tracing. We will center this
    // at a_oldTime, rendering the outcome first-order in time.
    CellToEdge(a_oldVel, a_advVel);
    m_levGeoPtr->multByJ(a_advVel);
    {
        EdgeVelBCHolder edgeVelBC(m_physBCPtr->advectingVelFuncBC(isViscous));
        edgeVelBC.setGhosts(a_advVel,
                            NULL,
                            dx,
                            &(m_levGeoPtr->getFCJgup()),
                            false, // inhomogeneous
                            a_oldTime);
    }

    // Compute predicted velocities
    {
        LevelData<FluxBox> predVel(grids, SpaceDim); // Changed from m_tracingGhosts to 0
        predictVelocities(predVel, a_oldVel, a_advVel, a_oldTime, a_dt);

        // Copy normal components to a_advVel. From here on, a_advVel will be
        // staggered in time.
        for (dit.reset(); dit.ok(); ++dit) {
            for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
                a_advVel[dit][dir].copy(predVel[dit][dir], dir, 0, 1);
            }
        }
    }
    nanCheck(a_advVel);

    // Scale a_advVel as a flux.
    m_levGeoPtr->multByJ(a_advVel);

    // Do the level MAC projection (This function exchanges a_advVel)
    if (s_isIncompressible) {
        const Real projDt = 0.5 * a_dt;     // a_advVel is staggered in time
        const Real projTime = a_oldTime + projDt;
        Vector<LevelData<FluxBox>*> amrVel(1, &a_advVel);

        pout() << "Level " << m_level << " MAC proj: " << flush;
        m_macProjector.levelProject(amrVel,
                                    m_levGeoPtr,
                                    projTime,
                                    projDt,
                                    true,        // a_advVel is a flux
                                    true,        // isLevelSolve
                                    false);      // forceHomogSolve
    }

    // Add volume discrepancy correction to advecting velocity
    if (s_etaLambda > 0.0 && s_applyFreestreamCorrection) {
        for (dit.reset(); dit.ok(); ++dit) {
            for (int dir = 0; dir < SpaceDim; ++dir) {
                a_advVel[dit][dir].plus(m_gradELambda[dit][dir], 0, 0, 1);
            }
        }
    }

    // Set BCs on traced, staggered, corrected velocities
    EdgeVelBCHolder edgeVelBC(m_physBCPtr->advectingVelFuncBC(isViscous));
    edgeVelBC.setGhosts(a_advVel,
                        NULL,
                        dx,
                        &(m_levGeoPtr->getFCJgup()),
                        false,  // inhomogeneous
                        halfTime);

    // This seems to help multi-level, mapped solves.
    a_advVel.exchange(m_oneGhostExCopier);
}


// -----------------------------------------------------------------------------
// Calculate d(lambda)/dt = RHS centered at the time of a_oldVel + a_dt/2
// -----------------------------------------------------------------------------
void AMRNavierStokes::getNewLambda (LevelData<FArrayBox>&       a_rhs,
                                    LevelData<FArrayBox>&       a_new_lambda,
                                    LevelData<FArrayBox>&       a_old_lambda,
                                    const LevelData<FArrayBox>& a_oldVel,
                                    const LevelData<FluxBox>&   a_advVel,
                                    const Real                  a_oldTime,
                                    const Real                  a_dt,
                                    const Real                  a_FRscale)
{
    CH_TIME("AMRNavierStokes::getNewLambda");

    CH_assert(a_rhs       .nComp() == 1);
    CH_assert(a_new_lambda.nComp() == 1);
    CH_assert(a_old_lambda.nComp() == 1);
    CH_assert(a_advVel    .nComp() == 1);
    CH_assert(a_oldVel    .nComp() == SpaceDim);

    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    DataIterator dit = grids.dataIterator();

    // Create/Reference data holders
    LevelData<FluxBox> lambda_flux(grids, 1); // Changed from m_tracingGhosts to 0
    if (s_set_bogus_values) {
        setValLevel(a_rhs, s_bogus_value);
        setValLevel(a_new_lambda, s_bogus_value);
        setValLevel(lambda_flux, s_bogus_value);
    }

    // Set all ghosts on lambda
    BCMethodHolder lambdaBC = m_physBCPtr->lambdaFuncBC();
    setGhostsLambda(a_old_lambda, lambdaBC, a_oldTime);

    // Get the BC holders
    Tuple<BCMethodHolder,SpaceDim> BCValues = m_physBCPtr->lambdaRiemannBC();
    Tuple<BCMethodHolder,SpaceDim> BCSlopes = m_physBCPtr->lambdaSlopeBC();

    // Do the tracing
    m_advectUtilLambda.predictScalar(lambda_flux,
                                     a_old_lambda,
                                     NULL,
                                     a_oldVel,
                                     a_advVel,
                                     a_dt,
                                     *m_levGeoPtr,
                                     BCValues,
                                     BCSlopes,
                                     a_oldTime,
                                     true); // return fluxes

    // Update lambda from fluxes
    Divergence::levelDivergenceMAC(a_rhs, lambda_flux, *m_levGeoPtr);
    for (dit.reset(); dit.ok(); ++dit) {
        a_rhs[dit] *= -1.0;

        a_new_lambda[dit].copy(a_old_lambda[dit]);
        a_new_lambda[dit].plus(a_rhs[dit], a_dt);
    }

    // Do flux register stuff...
    if (s_advective_lambda_reflux && a_FRscale != 0.0) {
        updateLambdaFluxRegister(lambda_flux, a_FRscale);
    }
}


// -----------------------------------------------------------------------------
// Calculate d(scalar)/dt = RHS centered at the time of a_oldVel + a_dt/2
// Since diffusive scalars are updated implicitly, this function actually
// calculates a_newScalar, then sets a_rhs = (a_newScalar - a_oldScalar) / dt
//
// WARNING: At some point, I just assumed that a_comp = 0. If you plan to use
// a scalar other than buoyancy, plan to debug!
// -----------------------------------------------------------------------------
void AMRNavierStokes::getNewScalar (LevelData<FArrayBox>&       a_rhs,
                                    LevelData<FArrayBox>&       a_newScalar,
                                    LevelData<FArrayBox>&       a_oldScalar,
                                    const LevelData<FArrayBox>& a_oldVel,
                                    const LevelData<FluxBox>&   a_advVel,
                                    const Real                  a_oldTime,
                                    const Real                  a_dt,
                                    const Real                  a_FRscale,
                                    const int                   a_comp)
{
    CH_TIME("AMRNavierStokes::getNewScalar");

    // Gather data
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    DataIterator dit = grids.dataIterator();
    const Interval& interv = a_oldScalar.interval();
    const bool isDiffusive = (s_scal_coeffs[a_comp] > 0.0);

    // Sanity checks
    CH_assert(a_newScalar.nComp() == 1);
    CH_assert(a_rhs      .nComp() == 1);
    CH_assert(a_oldScalar.nComp() == 1);
    CH_assert(a_advVel   .nComp() == 1);
    CH_assert(a_oldVel   .nComp() == SpaceDim);
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < s_num_scal_comps);

    CH_assert(grids == a_rhs.getBoxes());
    CH_assert(grids == a_newScalar.getBoxes());
    CH_assert(grids == a_oldScalar.getBoxes());
    CH_assert(grids == a_oldVel.getBoxes());
    CH_assert(grids == a_advVel.getBoxes());

    // Create/Reference data holders
    LevelData<FluxBox> scalarFlux(grids, 1); // Changed from m_tracingGhosts to 0
    if (s_set_bogus_values) {
        setValLevel(a_rhs, s_bogus_value);
        setValLevel(a_newScalar, s_bogus_value);
        setValLevel(scalarFlux, s_bogus_value);
    }

    // Set all ghosts on the scalar
    BCMethodHolder scalBC = m_physBCPtr->scalarTraceFuncBC(a_comp);
    setGhostsScalar(a_oldScalar, scalBC, a_oldTime, a_comp);

    // The background advective source term is -Div[Uad^a * backgroundScalar] where Uad is a flux.
    // TODO: This is a bit of a waste!!!
    LevelData<FArrayBox> bkgdSrc(grids, 1, IntVect::Unit);

    if (m_physBCPtr->useBackgroundScalar() && s_gravityMethod == ProblemContext::GravityMethod::EXPLICIT) {
        LevelData<FluxBox> bkgdFlux(grids, 1, IntVect::Unit);
#ifndef NDEBUG
        setValLevel(bkgdSrc, quietNAN);
        setValLevel(bkgdFlux, quietNAN);
#endif
        setValLevel(bkgdSrc, 0.0);
        setValLevel(bkgdFlux, 0.0);

        // Fill bkgdFlux with (Uad) * (background scalar)
        for (dit.reset(); dit.ok(); ++dit) {
            FluxBox& bkgdFluxFlub = bkgdFlux[dit];
            const FluxBox& UadFlub = a_advVel[dit];
            const Box CCRegion = bkgdFluxFlub.box();
            CH_assert(UadFlub.box().contains(CCRegion));

            // First fill background scalar
            for (int dir = 0; dir < SpaceDim; ++dir) {
                m_physBCPtr->setBackgroundScalar(bkgdFluxFlub[dir],
                                                 0,
                                                 *m_levGeoPtr,
                                                 dit(),
                                                 a_oldTime); // This shouldn't matter for the static background!
            }

            // Then multiply by Uad
            bkgdFluxFlub.mult(UadFlub, CCRegion, 0, 0, 1);
        }

        // Compute divergence
        Divergence::levelDivergenceMAC(bkgdSrc, bkgdFlux, *m_levGeoPtr);

        // Turn into -Div[...]
        for (dit.reset(); dit.ok(); ++dit) {
            bkgdSrc[dit].negate();
        }

        // Exchange data.       // TODO: Is this necessary?
        bkgdSrc.exchange(m_oneGhostExCopier);

    } else {
        // No background source. Just set to zero.
        setValLevel(bkgdSrc, 0.0);
    }

    if (isDiffusive) {
        // This scalar is diffusive...

        // Create data holders
        LevelData<FArrayBox> diffusiveSrc(grids, 1, m_tracingGhosts);
        LevelData<FArrayBox> tempStorage (grids, 1, m_tracingGhosts);     // TODO: Can this be a_rhs?

        // Set to bogus values
        if (s_set_bogus_values) {
            setValLevel(diffusiveSrc, s_bogus_value);
            setValLevel(tempStorage, s_bogus_value);
        }

        // Get coarse-level scalars for BC's if necessary
        LevelData<FArrayBox>* oldCrseScalPtr = NULL;
        LevelData<FArrayBox>* newCrseScalPtr = NULL;
        Real oldCrseTime = -1e8;
        Real newCrseTime = 1e8;

        if (m_level > 0) {
            newCrseScalPtr = &(crseNSPtr()->newScal(a_comp));
            oldCrseScalPtr = &(crseNSPtr()->oldScal(a_comp));
            newCrseTime = crseNSPtr()->time();
            oldCrseTime = newCrseTime - crseNSPtr()->dt();
        }

        // 1. Compute the diffusive source term
        {
            LevelData<FArrayBox>* crseDataPtr = NULL;

            // Copy scalar into temp storage so we don't overwrite C/F BC's.
            // Define our own copier here to prevent ghost cells
            // being filled (avoiding a potentially expensive exchange)
            // since they will be filled a bit later when we compute the Laplacian
            {
                Copier noGhostCopier(a_oldScalar.getBoxes(),
                                     tempStorage.getBoxes(),
                                     IntVect::Zero);
                a_oldScalar.copyTo(interv,
                                    tempStorage,
                                    tempStorage.interval(),
                                    noGhostCopier);
            }

            // Set up crse level BC.
            // Remember, coarse level may be at a more advanced time than this level.
            if (m_level > 0) {
                CH_assert(oldCrseScalPtr != NULL);
                CH_assert(newCrseScalPtr != NULL);

                const DisjointBoxLayout& crseGrids = crseNSPtr()->newVel().getBoxes();
                CH_assert(crseGrids == m_levGeoPtr->getCoarserPtr()->getBoxes());
                crseDataPtr = new LevelData<FArrayBox>(crseGrids, 1);

                timeInterp(*crseDataPtr, a_oldTime,
                           *oldCrseScalPtr, oldCrseTime,
                           *newCrseScalPtr, newCrseTime,
                           interv);
            }

            // Compute the diffusive source term, coeff * L[scalar].
            // (The op takes care of the exchanges.)
            this->computeLapScal(diffusiveSrc, tempStorage, crseDataPtr, NULL);
            for (dit.reset(); dit.ok(); ++dit) {
                diffusiveSrc[dit] *= s_scal_coeffs[a_comp];
            }

            // Free memory
            if (crseDataPtr != NULL) {
                delete crseDataPtr;
            }
        }

        // 2. Compute the advective update

        // The background advective source term is -w * d[backgroundScalar] / dz
        // TODO: This is a bit of a waste!!!
        if (m_physBCPtr->useBackgroundScalar()) {
            for (dit.reset(); dit.ok(); ++dit) {
                diffusiveSrc[dit].plus(bkgdSrc[dit], 0, 0, 1);
            }
        }

        // Get the BC holders
        // const RealVect fluxScale = m_levGeoPtr->getFluxScale();
        Tuple<BCMethodHolder,SpaceDim> BCValues = m_physBCPtr->scalarRiemannBC(a_comp);
        Tuple<BCMethodHolder,SpaceDim> BCSlopes = m_physBCPtr->scalarSlopeBC(a_comp);

        // Do the tracing
        m_advectUtilScal.predictScalar(scalarFlux,
                                       a_oldScalar,
                                       (s_useScalAdvectiveSource? &diffusiveSrc: NULL),
                                       a_oldVel,
                                       a_advVel,
                                       a_dt,
                                       *m_levGeoPtr,
                                       BCValues,
                                       BCSlopes,
                                       a_oldTime,
                                       true); // return fluxes

        // We need to create storage for -Div[u*s]. We can just
        // use tmpStorage since it is already lying around.
        LevelData<FArrayBox>& adv_src = tempStorage;
        adv_src.define(grids, 1);
        if (s_set_bogus_values) {
            setValLevel(adv_src, s_bogus_value);
        }

        // Calculate the advective source, -Div[scalarFlux].
        // DO NOT multiply by dt, LevelTGA does that for us.
        Divergence::levelDivergenceMAC(adv_src, scalarFlux, *m_levGeoPtr);
        for (dit.reset(); dit.ok(); ++dit) {
            adv_src[dit] *= -1.0;
        }

        // Add the contribution from bkgdSrc
        if (m_physBCPtr->useBackgroundScalar()) {
            for (dit.reset(); dit.ok(); ++dit) {
                adv_src[dit].plus(bkgdSrc[dit], 0, 0, 1);
            }
        }

        // If we are using a sponge layer, add its contribution to adv_src.
        if (m_physBCPtr->useSpongeLayer()) {
            LevelData<FArrayBox> spongeLayerSource(grids, 1);
            m_physBCPtr->fillSpongeLayerSrcTerm(spongeLayerSource,
                                                a_oldScalar,
                                                a_oldTime,
                                                a_dt,
                                                *m_levGeoPtr,
                                                a_comp);

            for (dit.reset(); dit.ok(); ++dit) {
                adv_src[dit].plus(spongeLayerSource[dit], 0, 0, 1);
            }
        }


        // Obtain a forward Euler estimate a_newScalar.
        for (dit.reset(); dit.ok(); ++dit) {
            // NOTE: Chombo examples never use the diffusive src term.
            a_newScalar[dit].copy(a_oldScalar[dit]);
            a_newScalar[dit].plus(adv_src[dit], a_dt);
            // a_newScalar[dit].plus(diffusiveSrc[dit], a_dt);
        }

        // 3. Increment the flux registers with advective flux only.
        // LevelTGA will increment the registers with the diffusive flux.
        if (s_advective_scalar_reflux && a_FRscale != 0.0) {
            updateScalarFluxRegister(scalarFlux, a_FRscale, a_comp);
        }

        // Clean up memory
        scalarFlux.clear();

        // 4. Perform the diffusive solve, if necessary.
        if (s_diffSolverScheme == ProblemContext::HeatSolverScheme::EXPLICIT) {
            // Just add the diffusive source term explicitly...

            // Compute the diffusive source term at half-time.
            const Real halfTime = a_oldTime + 0.5 * a_dt;
            const Real newTime = a_oldTime + a_dt;
            {
                LevelData<FArrayBox> halfScalar(grids, 1, IntVect::Unit);
                timeInterp(halfScalar, halfTime,
                           a_oldScalar, a_oldTime,
                           a_newScalar, newTime,
                           Interval(0,0));

                BCMethodHolder scalBC = m_physBCPtr->diffusiveSourceFuncBC();
                setGhostsScalar(halfScalar, scalBC, halfTime, a_comp);

                // Set up crse level BC.
                // Remember, coarse level may be at a more advanced time than this level.
                LevelData<FArrayBox>* crseDataPtr = NULL;
                if (m_level > 0) {
                    CH_assert(oldCrseScalPtr != NULL);
                    CH_assert(newCrseScalPtr != NULL);

                    const DisjointBoxLayout& crseGrids = crseNSPtr()->newVel().getBoxes();
                    CH_assert(crseGrids == m_levGeoPtr->getCoarserPtr()->getBoxes());
                    crseDataPtr = new LevelData<FArrayBox>(crseGrids, 1);

                    timeInterp(*crseDataPtr, halfTime,
                               *oldCrseScalPtr, oldCrseTime,
                               *newCrseScalPtr, newCrseTime,
                               Interval(0,0));
                }

                // Compute the diffusive source term, coeff * L[scalar].
                // (The op takes care of the exchanges.)
                this->computeLapScal(diffusiveSrc, halfScalar, crseDataPtr, NULL);
                for (dit.reset(); dit.ok(); ++dit) {
                    diffusiveSrc[dit] *= s_scal_coeffs[a_comp];
                }

                // Free memory
                delete crseDataPtr;
                crseDataPtr = NULL;
            }

            // Add diffusive source term.
            for (dit.reset(); dit.ok(); ++dit) {
                a_newScalar[dit].plus(diffusiveSrc[dit], a_dt);
            }

        } else {
            // More complicated, implicit diffusive solve...

            // Set for a 2-level solve
            int numberMGlevels = (m_level == 0) ? 0 : 1;

            // The LevelTGA solver requires this pointer
            MappedLevelFluxRegister* fineScalFluxRegPtr = NULL;
            MappedLevelFluxRegister* crseScalFluxRegPtr = NULL;
            if (s_diffusive_scalar_reflux) {
                if (!finestLevel()) {
                    fineScalFluxRegPtr = m_scal_fluxreg_ptrs[a_comp];
                    CH_assert(fineScalFluxRegPtr->isDefined());
                }
                if (m_level > 0) {
                    crseScalFluxRegPtr = (crseNSPtr()->m_scal_fluxreg_ptrs[a_comp]);
                    CH_assert(crseScalFluxRegPtr->isDefined());
                }
            }

            // The LevelTGA solver may emit convergence data.
            // Give the user some context for this data.
            if (s_verbosity >= 1) {
                pout() << "Level " << m_level << " diff solve on comp " << a_comp << ": " << flush;
            }

            // Set new time in all operators that will need to set BCs during solve.
            const Real newTime = a_oldTime + a_dt;
            Vector<MappedMGLevelOp<LevelData<FArrayBox> > * > allOps = m_diffAMRMGPtrs[a_comp]->getAllOperators();
            for (int idx = 0; idx < allOps.size(); ++idx) {
                MappedAMRPoissonOp* thisOp = dynamic_cast<MappedAMRPoissonOp*>(allOps[idx]);
                thisOp->setTime(newTime);
            }

            { // Do the diffusive update
                int ncomp = a_oldScalar.nComp();
                LevelData<FluxBox> diffFlux(grids, ncomp, IntVect::Zero); // The IntVect::Zero is important here!
                m_diffSolverPtrs[a_comp]->updateSoln(a_newScalar,
                                                     a_oldScalar,
                                                     adv_src,
                                                     diffFlux,
                                                     fineScalFluxRegPtr,
                                                     crseScalFluxRegPtr,
                                                     oldCrseScalPtr,
                                                     newCrseScalPtr,
                                                     a_oldTime,
                                                     oldCrseTime,
                                                     newCrseTime,
                                                     a_dt,
                                                     numberMGlevels,
                                                     false,  // zero-out a_newScalar?
                                                     true,   // already kappa weighted?
                                                     0);     // flux reg start comp
            }

            // Finally, compute the RHS via (a_newScalar - a_oldScalar) / a_dt
            for (dit.reset(); dit.ok(); ++dit) {
                a_rhs[dit].copy(a_newScalar[dit]);
                a_rhs[dit].plus(a_oldScalar[dit], -1.0);
                a_rhs[dit] /= a_dt;
            }
        } // end if not / if diffusive

    } else {
        // Explicit update...

        // Get the BC holders
        Tuple<BCMethodHolder,SpaceDim> BCValues = m_physBCPtr->scalarRiemannBC(a_comp);
        Tuple<BCMethodHolder,SpaceDim> BCSlopes = m_physBCPtr->scalarSlopeBC(a_comp);

        // Just perform tracing with diffusiveSrc = bkgdSrc.
        m_advectUtilScal.predictScalar(scalarFlux,
                                       a_oldScalar,
                                       (s_useScalAdvectiveSource? &bkgdSrc: NULL),
                                       a_oldVel,
                                       a_advVel,
                                       a_dt,
                                       *m_levGeoPtr,
                                       BCValues,
                                       BCSlopes,
                                       a_oldTime,
                                       true); // return fluxes

        // Calculate -div(u*s) from the fluxes
        Divergence::levelDivergenceMAC(a_rhs, scalarFlux, *m_levGeoPtr);
        for (dit.reset(); dit.ok(); ++dit) {
            a_rhs[dit] *= -1.0;
        }

        // If we are using a sponge layer, add its contribution to a_rhs.
        if (m_physBCPtr->useSpongeLayer()) {
            LevelData<FArrayBox> spongeLayerSource(grids, 1);
            m_physBCPtr->fillSpongeLayerSrcTerm(spongeLayerSource,
                                                a_oldScalar,
                                                a_oldTime,
                                                a_dt,
                                                *m_levGeoPtr,
                                                a_comp);

            for (dit.reset(); dit.ok(); ++dit) {
                a_rhs[dit].plus(spongeLayerSource[dit], 0, 0, 1);
            }
        }

        // Add the background advective source term
        if (m_physBCPtr->useBackgroundScalar()) {
            for (dit.reset(); dit.ok(); ++dit) {
                a_rhs[dit].plus(bkgdSrc[dit]);
            }
        }

        // Calculation of a_rhs is complete. Update the scalar.
        for (dit.reset(); dit.ok(); ++dit) {
            a_newScalar[dit].copy(a_oldScalar[dit]);
            a_newScalar[dit].plus(a_rhs[dit], a_dt);
        }

        // Do flux register stuff.
        if (s_advective_scalar_reflux) {
            updateScalarFluxRegister(scalarFlux, a_dt, a_comp);
        }
    } // end if not diffusive
}


// -----------------------------------------------------------------------------
// Calculate d(velocity)/dt = RHS centered at the time of a_oldVel + a_dt/2
// -----------------------------------------------------------------------------
void AMRNavierStokes::getNewVelocity (LevelData<FArrayBox>&       a_rhs,
                                      LevelData<FArrayBox>&       a_newVel,
                                      LevelData<FArrayBox>&       a_oldVel,
                                      const LevelData<FluxBox>&   a_advVel,
                                      const Real                  a_oldTime,
                                      const Real                  a_dt,
                                      const Real                  a_FRscale)
{
    CH_TIME("AMRNavierStokes::getNewVelocity");

    // Sanity checks
    CH_assert(a_rhs   .nComp() == SpaceDim);
    CH_assert(a_newVel.nComp() == SpaceDim);
    CH_assert(a_oldVel.nComp() == SpaceDim);
    CH_assert(a_advVel.nComp() == 1);

    // Gather some needed info
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    DataIterator dit = grids.dataIterator();
    const RealVect& dx = m_levGeoPtr->getDx();
    const bool isViscous = (s_nu > 0.0);
    const Real half_time = a_oldTime + 0.5 * a_dt;
    const Real newTime = a_oldTime + a_dt;

    // Initialize the outputs
    if (s_set_bogus_values) {
        setValLevel(a_rhs, s_bogus_value);
        setValLevel(a_newVel, s_bogus_value);
    }

    // The nonlinear advection term will be placed in a_rhs.
    LevelData<FArrayBox>& adv_term = a_rhs;

    if (s_nonlinearDifferencingForm != ProblemContext::NonlinearDifferencingForm::NONE) {
        // 1. Predict velocities ---------------------------------------------------
        // This is the FC, half-time, vector scaled version of the advecting
        // velocity with the normal and transverse components corrected with phi and
        // NOT corrected for freestream preservation.

        // Initialize the data holder
        LevelData<FluxBox> pred_vel(grids, SpaceDim); // Changed from m_tracingGhosts to 0
        if (s_set_bogus_values) {
            setValLevel(pred_vel, s_bogus_value);
        }

        // Do the prediction. All inputs and outputs of this function
        // are in the mapped basis.
        predictVelocities(pred_vel, a_oldVel, a_advVel, a_oldTime, a_dt);

        // It seems copying the normal components of a_advVel into pred_vel
        // (appropriately scaled, etc) is more stable than using the normal
        // components that come out of predictVelocities. The differences are...
        // Using pred_vel:
        // 1. Was generated with the projected, VD-corrected, half time-centered
        //    a_advVel.
        // 2. Is not projected, but is corrected using gradPhi from the a_advVel
        //    projection.
        // Using a_advVel:
        // 1. Was generated with the approximately projected, non VD-corrected,
        //    a_oldTime centered, Av[a_oldVel].
        // 2. Was explicitly MAC projected.
        for (dit.reset(); dit.ok(); ++dit) {
            for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
                FArrayBox& predVelFAB = pred_vel[dit][dir];
                const FArrayBox& advVelFAB = a_advVel[dit][dir];

                // Copy the normal component of the advecting velocity.
                predVelFAB.copy(advVelFAB, 0, dir, 1);

                // Remove the VD correction if necessary.
                if (s_etaLambda > 0.0 && s_applyFreestreamCorrection) {
                    const FArrayBox& corrFAB = m_gradELambda[dit][dir];
                    predVelFAB.minus(corrFAB, 0, dir, 1);
                }

                // Scale as a vector.
                m_levGeoPtr->divByJ(predVelFAB, dit(), dir);
            }
        }

        // Apply the projection correction (subtract Grad[phi]).
        // We are using phi because this is an FC velocity field!
        // NOTE: phi = dt * macPressure
        if (s_isIncompressible) {
            // Fill with grad[phi]
            LevelData<FluxBox> gradPhi(grids, SpaceDim);
            gradMACPressure(gradPhi, 0.5*a_dt);
            // m_projection.gradPhi(gradPhi);

            // Loop over grids and apply lagged correction
            for (dit.reset(); dit.ok(); ++dit) {
                // Create references for convenience
                FluxBox& predVelFB = pred_vel[dit];
                FluxBox& gradPhiFB = gradPhi[dit];

                // pred_vel is a vector, but gradPhi was given to us as a flux.
                // So, convert gradPhi to a vector.
                m_levGeoPtr->divByJ(gradPhiFB, dit());

                // Loop over FC directions and apply the correction to (u,v,w).
                // NOTE: This does not need to be done to the normal component!
                for (int FCdir = 0; FCdir < CH_SPACEDIM; ++FCdir) {
                    for (int comp = 0; comp < CH_SPACEDIM; ++comp) {
                        if (comp == FCdir) continue;
                        predVelFB[FCdir].minus(gradPhiFB[FCdir], comp, comp, 1);
                    }
                }
            }
        }

        // Set physical BCs on predicted velocity.
        if (isViscous) {
            EdgeVelBCHolder velBC(m_physBCPtr->viscousVelFuncBC());
            velBC.setGhosts(pred_vel,
                            NULL,
                            dx,
                            &(m_levGeoPtr->getFCJgup()),
                            false,  // not homogeneous
                            half_time);
        } else {
            EdgeVelBCHolder velBC(m_physBCPtr->tracingVelFuncBC());
            velBC.setGhosts(pred_vel,
                            NULL,
                            dx,
                            &(m_levGeoPtr->getFCJgup()),
                            false,  // not homogeneous
                            half_time);
        }

        // We are done calculating the predicted velocity.
        // Convert pred_vel back to the Cartesian basis.
        m_levGeoPtr->sendToCartesianBasis(pred_vel, false); // Changed from true on May 9, 2013.

        // The pred_vel computation is complete.


        // 2. Compute -u.Grad[u] ---------------------------------------------------
        // To do this, we need a vector-scaled, mapped basis, CC version of adv_vel.
        // The vector scaling is important because it will bring d to a physical
        // scale. Also, pred_vel should be in the Cartesian basis here.
        switch (s_nonlinearDifferencingForm) {
        case ProblemContext::NonlinearDifferencingForm::CONSERVATIVE:
            {
                // Loop over grids and multiply the Cartesian based pred_vel
                // by the mapping based a_advVel. This will give us the
                // fluxes adv_vel * pred_vel.
                for (dit.reset(); dit.ok(); ++dit) {
                    for (int dir = 0; dir < SpaceDim; ++dir) {
                        for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                            pred_vel[dit][dir].mult(a_advVel[dit][dir], 0, velComp, 1);
                        }
                    }
                }
                // Now, pred_vel contains the momentum fluxes
                //                   comp 0                comp 1                comp 2
                // FC dir 0  (J*Uad[0]*pred_vel[0], J*Uad[0]*pred_vel[1], J*Uad[0]*pred_vel[2])
                // FC dir 1  (J*Uad[1]*pred_vel[0], J*Uad[1]*pred_vel[1], J*Uad[1]*pred_vel[2])
                // FC dir 2  (J*Uad[2]*pred_vel[0], J*Uad[2]*pred_vel[1], J*Uad[2]*pred_vel[2])

                // Compute adv_term = -Div[fluxes]
                if (s_set_bogus_values) {
                    setValLevel(adv_term, s_bogus_value);
                }

                Divergence::levelDivergenceMAC(adv_term, pred_vel, *m_levGeoPtr, half_time, NULL);

                for (dit.reset(); dit.ok(); ++dit) {
                    adv_term[dit] *= -1.0;
                }

                if (s_advective_momentum_reflux) {
                    // Use these fluxes to update the momentum flux registers.
                    updateVelFluxRegister(pred_vel, a_FRscale);
                }
            }
            break;
        case ProblemContext::NonlinearDifferencingForm::ADVECTIVE:
            {
                // Initialize adv_term
                if (s_set_bogus_values) {
                    setValLevel(adv_term, s_bogus_value);
                }

                // Compute Av[adv_vel / J]. Use a_newVel as a temp holder.
                LevelData<FArrayBox>& half_vel = a_newVel;
                EdgeToCell(a_advVel, half_vel);
                m_levGeoPtr->divByJ(half_vel);

                // Compute the negative of the advective term, half_vel.Grad[pred_vel]
                Gradient::levelCCDotGradFC(adv_term, half_vel, pred_vel, *m_levGeoPtr);

                // Set adv_vel to -u.Grad[u] to mimic other source terms.
                for (dit.reset(); dit.ok(); ++dit) {
                    adv_term[dit] *= -1.0;
                }

                if (s_advective_momentum_reflux) {
                    // Loop over grids and multiply the Cartesian based pred_vel
                    // by the mapping based a_advVel. This will give us the
                    // fluxes adv_vel * pred_vel.
                    for (dit.reset(); dit.ok(); ++dit) {
                        for (int dir = 0; dir < SpaceDim; ++dir) {
                            for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                                pred_vel[dit][dir].mult(a_advVel[dit][dir], 0, velComp, 1);
                            }
                        }
                    }

                    // Now, pred_vel contains the momentum fluxes
                    //                   comp 0                comp 1                comp 2
                    // FC dir 0  (J*Uad[0]*pred_vel[0], J*Uad[0]*pred_vel[1], J*Uad[0]*pred_vel[2])
                    // FC dir 1  (J*Uad[1]*pred_vel[0], J*Uad[1]*pred_vel[1], J*Uad[1]*pred_vel[2])
                    // FC dir 2  (J*Uad[2]*pred_vel[0], J*Uad[2]*pred_vel[1], J*Uad[2]*pred_vel[2])

                    // Use these fluxes to update the momentum flux registers.
                    updateVelFluxRegister(pred_vel, a_FRscale);
                }
            }
            break;
        default:
            MayDay::Error("Unknown nonlinear differencing form");
        };

    } else {
        // Don't use nonlinear advection -- set to zero.
        setValLevel(adv_term, 0.0);
    }

    // At this point, we've computed the advective source term, updated the
    // flux registers with the momentum flux, and cleaned up after ourselves.
    // Don't forget that adv_term is using the a_rhs holder!


    // 3. Compute uStar --------------------------------------------------------
    // This is the approximation of the new time CC velocity field.

    // We are about to perform a forward time step. From here on,
    // sources should be at half_time.

    // If there is gravity, we need to add its contribution to the forcing.
    // To do this, we can simply add the grav source term to adv_term.
    if (s_gravityMethod == ProblemContext::GravityMethod::EXPLICIT) {
        // Sanity check
        CH_assert(adv_term.nComp() == SpaceDim);

        // Compute the gravitational source term.
        LevelData<FArrayBox> gravSource(grids, SpaceDim);
        this->fillGravSource(gravSource, half_time,
                             false);  // add background?

        // Combine the advective and gravitational source terms.
        for (dit.reset(); dit.ok(); ++dit) {
            adv_term[dit].plus(gravSource[dit], 0, 0, SpaceDim);
        }
    }

    // Compute the tidal forcing term
    if (s_tidalU0 * s_tidalOmega != 0.0) {
        // Compute the gravitational source term.
        LevelData<FArrayBox> tidalSource(grids, SpaceDim);
        this->fillTidalSource(tidalSource, a_oldTime, a_dt);

        // Combine the advective and gravitational source terms.
        for (dit.reset(); dit.ok(); ++dit) {
            adv_term[dit].plus(tidalSource[dit], 0, 0, SpaceDim);
        }
    }

    // If we are using a sponge layer, add its contribution.
    if (m_physBCPtr->useSpongeLayer()) {
        LevelData<FArrayBox> cartOldVel(grids, SpaceDim);  // TODO: I think this can use the same space as spongeLayerSource.
        a_oldVel.copyTo(cartOldVel);
        m_levGeoPtr->sendToCartesianBasis(cartOldVel);

        LevelData<FArrayBox> spongeLayerSource(grids, SpaceDim);
        m_physBCPtr->fillSpongeLayerSrcTerm(spongeLayerSource,
                                            cartOldVel,
                                            a_oldTime,
                                            a_dt,
                                            *m_levGeoPtr);

        for (dit.reset(); dit.ok(); ++dit) {
            adv_term[dit].plus(spongeLayerSource[dit], 0, 0, SpaceDim);
        }
    }

    // Inviscid flow: Simple, explicit update...
    // Right now, a_rhs is du/dt = -(Uad.Grad)U + gravSource, and a_newVel
    // contains garbage. Just overwrite a_newVel with a forward Euler estimate.
    // NOTE: If we are dealing with viscous flow, this will serve as the new
    // velocity guess handed to the LevelTGA solver.

    // Initialize a_newVel.
    if (s_set_bogus_values) {
        setValLevel(a_newVel, s_bogus_value);
    }

    // Convert adv_term to a mapped basis to match a_oldVel.
    m_levGeoPtr->sendToMappedBasis(adv_term, true);

    // Perform the forward Euler timestep.
    for (dit.reset(); dit.ok(); ++dit) {
        a_newVel[dit].copy(a_oldVel[dit]);
        a_newVel[dit].plus(adv_term[dit], a_dt);
    }

    if (isViscous && s_viscSolverScheme == ProblemContext::HeatSolverScheme::EXPLICIT) {
        // Viscous flow: Simple, explicit update...

        // Compute half-time, Cartesian velocity
        LevelData<FArrayBox> cartVel(grids, SpaceDim, IntVect::Unit);
        {
            timeInterp(cartVel, half_time,
                       a_oldVel, a_oldTime,
                       a_newVel, newTime,
                       Interval(0, SpaceDim-1));

            VelBCHolder velBC(m_physBCPtr->viscousVelFuncBC());   // Changed from uDelUFuncBC
            velBC.setGhosts(cartVel,
                            NULL,
                            dx,
                            &(m_levGeoPtr->getFCJgup()),
                            false,  // not homogeneous
                            half_time);

            m_levGeoPtr->sendToCartesianBasis(cartVel, true);
        }

        // Compute half-time viscous source.
        LevelData<FArrayBox> viscSource(grids, SpaceDim);
        fillViscousSource(viscSource, cartVel, half_time);
        m_levGeoPtr->sendToMappedBasis(viscSource, false);

        // Add viscous source term.
        for (dit.reset(); dit.ok(); ++dit) {
            a_newVel[dit].plus(viscSource[dit], a_dt);
        }

    } else if (isViscous && s_viscSolverScheme != ProblemContext::HeatSolverScheme::EXPLICIT) {
        // Viscous flow: More complicated, semi-implicit update...
        // Right now, adv_term = a_rhs = du/dt = -(Uad.Grad)U + gravSource,
        // and a_newVel contains a guess to the new velocity solution.

        {// Reset BCs on old-time velocity
            VelBCHolder velBC(m_physBCPtr->viscousVelFuncBC());   // Changed from uDelUFuncBC
            velBC.setGhosts(a_oldVel,
                            NULL,
                            dx,
                            &(m_levGeoPtr->getFCJgup()),
                            false,  // not homogeneous
                            a_oldTime);
        }

        // We need to compute new_vel = old_vel + dt * [ adv_term - Grad(Pi) + viscousSource ]
        // where Grad[Pi] is the lagged pressure correction at t^{n-1/2}, but LevelTGA will
        // add the viscousSource for us. So we really just need to lump the adv_term and
        // Grad(Pi) terms into one holder to send the LevelTGA. To learn why we use the lagged
        // pressure correction, see JCOMP 168, 464–499 (2001), doi:10.1006/jcph.2001.6715.

        // Start by getting gradPi and converting it to a vector.
        LevelData<FArrayBox> gradPi;
        if (s_useLaggedPressure) {
            gradPi.define(grids, SpaceDim);
            gradCCPressure(gradPi, 1.0);
            m_levGeoPtr->divByJ(gradPi);

            DataIterator localDit = adv_term.dataIterator();
            for (localDit.reset(); localDit.ok(); ++localDit) {
                adv_term[localDit].plus(gradPi[localDit], -1.0);

                // NOTE: Chombo examples never use the viscous src or gradPi terms,
                // so I am commenting this out.
                // a_newVel[dit].plus(gradPi[dit], -a_dt);
            }
        }

        // Now, let's set up the viscous solve...
        const int numberMGlevels = (m_level == 0) ? 0 : 1;

        // Update fine flux registers if and only if NOT finest level
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


        // Convert everything to a Cartesian basis.
        // TODO: We should re-use the Jacobian matrix here.

        // First, the old velocity
        LevelData<FArrayBox> tmpOldVel(grids, SpaceDim, a_oldVel.ghostVect()); // TODO: Try with less ghosts.
        if (s_set_bogus_values) {
            setValLevel(tmpOldVel, s_bogus_value);
        }
        for (DataIterator localDit = tmpOldVel.dataIterator(); localDit.ok(); ++localDit) {
            tmpOldVel[localDit].copy(a_oldVel[localDit]);
        }
        m_levGeoPtr->sendToCartesianBasis(tmpOldVel, true);     // TODO: Try false.

        // Then the coarse-level velocities
        LevelData<FArrayBox>* tmpCrseOldVelPtr = NULL;
        LevelData<FArrayBox>* tmpCrseNewVelPtr = NULL;

        if (numberMGlevels > 0) {
            const AMRNavierStokes* crseLevelPtr = this->crseNSPtr();
            const LevelData<FArrayBox>& crseOldVel = *(crseLevelPtr->m_vel_old_ptr);
            const LevelData<FArrayBox>& crseNewVel = *(crseLevelPtr->m_vel_new_ptr);
            const DisjointBoxLayout& crseGrids = crseNewVel.getBoxes();

            tmpCrseOldVelPtr = new LevelData<FArrayBox>(crseGrids, SpaceDim, a_oldVel.ghostVect()); // TODO: Try with less ghosts.
            tmpCrseNewVelPtr = new LevelData<FArrayBox>(crseGrids, SpaceDim, a_newVel.ghostVect()); // TODO: Try with less ghosts.

            if (s_set_bogus_values) {
                setValLevel(*tmpCrseOldVelPtr, s_bogus_value);
                setValLevel(*tmpCrseNewVelPtr, s_bogus_value);
            }

            DataIterator crseDit = crseNewVel.dataIterator();
            for (crseDit.reset(); crseDit.ok(); ++crseDit) {
                (*tmpCrseOldVelPtr)[crseDit].copy(crseOldVel[crseDit]);
                (*tmpCrseNewVelPtr)[crseDit].copy(crseNewVel[crseDit]);
            }

            m_levGeoPtr->sendToCartesianBasis(*tmpCrseOldVelPtr, true);     // TODO: Try false.
            m_levGeoPtr->sendToCartesianBasis(*tmpCrseNewVelPtr, true);     // TODO: Try false.
        } // end convert coarse-level velocities

        // Then everything else
        m_levGeoPtr->sendToCartesianBasis(a_newVel, true);      // TODO: Try false.
        m_levGeoPtr->sendToCartesianBasis(adv_term, true);      // TODO: Try false.

        for (int comp = 0; comp < SpaceDim; comp++) {
            // Set new time in all operators that will need to set BCs during solve.
            Vector<MappedMGLevelOp<LevelData<FArrayBox> > * > allOps = m_viscAMRMGPtrs[comp]->getAllOperators();
            for (int idx = 0; idx < allOps.size(); ++idx) {
                MappedAMRPoissonOp* thisOp = dynamic_cast<MappedAMRPoissonOp*>(allOps[idx]);
                thisOp->setTime(newTime);
            }

            // Alias the needed velocity and source term components.
            const Interval intvl(comp, comp);

            LevelData<FArrayBox> compNewVelocity;
            aliasLevelData(compNewVelocity, &a_newVel, intvl);

            LevelData<FArrayBox> compOldVelocity;
            aliasLevelData(compOldVelocity, &tmpOldVel, intvl);

            LevelData<FArrayBox> compSrc;
            aliasLevelData(compSrc, &adv_term, intvl);

            if (s_verbosity >= 1) {
                pout() << "Level " << m_level << " visc solve on comp " << comp << ": " << flush;
            }

            if (numberMGlevels == 0) {
                // Bottom level - no coarser.
                int ncomp = compOldVelocity.nComp();
                LevelData<FluxBox> viscFlux(grids, ncomp, IntVect::Zero); // The IntVect::Zero is important here!
                m_viscSolverPtrs[comp]->updateSoln(compNewVelocity,
                                                   compOldVelocity,
                                                   compSrc,
                                                   viscFlux,
                                                   fineFluxRegisterPtr,
                                                   NULL,  // crse flux reg
                                                   NULL,  // crse old phi ptr
                                                   NULL,  // crse new phi ptr
                                                   a_oldTime,
                                                   crseOldTime,
                                                   crseNewTime,
                                                   a_dt,
                                                   numberMGlevels, // level
                                                   false, // zero-out compNewVelocity?
                                                   true,  // already kappa weighted?
                                                   comp); // flux reg start comp

            } else {
                // Coarser level exists - needed for BCs.
                LevelData<FArrayBox> compCrseOldVel;
                aliasLevelData(compCrseOldVel,
                               tmpCrseOldVelPtr,
                               intvl);

                LevelData<FArrayBox> compCrseNewVel;
                aliasLevelData(compCrseNewVel,
                               tmpCrseNewVelPtr,
                               intvl);

                // Do the solve
                int ncomp = compOldVelocity.nComp();
                LevelData<FluxBox> viscFlux(grids, ncomp, IntVect::Zero); // The IntVect::Zero is important here!
                m_viscSolverPtrs[comp]->updateSoln(compNewVelocity,
                                                   compOldVelocity,
                                                   compSrc,
                                                   viscFlux,
                                                   fineFluxRegisterPtr,
                                                   crseFluxRegisterPtr,
                                                   &compCrseOldVel,
                                                   &compCrseNewVel,
                                                   a_oldTime,
                                                   crseOldTime,
                                                   crseNewTime,
                                                   a_dt,
                                                   numberMGlevels,
                                                   false, // zero-out compNewVelocity?
                                                   true,  // already kappa weighted?
                                                   comp); // flux reg start comp
            } // end if numberMGlevels == 1
        } // end loop over velocity components (comp)

        // Free memory
        tmpOldVel.clear();
        if (numberMGlevels > 0) {
            delete tmpCrseOldVelPtr;
            delete tmpCrseNewVelPtr;
        }

        // Convert back to the mapped basis.
        m_levGeoPtr->sendToMappedBasis(a_newVel, true);    // TODO: Try with false.

        for (DataIterator localDit = a_newVel.dataIterator(); localDit.ok(); ++localDit) {
            if (s_useLaggedPressure) {
                // Remember: We computed a_newVel - dt * Grad[Pi],
                // where Grad[Pi] is the lagged pressure correction at t^{n-1/2}.
                // We now need to remove the Grad[Pi] term.
                a_newVel[localDit].plus(gradPi[localDit], a_dt);
            }

            // And finally, set the update, a_rhs = du/dt = (new_vel - old_vel) / dt
            a_rhs[localDit].copy(a_newVel[localDit]);
            a_rhs[localDit].plus(a_oldVel[localDit], -1.0);
            a_rhs[localDit] /= a_dt;
        }

    } else {
        // not viscous...
        // since adv_term is just an alias for a_rhs, we have nothing to do.
    }
}


// -----------------------------------------------------------------------------
// Uses the CC velocity at a_oldTime to predict each component at all FCs at
// time a_oldTime + 0.5 * a_dt. The exchange on the result is up to the caller.
// This function also does not perform a projection or VD correction.
// the input and output velocities are all in the mapped basis.
// -----------------------------------------------------------------------------
void AMRNavierStokes::predictVelocities (LevelData<FluxBox>&       a_predVel,
                                         LevelData<FArrayBox>&     a_oldVel,
                                         const LevelData<FluxBox>& a_advVel,
                                         const Real                a_oldTime,
                                         const Real                a_dt)
{
    CH_TIME("AMRNavierStokes::predictVelocities");

    // Declare variables
    const bool isViscous = (s_nu > 0.0);

    const RealVect& dx = m_levGeoPtr->getDx();
    const DisjointBoxLayout& grids = newVel().getBoxes();
    DataIterator dit = grids.dataIterator();

    LevelData<FArrayBox> srcTerms(grids, SpaceDim, m_tracingGhosts);
    LevelData<FArrayBox> cartOldVel(grids, SpaceDim, m_tracingGhosts);

    // Sanity checks
    CH_assert(a_predVel.nComp() == SpaceDim);
    CH_assert(a_oldVel .nComp() == SpaceDim);
    CH_assert(a_advVel .nComp() == 1);
    CH_assert(a_predVel.getBoxes().compatible(grids));
    CH_assert(a_oldVel .getBoxes().compatible(grids));
    CH_assert(a_advVel .getBoxes().compatible(grids));

    // Initialize fields with bogus data
    if (s_set_bogus_values) {
        setValLevel(a_predVel, s_bogus_value);
        setValLevel(srcTerms, s_bogus_value);
        setValLevel(cartOldVel, s_bogus_value);
    }

    // Set all ghosts on old velocity
    VelBCHolder velBC(m_physBCPtr->tracingVelFuncBC());
    setGhostsVelocity(a_oldVel, velBC, a_oldTime);

    // Convert a_oldVel to the Cartesian basis.
    // NOTE: For some reason, using a Copier here causes infs/nans.
    for (dit.reset(); dit.ok(); ++dit) {
        cartOldVel[dit].copy(a_oldVel[dit]);
    }
    m_levGeoPtr->sendToCartesianBasis(cartOldVel, true);


    // du/dt = -(Uad.Grad)U + [viscous source] + [grav source]
    // We need to estimate the sum of the source terms here.

    // Compute the viscous source term. This initializes srcTerms.
    if (isViscous) {
        this->fillViscousSource(srcTerms, cartOldVel, a_oldTime);
    } else {
        setValLevel(srcTerms, 0.0);
    }

    // Compute the gravitational source term.
    if (s_gravityMethod == ProblemContext::GravityMethod::EXPLICIT) {
        const bool addBackground = false;
        LevelData<FArrayBox> gravSource(grids, SpaceDim, m_tracingGhosts);
        this->fillGravSource(gravSource, a_oldTime, addBackground);

        // Note: fillGravSource produces a vector, not a flux,
        // which allows us to simply add it to the sources.
        for (dit.reset(); dit.ok(); ++dit) {
            srcTerms[dit].plus(gravSource[dit], 0, 0, SpaceDim);
        }
    }

    // Compute the tidal forcing term.
    if (s_tidalU0 * s_tidalOmega != 0.0) {
        // The a_dt passed in is used to compute a derivative.
        // This only needs to be first order, so passing in 0.5*a_dt ot a_dt
        // shouldn't make much of a difference.
        LevelData<FArrayBox> tidalSource(grids, SpaceDim, m_tracingGhosts);
        this->fillTidalSource(tidalSource, a_oldTime, 0.5*a_dt);

        // Note: fillTidalSource produces a vector, not a flux,
        // which allows us to simply add it to viscousSource.
        DataIterator dit = tidalSource.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            srcTerms[dit()].plus(tidalSource[dit()], 0, 0, SpaceDim);
        }
    }

    // If we are using a sponge layer, add its contribution.
    if (m_physBCPtr->useSpongeLayer()) {
        LevelData<FArrayBox> spongeLayerSource(grids, SpaceDim);
        m_physBCPtr->fillSpongeLayerSrcTerm(spongeLayerSource,
                                            cartOldVel,
                                            a_oldTime,
                                            a_dt,
                                            *m_levGeoPtr);

        for (dit.reset(); dit.ok(); ++dit) {
            srcTerms[dit].plus(spongeLayerSource[dit], 0, 0, SpaceDim);
        }
    }

    // These aren't a must, but they do make a difference.
    srcTerms.exchange(m_tracingExCopier);
    srcTerms.exchange(m_tracingExCornerCopier);

    // Loop through velocity components and compute the fluxes for each.
    for (int veldir = 0; veldir < CH_SPACEDIM; ++veldir) {
        const Interval velInt(veldir,veldir);

        // predictScalar can only put data into a one-component holder.
        // Since predVel is a SpaceDim-component holder, we need a temp
        // data holder for predictScalar to write to, then we can copy the
        // result to the appropriate component of predVel.
        //
        // NOTE: This would be more efficient if aliasLevelData worked
        // with T = FluxBox or if predictScalar had a destComp parameter.
        LevelData<FluxBox> predVelDir(grids, 1, m_tracingGhosts);     // TODO: Try with less ghosts
        if (s_set_bogus_values) {
            setValLevel(predVelDir, s_bogus_value);
        }

        // Alias Cartesian comp of initial velocity
        LevelData<FArrayBox> cartOldVelDir;
        aliasLevelData(cartOldVelDir, &cartOldVel, velInt);

        // Alias source term comp
        LevelData<FArrayBox> thisViscSrcDir;
        aliasLevelData(thisViscSrcDir, &srcTerms, velInt);

        // Get the BC holders
        Tuple<BCMethodHolder,SpaceDim> BCValues = m_physBCPtr->velRiemannBC(veldir, isViscous);
        Tuple<BCMethodHolder,SpaceDim> BCSlopes = m_physBCPtr->velSlopeBC(veldir, isViscous);

        // Do the tracing
        m_advectUtilVel.predictScalar(predVelDir,
                                      cartOldVelDir,
                                      (s_useVelAdvectiveSource? &thisViscSrcDir: NULL),
                                      a_oldVel,
                                      a_advVel,
                                      a_dt,
                                      *m_levGeoPtr,
                                      BCValues,
                                      BCSlopes,
                                      a_oldTime,
                                      false);         // Not returning fluxes

        // Copy to permanent holder
        for (dit.reset(); dit.ok(); ++dit) {
            for (int FCdir = 0; FCdir < CH_SPACEDIM; ++FCdir) {
                a_predVel[dit][FCdir].copy(predVelDir[dit][FCdir], 0, veldir, 1);
            }
        }
    }

    // Convert to a mapped basis.
    CH_assert(a_predVel.ghostVect() == IntVect::Zero);
    m_levGeoPtr->sendToMappedBasis(a_predVel, false);
}


// -----------------------------------------------------------------------------
// Performs a CC (approximate) level projection of a_vel
// -----------------------------------------------------------------------------
void AMRNavierStokes::doCCProjection (LevelData<FArrayBox>& a_vel,
                                      const Real            a_new_time,
                                      const Real            a_dt,
                                      const bool            a_doProj)
{
    CH_TIME("AMRNavierStokes::doCCProjection");

    // Set physical BCs
    const RealVect& dx = m_levGeoPtr->getDx();
    const bool isViscous = (s_nu > 0.0);
    VelBCHolder velBC(m_physBCPtr->uStarFuncBC(isViscous));
    velBC.setGhosts(a_vel,                          // state
                    NULL,                           // extrap
                    dx,                             // dx
                    &(m_levGeoPtr->getFCJgup()),    // Jgup
                    false,                          // isHomogeneous
                    a_new_time);                    // time

    // Just for good measure
    a_vel.exchange();

    // Ghosts are set. If we aren't projecting, then we are done here.
    if (!a_doProj || !s_isIncompressible) return;
    if (s_level_projection_iters == 0) return;

    // If a coarser level exists, will need coarse-level data for proj.
    // No need to gather coarse-level pressure pointer -- the projector has it.
    LevelData<FArrayBox>* crseVelPtr = NULL;
    if (m_level > 0) {
        // Get coarse level ptr
        AMRNavierStokes* crseLevPtr = crseNSPtr();
        CH_assert(crseLevPtr != NULL);

        // Coarse velocity BC data is interpolated in time...
        const DisjointBoxLayout& crseGrids = crseLevPtr->newVel().getBoxes();
        crseVelPtr = new LevelData<FArrayBox>(crseGrids, SpaceDim, IntVect::Unit);  // TODO: do we need a ghost?
        crseLevPtr->fillVelocity(*crseVelPtr, a_new_time);

        // ...and "unprojected."
        LevelData<FArrayBox> crseGradPres(crseGrids, SpaceDim);
        crseLevPtr->gradCCPressure(crseGradPres, 1.0);

        DataIterator cdit = crseGrids.dataIterator();
        for (cdit.reset(); cdit.ok(); ++cdit) {
            FArrayBox& crseVelFAB = (*crseVelPtr)[cdit];
            const FArrayBox& crseGradPresFAB = crseGradPres[cdit];

            crseVelFAB.plus(crseGradPresFAB, a_dt);
        }
    }

    // Package data for projector.
    Vector<LevelData<FArrayBox>*> amrVel(0);
    if (m_level > 0) {
        amrVel.push_back(crseVelPtr);
    }
    amrVel.push_back(&a_vel);

    // Project!
    pout() << "Level " << m_level << " CC proj:  " << flush;
    m_ccProjector.levelProject(amrVel,
                               m_levGeoPtr,
                               a_new_time,
                               a_dt,
                               false,       // is vel multiplied by J?
                               true,        // initialize pressure to zero? I see better results with 'true'.
                               false);      // force homog solve?

    // This projection (assuming it converged) validates the pressure.
    m_ccPressureState = CCPressureState::VALID;

    // Free memory
    if (crseVelPtr != NULL) {
        delete crseVelPtr;
        crseVelPtr = NULL;
    }
}

