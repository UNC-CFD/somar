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
#include "SetValLevel.H"
#include "BoxIterator.H"
#include "AlteredMetric.H"
#include "ProblemContext.H"
#include <iomanip>


// -----------------------------------------------------------------------------
// This does the actual computation to update the state variables.
// I've included this function since the initializeGlobalPressure() and advance()
// functions are, for the most part, the same piece of code.
// -----------------------------------------------------------------------------
void AMRNavierStokes::PPMIGTimeStep (const Real a_oldTime,
                                     const Real a_dt,
                                     const bool a_updatePassiveScalars,
                                     const bool a_doLevelProj)
{
    CH_TIME("AMRNavierStokes::PPMIGTimeStep");
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

    LevelData<FluxBox> adv_vel(grids, 1, IntVect::Unit); // Changed from m_tracingGhosts to 1
    computeAdvectingVelocities(adv_vel, old_vel, a_oldTime, a_dt);

    if (a_updatePassiveScalars) {
        // Lambda update
        LevelData<FArrayBox> old_lambda;
        fillLambda(old_lambda, a_oldTime);
        old_lambda.exchange(m_tracingExCopier);

        LevelData<FArrayBox> dLdt(grids, 1);
        getNewLambda(dLdt, new_lambda, old_lambda, old_vel, adv_vel, a_oldTime, a_dt, a_dt);
    }

    LevelData<FArrayBox> old_b;
    if (s_num_scal_comps > 0) {
        // Scalar update
        fillScalars(old_b, a_oldTime, 0);
        old_b.exchange(m_tracingExCopier);

        LevelData<FArrayBox> dSdt(grids, 1);
        getNewScalar(dSdt, new_b, old_b, old_vel, adv_vel, a_oldTime, a_dt, a_dt, 0);
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
    }

    if (s_num_scal_comps > 0) {
        // Do implicit gravity update and CC projection.
        doCCIGProjection(new_vel, new_b, old_vel, old_b, adv_vel, a_oldTime, a_dt, a_doLevelProj);
    } else {
        // Do a standard CC projection.
        doCCProjection(new_vel, a_oldTime + a_dt, a_dt, a_doLevelProj);
    }
}


// -----------------------------------------------------------------------------
// Semi-implicitly handles the gravity forcing and projection.
// The old* inputs should be the values at t^n.
// The new* inputs should be the updated values from the TGA solver.
// -----------------------------------------------------------------------------
void AMRNavierStokes::doCCIGProjection (LevelData<FArrayBox>&       a_newVel,
                                        LevelData<FArrayBox>&       a_newB,
                                        const LevelData<FArrayBox>& a_oldVel,
                                        const LevelData<FArrayBox>& a_oldB,
                                        const LevelData<FluxBox>&   a_advVel,
                                        const Real                  a_oldTime,
                                        const Real                  a_dt,
                                        const bool                  a_doProj)
{
    CH_TIME("AMRNavierStokes::doCCIGProjection");

    const Real halfTime = a_oldTime + 0.5 * a_dt;
    const Real newTime  = a_oldTime + 1.0 * a_dt;
    const Real dummyTime = -1.0e300;
    const GeoSourceInterface& geoSource = *(m_levGeoPtr->getGeoSourcePtr());
    const RealVect& dx = m_levGeoPtr->getDx();
    const DisjointBoxLayout& grids = a_newVel.getBoxes();
    DataIterator dit = grids.dataIterator();


    // 1. Compute the background buoyancy, N, Dinv, etc...

    // Fill the FC background buoyancy field.
    LevelData<FluxBox> bbar(grids, 1, 2*IntVect::Unit);
    for (dit.reset(); dit.ok(); ++dit) {
        FluxBox& bbarFB = bbar[dit];

        D_TERM(m_physBCPtr->setBackgroundScalar(bbarFB[0], 0, *m_levGeoPtr, dit(), dummyTime);,
               m_physBCPtr->setBackgroundScalar(bbarFB[1], 0, *m_levGeoPtr, dit(), dummyTime);,
               m_physBCPtr->setBackgroundScalar(bbarFB[2], 0, *m_levGeoPtr, dit(), dummyTime);)
    }

    // Fill the CC dXi^i/dz field.
    LevelData<FArrayBox> dXidz(grids, SpaceDim, 2*IntVect::Unit);
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& dXidzFAB = dXidz[dit];

        D_TERM(geoSource.fill_dXidx(dXidzFAB, 0, 0, SpaceDim-1, dx);,
               geoSource.fill_dXidx(dXidzFAB, 1, 1, SpaceDim-1, dx);,
               geoSource.fill_dXidx(dXidzFAB, 2, 2, SpaceDim-1, dx);)
    }

    // Fill the CC dz/dXi^i field.
    LevelData<FArrayBox> dzdXi(grids, SpaceDim, 2*IntVect::Unit);
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& dzdXiFAB = dzdXi[dit];

        D_TERM(geoSource.fill_dxdXi(dzdXiFAB, 0, SpaceDim-1, 0, dx);,
               geoSource.fill_dxdXi(dzdXiFAB, 1, SpaceDim-1, 1, dx);,
               geoSource.fill_dxdXi(dzdXiFAB, 2, SpaceDim-1, 2, dx);)
    }

    // Compute CC Nsq and Dinv
    LevelData<FArrayBox> Nsq(grids, 1);
    LevelData<FArrayBox> Dinv(grids, 1);
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& NsqFAB = Nsq[dit];
        FArrayBox& DinvFAB = Dinv[dit];
        const FluxBox& bbarFB = bbar[dit];
        const FArrayBox& dXidzFAB = dXidz[dit];
        const Box& region = grids[dit];

#       if CH_SPACEDIM == 2
            FORT_COMPUTENSQANDDINV2D(
                CHF_FRA1(NsqFAB,0),
                CHF_FRA1(DinvFAB,0),
                CHF_CONST_FRA1(bbarFB[0],0),
                CHF_CONST_FRA1(bbarFB[1],0),
                CHF_CONST_FRA(dXidzFAB),
                CHF_CONST_REAL(a_dt),
                CHF_CONST_REAL(s_gravityTheta),
                CHF_CONST_REALVECT(dx),
                CHF_BOX(region));
#       else
            FORT_COMPUTENSQANDDINV3D(
                CHF_FRA1(NsqFAB,0),
                CHF_FRA1(DinvFAB,0),
                CHF_CONST_FRA1(bbarFB[0],0),
                CHF_CONST_FRA1(bbarFB[1],0),
                CHF_CONST_FRA1(bbarFB[2],0),
                CHF_CONST_FRA(dXidzFAB),
                CHF_CONST_REAL(a_dt),
                CHF_CONST_REAL(s_gravityTheta),
                CHF_CONST_REALVECT(dx),
                CHF_BOX(region));
#       endif
    }


    // 2. Perform as much of the update as we can without the pressure.

    // Compute thetaVel = theta*newVel + (1-theta)*oldVel
    LevelData<FArrayBox> thetaVel(grids, SpaceDim);
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& thetaVelFAB = thetaVel[dit];
        const FArrayBox& newVelFAB = a_newVel[dit];
        const FArrayBox& oldVelFAB = a_oldVel[dit];
        const Box& region = thetaVelFAB.box();

        // Make sure regions are sized correctly
        CH_assert(newVelFAB.box().contains(region));
        CH_assert(oldVelFAB.box().contains(region));

        // Compute thetaVel
        FORT_WEIGHTEDAVG(
            CHF_FRA(thetaVelFAB),
            CHF_CONST_FRA(newVelFAB),
            CHF_CONST_FRA(oldVelFAB),
            CHF_CONST_REAL(s_gravityTheta),
            CHF_BOX(region));
    }

    // Compute thetaB = theta*newB + (1-theta)*oldB
    LevelData<FArrayBox> thetaB(grids, 1);
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& thetaBFAB = thetaB[dit];
        const FArrayBox& newBFAB = a_newB[dit];
        const FArrayBox& oldBFAB = a_oldB[dit];
        const Box& region = thetaBFAB.box();

        // Make sure regions are sized correctly
        CH_assert(newBFAB.box().contains(region));
        CH_assert(oldBFAB.box().contains(region));

        // Compute thetaB
        FORT_WEIGHTEDAVG(
            CHF_FRA(thetaBFAB),
            CHF_CONST_FRA(newBFAB),
            CHF_CONST_FRA(oldBFAB),
            CHF_CONST_REAL(s_gravityTheta),
            CHF_BOX(region));
    }

    // Update the velocity.
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& newVelFAB = a_newVel[dit];
        const FArrayBox& thetaBFAB = thetaB[dit];
        const FArrayBox& thetaVelFAB = thetaVel[dit];
        const FArrayBox& dXidzFAB = dXidz[dit];
        const FArrayBox& dzdXiFAB = dzdXi[dit];
        const FArrayBox& DinvFAB = Dinv[dit];
        const FArrayBox& NsqFAB = Nsq[dit];
        const Box& region = grids[dit];
        const Real dtTheta = a_dt * s_gravityTheta;

        BoxIterator bit(region);
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& cc = bit();

            Real W = D_TERM(  dzdXiFAB(cc,0) * thetaVelFAB(cc,0),
                            + dzdXiFAB(cc,1) * thetaVelFAB(cc,1),
                            + dzdXiFAB(cc,2) * thetaVelFAB(cc,2));

            Real btilde = thetaBFAB(cc) + dtTheta * NsqFAB(cc) * W;

            D_TERM(newVelFAB(cc,0) -= a_dt * DinvFAB(cc) * dXidzFAB(cc,0) * btilde;,
                   newVelFAB(cc,1) -= a_dt * DinvFAB(cc) * dXidzFAB(cc,1) * btilde;,
                   newVelFAB(cc,2) -= a_dt * DinvFAB(cc) * dXidzFAB(cc,2) * btilde;)
        }
    }


    // 3. Project the velocity.
    if (s_isIncompressible) {
        // const DisjointBoxLayout* crseGridsPtr = (m_level == 0)? NULL: &(crseNSPtr()->newVel().getBoxes());
        const DisjointBoxLayout* crseGridsPtr = NULL;
        if (m_level > 0) {
            crseGridsPtr = &(crseNSPtr()->newVel().getBoxes());
        }

        // Define the metric.
        m_alteredMetric.define(&geoSource, m_physBCPtr, a_dt*s_gravityTheta);

        // Install the new metric into the pressure solver. Unfortunately, we
        // need to do this at each timestep since the metric keeps changing.
        const LevelData<FArrayBox>* crsePressurePtr = NULL;
        if (m_level > 0) {
            crsePressurePtr = &(crseNSPtr()->m_ccPressure);
        }
        m_ccProjector.define(&m_ccPressure,
                             crsePressurePtr,
                             *m_physBCPtr,
                             *m_levGeoPtr,
                             &m_alteredMetric);

        // Perform the elliptic solve.
        doCCProjection(a_newVel, newTime, a_dt, a_doProj);

        // Restore the projector
        // TODO: Is this necessary?
        m_ccProjector.define(&m_ccPressure,
                             crsePressurePtr,
                             *m_physBCPtr,
                             *m_levGeoPtr,
                             NULL);
    }


    // 4. Perform the semi-implicit buoyancy update.

    // newB -= dt * Div[bbar*(theta*newVel + (1-theta)*oldVel)]
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& newBFAB = a_newB[dit];
        FArrayBox& thetaVelFAB = thetaVel[dit];
        const FArrayBox& newVelFAB = a_newVel[dit];
        const FArrayBox& oldVelFAB = a_oldVel[dit];
        const FArrayBox& NsqFAB = Nsq[dit];
        const FArrayBox& dzdXiFAB = dzdXi[dit];
        const Box& region = grids[dit];

        // Compute the CC advecting velocity, theta*newVel + (1-theta)*oldVel.
        FORT_WEIGHTEDAVG(
            CHF_FRA(thetaVelFAB),
            CHF_CONST_FRA(newVelFAB),
            CHF_CONST_FRA(oldVelFAB),
            CHF_CONST_REAL(s_gravityTheta),
            CHF_BOX(region));

        BoxIterator bit(region);
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& cc = bit();

            Real W = D_TERM(  dzdXiFAB(cc,0) * thetaVelFAB(cc,0),
                            + dzdXiFAB(cc,1) * thetaVelFAB(cc,1),
                            + dzdXiFAB(cc,2) * thetaVelFAB(cc,2));

            newBFAB(cc) += a_dt * NsqFAB(cc) * W;
        }
    }
}

