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
#include "timeInterp.H"
#include "MappedFineInterp.H"
#include "MappedPiecewiseLinearFillPatch.H"
#include "SetValLevel.H"
#include "ProblemContext.H"


// -----------------------------------------------------------------------------
// Sets all ghosts on velocity. (Does not exchange.)
// -----------------------------------------------------------------------------
void AMRNavierStokes::setGhostsVelocity (LevelData<FArrayBox>& a_vel,
                                         VelBCHolder&          a_velBC,
                                         const Real            a_time) const
{
    CH_TIME("AMRNavierStokes::setGhostsVelocity");

    // This function will NOT handle allocation
    CH_assert(a_vel.isDefined());
    CH_assert(a_vel.nComp() == SpaceDim);

    // Get a_vel's grids. These do not need to be the same as
    // the grids in m_levGeoPtr.
    const DisjointBoxLayout& grids = a_vel.getBoxes();

    // Fill ghost cells by linear interpolation of coarse data
    if (m_level > 0) {
        // Gather coarse level data
        AMRNavierStokes& crseLevel = *crseNSPtr();
        LevelData<FArrayBox>& oldCrseVel = crseLevel.oldVel();
        LevelData<FArrayBox>& newCrseVel = crseLevel.newVel();
        const DisjointBoxLayout& crseGrids = oldCrseVel.getBoxes();
        const ProblemDomain& crseDomain = crseLevel.problemDomain();
        const IntVect& nRefCrse = crseLevel.refRatio();

        const Real crse_new_time = crseLevel.m_time;
        const Real crse_dt = crseLevel.dt();
        const Real crse_old_time = crse_new_time - crse_dt;
        Real crse_time_interp_coeff;

        // Scale crse_time_interp_coeff as 0 for old time, 1 for new time.
        if (abs(a_time - crse_old_time) < TIME_EPS) {
            crse_time_interp_coeff = 0.0;
        } else if (abs(a_time - crse_new_time) < TIME_EPS) {
            crse_time_interp_coeff = 1.0;
        } else {
            crse_time_interp_coeff = (a_time - crse_old_time) / crse_dt;
        }

#       ifndef NDEBUG
            // This is useful for debugging erroneous time_interp_coeffs
            // which has been a difficult bug for me in the past.
            if (crse_time_interp_coeff < 0. || crse_time_interp_coeff > 1.) {
                pout() << "\ncrse_time_interp_coeff = " << crse_time_interp_coeff
                       << "\n\tcrse_old_time = " << crse_old_time
                       << "\n\tm_time = " << m_time
                       << "\n\ta_time = " << a_time
                       << "\n\tcrse_dt = " << crse_dt
                       << std::endl;
            }

            CH_assert(crse_time_interp_coeff >= 0.);
            CH_assert(crse_time_interp_coeff <= 1.);
#       endif

        // Define the ghost interpolator
        MappedPiecewiseLinearFillPatch filpatcher(grids, crseGrids,
                                                  SpaceDim, crseDomain,
                                                  nRefCrse, a_vel.ghostVect());

        // Fill in CF BC data by conservative linear interpolation.
        filpatcher.fillInterp(a_vel, oldCrseVel, newCrseVel,
                              crse_time_interp_coeff, 0, 0, SpaceDim);
    }

    // Apply physical BCs to our new data..

    // Start by collecting an appropriate Jgup.
    // If Jgup's grids are not the same as a_vel's grids, then we cannot use
    // the cached metric data. We must allocate a new data holder and fill
    // it over a_vel's grids.
    RefCountedPtr<LevelData<FluxBox> > JgupPtr = m_levGeoPtr->getFCJgupPtr();
    CH_assert(JgupPtr->getBoxes() == grids);

    // We can now use JgupPtr to apply the BCs.
    a_velBC.setGhosts(a_vel,          // state
                      NULL,           // extrap
                      m_levGeoPtr->getDx(),  // dx
                      &(*JgupPtr),    // JgupPtr
                      false,          // inhomogeneous
                      a_time);        // time
}


// -----------------------------------------------------------------------------
// Sets all ghosts on lambda. (Does not exchange.)
// -----------------------------------------------------------------------------
void AMRNavierStokes::setGhostsLambda (LevelData<FArrayBox>& a_lambda,
                                       BCMethodHolder&       a_lambdaBC,
                                       const Real            a_time) const
{
    CH_TIME("AMRNavierStokes::setGhostsLambda");

    // Collect some needed data
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    DataIterator dit = grids.dataIterator();
    const ProblemDomain& domain = this->problemDomain();
    const IntVect tracingGhosts(D_DECL(ADVECT_GROW, ADVECT_GROW, ADVECT_GROW));

    const Real old_time = m_time - m_dt;
    const Real new_time = m_time;

    // This function will NOT handle allocation
    CH_assert(a_lambda.isDefined());

    CH_assert(grids.compatible(m_levGeoPtr->getBoxes()));
    CH_assert(grids == m_levGeoPtr->getBoxes());
    CH_assert(grids.compatible(a_lambda.getBoxes()));
    CH_assert(grids == a_lambda.getBoxes());

    CH_assert(a_lambda.ghostVect() == tracingGhosts);

    CH_assert(old_time - TIME_EPS < a_time);
    CH_assert(a_time < new_time + TIME_EPS);

    // Fill ghost cells by linear interpolation of coarse data
    if (m_level > 0) {
        // Gather coarse level data
        AMRNavierStokes& crseLevel = *crseNSPtr();
        LevelData<FArrayBox>& oldCrseLambda = crseLevel.oldLambda();
        LevelData<FArrayBox>& newCrseLambda = crseLevel.newLambda();
        const DisjointBoxLayout& crseGrids = oldCrseLambda.getBoxes();
        const ProblemDomain& crseDomain = crseLevel.problemDomain();
        const IntVect& nRefCrse = crseLevel.refRatio();

        const Real crse_new_time = crseLevel.m_time;
        const Real crse_dt = crseLevel.dt();
        const Real crse_old_time = crse_new_time - crse_dt;
        Real crse_time_interp_coeff;

        // Scale time_interp_coeff as 0 for old time, 1 for new time.
        if (abs(a_time - crse_old_time) < TIME_EPS) {
            crse_time_interp_coeff = 0.0;
        } else if (abs(a_time - crse_new_time) < TIME_EPS) {
            crse_time_interp_coeff = 1.0;
        } else {
            crse_time_interp_coeff = (a_time - crse_old_time) / crse_dt;
        }

#       ifndef NDEBUG
            // This is useful for debugging erroneous time_interp_coeffs
            // which has been a difficult bug for me in the past.
            if (crse_time_interp_coeff < 0. || crse_time_interp_coeff > 1.) {
                pout() << "\ncrse_time_interp_coeff = " << crse_time_interp_coeff
                       << "\n\tcrse_old_time = " << crse_old_time
                       << "\n\tm_time (new time) = " << new_time
                       << "\n\ta_time = " << a_time
                       << "\n\tcrse_dt = " << crse_dt
                       << std::endl;
            }

            CH_assert(crse_time_interp_coeff >= 0.);
            CH_assert(crse_time_interp_coeff <= 1.);
#       endif

        // Define the ghost interpolator
        MappedPiecewiseLinearFillPatch filpatcher(grids, crseGrids,
                                                  a_lambda.nComp(), crseDomain,
                                                  nRefCrse, a_lambda.ghostVect());

        // Fill in CF BC data by conservative linear interpolation.
        filpatcher.fillInterp(a_lambda, oldCrseLambda, newCrseLambda,
                              crse_time_interp_coeff, 0, 0, a_lambda.nComp());
    }

    // Apply physical BCs to our new data
    const int nJgupComp = m_levGeoPtr->getFCJgup().nComp();
    for (dit.reset(); dit.ok(); ++dit) {
        // const FluxBox& JgupFB = m_levGeoPtr->getFCJgup()[dit];

        Box region = grow(grids[dit], ADVECT_GROW) & a_lambda[dit].box();
        FluxBox JgupFB(region, nJgupComp);
        m_levGeoPtr->fill_Jgup(JgupFB);

        a_lambdaBC.setGhosts(a_lambda[dit], // stateFAB
                             NULL,          // extrapFABPtr
                             a_lambda[dit].box() & domain.domainBox(),    // valid
                             domain,        // domain
                             m_levGeoPtr->getDx(),    // dx
                             dit(),         // DataIndex
                             &JgupFB,       // JgupFBPtr
                             false,         // isHomogeneous
                             a_time);       // time
    }
}


// -----------------------------------------------------------------------------
// Sets all ghosts on a scalar. (Does not exchange.)
// -----------------------------------------------------------------------------
void AMRNavierStokes::setGhostsScalar (LevelData<FArrayBox>& a_scal,
                                       BCMethodHolder&       a_scalBC,
                                       const Real            a_time,
                                       const int             a_comp) const
{
    CH_TIME("AMRNavierStokes::setGhostsScalar");

    // Collect some needed data
    const ProblemDomain& domain = m_levGeoPtr->getDomain();
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    DataIterator dit = grids.dataIterator();

    const Real old_time = m_time - m_dt;
    const Real new_time = m_time;

    // This function will NOT handle allocation
    CH_assert(a_scal.isDefined());

    CH_assert(grids.compatible(m_levGeoPtr->getBoxes()));
    CH_assert(grids == m_levGeoPtr->getBoxes());
    CH_assert(grids.compatible(a_scal.getBoxes()));
    CH_assert(grids == a_scal.getBoxes());

    CH_assert(old_time - TIME_EPS < a_time);
    CH_assert(a_time < new_time + TIME_EPS);

    // Fill ghost cells by linear interpolation of coarse data
    if (m_level > 0) {
        // Gather coarse level data
        AMRNavierStokes& crseLevel = *crseNSPtr();
        LevelData<FArrayBox>& oldCrseScal = crseLevel.oldScal(a_comp);
        LevelData<FArrayBox>& newCrseScal = crseLevel.newScal(a_comp);
        const DisjointBoxLayout& crseGrids = oldCrseScal.getBoxes();
        const ProblemDomain& crseDomain = crseLevel.problemDomain();
        const IntVect& nRefCrse = crseLevel.refRatio();

        const Real crse_new_time = crseLevel.m_time;
        const Real crse_dt = crseLevel.dt();
        const Real crse_old_time = crse_new_time - crse_dt;
        Real crse_time_interp_coeff;

        // Scale time_interp_coeff as 0 for old time, 1 for new time.
        if (abs(a_time - crse_old_time) < TIME_EPS) {
            crse_time_interp_coeff = 0.0;
        } else if (abs(a_time - crse_new_time) < TIME_EPS) {
            crse_time_interp_coeff = 1.0;
        } else {
            crse_time_interp_coeff = (a_time - crse_old_time) / crse_dt;
        }

#       ifndef NDEBUG
            // This is useful for debugging erroneous time_interp_coeffs
            // which has been a difficult bug for me in the past.
            if (crse_time_interp_coeff < 0. || crse_time_interp_coeff > 1.) {
                pout() << "\ncrse_time_interp_coeff = " << crse_time_interp_coeff
                       << "\n\tcrse_old_time = " << crse_old_time
                       << "\n\tm_time (new time) = " << new_time
                       << "\n\ta_time = " << a_time
                       << "\n\tcrse_dt = " << crse_dt
                       << std::endl;
            }

            CH_assert(crse_time_interp_coeff >= 0.);
            CH_assert(crse_time_interp_coeff <= 1.);
#       endif

        // Define the ghost interpolator
        MappedPiecewiseLinearFillPatch filpatcher(grids, crseGrids,
                                                  a_scal.nComp(), crseDomain,
                                                  nRefCrse, a_scal.ghostVect());

        // Fill in CF BC data by conservative linear interpolation.
        filpatcher.fillInterp(a_scal, oldCrseScal, newCrseScal,
                              crse_time_interp_coeff, 0, 0, a_scal.nComp());
    }

    const int nJgupComp = m_levGeoPtr->getFCJgup().nComp();
    for (dit.reset(); dit.ok(); ++dit) {
        // const FluxBox& JgupFB = m_levGeoPtr->getFCJgup()[dit];

        Box region = grow(grids[dit], a_scal.ghostVect()) & a_scal[dit].box();
        FluxBox JgupFB(region, nJgupComp);
        m_levGeoPtr->fill_Jgup(JgupFB);

        a_scalBC.setGhosts(a_scal[dit],   // stateFAB
                           NULL,          // extrapFABPtr
                           a_scal[dit].box() & domain.domainBox(),    // valid
                           domain,        // domain
                           m_levGeoPtr->getDx(),    // dx
                           dit(),         // DataIndex
                           &JgupFB,       // JgupFBPtr
                           false,         // isHomogeneous
                           a_time);       // time
    }
}


// -----------------------------------------------------------------------------
// Returns velocity at time a_time.
// This does a linear interpolation in time between old and new time velocities,
// but does not fill ghosts or interpolate in space from a coarser level.
// -----------------------------------------------------------------------------
void AMRNavierStokes::velocity (LevelData<FArrayBox>& a_vel,
                                Real                  a_time) const
{

    const Interval velComps(0, SpaceDim-1);
    const Real old_time = m_time - m_dt;

    CH_assert(old_time - TIME_EPS < a_time);
    CH_assert(a_time < m_time + TIME_EPS);

    // Set entire level data to a bogus value, if requested
    if (s_set_bogus_values) {
        setValLevel(a_vel, s_bogus_value);
    }

    if (abs(a_time - m_time) < TIME_EPS) {
        // Copy from new values
        m_vel_new_ptr->copyTo(velComps, a_vel, velComps);
    } else if (abs(a_time - old_time) < TIME_EPS) {
        // Copy from old values
        m_vel_old_ptr->copyTo(velComps, a_vel, velComps);
    } else {
        // Interpolate in time
        timeInterp(a_vel, a_time, *m_vel_old_ptr, old_time,
                   *m_vel_new_ptr, m_time, velComps);
    }
}


// -----------------------------------------------------------------------------
// Returns lambda at time a_time.
// This does a linear interpolation in time between old and new time lambdas,
// but does not fill ghosts or interpolate in space from a coarser level.
// -----------------------------------------------------------------------------
void AMRNavierStokes::lambda (LevelData<FArrayBox>& a_lambda,
                              Real                  a_time) const
{

    const Interval lambdaComps(0,0);
    const Real old_time = m_time - m_dt;

    CH_assert(old_time - TIME_EPS < a_time);
    CH_assert(a_time < m_time + TIME_EPS);

    // Set entire level data to a bogus value, if requested
    if (s_set_bogus_values) {
        setValLevel(a_lambda, s_bogus_value);
    }

    if (abs(a_time - m_time) < TIME_EPS) {
        // Copy from new values
        m_lambda_new_ptr->copyTo(lambdaComps, a_lambda, lambdaComps);
    } else if (abs(a_time - old_time) < TIME_EPS) {
        // Copy from old values
        m_lambda_old_ptr->copyTo(lambdaComps, a_lambda, lambdaComps);
    } else {
        // Interpolate in time
        timeInterp(a_lambda, a_time, *m_lambda_old_ptr, old_time,
                   *m_lambda_new_ptr, m_time, lambdaComps);
    }
}


// -----------------------------------------------------------------------------
// Returns scalar at time a_time.
// This does a linear interpolation in time between old and new time scalars,
// but does not fill ghosts or interpolate in space from a coarser level.
// -----------------------------------------------------------------------------
void AMRNavierStokes::scalar (LevelData<FArrayBox>& a_scal,
                              const Real            a_time,
                              const int             a_comp) const
{

    CH_assert(0 <= a_comp);
    CH_assert(a_comp < s_num_scal_comps);

    const Interval& scalComps = m_scal_new[a_comp]->interval();
    const Real old_time = m_time - m_dt;

    CH_assert(old_time - TIME_EPS < a_time);
    CH_assert(a_time < m_time + TIME_EPS);

    // Set entire level data to a bogus value, if requested
    if (s_set_bogus_values) {
        setValLevel(a_scal, s_bogus_value);
    }

    if (abs(a_time - m_time) < TIME_EPS) {
        // Copy from new values
        m_scal_new[a_comp]->copyTo(scalComps, a_scal, scalComps);
    } else if (abs(a_time - old_time) < TIME_EPS) {
        // Copy from old values
        m_scal_old[a_comp]->copyTo(scalComps, a_scal, scalComps);
    } else {
        // Interpolate in time
        timeInterp(a_scal, a_time, *(m_scal_old[a_comp]), old_time,
                   *(m_scal_new[a_comp]), m_time, scalComps);
    }
}


// -----------------------------------------------------------------------------
// Fill grown velocity field using piecewise-constant interpolation.
// This does not require a_vel's grids to be the same as this level's grids.
// -----------------------------------------------------------------------------
void AMRNavierStokes::fillVelocity (LevelData<FArrayBox>& a_vel,
                                    const Real            a_time) const
{
    CH_TIME("AMRNavierStokes::fillVelocity 1");

    // Interpolate the velocity in time
    this->velocity(a_vel, a_time);

    // Set all ghosts on velocity
    VelBCHolder velBC(m_physBCPtr->tracingVelFuncBC()); // Set to tracingBCs. User can override if needed.
    this->setGhostsVelocity(a_vel, velBC, a_time);

    // Summary:
    // 1. Valid data of a_vel should now be completely filled with
    //    interpolated data at a_time, with no bogus values left.
    // 2. Ghost values should be filled at the CF interfaces and cells
    //    adjacent to the physical boundaries.
    // 3. The only bogus values that remain are at the edges and vertices
    //    of each grid and where exchanges are needed. The exchange is
    //    not performed here to save expenses.
}


// -----------------------------------------------------------------------------
// Resize and fill temporary lambda for advection tracing computation.
// This fills both interior and boundary cells.
// This function WILL handle allocation
// -----------------------------------------------------------------------------
void AMRNavierStokes::fillLambda (LevelData<FArrayBox>& a_lambda,
                                  Real                  a_time) const
{
    CH_TIME("AMRNavierStokes::fillLambda");

    // Collect some needed data
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    const IntVect ghostVect(D_DECL(ADVECT_GROW, ADVECT_GROW, ADVECT_GROW));
    const int ncomp = 1;

    // Define the level data holder
    CH_assert(!a_lambda.isDefined());
    a_lambda.define(grids, ncomp, ghostVect);

    // Interpolate lambda in time
    this->lambda(a_lambda, a_time);

    // Set all ghosts on lambda
    BCMethodHolder lambdaBC = m_physBCPtr->lambdaFuncBC();
    setGhostsLambda(a_lambda, lambdaBC, a_time);

    // Summary:
    // 1. Valid data of a_lambda should now be completely filled with
    //    interpolated data at a_time, with no bogus values left.
    // 2. Ghost values should be filled at the CF interfaces and cells
    //    adjacent to the physical boundaries.
    // 3. The only bogus values that remain are at the edges and vertices
    //    of each grid and where exchanges are needed. The exchange is
    //    not performed here to save expenses.
}


// -----------------------------------------------------------------------------
// Resize and fill temporary scalar holder for advection tracing computation.
// This fills both interior and boundary cells.
// This function WILL handle allocation
// -----------------------------------------------------------------------------
void AMRNavierStokes::fillScalars (LevelData<FArrayBox>& a_scal,
                                   Real                  a_time,
                                   const int             a_comp,
                                   const bool            a_addBackground) const
{
    // Collect some needed data
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    const IntVect ghostVect(D_DECL(ADVECT_GROW, ADVECT_GROW, ADVECT_GROW));

    // Define the level data holder
    CH_assert(!a_scal.isDefined());
    a_scal.define(grids, m_scal_new[a_comp]->nComp(), ghostVect);

    // Interpolate the scalar in time
    this->scalar(a_scal, a_time, a_comp);

    // Set all ghosts on the scalar
    BCMethodHolder scalBC = m_physBCPtr->scalarTraceFuncBC(a_comp);
    this->setGhostsScalar(a_scal, scalBC, a_time, a_comp);

    if (a_addBackground) {
        m_physBCPtr->addBackgroundScalar(a_scal, 0, a_time, *m_levGeoPtr);
    }

    // Summary:
    // 1. Valid data of a_scal should now be completely filled with
    //    interpolated data at a_time, with no bogus values left.
    // 2. Ghost values should be filled at the CF interfaces and cells
    //    adjacent to the physical boundaries.
    // 3. The only bogus values that remain are at the edges and vertices
    //    of each grid and where exchanges are needed. The exchange is
    //    not performed here to save expenses.
}


// -----------------------------------------------------------------------------
// Computes the viscous source term, nu.L[u].
// All inputs and outputs are in the Cartesian basis.
// This assumes the ghosts of a_cartVel have already been set.
// -----------------------------------------------------------------------------
void AMRNavierStokes::fillViscousSource (LevelData<FArrayBox>&       a_viscSource,
                                         const LevelData<FArrayBox>& a_cartVel,
                                         const Real                  a_time)
{
    const IntVect ghostVect = a_viscSource.ghostVect();
    const DisjointBoxLayout& grids = a_cartVel.getBoxes();
    DataIterator dit = a_cartVel.dataIterator();
    const bool isViscous = (s_nu > 0.0);

    // Sanity checks
    CH_assert(a_viscSource.nComp() == SpaceDim);
    CH_assert(a_cartVel   .nComp() == SpaceDim);
    CH_assert(a_viscSource.getBoxes() == grids);
    CH_assert(a_cartVel   .getBoxes() == grids);

    if (s_set_bogus_values) {
        setValLevel(a_viscSource, s_bogus_value);
    }

    // Compute the viscous source term
    if (isViscous) {
        // If crse level exists, fill coarse velocity BC and convert to a
        // Cartesian basis as well.
        LevelData<FArrayBox>* crseVelPtr = NULL;
        if (m_level > 0) {
            AMRNavierStokes* crseLevelPtr = crseNSPtr();
            const DisjointBoxLayout& crseGrids = crseLevelPtr->newVel().getBoxes();

            crseVelPtr = new LevelData<FArrayBox>(crseGrids, SpaceDim);
            crseLevelPtr->fillVelocity(*crseVelPtr, a_time);
            crseLevelPtr->m_levGeoPtr->sendToCartesianBasis(*crseVelPtr, false);
        }

        // Compute viscous source term
        this->computeLapVel(a_viscSource, a_cartVel, crseVelPtr, NULL, a_time);

        // Multiply by the viscous coeff
        DataIterator localDit = a_viscSource.dataIterator();
        for (localDit.reset(); localDit.ok(); ++localDit) {
            a_viscSource[localDit].mult(s_nu);
        }

        // Clean up temp storage
        if (crseVelPtr != NULL) {
            delete crseVelPtr;
            crseVelPtr = NULL;
        }

    } else {
        // The fluid is not viscous - just set viscous source to zero.
        setValLevel(a_viscSource, 0.0);
    }
}


// -----------------------------------------------------------------------------
// Returns the gravitational source term.
// This is -g*zhat in the Cartesian basis. a_gravSource needs SpaceDim comps.
// -----------------------------------------------------------------------------
void AMRNavierStokes::fillGravSource (LevelData<FArrayBox>& a_gravSource,
                                      const Real            a_time,
                                      const bool            a_addBackground) const
{
    CH_TIME("AMRNavierStokes::fillGravSource");

    CH_assert(a_gravSource.nComp() == SpaceDim);
    DataIterator dit = a_gravSource.dataIterator();

    // Set output to bogus values, if requested
    if (s_set_bogus_values) {
        setValLevel(a_gravSource, s_bogus_value);
    }

    if (s_gravityMethod != ProblemContext::GravityMethod::NONE) {
        // Fill a z comp of gravSource with density
        LevelData<FArrayBox> density;
        this->fillScalars(density, a_time, 0, a_addBackground);

        // Set all other comps to zero and scale z-comp to -b.
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& thisFAB = a_gravSource[dit];
            thisFAB.setVal(0.0, thisFAB.box(), 0, SpaceDim-1);
            thisFAB.copy(density[dit], 0, SpaceDim-1);
            thisFAB.negate(SpaceDim-1);
        }
    } else {
        // Just set a_gravSource to zero.
        for (dit.reset(); dit.ok(); ++dit) {
            a_gravSource[dit].setVal(0.0);
        }
    }
}


// -----------------------------------------------------------------------------
// Returns the tidal forcing source term.
// This is U0*omega*cos(omega*t) in the Cartesian basis.
// a_tidalSource needs SpaceDim comps.
// -----------------------------------------------------------------------------
void AMRNavierStokes::fillTidalSource (LevelData<FArrayBox>& a_tidalSource,
                                       const Real            a_oldTime,
                                       const Real            a_dt) const
{
    CH_TIME("AMRNavierStokes::fillTidalSource");

    CH_assert(s_tidalOmega * s_tidalU0 != 0.0);
    CH_assert(a_tidalSource.nComp() == SpaceDim);
    DataIterator dit = a_tidalSource.dataIterator();

    // Set output to bogus values, if requested
    if (s_set_bogus_values) {
        setValLevel(a_tidalSource, s_bogus_value);
    }

    if (s_tidalOmega * s_tidalU0 != 0.0) {
        // const Real tidalForce = s_tidalU0 * s_tidalOmega * cos(s_tidalOmega * a_time);

        const Real oldArg = s_tidalOmega * a_oldTime;
        const Real newArg = s_tidalOmega * (a_oldTime + a_dt);
        const Real tidalForce = s_tidalU0 * (sin(newArg) - sin(oldArg)) / a_dt;

        for (dit.reset(); dit.ok(); ++dit) {
            a_tidalSource[dit].setVal(tidalForce, 0);
            a_tidalSource[dit].setVal(0.0, a_tidalSource[dit].box(), 1, SpaceDim-1);
        }

    } else {
        // Just set a_tidalSource to zero.
        for (dit.reset(); dit.ok(); ++dit) {
            a_tidalSource[dit].setVal(0.0);
        }
    }
}


// -----------------------------------------------------------------------------
// Returns the internal wave speed
// -----------------------------------------------------------------------------
void AMRNavierStokes::fillInternalWaveSpeed (LevelData<FArrayBox>& a_c0i) const
{
    CH_assert(a_c0i.nComp() == SpaceDim);
    CH_assert(a_c0i.getBoxes().compatible(m_levGeoPtr->getBoxes()));

    if (m_level == 0) {
        // Just perform a copy.
        m_c0iPtr->copyTo(a_c0i);

    } else if (m_level == 1) {
        // Direct interpolation from level 0.
        const LevelData<FArrayBox>* crseDataPtr = crseNSPtr()->m_c0iPtr;

        // Set up fine level structures
        const DisjointBoxLayout& fineGrids = m_levGeoPtr->getBoxes();
        const ProblemDomain& fineDomain = m_levGeoPtr->getDomain();
        const IntVect& refRatio = m_levGeoPtr->getCrseRefRatio();

        // Interpolate
        MappedFineInterp interpObj(fineGrids,
                                   crseDataPtr->nComp(),
                                   refRatio,
                                   fineDomain,
                                   m_levGeoPtr,
                                   true);  // considerCellVols
        interpObj.interpToFine(a_c0i, *crseDataPtr);

    } else {
        const AMRNavierStokes* crseLevelPtr = coarsestNSPtr();
        const AMRNavierStokes* fineLevelPtr = crseLevelPtr->fineNSPtr();
        CH_assert(fineLevelPtr->m_level < m_level);

        LevelData<FArrayBox>* crseDataPtr = crseLevelPtr->m_c0iPtr;
        LevelData<FArrayBox>* fineDataPtr = NULL;

        const int numComps = crseDataPtr->nComp();

        // Interpolate up to level 1
        {
            CH_assert(crseLevelPtr->m_level == 0);
            CH_assert(fineLevelPtr->m_level == 1);

            // Set up fine level structures
            const LevelGeometry* fineLevGeoPtr = fineLevelPtr->m_levGeoPtr;
            const DisjointBoxLayout& fineGrids = fineLevGeoPtr->getBoxes();
            const ProblemDomain& fineDomain = fineLevGeoPtr->getDomain();
            const IntVect& refRatio = fineLevGeoPtr->getCrseRefRatio();
            fineDataPtr = new LevelData<FArrayBox>(fineGrids, numComps);

            // Interpolate
            MappedFineInterp interpObj(fineGrids,
                                       numComps,
                                       refRatio,
                                       fineDomain,
                                       fineLevGeoPtr,
                                       true);  // considerCellVols
            interpObj.interpToFine(*fineDataPtr, *crseDataPtr);

            // Move up one level
            crseDataPtr = fineDataPtr;
            fineDataPtr = NULL;

            crseLevelPtr = fineLevelPtr;
            fineLevelPtr = fineLevelPtr->fineNSPtr();
        }

        // Interpolate data up to m_level-1
        if (fineLevelPtr->m_level < m_level) {
            fineDataPtr = new LevelData<FArrayBox>;

            while (fineLevelPtr->m_level < m_level) {
                CH_assert(crseLevelPtr->m_level == fineLevelPtr->m_level - 1);
                CH_assert(fineLevelPtr->m_level < m_level);

                // Set up fine level structures
                const LevelGeometry* fineLevGeoPtr = fineLevelPtr->m_levGeoPtr;
                const DisjointBoxLayout& fineGrids = fineLevGeoPtr->getBoxes();
                const ProblemDomain& fineDomain = fineLevGeoPtr->getDomain();
                const IntVect& refRatio = fineLevGeoPtr->getCrseRefRatio();
                fineDataPtr->define(fineGrids, numComps);

                // Interpolate
                MappedFineInterp interpObj(fineGrids,
                                           numComps,
                                           refRatio,
                                           fineDomain,
                                           fineLevGeoPtr,
                                           true);  // considerCellVols
                interpObj.interpToFine(*fineDataPtr, *crseDataPtr);

                // Move up one level
                std::swap(crseDataPtr, fineDataPtr);
                crseLevelPtr = fineLevelPtr;
                fineLevelPtr = fineLevelPtr->fineNSPtr();
            }

            delete fineDataPtr;
            fineDataPtr = NULL;
        }

        // Interpolate up to this level.
        {
            CH_assert(crseLevelPtr->m_level == fineLevelPtr->m_level - 1);
            CH_assert(fineLevelPtr->m_level == m_level);

            // Set up fine level structures
            const LevelGeometry* fineLevGeoPtr = fineLevelPtr->m_levGeoPtr;
            const DisjointBoxLayout& fineGrids = fineLevGeoPtr->getBoxes();
            const ProblemDomain& fineDomain = fineLevGeoPtr->getDomain();
            const IntVect& refRatio = fineLevGeoPtr->getCrseRefRatio();
            fineDataPtr = &a_c0i;

            // Interpolate
            MappedFineInterp interpObj(fineGrids,
                                       numComps,
                                       refRatio,
                                       fineDomain,
                                       fineLevGeoPtr,
                                       true);  // considerCellVols
            interpObj.interpToFine(*fineDataPtr, *crseDataPtr);
        }

        // Final cleanup.
        delete crseDataPtr;
        crseDataPtr = NULL;
    }
}


// -----------------------------------------------------------------------------
// Returns the vertical structure function.
// -----------------------------------------------------------------------------
void AMRNavierStokes::fillVerticalStructure (LevelData<FArrayBox>& a_phi0) const
{
    CH_assert(a_phi0.nComp() == 1);
    CH_assert(a_phi0.getBoxes().compatible(m_levGeoPtr->getBoxes()));

    if (m_level == 0) {
        // Just perform a copy.
        m_phi0Ptr->copyTo(a_phi0);

    } else if (m_level == 1) {
        // Direct interpolation from level 0.
        const LevelData<FArrayBox>* crseDataPtr = crseNSPtr()->m_phi0Ptr;

        // Set up fine level structures
        const DisjointBoxLayout& fineGrids = m_levGeoPtr->getBoxes();
        const ProblemDomain& fineDomain = m_levGeoPtr->getDomain();
        const IntVect& refRatio = m_levGeoPtr->getCrseRefRatio();

        // Interpolate
        MappedFineInterp interpObj(fineGrids,
                                   crseDataPtr->nComp(),
                                   refRatio,
                                   fineDomain,
                                   m_levGeoPtr,
                                   true);  // considerCellVols
        interpObj.interpToFine(a_phi0, *crseDataPtr);

    } else {
        const AMRNavierStokes* crseLevelPtr = coarsestNSPtr();
        const AMRNavierStokes* fineLevelPtr = crseLevelPtr->fineNSPtr();
        CH_assert(fineLevelPtr->m_level < m_level);

        LevelData<FArrayBox>* crseDataPtr = crseLevelPtr->m_phi0Ptr;
        LevelData<FArrayBox>* fineDataPtr = NULL;

        const int numComps = crseDataPtr->nComp();

        // Interpolate up to level 1
        {
            CH_assert(crseLevelPtr->m_level == 0);
            CH_assert(fineLevelPtr->m_level == 1);

            // Set up fine level structures
            const LevelGeometry* fineLevGeoPtr = fineLevelPtr->m_levGeoPtr;
            const DisjointBoxLayout& fineGrids = fineLevGeoPtr->getBoxes();
            const ProblemDomain& fineDomain = fineLevGeoPtr->getDomain();
            const IntVect& refRatio = fineLevGeoPtr->getCrseRefRatio();
            fineDataPtr = new LevelData<FArrayBox>(fineGrids, numComps);

            // Interpolate
            MappedFineInterp interpObj(fineGrids,
                                       numComps,
                                       refRatio,
                                       fineDomain,
                                       fineLevGeoPtr,
                                       true);  // considerCellVols
            interpObj.interpToFine(*fineDataPtr, *crseDataPtr);

            // Move up one level
            crseDataPtr = fineDataPtr;
            fineDataPtr = NULL;

            crseLevelPtr = fineLevelPtr;
            fineLevelPtr = fineLevelPtr->fineNSPtr();
        }

        // Interpolate data up to m_level-1
        if (fineLevelPtr->m_level < m_level) {
            fineDataPtr = new LevelData<FArrayBox>;

            while (fineLevelPtr->m_level < m_level) {
                CH_assert(crseLevelPtr->m_level == fineLevelPtr->m_level - 1);
                CH_assert(fineLevelPtr->m_level < m_level);

                // Set up fine level structures
                const LevelGeometry* fineLevGeoPtr = fineLevelPtr->m_levGeoPtr;
                const DisjointBoxLayout& fineGrids = fineLevGeoPtr->getBoxes();
                const ProblemDomain& fineDomain = fineLevGeoPtr->getDomain();
                const IntVect& refRatio = fineLevGeoPtr->getCrseRefRatio();
                fineDataPtr->define(fineGrids, numComps);

                // Interpolate
                MappedFineInterp interpObj(fineGrids,
                                           numComps,
                                           refRatio,
                                           fineDomain,
                                           fineLevGeoPtr,
                                           true);  // considerCellVols
                interpObj.interpToFine(*fineDataPtr, *crseDataPtr);

                // Move up one level
                std::swap(crseDataPtr, fineDataPtr);
                crseLevelPtr = fineLevelPtr;
                fineLevelPtr = fineLevelPtr->fineNSPtr();
            }

            delete fineDataPtr;
            fineDataPtr = NULL;
        }

        // Interpolate up to this level.
        {
            CH_assert(crseLevelPtr->m_level == fineLevelPtr->m_level - 1);
            CH_assert(fineLevelPtr->m_level == m_level);

            // Set up fine level structures
            const LevelGeometry* fineLevGeoPtr = fineLevelPtr->m_levGeoPtr;
            const DisjointBoxLayout& fineGrids = fineLevGeoPtr->getBoxes();
            const ProblemDomain& fineDomain = fineLevGeoPtr->getDomain();
            const IntVect& refRatio = fineLevGeoPtr->getCrseRefRatio();
            fineDataPtr = &a_phi0;

            // Interpolate
            MappedFineInterp interpObj(fineGrids,
                                       numComps,
                                       refRatio,
                                       fineDomain,
                                       fineLevGeoPtr,
                                       true);  // considerCellVols
            interpObj.interpToFine(*fineDataPtr, *crseDataPtr);
        }

        // Final cleanup.
        delete crseDataPtr;
        crseDataPtr = NULL;
    }
}
