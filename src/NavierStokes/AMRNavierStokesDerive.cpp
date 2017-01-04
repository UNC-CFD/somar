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
#include "Curl.H"
#include "ExtrapolationUtils.H"
#include "computeMappedSum.H"
#include "MappedAMRPoissonOpFactory.H"
#include "StratUtils.H"
#include "EllipticBCUtils.H"
#include <iomanip>


// -----------------------------------------------------------------------------
// Compute Grad[macPressure]. This is really just a wrapper around
// m_macProjector.computeGrad to help keep the code neat. This also avoids
// using a possibly undefined projector in the compressible case.
// -----------------------------------------------------------------------------
void AMRNavierStokes::gradMACPressure (LevelData<FluxBox>& a_gradPressure,
                                       const Real          a_scale)
{
    // Note that we do not query isPressureAvail(). If the flow is
    // incompressible and we are trying to use Grad[pressure] irresponsibly,
    // I want an error to be thrown!
    if (s_isIncompressible) {
        m_macProjector.computeLevelGradPressure(a_gradPressure, a_scale);
    } else {
        DataIterator dit = a_gradPressure.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            a_gradPressure[dit].setVal(0.0);
        }
    }
}


// -----------------------------------------------------------------------------
// Compute Grad[ccPressure]. This is really just a wrapper around
// m_ccProjector.computeGrad to help keep the code neat. This also avoids
// using a possibly undefined projector in the compressible case.
// -----------------------------------------------------------------------------
void AMRNavierStokes::gradCCPressure (LevelData<FArrayBox>& a_gradPressure,
                                      const Real            a_scale)
{
    // Note that we do not query isPressureAvail(). If the flow is
    // incompressible and we are trying to use Grad[pressure] irresponsibly,
    // I want an error to be thrown!
    if (s_isIncompressible && m_ccPressureState == CCPressureState::VALID) {
        m_ccProjector.computeLevelGradPressure(a_gradPressure, a_scale);
    } else {
        DataIterator dit = a_gradPressure.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            a_gradPressure[dit].setVal(0.0);
        }
    }
}


// -----------------------------------------------------------------------------
// BCs must already be set.
// -----------------------------------------------------------------------------
void AMRNavierStokes::computeLapVel (LevelData<FArrayBox>&       a_lapVel,
                                     const LevelData<FArrayBox>& a_vel,
                                     const LevelData<FArrayBox>* a_crseVelPtr,
                                     const CFRegion*             a_CFRegionPtr,
                                     const Real                  a_time)
{
    CH_assert(a_lapVel.nComp() == SpaceDim);
    CH_assert(a_vel   .nComp() == SpaceDim);
    CH_assert(a_lapVel.getBoxes() == a_vel.getBoxes());

    const DisjointBoxLayout grids = a_vel.getBoxes();
    DataIterator dit = grids.dataIterator();
    LevelData<FArrayBox> tmpVel(a_vel.getBoxes(), a_vel.nComp(), a_vel.ghostVect());
    for (dit.reset(); dit.ok(); ++dit) {
        tmpVel[dit].copy(a_vel[dit]);
    }

    // Apply the op
    bool isHomogeneous = false;
    VelBCHolder velBC(m_physBCPtr->viscousSourceFuncBC());

    m_velocityAMRPoissonOp.setTime(a_time);
    m_velocityAMRPoissonOp.applyOp(a_lapVel,
                                   tmpVel,
                                   a_crseVelPtr,
                                   isHomogeneous,
                                   &velBC);

    // Extrap all ghosts
    extrapAllGhosts(a_lapVel,0);

    // And exchange
    if (a_lapVel.ghostVect() == m_copierCache.getTracingGhosts()) {
        a_lapVel.exchange(m_copierCache.getTracingExCopier(a_lapVel.getBoxes()));
    } else if(a_lapVel.ghostVect() == IntVect::Unit) {
        a_lapVel.exchange(m_copierCache.getOneGhostExCopier(a_lapVel.getBoxes()));
    } else {
        a_lapVel.exchange();
    }
}


// -----------------------------------------------------------------------------
// Compute D[nu G[u]].
// -----------------------------------------------------------------------------
void AMRNavierStokes::computeViscousSrc (LevelData<FArrayBox>&       a_viscSrc,
                                         const LevelData<FArrayBox>& a_cartVel,
                                         const Real                  a_time) const
{
    ((AMRNavierStokes*)this)->fillViscousSource(a_viscSrc, a_cartVel, a_time);
    return;
}


// -----------------------------------------------------------------------------
// Compute D[kappa G[b]].
// -----------------------------------------------------------------------------
void AMRNavierStokes::computeDiffusiveSrc (LevelData<FArrayBox>&       a_diffSrc,
                                           const LevelData<FArrayBox>& a_scal,
                                           const LevelData<FArrayBox>* a_crseScalPtr,
                                           const int                   a_comp,
                                           const Real                  a_time) const
{
    // Sanity checks
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < PhysBCUtil::getNumScalars());
    CH_assert(m_levGeoPtr->getBoxes() == a_diffSrc.getBoxes());
    CH_assert(m_levGeoPtr->getBoxes() == a_scal.getBoxes());

    // Gather data structures
    const ProblemDomain& domain = m_levGeoPtr->getDomain();
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    DataIterator dit = a_diffSrc.dataIterator();

    // Send scalar to a new container.
    // This way, we can change the ghost values.
    LevelData<FArrayBox> localScal(grids, 1, a_scal.ghostVect());
    {
        const Interval srcIvl(a_comp, a_comp);
        const Interval& destIvl = localScal.interval();
        Copier noGhostCopier(grids, grids, domain, IntVect::Zero, false);
        a_scal.copyTo(srcIvl, localScal, destIvl, noGhostCopier);
    }

    MappedAMRPoissonOp& lapOp = (MappedAMRPoissonOp&)m_scalarsAMRPoissonOp;

    // Apply the op
    BCMethodHolder scalPhysBCs = m_physBCPtr->diffusiveSourceFuncBC();
    lapOp.setBC(scalPhysBCs);

    bool isHomogeneous = false;
    if (a_crseScalPtr != NULL) {
        lapOp.AMROperatorNF(a_diffSrc,
                            localScal,
                            *a_crseScalPtr,
                            isHomogeneous);
    } else {
        lapOp.applyOpI(a_diffSrc,
                       localScal,
                       isHomogeneous);
    }

    for (dit.reset(); dit.ok(); ++dit) {
        a_diffSrc[dit] *= s_scal_coeffs[a_comp];
    }

    // Extrapolate all ghosts.
    extrapAllGhosts(a_diffSrc,0);

    // Do exchanges.
    const IntVect& ghostVect = a_diffSrc.ghostVect();
    if (ghostVect == m_copierCache.getTracingGhosts()) {
        a_diffSrc.exchange(m_copierCache.getTracingExCopier(a_diffSrc.getBoxes()));
    } else if (ghostVect == IntVect::Unit) {
        a_diffSrc.exchange(m_copierCache.getOneGhostExCopier(a_diffSrc.getBoxes()));
    } else if (ghostVect != IntVect::Zero) {
        a_diffSrc.exchange();
    }

    return;
}


// -----------------------------------------------------------------------------
// Compute the gradient Richardson number.
// -----------------------------------------------------------------------------
void AMRNavierStokes::computeRiNumber (LevelData<FArrayBox>& a_Ri,
                                       const int             a_RiComp,
                                       const Real            a_time) const
{
    CH_TIME("AMRNavierStokes::computeRiNumber");

    // Collect references
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    const Copier& cp = m_copierCache.getOneGhostExCopier(grids);

    // If useExtrapBCs is true, cartVel and buoyancy will simply be extrapolated
    // to the ghost cells. This is an attempt to avoid excessive tagging at
    // the physical boundaries.
    const bool useExtrapBCs = true;
    const int extrapOrder = 1;

    // Set up extrapolation BCs if needed.
    BCMethodHolder extrapBCHolder;
    if (useExtrapBCs) {
        RefCountedPtr<BCGhostClass> extrapBCPtr(
            new EllipticExtrapBCGhostClass(extrapOrder,
                                           IntVect::Unit,
                                           IntVect::Unit)
        );
        extrapBCHolder.addBCMethod(extrapBCPtr);
    }

    // Sanity checks
    CH_assert(m_problem_domain == grids.physDomain());
    CH_assert(newVel().getBoxes() == grids);
    CH_assert(oldVel().getBoxes() == grids);
    CH_assert(newScal(0).getBoxes() == grids);
    CH_assert(oldScal(0).getBoxes() == grids);

    // Compute Cartesian-based velocity with one ghost.
    LevelData<FArrayBox> cartVel(grids, SpaceDim, IntVect::Unit);
    {
        // Interpolate the velocity in time
        this->velocity(cartVel, a_time);

        // Fill physical and CF BC ghosts.
        if (useExtrapBCs) {
            Tuple<BCMethodHolder, SpaceDim> extrapTuple;
            D_TERM(extrapTuple[0] = extrapBCHolder;,
                   extrapTuple[1] = extrapBCHolder;,
                   extrapTuple[2] = extrapBCHolder;)
            VelBCHolder velBCHolder(extrapTuple);
            this->setGhostsVelocity(cartVel, velBCHolder, a_time);
        } else {
            if (s_nu > 0.0) {
                // Use viscous BCs
                VelBCHolder velBCHolder(m_physBCPtr->viscousVelFuncBC());
                this->setGhostsVelocity(cartVel, velBCHolder, a_time);
            } else {
                // Use inviscid BCs for tracing
                VelBCHolder velBCHolder(m_physBCPtr->tracingVelFuncBC());
                this->setGhostsVelocity(cartVel, velBCHolder, a_time);
            }
        }

        // Exchange ghosts
        cartVel.exchange(cp);

        // Convert to Cartesian basis.
        m_levGeoPtr->sendToCartesianBasis(cartVel, true);
    }

    // Compute total buoyancy
    LevelData<FArrayBox> buoyancy(grids, 1, IntVect::Unit);
    {
        const int comp = 0;
        const bool onlyAtValid = true;

        // Interpolate the scalar in time
        this->scalar(buoyancy, a_time, comp);

        // Set all ghosts on the scalar
        if (useExtrapBCs) {
            this->setGhostsScalar(buoyancy, extrapBCHolder, a_time, comp, onlyAtValid);
        } else {
            if (s_scal_coeffs[0] > 0.0) {
                BCMethodHolder scalBCHolder = m_physBCPtr->diffusiveSourceFuncBC();
                this->setGhostsScalar(buoyancy, scalBCHolder, a_time, comp, onlyAtValid);
            } else {
                BCMethodHolder scalBCHolder = m_physBCPtr->scalarTraceFuncBC(comp);
                this->setGhostsScalar(buoyancy, scalBCHolder, a_time, comp, onlyAtValid);
            }
        }

        // Compute total buoyancy by adding the background to the deviation.
        m_physBCPtr->addBackgroundScalar(buoyancy, comp, a_time, *m_levGeoPtr);

        // Exchange ghosts
        buoyancy.exchange(cp);
    }

    // Compute the gradient Richardson number via the StratUtils function.
    const bool useHorizSsq = true;
    computeGradRiNumber(a_Ri, a_RiComp, cartVel,
                        buoyancy, *m_levGeoPtr,
                        a_time, useHorizSsq);
}


// -----------------------------------------------------------------------------
// computeVorticity
// NOTE: This returns two-form components, not a vector!
// -----------------------------------------------------------------------------
void AMRNavierStokes::computeVorticity (LevelData<FArrayBox>& a_vorticity) const
{
    CH_TIME("AMRNavierStokes::computeVorticity");

    // This breaks the "const"-ness of the function,
    // but is necessary to ensure that boundary
    // conditions are properly set.
    const DisjointBoxLayout& grids = a_vorticity.getBoxes();
    DataIterator dit = a_vorticity.dataIterator();
    const RealVect& dx = m_levGeoPtr->getDx();
    const bool isViscous = (s_nu > 0.0);

    // Copy velocity to new holder so we can set BCs.
    LevelData<FArrayBox> vel(grids, SpaceDim, IntVect::Unit);
    newVel().copyTo(vel);

    // CF BCs
    if (m_level > 0) {
        // do quadratic C/F BC's
        // for now, assume that BC's are with new velocity
        // (may need to be changed for fine-level regridding)
        LevelData<FArrayBox>& crseVel = crseNSPtr()->newVel();
        const DisjointBoxLayout& crseGrids = crseVel.getBoxes();
        const IntVect& nRefCrse = crseNSPtr()->refRatio();

        MappedQuadCFInterp interpolator(grids, &crseGrids, dx, nRefCrse,
                                        SpaceDim, m_problem_domain);

        interpolator.coarseFineInterp(vel, crseVel);
    }

    // Set physical BCs
    Tuple<BCMethodHolder, SpaceDim> vortVelBC = m_physBCPtr->vortFuncBC(isViscous);

    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& velFAB = vel[dit];
        const FluxBox& JgupFB = m_levGeoPtr->getFCJgup()[dit];
        const Box& valid = grids[dit];

        // Apply physical BC's
        for (int velcomp = 0; velcomp < SpaceDim; ++velcomp) {
            vortVelBC[velcomp].setGhosts(velFAB,             // stateFAB
                                         NULL,               // extrapFABPtr
                                         valid,              // valid
                                         m_problem_domain,   // domain
                                         dx,                 // dx
                                         dit(),              // DataIndex
                                         &JgupFB,            // JgupFBPtr
                                         false);             // isHomogeneous
        }
    }

    // Do exchange.
    vel.exchange(m_copierCache.getOneGhostExCopier(vel.getBoxes()));


    // Compute vorticity
    for (dit.reset(); dit.ok(); ++dit) {
        Curl::simpleCurlCC(a_vorticity[dit],
                           vel[dit],
                           dit(),
                           *m_levGeoPtr);

    } // end loop over grids
}


// -----------------------------------------------------------------------------
// Compute the stream function on this level.
// -----------------------------------------------------------------------------
void AMRNavierStokes::computeStreamFunction (LevelData<FArrayBox>& a_stream) const
{
    CH_TIME("AMRNavierStokes::computeStreamFunction");

    // Print warnings
    if (SpaceDim > 2) {
        MayDay::Warning("computeStreamFunction needs testing in 3D");
    }
    if (m_level > 0) {
        MayDay::Warning("computeStreamFunction needs testing when m_level > 0");
    }

    // Gather needed info
    const int streamComps = D_TERM(0,+1,+2);
    const Vector<const LevelGeometry*> amrLevGeos = ((const LevelGeometry*)m_levGeoPtr)->getAMRLevGeos();
    const int numAMRLevels = m_level + 1;
    CH_assert(amrLevGeos.size() >= numAMRLevels);

    // Set up the Poisson op factory.
    const Real alpha = 0.0;
    const Real beta = 1.0;
    BCMethodHolder bcHolder(m_physBCPtr->streamSolverBC(0));
    MappedAMRPoissonOpFactory localPoissonOpFactory;
    localPoissonOpFactory.define(m_levGeoPtr,
                                 alpha,
                                 beta,
                                 bcHolder,
                                 s_AMRMG_maxDepth,
                                 s_AMRMG_num_smooth_precond,
                                 s_AMRMG_precondMode,
                                 s_AMRMG_relaxMode);

    // Set up bottom solver for AMRMultigrid.
    BiCGStabSolver<LevelData<FArrayBox> > bottomSolver;
    bottomSolver.m_eps = s_bottom_eps;
    bottomSolver.m_reps = s_bottom_reps;
    bottomSolver.m_imax = s_bottom_imax;
    bottomSolver.m_numRestarts = s_bottom_numRestarts;
    bottomSolver.m_hang = s_bottom_hang;
    bottomSolver.m_small = s_bottom_small;
    bottomSolver.m_verbosity = s_bottom_verbosity;
    bottomSolver.m_normType = s_bottom_normType;

    // Set up the AMRMultigrid solver.
    MappedAMRMultiGrid<LevelData<FArrayBox> > streamSolver;
    streamSolver.define(amrLevGeos[0]->getDomain(),
                        localPoissonOpFactory,
                        &bottomSolver,
                        m_level + 1);  // num amr levels
    streamSolver.m_verbosity = s_AMRMG_verbosity;
    streamSolver.m_imin = s_AMRMG_imin;
    streamSolver.setSolverParameters(s_AMRMG_num_smooth_down,
                                     s_AMRMG_num_smooth_up,
                                     s_AMRMG_num_smooth_bottom,
                                     s_AMRMG_numMG,
                                     s_AMRMG_imax,
                                     s_AMRMG_eps,
                                     s_AMRMG_hang,
                                     s_AMRMG_normThresh);

    // Set up data holders
    Vector<LevelData<FArrayBox>*> amrRHS(numAMRLevels, NULL);
    Vector<LevelData<FArrayBox>*> amrSol(numAMRLevels, NULL);

    const int lev = m_level;
    const AMRNavierStokes* levPtr = this;

    amrRHS[lev] = new LevelData<FArrayBox>(amrLevGeos[lev]->getBoxes(), streamComps, IntVect::Unit);
    amrSol[lev] = &a_stream;

    levPtr->computeVorticity(*amrRHS[lev]);
    levPtr->m_levGeoPtr->divByJ(*amrRHS[lev]);

    for (DataIterator dit = amrRHS[lev]->dataIterator(); dit.ok(); ++dit) {
        FArrayBox& rhsFAB = (*amrRHS[lev])[dit];
        rhsFAB.negate();
    }

    // Loop over streamComps and solve.
    for (int comp = 0; comp < streamComps; ++comp) {
        const Interval rhsComps(comp,comp);

        // Alias the current comps.
        Vector<LevelData<FArrayBox>*> amrRHSComp(numAMRLevels, NULL);
        Vector<LevelData<FArrayBox>*> amrSolComp(numAMRLevels, NULL);

        amrRHSComp[lev] = new LevelData<FArrayBox>;
        amrSolComp[lev] = new LevelData<FArrayBox>;

        aliasLevelData(*amrRHSComp[lev], amrRHS[lev], rhsComps);
        aliasLevelData(*amrSolComp[lev], amrSol[lev], rhsComps);

        // Create notice of the solve.
        if (s_verbosity >= 2) {
            pout() << "\nSolving for comp " << comp << " of the streamfunction"
                   << setiosflags(ios::scientific) << setprecision(8)
                   << endl;
        }

        // Solve!
        streamSolver.solve(amrSolComp,
                           amrRHSComp,
                           m_level, // finest level
                           m_level, // base level
                           true,    // Initialize solution to zero
                           true);   // force homogeneous

        // Free memory
        delete amrRHSComp[lev];
        delete amrSolComp[lev];
    } // end loop over streamfunction components

    // Free memory
    for (int idx = 0; idx < amrRHS.size(); ++idx) {
        if (amrRHS[idx] != NULL) {
            delete amrRHS[idx];
            amrRHS[idx] = NULL;
        }
        amrSol[idx] = NULL;
    }
}



// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
Real AMRNavierStokes::totalEnergy () const
{
    CH_TIME("AMRNavierStokes::totalEnergy");

    // Gather composite data
    Vector<LevelData<FArrayBox>*> amrVel(0), amrB(0);
    const AMRNavierStokes* levelNSPtr = this;
    while(levelNSPtr != NULL) {
        amrVel.push_back(levelNSPtr->m_vel_new_ptr);
        if (s_num_scal_comps > 0) amrB.push_back(levelNSPtr->m_scal_new[0]);

        levelNSPtr = levelNSPtr->fineNSPtr();
    }

    Real globalEnergy;
    Vector<LevelData<FArrayBox>*> amrEnergy(0);
    const int numLevels = amrVel.size();
    const LevelGeometry* levGeoPtr = m_levGeoPtr;

    for (int ilev = m_level; ilev < numLevels; ++ilev) {
        // Check if this level is in use
        if (amrVel[ilev] == NULL) break;
        if (s_num_scal_comps > 0) {
            if (amrB[ilev] == NULL) break;
        }

        // Gather geometric data
        CH_assert(levGeoPtr != NULL);
        const RealVect levelDx = levGeoPtr->getDx();
        const DisjointBoxLayout& levelGrids = levGeoPtr->getBoxes();
        CH_assert(amrVel[ilev]->getBoxes().compatible(levelGrids));

        // Allocate and define holder for energy computation
        amrEnergy.push_back(new LevelData<FArrayBox>);
        amrEnergy[ilev]->define(levelGrids, 1);

        // Loop over grids and compute the energy at each point
        DataIterator dit = amrEnergy[ilev]->dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            // Create references on this grid for convenience
            FArrayBox& eFAB = (*amrEnergy[ilev])[dit];
            const FArrayBox& velFAB = (*amrVel[ilev])[dit];
            const FArrayBox& gdn = levGeoPtr->getCCgdn()[dit];
            const Box& region = levelGrids[dit];

            CH_assert(levGeoPtr->getBoxes().check(dit()));

            if (s_num_scal_comps > 0) {
                const FArrayBox& bFAB = (*amrB[ilev])[dit];
                FORT_COMPUTEENERGY(CHF_FRA1(eFAB,0),
                                   CHF_CONST_FRA(velFAB),
                                   CHF_CONST_FRA1(bFAB,0),
                                   CHF_BOX(region),
                                   CHF_CONST_FRA(gdn),
                                   CHF_CONST_REALVECT(levelDx));
            } else {
                FORT_COMPUTEKINETICENERGY(CHF_FRA1(eFAB,0),
                                          CHF_CONST_FRA(velFAB),
                                          CHF_BOX(region),
                                          CHF_CONST_FRA(gdn),
                                          CHF_CONST_REALVECT(levelDx));
            }
        } // end loop over grids

        levGeoPtr = levGeoPtr->getFinerPtr();
    } // end loop over levels

    // Now, integrate the energy
    globalEnergy = computeMappedSum(amrEnergy, *m_levGeoPtr);

    // Free memory (we no longer need the AMR energy data holder)
    for (int ilev = 0; ilev < amrEnergy.size(); ++ilev) {
        if (amrEnergy[ilev] != NULL) {
            delete amrEnergy[ilev];
            amrEnergy[ilev] = NULL;
        }
    }
    amrEnergy.resize(0);

    return globalEnergy;
}
