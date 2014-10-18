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
#include "Curl.H"
#include "ExtrapolationUtils.H"
#include "computeMappedSum.H"
#include "MappedAMRPoissonOpFactory.H"
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
    m_velocityAMRPoissonOp.setTime(a_time);
    m_velocityAMRPoissonOp.applyOp(a_lapVel,
                                   tmpVel,
                                   a_crseVelPtr,
                                   isHomogeneous,
                                   NULL);

    // Extrap all ghosts
    extrapAllGhosts(a_lapVel,0);

    // And exchange
    if (a_lapVel.ghostVect() == m_tracingGhosts) {
        a_lapVel.exchange(m_tracingExCopier);
    } else if(a_lapVel.ghostVect() == IntVect::Unit) {
        a_lapVel.exchange(m_oneGhostExCopier);
    } else {
        a_lapVel.exchange();
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void AMRNavierStokes::computeLapScal (LevelData<FArrayBox>&       a_lapScal,
                                      LevelData<FArrayBox>&       a_scal,               // TODO: Make const!
                                      const LevelData<FArrayBox>* a_crseScalPtr,
                                      const CFRegion*             a_CFRegionPtr)
{
    // Apply the op
    BCMethodHolder scalPhysBCs = m_physBCPtr->diffusiveSourceFuncBC();
    m_scalarsAMRPoissonOp.setBC(scalPhysBCs);

    bool isHomogeneous = false;
    if (a_crseScalPtr != NULL) {
        m_scalarsAMRPoissonOp.AMROperatorNF(a_lapScal,
                                            a_scal,
                                            *a_crseScalPtr,
                                            isHomogeneous);
    } else {
        m_scalarsAMRPoissonOp.applyOpI(a_lapScal,
                                       a_scal,
                                       isHomogeneous);
    }

    // Extrap all ghosts
    extrapAllGhosts(a_lapScal,0);

    // And exchange
    if (a_lapScal.ghostVect() == m_tracingGhosts) {
        a_lapScal.exchange(m_tracingExCopier);
    } else if(a_lapScal.ghostVect() == IntVect::Unit) {
        a_lapScal.exchange(m_oneGhostExCopier);
    } else {
        a_lapScal.exchange();
    }
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
    LevelData<FArrayBox>& vel = *(const_cast<LevelData<FArrayBox>*>(m_vel_new_ptr));
    const DisjointBoxLayout& grids = vel.getBoxes();
    DataIterator dit = vel.dataIterator();
    const Interval& velComps = vel.interval();
    const RealVect& dx = m_levGeoPtr->getDx();
    const bool isViscous = (s_nu > 0.0);

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
    vel.exchange(m_oneGhostExCopier);

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
            // Create references on this grid for convienience
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
