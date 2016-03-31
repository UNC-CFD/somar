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

// This is a rework of code found in the Chombo library. While these functions
// have been reproduced for ease of modification, some of these functions may
// not have been altered at all. To that end, you can find Chombo's copyright
// file at somar/Chombo/Copyright.txt.

#include "VelocityAMRPoissonOp.H"
#include "MappedAMRPoissonOpF_F.H"
#include "DivCurlGradF_F.H"
#include "ExtrapolationUtils.H"
#include "Constants.H"

//#define _USE_EXTRAP_AT_CF_
#define _CF_EXTRAP_ORDER_ 2

#if CH_SPACEDIM == 2
#   define D_HTERM(a,b) a
#else
#   define D_HTERM(a,b) a b
#endif


// -----------------------------------------------------------------------------
// Defines this op to simply calculate Laplacians and residuals.
// -----------------------------------------------------------------------------
void VelocityAMRPoissonOp::define(const LevelGeometry&                        a_levGeo,
                                  const Real                                  a_alpha,
                                  const Real                                  a_beta,
                                  const RefCountedPtr<LevelData<FArrayBox> >& a_lapDiagsPtr)
{
    // Just use the original version of this define with bogus BCs
    BCMethodHolder bogusBCs;
    MappedAMRPoissonOp::define(bogusBCs, a_levGeo, a_alpha, a_beta, a_lapDiagsPtr);
}


// -----------------------------------------------------------------------------
// Evaluate operator, inhomogeneous C/F boundary conditions
// note that operator applies C/F interpolation and copy bc's
// -----------------------------------------------------------------------------
void VelocityAMRPoissonOp::applyOpI (LevelData<FArrayBox>&        a_lhs,
                                     const LevelData<FArrayBox>&  a_phi,
                                     bool                         a_homogeneous,
                                     VelBCHolder*                 a_velBCPtr)
{
    CH_TIME("VelocityAMRPoissonOp::applyOpI");

    const DisjointBoxLayout& grids = a_lhs.disjointBoxLayout();
    DataIterator dit = a_phi.dataIterator();
    const int ncomps = a_phi.nComp();

    // Sanity check
    CH_assert(m_FCJgup->getBoxes().compatible(grids));

    // This is OK if we use it to only change ghost cells
    LevelData<FArrayBox>& phiRef = const_cast< LevelData<FArrayBox>& >(a_phi);

    // Begin exchange of valid data
    this->exchangeComplete(phiRef);

    // Loop through grids and perform calculations that do not need ghosts.
    for (dit.reset(); dit.ok(); ++dit) {
        const FArrayBox& phiFAB    = a_phi[dit];
        FArrayBox&       phiRefFAB = phiRef[dit];
        FArrayBox&       lhsFAB    = a_lhs[dit];
        const FluxBox&   JgupFB    = (*m_FCJgup)[dit];
        const FArrayBox& JinvFAB   = (*m_CCJinv)[dit];
        Box              valid     = grids[dit];

        // Create temp holders
        FluxBox fluxFB(valid, ncomps);
#ifndef NDEBUG
        fluxFB.setVal(quietNAN);
#endif

        // Extrapolate ghosts for non-diagonal derivatives
        FArrayBox extrapFAB(phiFAB.box(), ncomps);
        this->fillExtrap(extrapFAB, phiFAB, 2);

        // Set physical BCs
        if (a_velBCPtr != NULL) {
            a_velBCPtr->setGhosts(phiRefFAB, &extrapFAB, valid, m_domain, m_dx, dit(), &JgupFB, a_homogeneous, m_time);
        }

        // Compute fluxes
        if (!m_horizontalOp) {
            D_TERM(this->getFluxComplete(fluxFB[0], phiFAB, extrapFAB, surroundingNodes(valid,0), dit(), 0);,
                   this->getFluxComplete(fluxFB[1], phiFAB, extrapFAB, surroundingNodes(valid,1), dit(), 1);,
                   this->getFluxComplete(fluxFB[2], phiFAB, extrapFAB, surroundingNodes(valid,2), dit(), 2);)
        } else {
            D_HTERM(this->getFluxComplete(fluxFB[0], phiFAB, extrapFAB, surroundingNodes(valid,0), dit(), 0);,
                    this->getFluxComplete(fluxFB[1], phiFAB, extrapFAB, surroundingNodes(valid,1), dit(), 1);)
        }

        // Set boundary fluxes
        if (a_velBCPtr != NULL) {
            a_velBCPtr->setFluxes(fluxFB, NULL, valid, m_domain, m_dx, dit(), &JgupFB, a_homogeneous, m_time);
        }

        // Compute divergence
#if CH_SPACEDIM == 2
        if (!m_horizontalOp) {
            FORT_MAPPEDFLUXDIVERGENCE2D(CHF_FRA(lhsFAB),
                                        CHF_CONST_FRA(fluxFB[0]),
                                        CHF_CONST_FRA(fluxFB[1]),
                                        CHF_CONST_FRA1(JinvFAB,0),
                                        CHF_BOX(valid),
                                        CHF_CONST_REALVECT(m_dx));
        } else {
            MayDay::Error("Need to write FORT_MAPPEDFLUXDIVERGENCE1D");
        }
#elif CH_SPACEDIM == 3
        if (!m_horizontalOp) {
            FORT_MAPPEDFLUXDIVERGENCE3D(CHF_FRA(lhsFAB),
                                        CHF_CONST_FRA(fluxFB[0]),
                                        CHF_CONST_FRA(fluxFB[1]),
                                        CHF_CONST_FRA(fluxFB[2]),
                                        CHF_CONST_FRA1(JinvFAB,0),
                                        CHF_BOX(valid),
                                        CHF_CONST_REALVECT(m_dx));
        } else {
            FORT_MAPPEDFLUXDIVERGENCE2D(CHF_FRA(lhsFAB),
                                        CHF_CONST_FRA(fluxFB[0]),
                                        CHF_CONST_FRA(fluxFB[1]),
                                        CHF_CONST_FRA1(JinvFAB,0),
                                        CHF_BOX(valid),
                                        CHF_CONST_REALVECT(m_dx));
        }
#else
#error Bad CH_SPACEDIM
#endif

        // Right now, lhsFAB = beta*L[phiFAB].
        // We want lhsFAB = (alpha + beta*L)[phiFAB].
        // This utility function does just that for us.
        if (m_alpha != 0.0) {
            const Real betaDummy = 1.0;
            FORT_AXBYIP(CHF_FRA(lhsFAB),
                        CHF_CONST_FRA(phiFAB),
                        CHF_CONST_REAL(m_alpha),
                        CHF_CONST_REAL(betaDummy),
                        CHF_BOX(valid));
        }
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void VelocityAMRPoissonOp::AMROperatorNF(LevelData<FArrayBox>&        a_LofPhi,
                                         const LevelData<FArrayBox>&  a_phi,
                                         const LevelData<FArrayBox>&  a_phiCoarse,
                                         bool                         a_homogeneousPhysBC,
                                         VelBCHolder*                 a_velBCPtr)
{
    CH_TIME("VelocityAMRPoissonOp::AMROperatorNF");

    CH_assert(a_phi.isDefined());
    CH_assert(a_phiCoarse.isDefined());

    LevelData<FArrayBox>& phi = (LevelData<FArrayBox>&)a_phi;

    if (a_phiCoarse.isDefined()) {
        for (int comp = 0; comp < SpaceDim; comp++) {
            Interval intvl(comp, comp);

            LevelData<FArrayBox> phiComp;
            aliasLevelData(phiComp, &phi, intvl);

            LevelData<FArrayBox> phiCoarseComp;
            aliasLevelData(phiCoarseComp, (LevelData<FArrayBox>*) &a_phiCoarse, intvl);

            // m_interpWithCoarser.coarseFineInterp(phiComp, phiCoarseComp);
            this->interpCFGhosts(phiComp, &phiCoarseComp, false);
        }
    }

    // Physical BCs are applyed in applyOp
    this->applyOpI(a_LofPhi, a_phi, a_homogeneousPhysBC, a_velBCPtr);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void VelocityAMRPoissonOp::applyOp(LevelData<FArrayBox>&        a_lhs,
                                   const LevelData<FArrayBox>&  a_phi,
                                   const LevelData<FArrayBox>*  a_phiCoarsePtr,
                                   bool                         a_homogeneous,
                                   VelBCHolder*                 a_velBCPtr)
{
    if (a_phiCoarsePtr == NULL) {
        applyOpI(a_lhs, a_phi, a_homogeneous, a_velBCPtr);
    } else {
        AMROperatorNF(a_lhs, a_phi, *a_phiCoarsePtr, a_homogeneous, a_velBCPtr);
    }
}

