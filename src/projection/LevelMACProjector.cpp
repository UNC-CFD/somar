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
#include "LevelMACProjector.H"
#include "Debug.H"


// -----------------------------------------------------------------------------
// Default constructor
// This sets the solver parameters, but leaves object unusable.
// -----------------------------------------------------------------------------
LevelMACProjector::LevelMACProjector ()
: m_isDefined(false),
  m_levGeoPtr(NULL)
{
    const ProblemContext* ctx = ProblemContext::getInstance();

    m_solver.setAMRMGParameters(ctx->MACprojection_AMRMG_imin,
                                ctx->MACprojection_AMRMG_imax,
                                ctx->MACprojection_AMRMG_eps,
                                ctx->MACprojection_AMRMG_maxDepth,
                                ctx->MACprojection_AMRMG_num_smooth_precond,
                                ctx->MACprojection_AMRMG_num_smooth_down,
                                ctx->MACprojection_AMRMG_num_smooth_up,
                                ctx->MACprojection_AMRMG_num_smooth_bottom,
                                ctx->MACprojection_AMRMG_precondMode,
                                ctx->MACprojection_AMRMG_relaxMode,
                                ctx->MACprojection_AMRMG_numMG,
                                ctx->MACprojection_AMRMG_hang,
                                ctx->MACprojection_AMRMG_normThresh,
                                ctx->MACprojection_AMRMG_verbosity);

    m_solver.setBottomParameters(ctx->MACprojection_bottom_imax,
                                 ctx->MACprojection_bottom_numRestarts,
                                 ctx->MACprojection_bottom_eps,
                                 ctx->MACprojection_bottom_reps,
                                 ctx->MACprojection_bottom_hang,
                                 ctx->MACprojection_bottom_small,
                                 ctx->MACprojection_bottom_normType,
                                 ctx->MACprojection_bottom_verbosity);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
LevelMACProjector::~LevelMACProjector ()
{
    undefine();
}


// -----------------------------------------------------------------------------
// Allocates memory and leaves object useable.
// -----------------------------------------------------------------------------
void LevelMACProjector::define (LevelData<FArrayBox>*       a_phiPtr,
                                const LevelData<FArrayBox>* a_crsePhiPtr,
                                const PhysBCUtil&           a_physBCUtil,
                                const LevelGeometry&        a_levGeo,
                                const FillJgupInterface*    a_customFillJgupPtr)
{
    // Clear the old data.
    undefine();

    // Collect geometric structures
    m_levGeoPtr = &a_levGeo;
    const DisjointBoxLayout& grids = a_levGeo.getBoxes();

    const LevelGeometry* crseLevGeoPtr = a_levGeo.getCoarserPtr();
    const DisjointBoxLayout* crseGridsPtr = NULL;
    if (crseLevGeoPtr != NULL) {
        crseGridsPtr = &(crseLevGeoPtr->getBoxes());
        CH_assert(a_crsePhiPtr != NULL);
        CH_assert(a_crsePhiPtr->getBoxes() == *crseGridsPtr);
    }

    // Collect BCs.
    const ProblemContext* ctx = ProblemContext::getInstance();
    const bool isViscous = (ctx->nu > 0.0);

    m_solverBC = a_physBCUtil.MacPressureFuncBC();
    m_divBC = a_physBCUtil.uStarFuncBC(isViscous);
    m_gradBC = a_physBCUtil.gradMacPressureFuncBC();

    // Define CF-BC interpolator.
    if (crseGridsPtr != NULL) {
        m_cfInterp.define(grids,
                          crseGridsPtr,
                          a_levGeo.getDx(),
                          a_levGeo.getCrseRefRatio(),
                          1, // ncomp
                          a_levGeo.getDomain());
    }

    // Collect pressure pointers.
    CH_assert(a_phiPtr != NULL);
    CH_assert(a_phiPtr->isDefined());
    CH_assert(a_phiPtr->getBoxes() == grids);
    setLevelPressure(a_phiPtr, a_crsePhiPtr);

    // Define the pressure solver.
    int numLevels = ((a_crsePhiPtr == NULL)? 1: 2);
    m_solver.levelDefine(m_solverBC, a_levGeo, numLevels, a_customFillJgupPtr);

    // This object is ready to be used.
    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Frees memory and leaves object unuseable.
// -----------------------------------------------------------------------------
void LevelMACProjector::undefine ()
{
    if (m_isDefined) {
        // Base members
        m_time = BOGUS_TIME;
        m_pressure.clear();
        m_solver.undefine();

        // Our members
        m_levGeoPtr = NULL;
        m_cfInterp.clear();
        m_isDefined = false;
    }
}


// -----------------------------------------------------------------------------
// From BaseProjector:
//  Computes J^{-1}\partial_i(J u^i) over an AMR hierarchy.
//  This must be overriden or an error will be thrown.
//
// This must receive lmin = lmax because it is a single-level divergence op.
// -----------------------------------------------------------------------------
void LevelMACProjector::computeDiv (Vector<LevelData<FArrayBox>*>&     a_div,
                                    const Vector<LevelData<FluxBox>*>& a_flux,
                                    const int                          a_lmin,
                                    const int                          a_lmax) const
{
    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(a_lmin == a_lmax);
    CH_assert(a_div .size() >= a_lmax);
    CH_assert(a_flux.size() >= a_lmax);

    // Collect field references.
    LevelData<FArrayBox>& divRef = *(a_div[a_lmax]);
    LevelData<FluxBox>& fluxRef = *(a_flux[a_lmax]);

    // Compute the divergence
    Divergence::levelDivergenceMAC(divRef,
                                   fluxRef,
                                   *m_levGeoPtr,
                                   m_time,
                                   &m_divBC);
}


// -----------------------------------------------------------------------------
// From BaseProjector:
//  Computes Jg^{i,j}\partial_j(phi) over an AMR hierarchy.
//  This must be overriden or an error will be thrown.
//
// This must receive lmin = lmax because it is a single-level gradient op.
// Coarser phi data will be used for CF-BCs.
// -----------------------------------------------------------------------------
void LevelMACProjector::computeGrad (Vector<LevelData<FluxBox>*>&         a_flux,
                                     const Vector<LevelData<FArrayBox>*>& a_phi,
                                     const int                            a_lmin,
                                     const int                            a_lmax) const
{
    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(a_lmin == a_lmax);
    CH_assert(a_flux.size() >= a_lmax);
    CH_assert(a_phi .size() >= a_lmax);

    // Collect field references.
    LevelData<FluxBox>& fluxRef = *(a_flux[a_lmax]);
    LevelData<FArrayBox>& phiRef = *(a_phi[a_lmax]);
    const LevelData<FArrayBox>* crsePhiPtr = NULL;
    if (a_lmax > 0) {
        crsePhiPtr = a_phi[a_lmax-1];
        CH_assert(crsePhiPtr != NULL);
    }

    // Compute the gradient.
    Gradient::levelGradientMAC(fluxRef,
                               phiRef,
                               crsePhiPtr,
                               m_cfInterp,
                               *m_levGeoPtr,
                               m_time,
                               &m_gradBC);
}


// -----------------------------------------------------------------------------
// Applies vel = vel - corr. Dt is built into corr.
// -----------------------------------------------------------------------------
void LevelMACProjector::applyCorrection (Vector<LevelData<FluxBox>*>&       a_amrVel,
                                         const Vector<LevelData<FluxBox>*>& a_amrCorr,
                                         const Real                         a_dt,
                                         const int                          a_lmin,
                                         const int                          a_lmax) const
{
    const Real dtScale = (a_dt == 0.0)? -1.0: -a_dt;

    for (int lev = a_lmin; lev <= a_lmax; ++lev) {
        DataIterator dit = a_amrVel[lev]->dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            FluxBox& velFB = (*a_amrVel[lev])[dit];
            const FluxBox& corrFB = (*a_amrCorr[lev])[dit];

            D_TERM(velFB[0].plus(corrFB[0], dtScale);,
                   velFB[1].plus(corrFB[1], dtScale);,
                   velFB[2].plus(corrFB[2], dtScale);)
        }
    }
}
