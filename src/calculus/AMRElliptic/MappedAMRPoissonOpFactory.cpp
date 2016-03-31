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
#include "MappedAMRPoissonOpFactory.H"
#include "MappedAMRPoissonOpF_F.H"
#include "AnisotropicRefinementTools.H"
#include "MappedCoarseAverage.H"
#include "Constants.H"
#include "ProblemContext.H"

#include "Jacobi.H"
#include "GSRB.H"


// Setting this to true will force full derivative and contraction
// calculations using all metric components. If set to false, these
// calculations will be optimized to only use non-zero entries.
// This only makes a difference when using an orthogonal system.
static const bool forceFullCalcs = false;


// -----------------------------------------------------------------------------
// Full AMR factory define
// Assembles the hierarchy of data via a_baseLevGeoPtr->getAMRLevGeos().
// Really, you can send in any levGeoPtr and this function will find the base.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOpFactory::define(const LevelGeometry*     a_baseLevGeoPtr,
                                       Real                     a_alpha,
                                       Real                     a_beta,
                                       BCMethodHolder&          a_bc,
                                       int                      a_maxDepth,
                                       int                      a_preCondSmoothIters,
                                       int                      a_precondMode,
                                       int                      a_relaxMode,
                                       bool                     a_horizontalFactory,
                                       const FillJgupInterface* a_customFillJgupPtr)
{
    CH_TIME("MappedAMRPoissonOpFactory::define 1");

    // Start clean.
    clear();

    // Collect the AMR levGeos
    CH_assert(a_baseLevGeoPtr != NULL);
    Vector<const LevelGeometry*> vLevGeoPtrs = a_baseLevGeoPtr->getAMRLevGeos();
    CH_assert(vLevGeoPtrs.size() > 0);

    // Collect the coarsest ProblemDomain
    const ProblemDomain& coarseDomain = vLevGeoPtrs[0]->getDomain();

    // Collect the coarsest derivScale
    RealVect coarseDx = vLevGeoPtrs[0]->getDx();

    // Collect the AMR data
    Vector<DisjointBoxLayout> vGrids(0);
    Vector<IntVect> vRefRatios(0);
    for (int ilev = 0; ilev < vLevGeoPtrs.size(); ++ilev) {
        CH_assert(vLevGeoPtrs[ilev] != NULL);
        vGrids.push_back(vLevGeoPtrs[ilev]->getBoxes());
        vRefRatios.push_back(vLevGeoPtrs[ilev]->getFineRefRatio());
    }

    // Call the full define function
    this->define(coarseDomain, vGrids, vRefRatios, coarseDx, a_bc,
                 a_maxDepth, a_preCondSmoothIters, a_precondMode,
                 a_alpha, a_beta, vLevGeoPtrs, a_relaxMode,
                 a_horizontalFactory, a_customFillJgupPtr);
}


// -----------------------------------------------------------------------------
// Full AMRMG factory define
// -----------------------------------------------------------------------------
void MappedAMRPoissonOpFactory::define(const ProblemDomain&                a_coarseDomain,
                                       const Vector<DisjointBoxLayout>&    a_grids,
                                       const Vector<IntVect>&              a_refRatios,
                                       const RealVect&                     a_coarsedx,
                                       BCMethodHolder&                     a_bc,
                                       int                                 a_maxDepth,
                                       int                                 a_preCondSmoothIters,
                                       int                                 a_precondMode,
                                       Real                                a_alpha,
                                       Real                                a_beta,
                                       const Vector<const LevelGeometry*>& a_vlevGeoPtr,
                                       int                                 a_relaxMode,
                                       bool                                a_horizontalFactory,
                                       const FillJgupInterface*            a_customFillJgupPtr)
{
    CH_TIME("MappedAMRPoissonOpFactory::define 2");

    // Start clean.
    clear();

    // Set all member data that is not a vector over AMR levels
    CH_assert(a_maxDepth >= -1);
    m_maxDepth = a_maxDepth;

    CH_assert(a_preCondSmoothIters >= 0);
    m_preCondSmoothIters = a_preCondSmoothIters;
    m_precondMode = a_precondMode;

    CH_assert(-1 <= a_relaxMode && a_relaxMode < ProblemContext::RelaxMode::NUM_RELAX_MODES);
    m_relaxMode = a_relaxMode;

    m_bc = a_bc;
    m_alpha = a_alpha;
    m_beta = a_beta;

    m_horizontalFactory = a_horizontalFactory;
    m_maskedMaxCoarse = MappedAMRPoissonOp::s_maxCoarse * IntVect::Unit; // The smallest allowable grid
    if (a_horizontalFactory) m_maskedMaxCoarse[CH_SPACEDIM-1] = 1;

    // Calculate the number of extant levels
    int numlevels = 0;
    for (int i = 0; i < a_grids.size(); ++i) {
        if (a_grids[i].size() > 0)
            ++numlevels;
    }

    // Resize vector data holders to include only extant levels
    m_boxes.resize(numlevels);
    m_refRatios.resize(numlevels);
    m_dx.resize(numlevels);
    m_domains.resize(numlevels);
    m_exchangeCopiers.resize(numlevels);
    m_cfregion.resize(numlevels);
    m_vlevGeoPtr.resize(numlevels);
    m_vvJgup.resize(numlevels);
    m_vvJinv.resize(numlevels);
    m_vvlapDiag.resize(numlevels);

    m_customFillJgupPtr = a_customFillJgupPtr;

    // Set level 0 data
    m_boxes[0] = a_grids[0];
    m_refRatios[0] = a_refRatios[0];
    m_dx[0] = a_coarsedx;
    m_domains[0] = a_coarseDomain;

    IntVect ghostVect = IntVect::Unit;
    if (a_horizontalFactory) ghostVect[CH_SPACEDIM-1] = 0;

    m_exchangeCopiers[0].define(a_grids[0], a_grids[0], m_domains[0], ghostVect, true);
    // m_exchangeCopiers[0].exchangeDefine(a_grids[0], ghostVect);
    // if (LevelGeometry::isDiagonal()) {
    //     m_exchangeCopiers[0].trimEdges(a_grids[0], ghostVect);
    // }

    m_cfregion[0].define(a_grids[0], m_domains[0]);
    m_vlevGeoPtr[0] = a_vlevGeoPtr[0];
    m_vvJgup[0].resize(0);
    m_vvJinv[0].resize(0);
    m_vvlapDiag[0].resize(0);

    // Set finer level data
    for (int i = 1; i < numlevels; ++i) {
        m_boxes[i] = a_grids[i];

        m_refRatios[i] = a_refRatios[i];
        m_dx[i] = m_dx[i-1] / RealVect(m_refRatios[i-1]);

        refine(m_domains[i], m_domains[i-1], m_refRatios[i-1]);

        if (a_grids[i].isClosed()) {
            m_exchangeCopiers[i].define(a_grids[i], a_grids[i], m_domains[i], ghostVect, true);
            // m_exchangeCopiers[i].exchangeDefine(a_grids[i], ghostVect);
            // if (LevelGeometry::isDiagonal()) {
            //     m_exchangeCopiers[i].trimEdges(a_grids[i], ghostVect);
            // }

            m_cfregion[i].define(a_grids[i], m_domains[i]);
        }

        m_vlevGeoPtr[i] = a_vlevGeoPtr[i];

        m_vvJgup[i].resize(0);
        m_vvJinv[i].resize(0);
        m_vvlapDiag[i].resize(0);
    }

    // By default, do not hard-code the coarse-level dx.
    m_useForceDxCrse = false;
    m_forceDxCrse = RealVect(D_DECL(quietNAN, quietNAN, quietNAN));
}


// -----------------------------------------------------------------------------
// Single level AMR factory define -- used by the Leptic solver.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOpFactory::define(const RefCountedPtr<LevelData<FluxBox> >&   a_JgupPtr,
                                       const RefCountedPtr<LevelData<FArrayBox> >& a_JinvPtr,
                                       const Copier&                               a_copier,
                                       const CFRegion&                             a_cfregion,
                                       const RealVect&                             a_dx,
                                       Real                                        a_alpha,
                                       Real                                        a_beta,
                                       BCMethodHolder&                             a_bc,
                                       int                                         a_maxDepth,
                                       int                                         a_preCondSmoothIters,
                                       int                                         a_precondMode,
                                       int                                         a_relaxMode,
                                       bool                                        a_horizontalFactory)
{
    // Start clean.
    clear();

    // Set all member data that is not a vector over AMR levels
    CH_assert(a_maxDepth >= -1);
    m_maxDepth = a_maxDepth;

    CH_assert(a_preCondSmoothIters >= 0);
    m_preCondSmoothIters = a_preCondSmoothIters;
    m_precondMode = a_precondMode;

    CH_assert(-1 <= a_relaxMode && a_relaxMode < ProblemContext::RelaxMode::NUM_RELAX_MODES);
    m_relaxMode = a_relaxMode;

    m_bc = a_bc;
    m_alpha = a_alpha;
    m_beta = a_beta;

    m_horizontalFactory = a_horizontalFactory;
    m_maskedMaxCoarse = MappedAMRPoissonOp::s_maxCoarse * IntVect::Unit;    // The smallest allowable grid
    if (a_horizontalFactory) m_maskedMaxCoarse[CH_SPACEDIM-1] = 1;

    const int numlevels = 1; // This is a single level define

    // Resize vector data holders to include only extant levels
    m_boxes.resize(numlevels);
    m_refRatios.resize(numlevels);
    m_dx.resize(numlevels);
    m_domains.resize(numlevels);
    m_exchangeCopiers.resize(numlevels);
    m_cfregion.resize(numlevels);
    m_vlevGeoPtr.resize(numlevels);
    m_vvJgup.resize(numlevels);
    m_vvJinv.resize(numlevels);
    m_vvlapDiag.resize(numlevels);

    // Set level 0 data
    const DisjointBoxLayout& grids = a_JgupPtr->getBoxes();
    m_boxes[0] = grids;
    m_refRatios[0] = IntVect::Unit;
    m_dx[0] = a_dx;
    m_domains[0] = grids.physDomain();

    IntVect ghostVect = IntVect::Unit;
    if (a_horizontalFactory) ghostVect[CH_SPACEDIM-1] = 0;

    m_exchangeCopiers[0] = a_copier;
    m_cfregion[0] = a_cfregion;

    m_vlevGeoPtr[0] = NULL;
    m_vvJgup[0].resize(1);
    m_vvJinv[0].resize(1);
    m_vvlapDiag[0].resize(1);

    m_vvJgup[0][0] = a_JgupPtr;
    m_vvJinv[0][0] = a_JinvPtr;
    m_vvlapDiag[0][0] = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(grids, 1));

    DataIterator dit = grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& lapDiagFAB = (*m_vvlapDiag[0][0])[dit];
        const FluxBox& JgupFAB = (*a_JgupPtr)[dit];
        const FArrayBox& JinvFAB = (*a_JinvPtr)[dit];

#if CH_SPACEDIM == 2
        FORT_FILLMAPPEDLAPDIAG2D(CHF_FRA1(lapDiagFAB,0),
                                 CHF_CONST_FRA(JgupFAB[0]),
                                 CHF_CONST_FRA(JgupFAB[1]),
                                 CHF_CONST_FRA1(JinvFAB,0),
                                 CHF_BOX(lapDiagFAB.box()),
                                 CHF_CONST_REALVECT(a_dx));
#else
        FORT_FILLMAPPEDLAPDIAG3D(CHF_FRA1(lapDiagFAB,0),
                                 CHF_CONST_FRA(JgupFAB[0]),
                                 CHF_CONST_FRA(JgupFAB[1]),
                                 CHF_CONST_FRA(JgupFAB[2]),
                                 CHF_CONST_FRA1(JinvFAB,0),
                                 CHF_BOX(lapDiagFAB.box()),
                                 CHF_CONST_REALVECT(a_dx));
#endif
    }

    m_bypassValidation = true;

    // By default, do not hard-code the coarse-level dx.
    m_useForceDxCrse = false;
    m_forceDxCrse = RealVect(D_DECL(quietNAN, quietNAN, quietNAN));
}


// -----------------------------------------------------------------------------
// The default destructor
// -----------------------------------------------------------------------------
MappedAMRPoissonOpFactory::~MappedAMRPoissonOpFactory ()
{
    clear();
}


// -----------------------------------------------------------------------------
// Clears all memory pools and undefines this object.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOpFactory::clear ()
{
    m_domains.clear();
    m_boxes.clear();
    m_exchangeCopiers.clear();
    m_cfregion.clear();

    m_vlevGeoPtr.clear();
    m_customFillJgupPtr = NULL;

    m_vvJgup.clear();
    m_vvJinv.clear();
    m_vvlapDiag.clear();
}


// -----------------------------------------------------------------------------
// When setting up horizontal solves, this level's LevGeo will not be tied
// to a coarse-level LevGeo. Therefore, we need to be told what the coarse
// level dx is.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOpFactory::forceDxCrse (const RealVect& a_dxCrse)
{
    m_useForceDxCrse = true;
    m_forceDxCrse = a_dxCrse;
}



// -----------------------------------------------------------------------------
// Create an operator at a coarsened index space.
// Return NULL if no such Multigrid level can be created at this a_depth.
// If a_homoOnly = true, then only homogeneous BCs will be needed.
// Caller is responsible for deletion.
// -----------------------------------------------------------------------------
MappedMGLevelOp<LevelData<FArrayBox> >*
MappedAMRPoissonOpFactory::MGnewOp(const ProblemDomain&   a_indexSpace,
                                   const int              a_depth,
                                   const bool             a_homoOnly,
                                   Vector<IntVect>*       a_allMGRefRatiosPtr,
                                   const Vector<IntVect>* a_forceAllMGRefRatiosPtr)
{
    CH_TIME("MappedAMRPoissonOpFactory::MGnewOp");
    // pout() << "MappedAMRPoissonOpFactory::MGnewOp" << endl;

    // Are we past the requested max depth?
    if ((m_maxDepth >= 0) && (a_depth > m_maxDepth)) {
        if (a_forceAllMGRefRatiosPtr != NULL) {
            if (a_depth <= a_forceAllMGRefRatiosPtr->size()) {
                MayDay::Error("You must make the maxDepth large enough to accomodate the mini V-cycles");
            }
        }
        return NULL;
    }

    // Find the indexSpace's level
    int ref;
    for (ref = 0; ref < m_domains.size(); ++ref) {
        if (a_indexSpace.domainBox() == m_domains[ref].domainBox()) {
            break;
        }
    }
    CH_assert(ref != m_domains.size()); // didn't find domain

    // Get dx at one AMR level coarser than ref level.
    // This will be used by homogeneousCFInterp
    RealVect dxCrse(D_DECL(quietNAN, quietNAN, quietNAN));
    if (ref > 0) {
        dxCrse = m_dx[ref-1];
    }

    // Were we given an explicit dxCrse? If so, there had better not be a coarser level!
    if (m_useForceDxCrse) {
        CH_assert(ref == 0);
        dxCrse = m_forceDxCrse;
    }

    // Time to refine...

    ProblemDomain domain(m_domains[ref]);
    RealVect dx = m_dx[ref];
    IntVect coarsening = IntVect::Unit;
    IntVect mgRefRatio;

    int forcedDepth = -1;
    if (a_forceAllMGRefRatiosPtr != NULL) {
        forcedDepth = a_forceAllMGRefRatiosPtr->size();
    }

    if (a_forceAllMGRefRatiosPtr != NULL && a_depth <= forcedDepth) {
        // Force the mgRefRatios to be whatever the user wants.
        CH_assert(a_forceAllMGRefRatiosPtr->size() >= a_depth-1);
        for (int i = 0; i < a_depth; ++i) {
            mgRefRatio = (*a_forceAllMGRefRatiosPtr)[i];

            ProblemDomain fineMGDomain = domain;
            coarsen(domain, fineMGDomain, mgRefRatio);
            dx *= RealVect(mgRefRatio);
            coarsening *= mgRefRatio;
        }

        // Is this amount of coarsening allowed?
        if (coarsening.product() > 1 && !coarsenable(m_boxes[ref], coarsening * m_maskedMaxCoarse)) {
            D_TERM(int minBlockFactor = coarsening[0] * m_maskedMaxCoarse[0];,
                   minBlockFactor = max(minBlockFactor, coarsening[1] * m_maskedMaxCoarse[1]);,
                   minBlockFactor = max(minBlockFactor, coarsening[2] * m_maskedMaxCoarse[2]);)
            ostringstream sstr;
            sstr << "Could not coarsen grids for mini V-cycle. Your block factor needs to be at least ";
            sstr << minBlockFactor;

            pout() << "m_boxes[" << ref << "] = " << m_boxes[ref] << endl;
            MayDay::Error(sstr.str().c_str());

            return NULL;
        }

    } else {
        // We are not forcing the mgRefRatio...

// #define ISOTROPIC_REFINEMENT_ONLY
#ifdef ISOTROPIC_REFINEMENT_ONLY
        // Coarsen the reference level's ProblemDomain and dx to the MG depth
        mgRefRatio = IntVect(D_DECL(2,2,2));
        if (m_horizontalFactory) mgRefRatio[CH_SPACEDIM-1] = 1;

        for (int i = 0; i < a_depth; ++i) {
            domain.coarsen(2);
            // domain.define(coarsen(domain.domainBox(), mgRefRatio));
            dx *= RealVect(mgRefRatio);
            coarsening *= mgRefRatio;
        }

        // Is this amount of coarsening allowed?
        if (coarsening.product() > 1 && !coarsenable(m_boxes[ref], coarsening * m_maskedMaxCoarse)) {
            return NULL;
        }

#else
        // Coarsen the reference level's ProblemDomain and dx to the MG depth
        for (int i = 0; i < a_depth; ++i) {
            if (a_allMGRefRatiosPtr != NULL && i < a_depth-1) {
                CH_assert((a_depth == 0 && a_allMGRefRatiosPtr->size() == 0) || (a_depth > 0 && a_allMGRefRatiosPtr->size() == a_depth-1));
                // CH_assert(a_allMGRefRatiosPtr->size() > i);
                mgRefRatio = (*a_allMGRefRatiosPtr)[i];
                D_TERM(CH_assert(mgRefRatio[0] == 1 || mgRefRatio[0] == 2);,
                       CH_assert(mgRefRatio[1] == 1 || mgRefRatio[1] == 2);,
                       CH_assert(mgRefRatio[2] == 1 || mgRefRatio[2] == 2);)
                CH_assert(mgRefRatio.product() > 1);
            } else {
                // Figure out which directions will benefit from coarsening.
                // This is kind of a dumb way of doing things, but it's good enough for now.
                const int dims = m_horizontalFactory? SpaceDim-1: SpaceDim;

                Real maxDx = 0.0;
                for (int dir = 0; dir < dims; ++dir) {
                    maxDx = max(maxDx, dx[dir]);
                }

                mgRefRatio = IntVect::Unit;
                for (int dir = 0; dir < dims; ++dir) {
                   if (dx[dir] <= maxDx / 2.0) mgRefRatio[dir] = 2;
                }

                // Is anisotropic coarsening worth it?
                if (mgRefRatio.product() == 1) {
                    mgRefRatio = IntVect(D_DECL(2,2,2));
                    if (m_horizontalFactory) mgRefRatio[CH_SPACEDIM-1] = 1;
                }
            }

            ProblemDomain crseMGDomain = domain;
            coarsen(domain, crseMGDomain, mgRefRatio);
            dx *= RealVect(mgRefRatio);
            coarsening *= mgRefRatio;
        }

        // Is this amount of coarsening allowed?
        if (coarsening.product() > 1 && !coarsenable(m_boxes[ref], coarsening * m_maskedMaxCoarse)) {
            // Try coarsening in less directions...
            // pout() << "m_boxes[" << ref << "] = " << m_boxes[ref] << endl;

            // 1. Refine back up one MG level.
            ProblemDomain fineMGDomain = domain;
            refine(fineMGDomain, domain, mgRefRatio);
            dx /= RealVect(mgRefRatio);
            coarsening /= mgRefRatio;

            // 2. Figure out which directions have enough cells to be coarsened.
            const int dims = m_horizontalFactory? SpaceDim-1: SpaceDim;
            mgRefRatio = IntVect::Unit;
            for (int dir = 0; dir < dims; ++dir) {
                mgRefRatio[dir] = 2;
                if (!coarsenable(m_boxes[ref], coarsening * m_maskedMaxCoarse * mgRefRatio)) {
                    mgRefRatio[dir] = 1;
                }
            }

            // Can anything be coarsened?
            if (mgRefRatio.product() == 1) return NULL;

            // 3. Figure out which of those directions will benefit from coarsening.
            // This is kind of a dumb way of doing things, but it's good enough for now.
            int refDir;
            for (refDir = 0; refDir < dims; ++refDir) {
               if (mgRefRatio[refDir] == 1) break;
            }
            RealVect aspectRatio(dx / dx[refDir]);
            for (int dir = 0; dir < dims; ++dir) {
                if (mgRefRatio[dir] > 1 && aspectRatio[dir] > 0.5) mgRefRatio[dir] = 1;
            }

            // Is it worth it?
            if (mgRefRatio.product() == 1) return NULL;

            // 4. We can coarsen!
            coarsen(domain, fineMGDomain, mgRefRatio);
            dx *= RealVect(mgRefRatio);
            coarsening *= mgRefRatio;

            // It never hurts to double check.
            if (coarsening.product() > 1 && !coarsenable(m_boxes[ref], coarsening * m_maskedMaxCoarse)) {
                return NULL;
            }
        }
#endif
    }

    // Inform the caller of the refRatio.
    if (a_allMGRefRatiosPtr != NULL && a_depth > 0) {
        a_allMGRefRatiosPtr->push_back(mgRefRatio);
    }

    // Coarsen the grids
    DisjointBoxLayout layout;
    coarsen(layout, m_boxes[ref], coarsening);

    // Coarsen the copier and cfregion
    Copier ex = m_exchangeCopiers[ref];
    CFRegion cfregion = m_cfregion[ref];
    if (coarsening.product() > 1) {
        coarsen(ex, coarsening);
        coarsen(cfregion, coarsening);
    }

    // Fill the metric fields if needed
    this->validateMetricPtrs(ref, a_depth, mgRefRatio, coarsening, layout);

    // We have created all of our coarse level data.
    // Time to build the coarse operator.
    MappedAMRPoissonOp* newOp = new MappedAMRPoissonOp;
    newOp->define(layout, dx, domain, m_bc, ex, cfregion);

    newOp->m_mgDepth = a_depth;
    newOp->m_mgCrseRefRatio = mgRefRatio;
    newOp->m_horizontalOp = m_horizontalFactory;

    newOp->m_alpha  = m_alpha;
    newOp->m_beta   = m_beta;
    newOp->m_aCoef  = m_alpha;
    newOp->m_bCoef  = m_beta;
    newOp->m_dxCrse = dxCrse;
    newOp->m_time   = BOGUS_TIME;
    newOp->m_preCondSmoothIters = m_preCondSmoothIters;

    newOp->m_FCJgup  = m_vvJgup[ref][a_depth];
    newOp->m_CCJinv  = m_vvJinv[ref][a_depth];
    newOp->m_lapDiag = m_vvlapDiag[ref][a_depth];

    newOp->m_isDiagonal = m_vlevGeoPtr[ref]->isDiagonal() && !forceFullCalcs;

    IntVect activeDirs = IntVect::Unit;
    if (m_horizontalFactory) activeDirs[SpaceDim-1] = 0;

    // Set the relaxation method.
    if (m_relaxMode == ProblemContext::RelaxMode::NORELAX) {
        // If you don't have your relaxation, you can't have any MG.
        MayDay::Error("How can you have any MG if you don't have relaxation?");
    } else if (m_relaxMode == ProblemContext::RelaxMode::JACOBI) {
        newOp->m_relaxPtr = new Jacobi(newOp, newOp->m_alpha, newOp->m_beta);
    } else if (m_relaxMode == ProblemContext::RelaxMode::LEVEL_GSRB) {
        newOp->m_relaxPtr = new LevelGSRB;
    } else if (m_relaxMode == ProblemContext::RelaxMode::LOOSE_GSRB) {
        newOp->m_relaxPtr = new LooseGSRB;
    } else if (m_relaxMode == ProblemContext::RelaxMode::LINE_GSRB) {
        newOp->m_relaxPtr = new LineGSRB;
    } else {
        MayDay::Warning("Relaxation method not yet converted...using LevelGSRB");
        newOp->m_relaxPtr = new LevelGSRB;
    }
    newOp->m_relaxPtr->define(newOp->m_alpha,
                              newOp->m_beta,
                              newOp->m_dx,
                              newOp->m_dxCrse,
                              newOp->m_FCJgup,
                              newOp->m_CCJinv,
                              newOp->m_lapDiag,
                              newOp->m_bc,
                              newOp->m_exchangeCopier,
                              newOp->m_cfregion,
                              newOp->m_isDiagonal,
                              activeDirs);

    // Set the relaxation method for the preconditioner.
    newOp->m_precondMode = m_precondMode;
    if (m_precondMode == ProblemContext::PrecondMode::None) {
        newOp->m_precondRelaxPtr = NULL;

    } else if (m_precondMode == ProblemContext::PrecondMode::DiagRelax) {
        newOp->m_precondRelaxPtr = newOp->m_relaxPtr;

    } else if (m_precondMode == ProblemContext::PrecondMode::DiagLineRelax) {
        newOp->m_precondRelaxPtr = new LineGSRB;
        newOp->m_precondRelaxPtr->define(newOp->m_alpha,
                                         newOp->m_beta,
                                         newOp->m_dx,
                                         newOp->m_dxCrse,
                                         newOp->m_FCJgup,
                                         newOp->m_CCJinv,
                                         newOp->m_lapDiag,
                                         newOp->m_bc,
                                         newOp->m_exchangeCopier,
                                         newOp->m_cfregion,
                                         newOp->m_isDiagonal,
                                         activeDirs);
    } else {
        MayDay::Error("Bad m_precondMode");
    }

    // Set the restriction strategy.
    newOp->m_restrictPtr = new FullWeightingPS(newOp->m_CCJinv);

    // Set the prolongation strategy.
    bool removeAvg = false;
    do {
        // If we survive each test, we will remove the average.
        const int dims = m_horizontalFactory? SpaceDim-1: SpaceDim;
        const DisjointBoxLayout& grids = newOp->m_FCJgup->getBoxes();
        DataIterator dit = grids.dataIterator();
        CFRegion& cfregion = newOp->m_cfregion;

        // 1. Check if this is a parabolic operator.
        if (m_alpha != 0.0) break;

        // 2. Check physical BCs
        bool hasNullSpace = m_bc.hasNullSpace(grids, cfregion, activeDirs);
        if (hasNullSpace == false) break;

        // 3. Check for CF interfaces
        bool hasCFInterface = false;
        for (int dir = 0; dir < dims; ++dir) {
            for (dit.reset(); dit.ok(); ++dit) {
                const CFIVS& loCFIVS = cfregion.loCFIVS(dit(), dir);
                const CFIVS& hiCFIVS = cfregion.hiCFIVS(dit(), dir);

                hasCFInterface |= !loCFIVS.isEmpty();
                hasCFInterface |= !hiCFIVS.isEmpty();
            }
            if (hasCFInterface) break;
        }
        if (hasCFInterface) break;

        // Well, it looks like we made it. Remove the average.
        removeAvg = true;
    } while(0);

    if (removeAvg) {
        const Real dxProduct = m_horizontalFactory?
                               (D_TERM(1.0, *dx[0], *dx[1])):
                               (D_TERM(dx[0], *dx[1], *dx[2]));
        newOp->m_prolongPtr = new ZeroAvgConstInterpPS(dxProduct, newOp->m_CCJinv);
        // pout() << "Using ZeroAvgConstInterpPS" << endl;
    } else {
        newOp->m_prolongPtr = new ConstInterpPS;
        // pout() << "Using ConstInterpPS" << endl;
    }

    // TEMPORARY!!!
    ostringstream ss;
    ss << "(" << ref << ", " << a_depth << ")";
    newOp->m_levid = ss.str();

    return (MappedMGLevelOp<LevelData<FArrayBox> >*)newOp;
}


// -----------------------------------------------------------------------------
// Return a new operator. This is done with a new call.
// Caller is responsible for deletion.
// -----------------------------------------------------------------------------
MappedAMRLevelOp<LevelData<FArrayBox> >*
MappedAMRPoissonOpFactory::AMRnewOp(const ProblemDomain& a_indexSpace)
{
    CH_TIME("MappedAMRPoissonOpFactory::AMRnewOp");
    // pout() << "MappedAMRPoissonOpFactory::AMRnewOp" << endl;

    // Allocate a new operator
    MappedAMRPoissonOp* newOp = new MappedAMRPoissonOp;

    // Find the indexSpace's level
    int ref;
    for (ref = 0; ref < m_domains.size(); ref++) {
        if (a_indexSpace.domainBox() == m_domains[ref].domainBox()) {
            break;
        }
    }
    CH_assert(ref < m_domains.size());

    // Fill the metric fields if needed
    if (!m_bypassValidation) {
        this->validateMetricPtrs(ref, 0, IntVect::Unit, IntVect::Unit, m_boxes[ref]);
    }

    // Call the appropriate define function
    RealVect dxCrse(D_DECL(quietNAN, quietNAN, quietNAN));
    if (ref == 0) {
        if (m_domains.size() == 1) {
            // No coarser and no finer
            newOp->define(m_boxes[0], m_dx[0],
                          a_indexSpace, m_bc,
                          m_exchangeCopiers[0],
                          m_cfregion[0]);
        } else {
            // Bottom AMR level
            const IntVect& dummyRat = IntVect::Unit;    // argument so compiler can find right function
            const IntVect& refToFiner = m_refRatios[0]; // actual refinement ratio
            newOp->define(m_boxes[0], m_boxes[1], m_dx[0],
                          dummyRat, refToFiner,
                          a_indexSpace, m_bc,
                          m_exchangeCopiers[0],
                          m_cfregion[0]);
        }
    } else if (ref == m_domains.size()-1) {
        dxCrse = m_dx[ref-1];

        // Top AMR level
        newOp->define(m_boxes[ref], m_boxes[ref-1], m_dx[ref],
                      m_refRatios[ref-1],
                      a_indexSpace, m_bc,
                      m_exchangeCopiers[ref],
                      m_cfregion[ref]);
    } else if (ref == m_domains.size()) {
        MayDay::Abort("Did not find a domain to match AMRnewOp(const ProblemDomain& a_indexSpace)");
    } else {
        // Intermediate AMR level, full define
        dxCrse = m_dx[ref-1];
        newOp->define(m_boxes[ref], m_boxes[ref+1], m_boxes[ref-1], m_dx[ref],
                      m_refRatios[ref-1], m_refRatios[ref],
                      a_indexSpace, m_bc,
                      m_exchangeCopiers[ref],
                      m_cfregion[ref]);
    }

    // Were we given an explicit dxCrse? If so, there had better not be a coarser level!
    if (m_useForceDxCrse) {
        CH_assert(ref == 0);
        dxCrse = m_forceDxCrse;
    }

    // Set public data
    newOp->m_mgCrseRefRatio = IntVect::Unit;
    newOp->m_horizontalOp = m_horizontalFactory;

    newOp->m_alpha  = m_alpha;
    newOp->m_beta   = m_beta;
    newOp->m_aCoef  = m_alpha;
    newOp->m_bCoef  = m_beta;
    newOp->m_dxCrse = dxCrse;
    newOp->m_time   = BOGUS_TIME;
    newOp->m_preCondSmoothIters = m_preCondSmoothIters;

    newOp->m_FCJgup  = m_vvJgup[ref][0];
    newOp->m_CCJinv  = m_vvJinv[ref][0];
    newOp->m_lapDiag = m_vvlapDiag[ref][0];

    newOp->m_isDiagonal = m_vlevGeoPtr[ref]->isDiagonal() && !forceFullCalcs;

    IntVect activeDirs = IntVect::Unit;
    if (m_horizontalFactory) activeDirs[SpaceDim-1] = 0;

    // Set the relaxation method.
    if (m_relaxMode == ProblemContext::RelaxMode::NORELAX) {
        newOp->m_relaxPtr = NULL;
    } else if (m_relaxMode == ProblemContext::RelaxMode::JACOBI) {
        newOp->m_relaxPtr = new Jacobi(newOp, newOp->m_alpha, newOp->m_beta);
    } else if (m_relaxMode == ProblemContext::RelaxMode::LEVEL_GSRB) {
        newOp->m_relaxPtr = new LevelGSRB;
    } else if (m_relaxMode == ProblemContext::RelaxMode::LOOSE_GSRB) {
        newOp->m_relaxPtr = new LooseGSRB;
    } else if (m_relaxMode == ProblemContext::RelaxMode::LINE_GSRB) {
        newOp->m_relaxPtr = new LineGSRB;
    } else {
        MayDay::Warning("Relaxation method not yet converted...using LevelGSRB");
        newOp->m_relaxPtr = new LevelGSRB;
    }
    if (newOp->m_relaxPtr != NULL) {
        newOp->m_relaxPtr->define(newOp->m_alpha,
                                  newOp->m_beta,
                                  newOp->m_dx,
                                  newOp->m_dxCrse,
                                  newOp->m_FCJgup,
                                  newOp->m_CCJinv,
                                  newOp->m_lapDiag,
                                  newOp->m_bc,
                                  newOp->m_exchangeCopier,
                                  newOp->m_cfregion,
                                  newOp->m_isDiagonal,
                                  activeDirs);
    }

    // Set the relaxation method for the preconditioner.
    newOp->m_precondMode = m_precondMode;
    if (m_precondMode == ProblemContext::PrecondMode::None) {
        newOp->m_precondRelaxPtr = NULL;

    } else if (m_precondMode == ProblemContext::PrecondMode::DiagRelax) {
        newOp->m_precondRelaxPtr = newOp->m_relaxPtr;

    } else if (m_precondMode == ProblemContext::PrecondMode::DiagLineRelax) {
        newOp->m_precondRelaxPtr = new LineGSRB;
        newOp->m_precondRelaxPtr->define(newOp->m_alpha,
                                         newOp->m_beta,
                                         newOp->m_dx,
                                         newOp->m_dxCrse,
                                         newOp->m_FCJgup,
                                         newOp->m_CCJinv,
                                         newOp->m_lapDiag,
                                         newOp->m_bc,
                                         newOp->m_exchangeCopier,
                                         newOp->m_cfregion,
                                         newOp->m_isDiagonal,
                                         activeDirs);
    } else {
        MayDay::Error("Bad m_precondMode");
    }


    // Set the restriction strategy.
    newOp->m_restrictPtr = new FullWeightingPS(newOp->m_CCJinv);

    // Set the prolongation strategy.
    bool removeAvg = false;
    do {
        // If we survive each test, we will remove the average.
        const int dims = m_horizontalFactory? SpaceDim-1: SpaceDim;
        const DisjointBoxLayout& grids = newOp->m_FCJgup->getBoxes();
        DataIterator dit = grids.dataIterator();
        CFRegion& cfregion = newOp->m_cfregion;

        // 1. Check if this is a parabolic operator.
        if (m_alpha != 0.0) break;

        // 2. Check physical BCs
        bool hasNullSpace = m_bc.hasNullSpace(grids, cfregion, activeDirs);
        if (hasNullSpace == false) break;

        // 3. Check for CF interfaces
        bool hasCFInterface = false;
        for (int dir = 0; dir < dims; ++dir) {
            for (dit.reset(); dit.ok(); ++dit) {
                const CFIVS& loCFIVS = cfregion.loCFIVS(dit(), dir);
                const CFIVS& hiCFIVS = cfregion.hiCFIVS(dit(), dir);

                hasCFInterface |= !loCFIVS.isEmpty();
                hasCFInterface |= !hiCFIVS.isEmpty();
            }
            if (hasCFInterface) break;
        }
        if (hasCFInterface) break;

        // Well, it looks like we made it. Remove the average.
        removeAvg = true;
    } while(0);

    if (removeAvg) {
        const RealVect& dx = newOp->m_dx;
        const Real dxProduct = m_horizontalFactory?
                               (D_TERM(1.0, *dx[0], *dx[1])):
                               (D_TERM(dx[0], *dx[1], *dx[2]));
        newOp->m_prolongPtr = new ZeroAvgConstInterpPS(dxProduct, newOp->m_CCJinv);
        // pout() << "Using ZeroAvgConstInterpPS" << endl;
    } else {
        newOp->m_prolongPtr = new ConstInterpPS;
        // pout() << "Using ConstInterpPS" << endl;
    }

    // TEMPORARY!!!
    ostringstream ss;
    ss << "(" << ref << ", 0)";
    newOp->m_levid = ss.str();

    return (MappedAMRLevelOp<LevelData<FArrayBox> >*)newOp;
}


// -----------------------------------------------------------------------------
// Returns the refinement ratio to the next finer level
// -----------------------------------------------------------------------------
IntVect MappedAMRPoissonOpFactory::getFineRefRatio(const ProblemDomain& a_domain) const
{
    IntVect retval(D_DECL(-1,-1,-1));
    bool found = false;

    for (int ilev = 0; ilev < m_domains.size(); ++ilev) {
        if (m_domains[ilev].domainBox() == a_domain.domainBox()) {
            retval = m_refRatios[ilev];
            found = true;
        }
    }

    if (!found) {
        MayDay::Abort("Domain not found in AMR hierarchy");
    }

    return retval;
}


// -----------------------------------------------------------------------------
// Upon completion, the metric field ptrs will be defined properly.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOpFactory::validateMetricPtrs(const int                a_AMRlevel,
                                                   const int                a_MGdepth,
                                                   const IntVect&           a_mgRefRatio,
                                                   const IntVect&           a_coarsening,
                                                   const DisjointBoxLayout& a_layout)
{
    CH_TIME("MappedAMRPoissonOpFactory::validateMetricPtrs");

    // Sanity checks
    CH_assert(a_AMRlevel < m_boxes.size());
    CH_assert((m_maxDepth < 0) || (a_MGdepth <= m_maxDepth));

    // Do we need to resize the vector? If so, do this without destroying the shallow data.
    for (int d = 0; d <= a_MGdepth; ++d) {
        if (m_vvJgup[a_AMRlevel].size() < d+1) {
            m_vvJgup[a_AMRlevel].push_back(RefCountedPtr<LevelData<FluxBox> >(NULL));
            m_vvJinv[a_AMRlevel].push_back(RefCountedPtr<LevelData<FArrayBox> >(NULL));
            m_vvlapDiag[a_AMRlevel].push_back(RefCountedPtr<LevelData<FArrayBox> >(NULL));
        }
    }

    // Create references for convenience.
    RefCountedPtr<LevelData<FluxBox> >& Jgup = m_vvJgup[a_AMRlevel][a_MGdepth];
    RefCountedPtr<LevelData<FArrayBox> >& Jinv = m_vvJinv[a_AMRlevel][a_MGdepth];
    RefCountedPtr<LevelData<FArrayBox> >& lapDiag = m_vvlapDiag[a_AMRlevel][a_MGdepth];

    // Check if some of the fields are already defined in levGeo.
    // If so, use them and only define what is necessary.
    bool useCachedData = false;
    if (m_vlevGeoPtr.size() > a_AMRlevel) {
        if (m_vlevGeoPtr[a_AMRlevel] != NULL) {
            if (m_vlevGeoPtr[a_AMRlevel]->getBoxes() == a_layout && m_customFillJgupPtr == NULL) {
                useCachedData = true;
            }
        }
    }
    if (useCachedData) {
        CH_assert(a_MGdepth == 0);

        // Grab Jgup
        Jgup = m_vlevGeoPtr[a_AMRlevel]->getFCJgupPtr();
        CH_assert(!Jgup.isNull());

        // Grab Jinv
        Jinv = m_vlevGeoPtr[a_AMRlevel]->getCCJinvPtr();
        CH_assert(!Jinv.isNull());

        // Define lapDiags
        lapDiag = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>);
        lapDiag->define(a_layout, 1, IntVect::Zero);

        const RealVect& levelDx = m_vlevGeoPtr[a_AMRlevel]->getDx();
        DataIterator dit = lapDiag->dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            CH_assert(m_vlevGeoPtr[a_AMRlevel]->getBoxes().check(dit()));
            CH_assert(Jgup->getBoxes().check(dit()));
            CH_assert(Jinv->getBoxes().check(dit()));

            FArrayBox& lapDiagFAB = (*lapDiag)[dit];
            const FluxBox& JgupFAB = (*Jgup)[dit];
            const FArrayBox& JinvFAB = (*Jinv)[dit];

#if CH_SPACEDIM == 2
            if (!m_horizontalFactory) {
                FORT_FILLMAPPEDLAPDIAG2D(
                    CHF_FRA1(lapDiagFAB,0),
                    CHF_CONST_FRA(JgupFAB[0]),
                    CHF_CONST_FRA(JgupFAB[1]),
                    CHF_CONST_FRA1(JinvFAB,0),
                    CHF_BOX(lapDiagFAB.box()),
                    CHF_CONST_REALVECT(levelDx));
            } else {
                FORT_FILLMAPPEDLAPDIAG1D(
                    CHF_FRA1(lapDiagFAB,0),
                    CHF_CONST_FRA1(JgupFAB[0],0),
                    CHF_CONST_FRA1(JinvFAB,0),
                    CHF_BOX(lapDiagFAB.box()),
                    CHF_CONST_REAL(levelDx[0]));
           }
#elif CH_SPACEDIM == 3
            if (!m_horizontalFactory) {
                FORT_FILLMAPPEDLAPDIAG3D(
                    CHF_FRA1(lapDiagFAB,0),
                    CHF_CONST_FRA(JgupFAB[0]),
                    CHF_CONST_FRA(JgupFAB[1]),
                    CHF_CONST_FRA(JgupFAB[2]),
                    CHF_CONST_FRA1(JinvFAB,0),
                    CHF_BOX(lapDiagFAB.box()),
                    CHF_CONST_REALVECT(levelDx));
            } else {
                FORT_FILLMAPPEDLAPDIAG2D(
                    CHF_FRA1(lapDiagFAB,0),
                    CHF_CONST_FRA(JgupFAB[0]),
                    CHF_CONST_FRA(JgupFAB[1]),
                    CHF_CONST_FRA1(JinvFAB,0),
                    CHF_BOX(lapDiagFAB.box()),
                    CHF_CONST_REALVECT(levelDx));
            }
#else
#error Bad CH_SPACEDIM
#endif
        }

    } else {
        // The levGeo has nothing for us. Check if memory has already been allocated.
        if (Jgup.isNull()) {
            CH_assert(Jinv.isNull());
            CH_assert(lapDiag.isNull());
            Jgup = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>);
            Jinv = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>);
            lapDiag = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>);
        }

        // Do we need to redefine the metric holders at this MG level?
        if(!(Jgup->getBoxes() == a_layout)) {
            // Call defines
            Jgup->define(a_layout, SpaceDim, IntVect::Zero);
            Jinv->define(a_layout, 1, IntVect::Zero);
            lapDiag->define(a_layout, 1, IntVect::Zero);

            // Fill the fields
            this->fill_MGfields(a_AMRlevel, a_MGdepth, a_mgRefRatio, a_coarsening);
        }
    }

    // Sanity checks
    CH_assert(Jgup->getBoxes() == a_layout);
    CH_assert(Jinv->getBoxes() == a_layout);
    CH_assert(lapDiag->getBoxes() == a_layout);
}


// -----------------------------------------------------------------------------
// Fill the metric fields with data.
// This function was made special for filling MG levels.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOpFactory::fill_MGfields(const int      a_AMRlevel,
                                              const int      a_MGdepth,
                                              const IntVect& a_mgRefRatio,
                                              const IntVect& a_coarsening)
{
    CH_TIME("MappedAMRPoissonOpFactory::fill_MGfields");

    // Sanity checks
    CH_assert(a_AMRlevel < m_dx.size());
    CH_assert((m_maxDepth < 0) || (a_MGdepth <= m_maxDepth));

    // Grab the fields
    RefCountedPtr<LevelData<FluxBox> >& Jgup = m_vvJgup[a_AMRlevel][a_MGdepth];
    RefCountedPtr<LevelData<FArrayBox> >& Jinv = m_vvJinv[a_AMRlevel][a_MGdepth];
    RefCountedPtr<LevelData<FArrayBox> >& lapDiag = m_vvlapDiag[a_AMRlevel][a_MGdepth];

    // Grab the geometry data
    const RealVect mgDx = m_dx[a_AMRlevel] * a_coarsening;
    LevelGeometry localLevGeo(m_dx[a_AMRlevel]);
    const GeoSourceInterface& geoSource = *localLevGeo.getGeoSourcePtr();

    if (a_MGdepth == 0) {
        // There is no finer data. Get levGeo to calculate fields.

        DataIterator dit = Jgup->dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& JinvFAB = (*Jinv)[dit];
            FluxBox& JgupFAB = (*Jgup)[dit];
            FArrayBox& lapDiagFAB = (*lapDiag)[dit];

            // FC Jgup
            if (m_customFillJgupPtr == NULL) {
                // Use the levGeo's version of Jgup.
                for (int adir = 0; adir < SpaceDim; ++adir) {
                    for (int bdir = 0; bdir < SpaceDim; ++bdir) {
                        geoSource.fill_Jgup(JgupFAB[adir],
                                            bdir, adir, bdir,
                                            mgDx);
                    }
                }
            } else {
                // Use a custom Jgup
                for (int adir = 0; adir < SpaceDim; ++adir) {
                    for (int bdir = 0; bdir < SpaceDim; ++bdir) {
                        m_customFillJgupPtr->fill_Jgup(JgupFAB[adir],
                                                       bdir, adir, bdir,
                                                       mgDx, 1.0);
                    }
                }
            }

            // CC 1/J
            geoSource.fill_Jinv(JinvFAB, 0, mgDx);

            // CC Laplacian diagonals
#if CH_SPACEDIM == 2
            if (!m_horizontalFactory) {
                FORT_FILLMAPPEDLAPDIAG2D(
                    CHF_FRA1(lapDiagFAB,0),
                    CHF_CONST_FRA(JgupFAB[0]),
                    CHF_CONST_FRA(JgupFAB[1]),
                    CHF_CONST_FRA1(JinvFAB,0),
                    CHF_BOX(lapDiagFAB.box()),
                    CHF_CONST_REALVECT(mgDx));
            } else {
                FORT_FILLMAPPEDLAPDIAG1D(
                    CHF_FRA1(lapDiagFAB,0),
                    CHF_CONST_FRA1(JgupFAB[0],0),
                    CHF_CONST_FRA1(JinvFAB,0),
                    CHF_BOX(lapDiagFAB.box()),
                    CHF_CONST_REAL(mgDx[0]));
            }
#elif CH_SPACEDIM == 3
            if (!m_horizontalFactory) {
                FORT_FILLMAPPEDLAPDIAG3D(
                    CHF_FRA1(lapDiagFAB,0),
                    CHF_CONST_FRA(JgupFAB[0]),
                    CHF_CONST_FRA(JgupFAB[1]),
                    CHF_CONST_FRA(JgupFAB[2]),
                    CHF_CONST_FRA1(JinvFAB,0),
                    CHF_BOX(lapDiagFAB.box()),
                    CHF_CONST_REALVECT(mgDx));
            } else {
                FORT_FILLMAPPEDLAPDIAG2D(
                    CHF_FRA1(lapDiagFAB,0),
                    CHF_CONST_FRA(JgupFAB[0]),
                    CHF_CONST_FRA(JgupFAB[1]),
                    CHF_CONST_FRA1(JinvFAB,0),
                    CHF_BOX(lapDiagFAB.box()),
                    CHF_CONST_REALVECT(mgDx));
            }
#else
#error Bad CH_SPACEDIM
#endif
        }
    } else {
        // Finer data exists. Average fine data down to coarse data.

        // Grab the fine fields
        const RefCountedPtr<LevelData<FluxBox> >& fineJgup = m_vvJgup[a_AMRlevel][a_MGdepth-1];
        const RefCountedPtr<LevelData<FArrayBox> >& fineJinv = m_vvJinv[a_AMRlevel][a_MGdepth-1];
        // const RefCountedPtr<LevelData<FArrayBox> >& fineLapDiag = m_vvlapDiag[a_AMRlevel][a_MGdepth-1];

        const DisjointBoxLayout& fineGrids = fineJgup->getBoxes();
        const DisjointBoxLayout& crseGrids = Jgup->getBoxes();
        const IntVect& ghostVect = Jinv->ghostVect();

        MappedCoarseAverageFace coarseAvgFace;
        MappedCoarseAverage coarseAvg;

        // FC Jgup
        coarseAvgFace.define(fineGrids, Jgup->nComp(), a_mgRefRatio);
        coarseAvgFace.averageToCoarse(*Jgup, *fineJgup);

        // CC 1/J
        coarseAvg.define(fineGrids, crseGrids, Jinv->nComp(), a_mgRefRatio, ghostVect);
        coarseAvg.averageToCoarseHarmonic(*Jinv, *fineJinv);

        // CC Laplacian diagonals
        // coarseAvg.averageToCoarse(*lapDiag, *fineLapDiag);  // Does not work. Regenerate diags.
        DataIterator dit = Jinv->dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& lapDiagFAB = (*lapDiag)[dit];
            const FluxBox& JgupFAB = (*Jgup)[dit];
            const FArrayBox& JinvFAB = (*Jinv)[dit];

#if CH_SPACEDIM == 2
            if (!m_horizontalFactory) {
                FORT_FILLMAPPEDLAPDIAG2D(
                    CHF_FRA1(lapDiagFAB,0),
                    CHF_CONST_FRA(JgupFAB[0]),
                    CHF_CONST_FRA(JgupFAB[1]),
                    CHF_CONST_FRA1(JinvFAB,0),
                    CHF_BOX(lapDiagFAB.box()),
                    CHF_CONST_REALVECT(mgDx));
            } else {
                FORT_FILLMAPPEDLAPDIAG1D(
                    CHF_FRA1(lapDiagFAB,0),
                    CHF_CONST_FRA1(JgupFAB[0],0),
                    CHF_CONST_FRA1(JinvFAB,0),
                    CHF_BOX(lapDiagFAB.box()),
                    CHF_CONST_REAL(mgDx[0]));
            }
#elif CH_SPACEDIM == 3
            if (!m_horizontalFactory) {
                FORT_FILLMAPPEDLAPDIAG3D(
                    CHF_FRA1(lapDiagFAB,0),
                    CHF_CONST_FRA(JgupFAB[0]),
                    CHF_CONST_FRA(JgupFAB[1]),
                    CHF_CONST_FRA(JgupFAB[2]),
                    CHF_CONST_FRA1(JinvFAB,0),
                    CHF_BOX(lapDiagFAB.box()),
                    CHF_CONST_REALVECT(mgDx));
            } else {
                FORT_FILLMAPPEDLAPDIAG2D(
                    CHF_FRA1(lapDiagFAB,0),
                    CHF_CONST_FRA(JgupFAB[0]),
                    CHF_CONST_FRA(JgupFAB[1]),
                    CHF_CONST_FRA1(JinvFAB,0),
                    CHF_BOX(lapDiagFAB.box()),
                    CHF_CONST_REALVECT(mgDx));
            }
#else
#error Bad CH_SPACEDIM
#endif
        }
    }
}

