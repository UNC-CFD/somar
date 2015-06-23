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
#include "LevelGeometry.H"
#include "ProblemContext.H"
#include "Constants.H"
#include "SubspaceF_F.H" // For the vertical integral
#include "CornerCopier.H"
#include "CubicSpline.H"
#include "MiscUtils.H"
#include "BoxIterator.H"
#include "NodeAMRIO.H"
#include "Debug.H"

// The advection scheme needs more ghost cells than any other piece of code.
// We must accomodate...
#include "AdvectUtil.H"

// Define the constants used to look up tensor indices.
#if CH_SPACEDIM == 1
    const int LevelGeometry::m_tensorCompCC2[SpaceDim][SpaceDim]
    = {{0}};
    const int LevelGeometry::m_symTensorCompCC2[SpaceDim][SpaceDim]
    = {{0}};
    const int LevelGeometry::m_symTensorCompCC3[SpaceDim][SpaceDim][SpaceDim]
    = {{{0}}};
#elif CH_SPACEDIM == 2
    const int LevelGeometry::m_tensorCompCC2[SpaceDim][SpaceDim]
    = {{0,1},{2,3}};
    const int LevelGeometry::m_symTensorCompCC2[SpaceDim][SpaceDim]
    = {{0,1},{1,2}};
    const int LevelGeometry::m_symTensorCompCC3[SpaceDim][SpaceDim][SpaceDim]
    = {{{0,1},{1,2}},{{3,4},{4,5}}};
#elif CH_SPACEDIM == 3
    const int LevelGeometry::m_tensorCompCC2[SpaceDim][SpaceDim]
    = {{0,1,2},{3,4,5},{6,7,8}};
    const int LevelGeometry::m_symTensorCompCC2[SpaceDim][SpaceDim]
    = {{0,1,2},{1,3,4},{2,4,5}};
    const int LevelGeometry::m_symTensorCompCC3[SpaceDim][SpaceDim][SpaceDim]
    = { {{ 0, 1, 2},{ 1, 3, 4},{ 2, 4, 5}},
        {{ 6, 7, 8},{ 7, 9,10},{ 8,10,11}},
        {{12,13,14},{13,15,16},{14,16,17}} };
#else
#   error SpaceDim must be one, two, or three.
#endif


// The metric source
RefCountedPtr<GeoSourceInterface> LevelGeometry::s_metricSourcePtr;

// The default coordinate mapping
int LevelGeometry::s_coordMap = ProblemContext::CoordMap::UNDEFINED;

// The domain's physical extents
RealVect LevelGeometry::s_domainLength = RealVect::Zero;

// Used to define the cached LevelDatas
const IntVect LevelGeometry::s_ghostVectFC = ADVECT_GROW * IntVect::Unit;
const IntVect LevelGeometry::s_ghostVectCC = ADVECT_GROW * IntVect::Unit;

RealVect                  LevelGeometry::s_lev0dXi = RealVect::Zero; // dXi at the base level.
LevelData<NodeFArrayBox>* LevelGeometry::s_lev0xPtr = NULL;          // physCoor on the base level.
LevelData<NodeFArrayBox>* LevelGeometry::s_lev0d2xPtr = NULL;        // The spline's 2nd derivatives.

// The fields
LevelGeometry::t_CCJMap       LevelGeometry::s_CCJMap;
LevelGeometry::t_CCJinvMap    LevelGeometry::s_CCJinvMap;
LevelGeometry::t_FCgupMap     LevelGeometry::s_FCgupMap;
LevelGeometry::t_FCJgupMap    LevelGeometry::s_FCJgupMap;
LevelGeometry::t_CCgdnMap     LevelGeometry::s_CCgdnMap;
LevelGeometry::t_bdryNormMaps LevelGeometry::s_bdryNormMaps;


// -----------------------------------------------------------------------------
// Sets the coordinate system.
// -----------------------------------------------------------------------------
void LevelGeometry::staticDefine ()
{
    // Once this function is called, it should not be called again.
    if (isStaticDefined()) return;

    const ProblemContext* ctx = ProblemContext::getInstance();

    s_domainLength = ctx->domainLength;
    s_lev0dXi = s_domainLength / RealVect(ctx->nx);
    s_coordMap = ctx->coordMap;
    s_metricSourcePtr = RefCountedPtr<GeoSourceInterface>(ctx->newGeoSourceInterface());

    // Sanity check
    CH_assert(isStaticDefined());
}


// -----------------------------------------------------------------------------
// Clear all memory used by the static maps. Calling this before the program
// ends makes for easier memory leak diagnosis.
// -----------------------------------------------------------------------------
void LevelGeometry::staticUndefine ()
{
    if (!isStaticDefined()) return;

    s_domainLength = RealVect::Zero;
    s_coordMap = ProblemContext::CoordMap::UNDEFINED;
    s_metricSourcePtr = RefCountedPtr<GeoSourceInterface>(NULL);

    s_lev0dXi = RealVect::Zero;

    delete s_lev0xPtr;
    s_lev0xPtr = NULL;

    delete s_lev0d2xPtr;
    s_lev0d2xPtr = NULL;

    s_CCJMap.clear();
    s_CCJinvMap.clear();
    s_FCgupMap.clear();
    s_FCJgupMap.clear();
    s_CCgdnMap.clear();
    s_bdryNormMaps.clear();
}


// -----------------------------------------------------------------------------
// Weak constructor
// -----------------------------------------------------------------------------
LevelGeometry::LevelGeometry ()
{
    if (!isStaticDefined()) {
        this->staticDefine();
    }

    m_finerPtr = NULL;
    m_coarserPtr = NULL;
    m_dXi = RealVect::Zero; // This will be used by define().
    m_refToLev0 = IntVect::Zero;
    m_fineRefRatio = IntVect::Unit;
    m_crseRefRatio = IntVect::Unit;
}


// -----------------------------------------------------------------------------
// Full constructor
// -----------------------------------------------------------------------------
LevelGeometry::LevelGeometry (const RealVect& a_dXi)
{
    if (!isStaticDefined()) {
        this->staticDefine();
    }

    m_finerPtr = NULL;
    m_coarserPtr = NULL;

    this->define(a_dXi);
}


// -----------------------------------------------------------------------------
// Puts this Geometry object in a usable state. Called by the full constructor.
// This is the single-level version
// -----------------------------------------------------------------------------
void LevelGeometry::define (const RealVect& a_dXi)
{
    // Sanity checks
    CH_assert(isStaticDefined());
    D_TERM(
        CH_assert(a_dXi[0] != 0.0);,
        CH_assert(a_dXi[1] != 0.0);,
        CH_assert(a_dXi[2] != 0.0);
    )

    // Collect the mapping parameters
    for (int dir = 0; dir < SpaceDim; ++dir) {
        m_dXi[dir] = a_dXi[dir];

        if (m_dXi[dir] <= s_lev0dXi[dir]) {
            // This level is NOT coarser than level 0.
            // Compute the refinement ratio from level 0 to this level and
            // make sure we don't have float to int errors.
            m_refToLev0[dir] = int(s_lev0dXi[dir] / m_dXi[dir]);
            CH_assert(m_refToLev0[dir] * m_dXi[dir] == s_lev0dXi[dir]);

        } else {
            // This level IS coarser than level 0.
            // For now, just set ref to something that will throw an error.
            m_refToLev0[dir] = quietNAN;
        }
    }

    // we will set up the level 0 spline data when we regrid.
}


// -----------------------------------------------------------------------------
// Puts this Geometry object in a usable state. Called by the full constructor.
// This version links this with the coarser level, which sets the refRatios.
// The finer pointer will still be NULL. It is expected that these defines will
// be called from coarsest to finest.
// -----------------------------------------------------------------------------
void LevelGeometry::define (const RealVect&      a_dXi,
                            const LevelGeometry* a_coarserPtr)
{
    this->define(a_dXi);

    this->setCoarserPtr(a_coarserPtr);
    this->setFinerPtr(NULL);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
LevelGeometry::~LevelGeometry ()
{;}


// -----------------------------------------------------------------------------
// Regrids the current geometry data
// -----------------------------------------------------------------------------
void LevelGeometry::regrid (const DisjointBoxLayout& a_newGrids)
{
    CH_TIME("LevelGeometry::regrid");
    CH_assert(this->isDefined());

    // Get the new ProblemDomain and make sure it hasn't changed, not
    // because it will cause problems, but because there is no reason
    // for that to happpen and usually means something is ill-defined.
    ProblemDomain domain(a_newGrids.physDomain());
    CH_assert(m_grids.physDomain().isEmpty() || domain == m_grids.physDomain());

    // Do nothing if the grids haven't changed
    if(m_grids == a_newGrids) return;

    // Clean up the map by removing all fields that still use oldGrids.
    // It is assumed that these RefCountedPtrs are unique. It may be nice to
    // check this by including a !isNonUnique() assert, but I'll save that
    // for the day I begin dereferencing NULL pointers.
    s_CCJMap.erase(m_grids);
    s_CCJinvMap.erase(m_grids);
    s_FCgupMap.erase(m_grids);
    s_FCJgupMap.erase(m_grids);
    s_CCgdnMap.erase(m_grids);

    // This flushes the entire cache. I could only erase elements that belong to the
    // regridded levels, but I tried that and it requires a lot of code refactoring.
    // Besides, this is efficient enough.
    s_bdryNormMaps.clear();

    // Redefine using the new grids
    m_grids = a_newGrids;

    // If this is level 0, then update the GeoSourceInterface's cache.
    if (m_dXi == s_lev0dXi) {
        s_metricSourcePtr->suggestLev0Grids(m_grids);
    }

    // Check if the fields are in the field map.
    // If not, create new ones.
    m_CCJPtr      = this->createCCJPtr(m_grids);
    m_CCJinvPtr   = this->createCCJinvPtr(m_grids);
    m_FCgupPtr    = this->createFCgupPtr(m_grids);
    m_FCJgupPtr   = this->createFCJgupPtr(m_grids);
    m_CCgdnPtr    = this->createCCgdnPtr(m_grids);
}


// -----------------------------------------------------------------------------
// Regrids the current geometry data with flat (vertically averaged) metrics.
// -----------------------------------------------------------------------------
void LevelGeometry::regridVertAvg (const DisjointBoxLayout& a_newGrids,
                                   const Box&               a_fullDomainBox)
{
    CH_TIME("LevelGeometry::regridVertAvg");
    CH_assert(this->isDefined());

    // Get the new ProblemDomain and make sure it hasn't changed, not
    // because it will cause problems, but because there is no reason
    // for that to happpen and usually means something is ill-defined.
    ProblemDomain domain(a_newGrids.physDomain());
    CH_assert(m_grids.physDomain().isEmpty() || domain == m_grids.physDomain());

    // Do nothing if the grids haven't changed
    if(m_grids == a_newGrids) return;

    // Clean up the map by removing all fields that still use oldGrids.
    // It is assumed that these RefCountedPtrs are unique. It may be nice to
    // check this by including a !isNonUnique() assert, but I'll save that
    // for the day I begin dereferencing NULL pointers.
    s_CCJMap.erase(m_grids);
    s_CCJinvMap.erase(m_grids);
    s_FCgupMap.erase(m_grids);
    s_FCJgupMap.erase(m_grids);
    s_CCgdnMap.erase(m_grids);

    // This flushes the entire cache. I could only erase elements that belong to the
    // regridded levels, but I tried that and it requires a lot of code refactoring.
    // Besides, this is efficient enough.
    s_bdryNormMaps.clear();

    // Redefine using the new grids
    m_grids = a_newGrids;

    // Check if the fields are in the field map.
    // If not, create new ones.
    m_CCJPtr      = t_CCJPtr(NULL);
    m_FCgupPtr    = t_FCgupPtr(NULL);
    m_CCgdnPtr    = t_CCgdnPtr(NULL);
    m_CCJinvPtr   = this->createVertAvgCCJinvPtr(m_grids, a_fullDomainBox);
    m_FCJgupPtr   = this->createVertAvgFCJgupPtr(m_grids, a_fullDomainBox);
}


// -----------------------------------------------------------------------------
// Clears all data holders and grids. Sets coarser/finer ptrs to NULL.
// -----------------------------------------------------------------------------
void LevelGeometry::reset ()
{
    m_grids       = DisjointBoxLayout();
    m_CCJPtr      = t_CCJPtr(NULL);
    m_CCJinvPtr   = t_CCJinvPtr(NULL);
    m_FCJgupPtr   = t_FCJgupPtr(NULL);
    m_CCgdnPtr    = t_CCgdnPtr(NULL);
    s_bdryNormMaps.clear();
}


// -----------------------------------------------------------------------------
// Find the J field in the map or create a new one.
// -----------------------------------------------------------------------------
LevelGeometry::t_CCJPtr LevelGeometry::createCCJPtr(const DisjointBoxLayout& a_grids)
{
    CH_TIME("LevelGeometry::createCCJPtr");

    // Search for the field int the map.
    t_CCJPtr& thisFieldPtr = s_CCJMap[a_grids];

    // If the field exists in the map, just return it's pointer.
    if (!thisFieldPtr.isNull()) return thisFieldPtr;

    Box dombox = a_grids.physDomain().domainBox();

    // The field did not exist. Allocate and define a new one.
    thisFieldPtr = t_CCJPtr(new LevelData<FArrayBox>);
    thisFieldPtr->define(a_grids, 1, s_ghostVectCC);

    // Fill the new field.
    DataIterator dit = thisFieldPtr->dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        this->fill_J((*thisFieldPtr)[dit]);
    }

    // Return the new field pointer.
    return thisFieldPtr;
}


// -----------------------------------------------------------------------------
// Find the 1/J field in the map or create a new one.
// -----------------------------------------------------------------------------
LevelGeometry::t_CCJinvPtr LevelGeometry::createCCJinvPtr(const DisjointBoxLayout& a_grids)
{
    CH_TIME("LevelGeometry::createCCJinvPtr");

    // Search for the field int the map.
    t_CCJinvPtr& thisFieldPtr = s_CCJinvMap[a_grids];

    // If the field exists in the map, just return it's pointer.
    if (!thisFieldPtr.isNull()) return thisFieldPtr;

    Box dombox = a_grids.physDomain().domainBox();

    // The field did not exist. Allocate and define a new one.
    thisFieldPtr = t_CCJinvPtr(new LevelData<FArrayBox>);
    thisFieldPtr->define(a_grids, 1, s_ghostVectCC);

    // Fill the new field.
    DataIterator dit = thisFieldPtr->dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        this->fill_Jinv((*thisFieldPtr)[dit]);
    }

    // Return the new field pointer.
    return thisFieldPtr;
}


// -----------------------------------------------------------------------------
// Find the 1/J field in the map or create a new one.
// -----------------------------------------------------------------------------
LevelGeometry::t_CCJinvPtr LevelGeometry::createVertAvgCCJinvPtr(const DisjointBoxLayout& a_grids,
                                                                 const Box&               a_fullDomainBox)
{
    CH_TIME("LevelGeometry::createVertAvgCCJinvPtr");

    // Search for the field int the map.
    t_CCJinvPtr& thisFieldPtr = s_CCJinvMap[a_grids];

    // If the field exists in the map, just return it's pointer.
    if (!thisFieldPtr.isNull()) return thisFieldPtr;

    Box dombox = a_grids.physDomain().domainBox();

    // The field did not exist. Allocate and define a new one.
    thisFieldPtr = t_CCJinvPtr(new LevelData<FArrayBox>);
    thisFieldPtr->define(a_grids, 1); // Changed from s_ghostVectCC on Mar 23, 2014

    // Fill the new field.
    DataIterator dit = thisFieldPtr->dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        (*thisFieldPtr)[dit].setVal(1.0);   // Because we are solving L[phi] = J*rhs.
    }

    // Return the new field pointer.
    return thisFieldPtr;
}


// -----------------------------------------------------------------------------
// Find the gup field in the map or create a new one.
// -----------------------------------------------------------------------------
LevelGeometry::t_FCgupPtr LevelGeometry::createFCgupPtr(const DisjointBoxLayout& a_grids)
{
    CH_TIME("LevelGeometry::createFCgupPtr");

    // Search for the field int the map.
    t_FCgupPtr& thisFieldPtr = s_FCgupMap[a_grids];

    // If the field exists in the map, just return it's pointer.
    if (!thisFieldPtr.isNull()) return thisFieldPtr;

    Box dombox = a_grids.physDomain().domainBox();

    // The field did not exist. Allocate and define a new one.
    thisFieldPtr = t_FCgupPtr(new LevelData<FluxBox>);
    thisFieldPtr->define(a_grids, SpaceDim, s_ghostVectFC);

    // Fill the new field.
    DataIterator dit = thisFieldPtr->dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        this->fill_gup((*thisFieldPtr)[dit]);
    }

    // Return the new field pointer.
    return thisFieldPtr;
}


// -----------------------------------------------------------------------------
// Find the Jgup field in the map or create a new one.
// -----------------------------------------------------------------------------
LevelGeometry::t_FCJgupPtr LevelGeometry::createFCJgupPtr(const DisjointBoxLayout& a_grids)
{
    CH_TIME("LevelGeometry::createFCJgupPtr");

    // Search for the field int the map.
    t_FCJgupPtr& thisFieldPtr = s_FCJgupMap[a_grids];

    // If the field exists in the map, just return it's pointer.
    if (!thisFieldPtr.isNull()) return thisFieldPtr;

    Box dombox = a_grids.physDomain().domainBox();

    // The field did not exist. Allocate and define a new one.
    thisFieldPtr = t_FCJgupPtr(new LevelData<FluxBox>);
    thisFieldPtr->define(a_grids, SpaceDim, s_ghostVectFC);

    // Fill the new field.
    DataIterator dit = thisFieldPtr->dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        this->fill_Jgup((*thisFieldPtr)[dit]);
    }

    // Return the new field pointer.
    return thisFieldPtr;
}


// -----------------------------------------------------------------------------
// Find the Jgup field in the map or create a new one.
// -----------------------------------------------------------------------------
LevelGeometry::t_FCJgupPtr LevelGeometry::createVertAvgFCJgupPtr(const DisjointBoxLayout& a_grids,
                                                                 const Box&               a_fullDomainBox)
{
    CH_TIME("LevelGeometry::createVertAvgFCJgupPtr");

    // Search for the field int the map.
    t_FCJgupPtr& thisFieldPtr = s_FCJgupMap[a_grids];

    // If the field exists in the map, just return it's pointer.
    if (!thisFieldPtr.isNull()) return thisFieldPtr;

    Box dombox = a_grids.physDomain().domainBox();
    D_TERM(CH_assert(dombox.size(SpaceDim-1) == 1);,
           CH_assert(dombox.type(0) == a_fullDomainBox.type(0));,
           CH_assert(dombox.type(1) == a_fullDomainBox.type(1));)
    D_TERM(CH_assert(a_fullDomainBox.size(SpaceDim-1) >= 1);,
           CH_assert(a_fullDomainBox.smallEnd(0) <= dombox.smallEnd(0) && dombox.bigEnd(0) <= a_fullDomainBox.bigEnd(0));,
           CH_assert(a_fullDomainBox.smallEnd(1) <= dombox.smallEnd(1) && dombox.bigEnd(1) <= a_fullDomainBox.bigEnd(1));)

    const IntVect vmask = BASISV(SpaceDim-1);
    const IntVect hmask = IntVect::Unit - vmask;

    // The field did not exist. Allocate and define a new one.
    thisFieldPtr = t_FCJgupPtr(new LevelData<FluxBox>);
    thisFieldPtr->define(a_grids, SpaceDim);// Changed from s_ghostVectFC on Mar 23, 2014

    // Fill the new field.
    DataIterator dit = thisFieldPtr->dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        FluxBox& JgupFlub = (*thisFieldPtr)[dit];

        const Box& flatValid = a_grids[dit];
        const IntVect flatShift = flatValid.smallEnd() * vmask;
        CH_assert(flatValid.size(SpaceDim-1) == 1);

        Box fullValid = flatValid;
        fullValid.shift(-flatShift);
        fullValid.shift(SpaceDim-1, a_fullDomainBox.smallEnd(SpaceDim-1));
        fullValid.setBig(SpaceDim-1, a_fullDomainBox.bigEnd(SpaceDim-1));
        CH_assert(a_fullDomainBox.contains(fullValid));

#ifndef NDEBUG
        JgupFlub.setVal(quietNAN);
#endif
        // This scale will result in a vertical average.
        // const Real scale = m_dXi[SpaceDim-1] / s_domainLength[SpaceDim-1];
        const Real scale = 1.0 / Real(fullValid.size(SpaceDim-1));

        for (int adir = 0; adir < SpaceDim-1; ++adir) {
            Box FCFullValid = surroundingNodes(fullValid, adir);
            FArrayBox fullFAB(FCFullValid, 1);

            for (int bdir = 0; bdir < SpaceDim-1; ++bdir) {
                // Fill the metric over the full domain.
                s_metricSourcePtr->fill_Jgup(fullFAB, 0,
                                             adir, bdir,
                                             m_dXi);

                // Average vertically onto the flat domain.
                Interval interv(bdir, bdir);
                FArrayBox flatJgupComp(interv, JgupFlub[adir]);
                flatJgupComp.setVal(0.0);

                FORT_UNMAPPEDVERTINTEGRAL(
                    CHF_FRA1_SHIFT(flatJgupComp, 0, flatShift),
                    CHF_CONST_FRA1(fullFAB, 0),
                    CHF_BOX(FCFullValid),
                    CHF_CONST_REAL(scale));
            }
        }
    }

    // Return the new field pointer.
    return thisFieldPtr;
}


// -----------------------------------------------------------------------------
// Find the gdn field in the map or create a new one.
// -----------------------------------------------------------------------------
LevelGeometry::t_CCgdnPtr LevelGeometry::createCCgdnPtr(const DisjointBoxLayout& a_grids)
{
    CH_TIME("LevelGeometry::createCCgdnPtr");

    // Search for the field int the map.
    t_CCgdnPtr& thisFieldPtr = s_CCgdnMap[a_grids];

    // If the field exists in the map, just return it's pointer.
    if (!thisFieldPtr.isNull()) return thisFieldPtr;

    Box dombox = a_grids.physDomain().domainBox();

    // The field did not exist. Allocate and define a new one.
    int symComps = (SpaceDim * (SpaceDim + 1)) / 2;
    thisFieldPtr = t_CCgdnPtr(new LevelData<FArrayBox>);
    thisFieldPtr->define(a_grids, symComps, s_ghostVectCC);

    // Fill the new field.
    DataIterator dit = thisFieldPtr->dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        this->fill_gdn((*thisFieldPtr)[dit]);
    }

    // Return the new field pointer.
    return thisFieldPtr;
}


// -----------------------------------------------------------------------------
// Sets the pointer to the coarser LevelGeometry
// -----------------------------------------------------------------------------
void LevelGeometry::setCoarserPtr (const LevelGeometry* a_coarserPtr)
{
    CH_assert(isStaticDefined());

    // If bond is not changing, just leave
    if (a_coarserPtr == m_coarserPtr) return;

    // Swap coarser level ptrs
    const LevelGeometry* oldCoarserPtr = m_coarserPtr;
    m_coarserPtr = a_coarserPtr;

    // Break bond with old coarser level
    if (oldCoarserPtr != NULL) {
        LevelGeometry& coarserRef = (LevelGeometry&)(*oldCoarserPtr);
        coarserRef.setFinerPtr(NULL);
    }

    // Create bond with new coarser level
    if (m_coarserPtr != NULL) {
        LevelGeometry& coarserRef = (LevelGeometry&)(*m_coarserPtr);
        coarserRef.setFinerPtr(this);
    }

    // Set new refRatio
    if (m_coarserPtr != NULL) {
        RealVect q = m_coarserPtr->m_dXi / m_dXi;
        D_TERM(m_crseRefRatio[0] = q[0];,
               m_crseRefRatio[1] = q[1];,
               m_crseRefRatio[2] = q[2];)
        D_TERM(CH_assert(Real(m_crseRefRatio[0]) == q[0]);,
               CH_assert(Real(m_crseRefRatio[1]) == q[1]);,
               CH_assert(Real(m_crseRefRatio[2]) == q[2]);)
    } else {
        m_crseRefRatio = IntVect::Unit;
    }
}


// -----------------------------------------------------------------------------
// Sets the pointer to the finer LevelGeometry
// (usually called by the finer level)
// -----------------------------------------------------------------------------
void LevelGeometry::setFinerPtr (const LevelGeometry* a_finerPtr)
{
    CH_assert(this->isDefined());

    // If bond is not changing, just leave
    if (a_finerPtr == m_finerPtr) return;

    // Swap finer ptrs
    const LevelGeometry* oldFinerPtr = m_finerPtr;
    m_finerPtr = a_finerPtr;

    // Break bond with old finer level
    if (oldFinerPtr != NULL) {
        LevelGeometry& finerRef = (LevelGeometry&)(*oldFinerPtr);
        finerRef.setCoarserPtr(NULL);
    }

    // Create bond with new finer level
    if (m_finerPtr != NULL) {
        LevelGeometry& finerRef = (LevelGeometry&)(*m_finerPtr);
        finerRef.setCoarserPtr(this);
    }

    // Set new refRatio
    if (m_finerPtr != NULL) {
        RealVect q = m_dXi / m_finerPtr->m_dXi;
        D_TERM(m_fineRefRatio[0] = q[0];,
               m_fineRefRatio[1] = q[1];,
               m_fineRefRatio[2] = q[2];)
        D_TERM(CH_assert(Real(m_fineRefRatio[0]) == q[0]);,
               CH_assert(Real(m_fineRefRatio[1]) == q[1]);,
               CH_assert(Real(m_fineRefRatio[2]) == q[2]);)
    } else {
        m_fineRefRatio = IntVect::Unit;
    }
}
