#ifdef CH_LANG_CC
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
#endif // CH_LANG_CC

#ifndef __LEVELGEOMETRY_HPP__INCLUDED__
#define __LEVELGEOMETRY_HPP__INCLUDED__

#ifdef CH_LANG_FORT
#   if CH_SPACEDIM == 2
#       define xxComp 0
#       define xyComp 1
#       define yxComp 2
#       define yyComp 3
#       define xxSymComp 0
#       define xySymComp 1
#       define yxSymComp 1
#       define yySymComp 2
#   elif CH_SPACEDIM == 3
#       define xxComp 0
#       define xyComp 1
#       define xzComp 2
#       define yxComp 3
#       define yyComp 4
#       define yzComp 5
#       define zxComp 6
#       define zyComp 7
#       define zzComp 8
#       define xxSymComp 0
#       define xySymComp 1
#       define xzSymComp 2
#       define yxSymComp 1
#       define yySymComp 3
#       define yzSymComp 4
#       define zxSymComp 2
#       define zySymComp 4
#       define zzSymComp 5
#   else
#       error Bad CH_SPACEDIM
#   endif
#endif //CH_LANG_FORT

#ifdef CH_LANG_CC
// The rest of this file contains C++ code

#include <map>

#include "LevelData.H"
#include "FArrayBox.H"
#include "NodeFArrayBox.H"
#include "FluxBox.H"
#include "RealVect.H"
#include "Vector.H"
#include "ProblemContext.H"
class DisjointBoxLayout;
class ProblemDomain;
class LevelFluxRegister;

#include "GeoSourceInterface.H"


// -----------------------------------------------------------------------------
// A container for the metric, connection, Jacobian, etc...
// This class completely describes the geometry of the problem's domain.
// -----------------------------------------------------------------------------
class LevelGeometry
{
public:
    // Construction / destruction ----------------------------------------------

    // Weak constructor
    LevelGeometry ();

    // Full constructor (calls define)
    LevelGeometry (const RealVect& a_dXi);

    // Full constructor
    void define (const RealVect& a_dXi);

    // Puts this Geometry object in a usable state. Called by the full constructor.
    // This version links this with the coarser level, which sets the refRatios.
    // The finer pointer will still be NULL. It is expected that these defines will
    // be called from coarsest to finest.
    void define (const RealVect&      a_dXi,
                 const LevelGeometry* a_coarserPtr);

    // Destructor
    ~LevelGeometry ();

    // Regrids the current geometry data
    void regrid (const DisjointBoxLayout& a_newGrids);

    // Regrids the current geometry data with flat (vertically averaged) metrics.
    void regridVertAvg (const DisjointBoxLayout& a_newGrids,
                        const Box&               a_fullDomainBox);

    // Clears all data holders and grids. Sets coarser/finer ptrs to NULL.
    void reset ();

    // Basic accessors ---------------------------------------------------------

    // Tests if this LevelGeometry is in a useable state
    inline bool isDefined () const;

    // Returns the ProblemDomain
    inline const ProblemDomain& getDomain () const;

    // Returns the current DisjointBoxLayout
    inline const DisjointBoxLayout& getBoxes () const;

    // Sets the pointer to the coarser LevelGeometry
    void setCoarserPtr (const LevelGeometry* a_coarserPtr);

    // Sets the pointer to the finer LevelGeometry (usually called by the finer level)
    void setFinerPtr (const LevelGeometry* a_finerPtr);

    // Return a pointer to the coarser LevelGeometry
    const LevelGeometry* getCoarserPtr () const;

    // Set/Return a pointer to the finer LevelGeometry
    const LevelGeometry* getFinerPtr () const;

    // Returns all levGeos over the AMR hierarchy
    Vector<const LevelGeometry*> getAMRLevGeos () const;

    // Returns all fineRefRatios over the AMR hierarchy
    Vector<IntVect> getAMRRefRatios () const;

    // Returns all DisjointBoxLayouts over the AMR hierarchy
    Vector<DisjointBoxLayout> getAMRGrids () const;

    // Returns the refinement ratio to the next coarser level
    inline const IntVect& getCrseRefRatio () const;

    // Returns the refinement ratio to the next coarser level
    inline int getCrseRefRatio (const int a_dir) const;

    // Returns the refinement ratio to the next finer level
    inline const IntVect& getFineRefRatio () const;

    // Returns the refinement ratio to the next finer level
    inline int getFineRefRatio (const int a_dir) const;

    // Returns this level's grid spacing
    inline const RealVect& getDx () const;

    // Returns this level's grid spacing
    inline Real getDx (const int a_dir) const;

    // Returns the current CC J field
    inline const LevelData<FArrayBox>& getCCJ () const;
    inline const RefCountedPtr<LevelData<FArrayBox> >& getCCJPtr () const;

    // Returns the current CC 1/J field
    inline const LevelData<FArrayBox>& getCCJinv () const;
    inline const RefCountedPtr<LevelData<FArrayBox> >& getCCJinvPtr () const;

    // Returns the current FC contravariant metric tensor
    inline const LevelData<FluxBox>& getFCgup () const;
    inline const RefCountedPtr<LevelData<FluxBox> >& getFCgupPtr () const;

    // Returns the current FC contravariant metric tensor * J
    inline const LevelData<FluxBox>& getFCJgup () const;
    inline const RefCountedPtr<LevelData<FluxBox> >& getFCJgupPtr () const;

    // Returns the current CC covariant metric tensor
    inline const LevelData<FArrayBox>& getCCgdn () const;
    inline const RefCountedPtr<LevelData<FArrayBox> >& getCCgdnPtr () const;


    // Fill functions ----------------------------------------------------------

    // Fills a FAB with displacements from Xi to physical locations.
    // For use with VisIt's displace operator.
    void fill_displacement (FArrayBox& a_dest) const;

    // Fills a mapped box with physical locations
    // a_dest must have SpaceDim comps.
    void fill_physCoor (FArrayBox& a_dest) const;

    // Fill a FAB with a Jacobian element (mapped to Cartesian)
    void fill_dxdXi (FArrayBox& a_dest) const;

    // Fill a FAB with a inverse Jacobian element (Cartesian to mapped)
    void fill_dXidx (FArrayBox& a_dest) const;

    // Fills a FAB with J
    void fill_J (FArrayBox& a_dest) const;

    // Fills a FluxBox with J
    void fill_J (FluxBox& a_dest) const;

    // Fills an FArrayBox with 1/J
    void fill_Jinv (FArrayBox& a_dest) const;

    // Fills an FArrayBox with the covariant metric elements (static version)
    static void fill_gdn (FArrayBox& a_dest, const RealVect& a_dXi);

    // Fills an FArrayBox with the covariant metric elements
    void fill_gdn (FArrayBox& a_dest) const;

    // Fills an FArrayBox with the contravariant metric elements (static version)
    static void fill_gup (FArrayBox& a_dest, const RealVect& a_dXi);

    // Fills an FArrayBox with the contravariant metric elements
    void fill_gup (FArrayBox& a_dest) const;

    // Fills a FluxBox with the contravariant metric elements
    void fill_gup (FluxBox& a_dest) const;

    // Fills an FArrayBox with just one contravariant metric element
    void fill_gup (FArrayBox& a_dest, int a_1, int a_2) const;

    // Fills a FluxBox with J * the contravariant metric elements
    void fill_Jgup (FluxBox& a_dest) const;

    // Fills an FArrayBox with Jgup^{a_mu,*}
    // This is useful for computing transverse derivatives.
    // NOTE: The FArrayBox can have any centering, but must have SpaceDim comps. The
    // SpaceDim comps of Jgup^{a_mu,*} will fill the SpaceDim comps of a_dest.
    void fill_Jgup (FArrayBox& a_dest, int a_mu) const;

    // Fills a single-component FArrayBox with Gamma^{a_1}_{a_2, a_3}.
    // This is scaled as h^{a_2} h^{a_3} / h^{a_1} where h = physDxCoarse.
    void fill_Gamma (FArrayBox& a_dest, int a_1, int a_2, int a_3) const;

    // Ensures that metric data is the average of finer data.
    void averageMetricsDown ();

#ifndef NDEBUG
    // Debugging functions -----------------------------------------------------

    // Dumps debugging info to pout()
    void dump () const;
#endif


    // Utilities ---------------------------------------------------------------

    // Sends a vector from the scaled, mapped basis to an unscaled, Cartesian basis.
    // uCart^{a} = [dx^{a} / dXi^{b}] * uMapped^{b}
    void sendToCartesianBasis (LevelData<FArrayBox>& a_CCvect, bool a_doGhosts = true) const;
    void sendToCartesianBasis (LevelData<FluxBox>& a_FCvect, bool a_doGhosts = true) const;

    // Sends a vector from the unscaled, Cartesian basis to an scaled, mapped basis.
    // uMapped^{a} = [dXi^{a} / dx^{b}] * uCart^{b}
    void sendToMappedBasis (LevelData<FArrayBox>& a_CCvect, bool a_doGhosts = true) const;
    void sendToMappedBasis (LevelData<FluxBox>& a_FCvect, bool a_doGhosts = true) const;

    // Multiplies every element by J
    void multByJ (Vector<LevelData<FArrayBox>*>& a_data) const;
    void multByJ (LevelData<FArrayBox>& a_data) const;
    void multByJ (FArrayBox& a_data, const DataIndex& a_di, int a_comp = -1) const;
    void multByJ (Vector<LevelData<FluxBox>*>& a_data) const;
    void multByJ (LevelData<FluxBox>& a_data) const;
    void multByJ (FluxBox& a_data, const DataIndex& a_di, int a_comp = -1) const;

    // Divides every element by J
    void divByJ (Vector<LevelData<FArrayBox>*>& a_data) const;
    void divByJ (LevelData<FArrayBox>& a_data) const;
    void divByJ (FArrayBox& a_data, const DataIndex& a_di, int a_comp = -1) const;
    void divByJ (Vector<LevelData<FluxBox>*>& a_data) const;
    void divByJ (LevelData<FluxBox>& a_data) const;
    void divByJ (FluxBox& a_data, const DataIndex& a_di, int a_comp = -1) const;

    // Contracts a CC contra-vector with gdn, making it covariant.
    // Single grid version.
    void makeCovariant (FArrayBox&       a_coVect,
                        const FArrayBox& a_contraVect,
                        const Box&       a_region,
                        const DataIndex& a_di) const;

    // Contracts a CC covariant one-form with gup, making it a contravariant vector.
    // Static, single grid version.
    static void makeContravariant (FArrayBox&       a_contraVect,
                                   const FArrayBox& a_coVect,
                                   const Box&       a_region,
                                   const DataIndex& a_di,
                                   const RealVect&  a_dXi);

    // Contracts a CC covariant one-form with gup, making it a contravariant vector.
    // Single grid version.
    void makeContravariant (FArrayBox&       a_contraVect,
                            const FArrayBox& a_coVect,
                            const Box&       a_region,
                            const DataIndex& a_di) const;

    // Compute the magnitude of a vector
    void contractVectors (FArrayBox&       a_mag,
                          const FArrayBox& a_vec1,
                          const FArrayBox& a_vec2,
                          const DataIndex& a_di);

    // Compute the mu-th comp of the vector-scaled Gamma^{mu}_{a b} vec1^{a} vec2^{b}
    void contractWithGamma (FArrayBox&       a_res,
                            const FArrayBox& a_vec1,
                            const FArrayBox& a_vec2,
                            int              a_mu,
                            const Box&       a_region,
                            const DataIndex& a_di) const;


    // Static utilities and accessors ------------------------------------------

    // Puts this geometry object in a usable state. Called by the full constructor.
    static void staticDefine ();

    // Tests if the static members are have been set
    static inline bool isStaticDefined ();

    // Clear all memory used by the static maps. Calling this before the program
    // ends makes for easier memory leak diagnosis.
    static void staticUndefine ();

    // Returns a pointer to the GeoSourceInterface. For use in MG functions that
    // need to fill data with a manually set dXi.
    static inline RefCountedPtr<GeoSourceInterface> getGeoSourcePtr();

    // Returns the coordinate mapping type
    static inline int getCoordMap ();

    // Returns current coordinate map name
    static inline const char* getCoorMapName ();

    // Returns if this metric is diagonal
    static inline bool isDiagonal ();

    // Returns if this metric is constant in space
    static inline bool isUniform ();

    // Returns the domain's physical extents
    static inline const RealVect& getDomainLength ();
    static inline Real getDomainLength (int a_dir);

    // Utility for converting CC tensor indices to a FAB component
    static inline int tensorCompCC (int, int);
    static inline int symTensorCompCC (int, int);
    static inline int symTensorCompCC (int, int, int);


protected:
    // Returns all levGeos over the AMR hierarchy
    Vector<LevelGeometry*> getAMRLevGeos ();

    // Static members ----------------------------------------------------------
    static RefCountedPtr<GeoSourceInterface>   s_metricSourcePtr;   // The metric source functions
    static int                                 s_coordMap;          // The coordinate mapping
    static RealVect                            s_domainLength;      // The domain's physical extents
    static const IntVect                       s_ghostVectFC;       // Used to define the cached LevelDatas
    static const IntVect                       s_ghostVectCC;       // Used to define the cached LevelDatas

    static RealVect                            s_lev0dXi;           // dXi at the base level.
    static LevelData<NodeFArrayBox>*           s_lev0xPtr;          // physCoor on the base level.
    static LevelData<NodeFArrayBox>*           s_lev0d2xPtr;        // The spline's 2nd derivatives.

    static const int m_tensorCompCC2[SpaceDim][SpaceDim];       // Used by the tensorComp functions
    static const int m_symTensorCompCC2[SpaceDim][SpaceDim];
    static const int m_symTensorCompCC3[SpaceDim][SpaceDim][SpaceDim];


    // Non-static members ------------------------------------------------------
    RealVect             m_dXi;          // This level's grid spacing
    IntVect              m_refToLev0;    // Coarsen by this amount to arrive at level 0's domain.
    DisjointBoxLayout    m_grids;        // This level's grids
    IntVect              m_crseRefRatio; // Refinement ratio to the next coarser level
    IntVect              m_fineRefRatio; // Refinement ratio to the next finer level
    const LevelGeometry* m_coarserPtr;   // Pointer to coarser LevelGeometry
    const LevelGeometry* m_finerPtr;     // Pointer to finer LevelGeometry

    // The CC J field
    typedef RefCountedPtr<LevelData<FArrayBox> >    t_CCJPtr;
    typedef map<const BoxLayout, t_CCJPtr>          t_CCJMap;
    static t_CCJMap                                 s_CCJMap;
    t_CCJPtr createCCJPtr(const DisjointBoxLayout& a_grids);

    // The CC 1/J field
    typedef RefCountedPtr<LevelData<FArrayBox> >    t_CCJinvPtr;
    typedef map<const BoxLayout, t_CCJinvPtr>       t_CCJinvMap;
    static t_CCJinvMap                              s_CCJinvMap;
    t_CCJinvPtr createCCJinvPtr(const DisjointBoxLayout& a_grids);
    t_CCJinvPtr createVertAvgCCJinvPtr(const DisjointBoxLayout& a_grids,
                                       const Box&               a_fullDomainBox);

    // The FC contravariant metric tensor
    typedef RefCountedPtr<LevelData<FluxBox> >      t_FCgupPtr;
    typedef map<const BoxLayout, t_FCgupPtr>        t_FCgupMap;
    static t_FCgupMap                               s_FCgupMap;
    t_FCgupPtr createFCgupPtr(const DisjointBoxLayout& a_grids);

    // The FC contravariant metric tensor * J
    typedef RefCountedPtr<LevelData<FluxBox> >      t_FCJgupPtr;
    typedef map<const BoxLayout, t_FCJgupPtr>       t_FCJgupMap;
    static t_FCJgupMap                              s_FCJgupMap;
    t_FCJgupPtr createFCJgupPtr(const DisjointBoxLayout& a_grids);
    t_FCJgupPtr createVertAvgFCJgupPtr(const DisjointBoxLayout& a_grids,
                                       const Box&               a_fullDomainBox);

    // The CC covariant metric tensor
    typedef RefCountedPtr<LevelData<FArrayBox> >    t_CCgdnPtr;
    typedef map<const BoxLayout, t_CCgdnPtr>        t_CCgdnMap;
    static t_CCgdnMap                               s_CCgdnMap;
    t_CCgdnPtr createCCgdnPtr(const DisjointBoxLayout& a_grids);

    // Keep pointers to this object's fields for easier access.
    RefCountedPtr<LevelData<FArrayBox> > m_CCJPtr;             // The CC J field
    RefCountedPtr<LevelData<FArrayBox> > m_CCJinvPtr;          // The CC 1/J field
    RefCountedPtr<LevelData<FluxBox> >   m_FCgupPtr;           // The FC contravariant metric tensor
    RefCountedPtr<LevelData<FluxBox> >   m_FCJgupPtr;          // The FC contravariant metric tensor * J
    RefCountedPtr<LevelData<FArrayBox> > m_CCgdnPtr;           // The CC covariant metric tensor

    // Caches the unit normal vectors at the boundary.
    // This is a set of maps: dXi -> faceBox -> data.
    // This cache works a bit differently than the metric caches.
    // It is completely flushed upon reset() or regrid().
    typedef RefCountedPtr<FArrayBox>                t_bdryNormData;
    typedef map<Box, t_bdryNormData>                t_bdryBox2DataMap;
    typedef map<RealVect, t_bdryBox2DataMap>        t_bdryNormMaps;
    static t_bdryNormMaps                           s_bdryNormMaps;

private:
    // Assignment operator is not allowed
    void operator= (const LevelGeometry&)
    {
        MayDay::Error ("LevelGeometry assignment operator is not allowed");
    }

    // Copy constructor not allowed
    LevelGeometry (const LevelGeometry&)
    {
        MayDay::Error ("LevelGeometry copy constructor is not allowed");
    }
};



// -----------------------------------------------------------------------------
// Implementation of inline functions
// -----------------------------------------------------------------------------

// Utility for converting tensor indices to a FAB component.
// Mimics TensorCFInterp::gradIndex(i,graddir) which returns i*SpaceDim+graddir.
// mode. The assert is good enough to tell you where things went wrong.
int LevelGeometry::tensorCompCC (int a_1, int a_2) {
    CH_assert((0 <= a_1) && (a_1 < SpaceDim));
    CH_assert((0 <= a_2) && (a_2 < SpaceDim));
    return m_tensorCompCC2[a_1][a_2];
}


// Utility for converting symmetric tensor indices to a FAB component.
// mode. The assert is good enough to tell you where things went wrong.
int LevelGeometry::symTensorCompCC (int a_1, int a_2) {
    CH_assert((0 <= a_1) && (a_1 < SpaceDim));
    CH_assert((0 <= a_2) && (a_2 < SpaceDim));
    return m_symTensorCompCC2[a_1][a_2];
}


// Utility for converting tensor indices to a FAB component.
// Tensor is symmetric in the last two indices.
// mode. The assert is good enough to tell you where things went wrong.
int LevelGeometry::symTensorCompCC (int a_1, int a_2, int a_3) {
    CH_assert((0 <= a_1) && (a_1 < SpaceDim));
    CH_assert((0 <= a_2) && (a_2 < SpaceDim));
    CH_assert((0 <= a_3) && (a_3 < SpaceDim));
    return m_symTensorCompCC3[a_1][a_2][a_3];
}


// Tests if the static members are have been set
bool LevelGeometry::isStaticDefined ()
{
    return ( (s_coordMap != ProblemContext::CoordMap::UNDEFINED) &&
             (!s_metricSourcePtr.isNull()) );
}


// Returns a pointer to the GeoSourceInterface. For use in MG functions that
// need to fill data with a manually set dXi.
RefCountedPtr<GeoSourceInterface> LevelGeometry::getGeoSourcePtr ()
{
    CH_assert(isStaticDefined());
    return s_metricSourcePtr;
}


// Tests if this LevelGeometry is in a useable state
bool LevelGeometry::isDefined () const
{
    CH_assert(isStaticDefined());
    return D_TERM(   m_dXi[0] > 0.0,
                  && m_dXi[1] > 0.0,
                  && m_dXi[2] > 0.0   );
}


// Returns the coordinate mapping type
int LevelGeometry::getCoordMap ()
{
    CH_assert(isStaticDefined());
    return s_coordMap;
}


// Returns the coordinate map's name
const char* LevelGeometry::getCoorMapName ()
{
    CH_assert(isStaticDefined());
    return s_metricSourcePtr->getCoorMapName();
}


// Returns whether or not this metric is diagonal
bool LevelGeometry::isDiagonal ()
{
    CH_assert(isStaticDefined());
    return s_metricSourcePtr->isDiagonal();
}


// Returns whether or not this metric is uniform
bool LevelGeometry::isUniform ()
{
    CH_assert(isStaticDefined());
    return s_metricSourcePtr->isUniform();
}


// Returns the current ProblemDomain
const ProblemDomain& LevelGeometry::getDomain () const
{
    CH_assert(this->isDefined());
    return m_grids.physDomain();
}


// Returns the current DisjointBoxLayout
const DisjointBoxLayout& LevelGeometry::getBoxes () const
{
    CH_assert(this->isDefined());
    return m_grids;
}


// Returns the refinement ratio to the next coarser level
const IntVect& LevelGeometry::getCrseRefRatio () const
{
    CH_assert(this->isDefined());
    return m_crseRefRatio;
}


// Returns the refinement ratio to the next coarser level
int LevelGeometry::getCrseRefRatio (const int a_dir) const
{
    CH_assert(this->isDefined());
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    return m_crseRefRatio[a_dir];
}


// Returns the refinement ratio to the next finer level
const IntVect& LevelGeometry::getFineRefRatio () const
{
    CH_assert(this->isDefined());
    return m_fineRefRatio;
}


// Returns the refinement ratio to the next finer level
int LevelGeometry::getFineRefRatio (const int a_dir) const
{
    CH_assert(this->isDefined());
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    return m_fineRefRatio[a_dir];
}


// Returns the domain's physical extents
const RealVect& LevelGeometry::getDomainLength ()
{
    CH_assert(isStaticDefined());
    return s_domainLength;
}


// Returns the domain's physical extents in one direction
Real LevelGeometry::getDomainLength (int a_dir)
{
    CH_assert(isStaticDefined());
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    return s_domainLength[a_dir];
}


// Returns this level's grid spacing
const RealVect& LevelGeometry::getDx () const
{
    CH_assert(this->isDefined());
    return m_dXi;
}


// Returns this level's grid spacing
Real LevelGeometry::getDx (const int a_dir) const
{
    CH_assert(this->isDefined());
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    return m_dXi[a_dir];
}


// Returns the current CC J field
const LevelData<FArrayBox>& LevelGeometry::getCCJ () const {
    CH_assert(!m_CCJPtr.isNull());
    return *m_CCJPtr;
}


// Returns the current CC J field pointer
const RefCountedPtr<LevelData<FArrayBox> >& LevelGeometry::getCCJPtr () const {
    return m_CCJPtr;
}


// Returns the current CC 1/J field
const LevelData<FArrayBox>& LevelGeometry::getCCJinv () const {
    CH_assert(!m_CCJinvPtr.isNull());
    return *m_CCJinvPtr;
}


// Returns the current CC 1/J field pointer
const RefCountedPtr<LevelData<FArrayBox> >& LevelGeometry::getCCJinvPtr () const {
    return m_CCJinvPtr;
}


// Returns the current FC contravariant metric tensor
const LevelData<FluxBox>& LevelGeometry::getFCgup () const {
    CH_assert(false); // is this ever called?
    CH_assert(!m_FCgupPtr.isNull());
    return *m_FCgupPtr;
}


// Returns the current FC contravariant metric tensor pointer
const RefCountedPtr<LevelData<FluxBox> >& LevelGeometry::getFCgupPtr () const {
    CH_assert(false); // is this ever called?
    return m_FCgupPtr;
}


// Returns the current FC contravariant metric tensor * J
const LevelData<FluxBox>& LevelGeometry::getFCJgup () const {
    CH_assert(!m_FCJgupPtr.isNull());
    return *m_FCJgupPtr;
}


// Returns the current FC contravariant metric tensor * J pointer
const RefCountedPtr<LevelData<FluxBox> >& LevelGeometry::getFCJgupPtr () const {
    return m_FCJgupPtr;
}


// Returns the current CC covariant metric tensor
const LevelData<FArrayBox>& LevelGeometry::getCCgdn () const {
    CH_assert(!m_CCgdnPtr.isNull());
    return *m_CCgdnPtr;
}


// Returns the current CC covariant metric tensor pointer
const RefCountedPtr<LevelData<FArrayBox> >& LevelGeometry::getCCgdnPtr () const {
    return m_CCgdnPtr;
}


#endif //CH_LANG_CC
#endif //!__LEVELGEOMETRY_HPP__INCLUDED__
