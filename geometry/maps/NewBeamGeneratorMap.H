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


// A total rewrite of the beam topogrophy. This uses interpolation
// to resolve the locations of points interior to cells.
// NOT FOR PUBLIC USE!!! This is a work in progress.


#ifndef __NEWBEAMGENERATORMAP_H__INCLUDED__
#define __NEWBEAMGENERATORMAP_H__INCLUDED__

#include "GeoSourceInterface.H"
#include "Tuple.H"
#include "DisjointBoxLayout.H"
#include "NodeFArrayBox.H"
#include "RefCountedPtr.H"


// -----------------------------------------------------------------------------
// The abstract base class for the metric component functions
// Remember that your children classes should return the physical, not scaled
// metric elements. It is up to LevelGeometry to fix the scaling.
// -----------------------------------------------------------------------------
class NewBeamGeneratorMap: public GeoSourceInterface
{
public:
    // Constructor
    NewBeamGeneratorMap ();

    // Default destructor
    virtual ~NewBeamGeneratorMap ();


    // These need to be overridden ---------------------------------------------

    // 1. Must return the name of the coordinate mapping
    virtual const char* getCoorMapName () const;

    // 2. Must return whether or not this metric is diagonal
    virtual bool isDiagonal () const;

    // 3. Must return whether or not this metric is uniform
    virtual bool isUniform () const;

    // 4. Fills a mapped box with Cartesian locations.
    virtual void fill_physCoor (FArrayBox&      a_dest,
                                const int       a_destComp,
                                const int       a_mu,
                                const RealVect& a_dXi) const;

    // // Fills an FArrayBox with the Jacobian matrix elements d[x^mu] / d[xi^nu].
    // virtual void fill_dxdXi (FArrayBox&      a_dest,
    //                          const int       a_destComp,
    //                          const int       a_mu,
    //                          const int       a_nu,
    //                          const RealVect& a_dXi,
    //                          const Real      a_scale = 1.0) const;


    // These can be overridden if the analytic functions are known -------------

    // // Fills a mapped box with Cartesian locations (a_dest must have SpaceDim comps)
    // virtual void fill_physCoor (FArrayBox&      a_dest,
    //                             const RealVect& a_dXi,
    //                             const RealVect  a_scale = RealVect::Unit) const;

    // // Fills an FArrayBox with the Jacobian matrix elements d[x^mu] / d[xi^nu].
    // // This is a speed bottleneck!!!
    // virtual void fill_dxdXi (FArrayBox&      a_dest,
    //                          const int       a_destComp,
    //                          const int       a_mu,
    //                          const int       a_nu,
    //                          const RealVect& a_dXi,
    //                          const Real      a_scale = 1.0) const;

    // // Fills an FArrayBox with J = det[Jacobian]
    // // This is a speed bottleneck!!!
    // // TODO: try calling dxdXi and using its results. This way, the speed
    // // bottleneck may be circumvented by rewriting one function instead of two.
    // virtual void fill_J (FArrayBox&      a_dest,
    //                      const int       a_destComp,
    //                      const RealVect& a_dXi,
    //                      const Real      a_scale = 1.0) const;

    // // Fills an FArrayBox with the inverse Jacobian matrix elements d[xi^mu] / d[x^nu].
    // virtual void fill_dXidx (FArrayBox&       a_dest,
    //                          const int        a_destComp,
    //                          const int        a_mu,
    //                          const int        a_nu,
    //                          const RealVect&  a_dXi,
    //                          const Real       a_scale = 1.0) const;

    // // Fills an FArrayBox with 1/J
    // virtual void fill_Jinv (FArrayBox&       a_dest,
    //                         const int        a_destComp,
    //                         const RealVect&  a_dXi,
    //                         const Real       a_scale = 1.0) const;

    // // Fills an FArrayBox with the covariant metric elements
    // virtual void fill_gdn (FArrayBox&      a_dest,
    //                        const int       a_destComp,
    //                        const int       a_mu,
    //                        const int       a_nu,
    //                        const RealVect& a_dXi,
    //                        const Real      a_scale = 1.0) const;

    // // Fills an FArrayBox with the contravariant metric elements
    // virtual void fill_gup (FArrayBox&       a_dest,
    //                        const int        a_destComp,
    //                        const int        a_mu,
    //                        const int        a_nu,
    //                        const RealVect&  a_dXi,
    //                        const Real       a_scale = 1.0) const;

    // // Fills an FArrayBox with J * gup
    // virtual void fill_Jgup (FArrayBox&       a_dest,
    //                         const int        a_destComp,
    //                         const int        a_mu,
    //                         const int        a_nu,
    //                         const RealVect&  a_dXi,
    //                         const Real       a_scale = 1.0) const;

    // // Fills an FArrayBox with the connection elements
    // virtual void fill_Gamma (FArrayBox&      a_dest,
    //                          const int       a_destComp,
    //                          const int       a_up,
    //                          const int       a_dn1,
    //                          const int       a_dn2,
    //                          const RealVect& a_dXi,
    //                          const Real      a_scale = 1.0) const;

    // This can be used to speed up cache access.
    virtual void suggestLev0Grids (const DisjointBoxLayout& a_grids);

protected:
    // Fills a mapped box with Cartesian locations.
    // Assumes we are on level 0.
    virtual void fill_physCoorLev0 (FArrayBox&      a_dest,
                                    const int       a_destComp,
                                    const int       a_mu) const;

    // Fills a mapped box with Cartesian locations.
    // Assumes we are on level 0 and a_dest is NC.
    // Stretching happens here.
    virtual void fill_physCoorLev0Nodes (FArrayBox& a_dest,
                                         const int  a_destComp,
                                         const int  a_mu) const;

    // Returns the topographic elevations at each Cartesian node.
    // Assumes we are on level 0 and a_dest is NC.
    virtual void fill_elevationLev0 (FArrayBox&       a_dest,
                                     const int        a_destComp,
                                     const FArrayBox& a_physCoor) const;

    // Static utility
    // Computes a refRatio. This will also check if crseDXi >= fineDXi.
    static IntVect computeRefToFiner (const RealVect& a_crseDXi,
                                      const RealVect& a_fineDXi);

    // Static utility
    // Creates lev 0 grids. Does not use a layout suggestion.
    static BoxLayout createLev0Grids ();

    // Static utility
    // Fills the cubic spline cache.
    static void computeSplineData (const NewBeamGeneratorMap* a_srcPtr);

    // Copies the Cartesian locations from the cache.
    // Assumes we are on level 0 and a_dest is NC.
    static void copyLev0NodeCacheData (FArrayBox&                          a_dest,
                                       const int                           a_destComp,
                                       const BoxLayoutData<NodeFArrayBox>& a_srcData,
                                       const int                           a_srcComp);

    // If the cache data is local, we reference it to avoid MPI+copy. Returns
    // true if a_dest is a reference, false if we needed to perform the MPI+copy.
    // a_dest must enter this function undefined.
    static bool refOrCopy (FArrayBox&                    a_dest,
                           const Box&                    a_destNodeBox,
                           BoxLayoutData<NodeFArrayBox>& a_srcData,
                           const int                     a_srcComp);


    static bool          s_staticParamsSet;

    static RealVect      s_L;
    static IntVect       s_lev0Nx;
    static RealVect      s_lev0DXi;
    static ProblemDomain s_lev0Domain;

    static RefCountedPtr<BoxLayoutData<NodeFArrayBox> > s_lev0xPtr;
    static RefCountedPtr<BoxLayoutData<NodeFArrayBox> > s_lev0d2xPtr;

    static Real          s_alpha;

    static const IntVect s_ghostVect;

    // Tuple<DisjointBoxLayout, CH_SPACEDIM> m_lev0Grids;
};


#endif //!__NEWBEAMGENERATORMAP_H__INCLUDED__

