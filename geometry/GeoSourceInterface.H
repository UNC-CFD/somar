/*******************************************************************************
 *    FILE: GeoSourceInterface.H
 *    DESC: The abstract base class for the metric component functions.
 *    DATE: Sat 12 Nov 2011 02:27:54 PM EST
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __GEOSOURCEINTERFACE_H__INCLUDED__
#define __GEOSOURCEINTERFACE_H__INCLUDED__

// Chombo includes
#include "GNUC_Extensions.H"
#include "LevelData.H"
#include "RealVect.H"
class FArrayBox;
class FluxBox;


// -----------------------------------------------------------------------------
// Provides a unified interface for functions that act as a substitute for
// GeoSourceInterface::fill_Jgup.
// -----------------------------------------------------------------------------
class FillJgupInterface
{
public:
    virtual void fill_Jgup (FArrayBox&       a_dest,
                            const int        a_destComp,
                            const int        a_mu,
                            const int        a_nu,
                            const RealVect&  a_dXi,
                            const Real       a_scale = 1.0) const = 0;
};


// -----------------------------------------------------------------------------
// The abstract base class for the metric component functions
// Remember that your children classes should return the physical, not scaled
// metric elements. It is up to LevelGeometry to fix the scaling.
// -----------------------------------------------------------------------------
class GeoSourceInterface: public FillJgupInterface
{
public:
    // Default destructor
    virtual ~GeoSourceInterface ();


    // These need to be overridden ---------------------------------------------

    // 1. Must return the name of the coordinate mapping
    virtual const char* getCoorMapName () const = 0;

    // 2. Must return whether or not this metric is diagonal
    virtual bool isDiagonal () const = 0;

    // 3. Must return whether or not this metric is uniform
    virtual bool isUniform () const = 0;

    // 4. Fills a mapped box with Cartesian locations.
    virtual void fill_physCoor (FArrayBox&      a_dest,
                                const int       a_destComp,
                                const int       a_mu,
                                const RealVect& a_dXi) const = 0;


    // These can be overridden if the analytic functions are known -------------

    // Fills a mapped box with Cartesian locations (a_dest must have SpaceDim comps)
    virtual void fill_physCoor (FArrayBox&      a_dest,
                                const RealVect& a_dXi,
                                const RealVect  a_scale = RealVect::Unit) const;

    // Fills an FArrayBox with the Jacobian matrix elements d[x^mu] / d[xi^nu].
    // This is a speed bottleneck!!!
    virtual void fill_dxdXi (FArrayBox&      a_dest,
                             const int       a_destComp,
                             const int       a_mu,
                             const int       a_nu,
                             const RealVect& a_dXi,
                             const Real      a_scale = 1.0) const;

    // Fills an FArrayBox with J = det[Jacobian]
    // This is a speed bottleneck!!!
    // TODO: try calling dxdXi and using its results. This way, the speed
    // bottleneck may be circumvented by rewriting one function instead of two.
    virtual void fill_J (FArrayBox&      a_dest,
                         const int       a_destComp,
                         const RealVect& a_dXi,
                         const Real      a_scale = 1.0) const;

    // Fills an FArrayBox with the inverse Jacobian matrix elements d[xi^mu] / d[x^nu].
    virtual void fill_dXidx (FArrayBox&       a_dest,
                             const int        a_destComp,
                             const int        a_mu,
                             const int        a_nu,
                             const RealVect&  a_dXi,
                             const Real       a_scale = 1.0) const __hot;

    // Fills an FArrayBox with 1/J
    virtual void fill_Jinv (FArrayBox&       a_dest,
                            const int        a_destComp,
                            const RealVect&  a_dXi,
                            const Real       a_scale = 1.0) const;

    // Fills an FArrayBox with the covariant metric elements
    virtual void fill_gdn (FArrayBox&      a_dest,
                           const int       a_destComp,
                           const int       a_mu,
                           const int       a_nu,
                           const RealVect& a_dXi,
                           const Real      a_scale = 1.0) const;

    // Fills an FArrayBox with the contravariant metric elements
    virtual void fill_gup (FArrayBox&       a_dest,
                           const int        a_destComp,
                           const int        a_mu,
                           const int        a_nu,
                           const RealVect&  a_dXi,
                           const Real       a_scale = 1.0) const;

    // Fills an FArrayBox with J * gup
    virtual void fill_Jgup (FArrayBox&       a_dest,
                            const int        a_destComp,
                            const int        a_mu,
                            const int        a_nu,
                            const RealVect&  a_dXi,
                            const Real       a_scale = 1.0) const;

    // Fills an FArrayBox with the connection elements.
    // Note that although this function is provided and works, it is never used
    // by our NS algorithm.
    virtual void fill_Gamma (FArrayBox&      a_dest,
                             const int       a_destComp,
                             const int       a_up,
                             const int       a_dn1,
                             const int       a_dn2,
                             const RealVect& a_dXi,
                             const Real      a_scale = 1.0) const;


    // This can be used to speed up cache access. (If you use a cache.)
    virtual void suggestLev0Grids (const DisjointBoxLayout& a_grids)
    {;}
};


#endif //!__GEOSOURCEINTERFACE_H__INCLUDED__

