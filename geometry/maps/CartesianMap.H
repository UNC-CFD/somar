/*******************************************************************************
 *    FILE: CartesianMap.hpp
 *    DESC: The Cartestian coordinate system's metric source functions.
 *    DATE: Sat 12 Nov 2011 03:33:00 PM EST (Overhauled on Sun 13 Jan 2013)
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __CARTESIANMAP_HPP__INCLUDED__
#define __CARTESIANMAP_HPP__INCLUDED__

// Project includes
#include "GeoSourceInterface.H"

#include "UsingNamespace.H"


// -----------------------------------------------------------------------------
// The Cartestian coordinate system's metric source functions.
// -----------------------------------------------------------------------------
class CartesianMap: public GeoSourceInterface
{
public:
    // Destructor
    virtual ~CartesianMap ();


    // Mandatory overrides...

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


    // Optional overrides...

    // Fills a mapped box with Cartesian locations (a_dest must have SpaceDim comps)
    virtual void fill_physCoor (FArrayBox&      a_dest,
                                const RealVect& a_dXi,
                                const RealVect  a_scale = RealVect::Unit) const;

    // Fills an FArrayBox with the Jacobian matrix elements d[x^mu] / d[xi^nu].
    // Major speedup!!!
    virtual void fill_dxdXi (FArrayBox&      a_dest,
                             const int       a_destComp,
                             const int       a_mu,
                             const int       a_nu,
                             const RealVect& a_dXi,
                             const Real      a_scale = 1.0) const;

    // Fills an FArrayBox with J = det[Jacobian]
    // Major speedup!!!
    // Without the analytic version of this function, machine-epsilon sized
    // noise is introduced.
    virtual void fill_J (FArrayBox&      a_dest,
                         const int       a_destComp,
                         const RealVect& a_dXi,
                         const Real      a_scale = 1.0) const;

    // Fills an FArrayBox with the inverse Jacobian matrix elements d[xi^mu] / d[x^nu].
    virtual void fill_dXidx (FArrayBox&      a_dest,
                             const int       a_destComp,
                             const int       a_mu,
                             const int       a_nu,
                             const RealVect& a_dXi,
                             const Real      a_scale = 1.0) const;

    // Fills an FArrayBox with 1/J
    virtual void fill_Jinv (FArrayBox&      a_dest,
                            const int       a_destComp,
                            const RealVect& a_dXi,
                            const Real      a_scale = 1.0) const;

    // Fills an FArrayBox with the covariant metric elements
    virtual void fill_gdn (FArrayBox&      a_dest,
                           const int       a_destComp,
                           const int       a_mu,
                           const int       a_nu,
                           const RealVect& a_dXi,
                           const Real      a_scale = 1.0) const;

    // Fills an FArrayBox with the contravariant metric elements
    virtual void fill_gup (FArrayBox&      a_dest,
                           const int       a_destComp,
                           const int       a_mu,
                           const int       a_nu,
                           const RealVect& a_dXi,
                           const Real      a_scale = 1.0) const;

    // Fills an FArrayBox with J * gup
    virtual void fill_Jgup (FArrayBox&      a_dest,
                            const int       a_destComp,
                            const int       a_mu,
                            const int       a_nu,
                            const RealVect& a_dXi,
                            const Real      a_scale = 1.0) const;

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
};


#endif //!__CARTESIANMAP_HPP__INCLUDED__

