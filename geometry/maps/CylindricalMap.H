/*******************************************************************************
 *    FILE: CylindricalMap.H
 *    DESC: A simple cylindrical coordinate system.
 *    DATE: Sun 03 Nov 2013 01:07:06 PM EST
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __CYLINDRICALMAP_HPP__INCLUDED__
#define __CYLINDRICALMAP_HPP__INCLUDED__

// Project includes
#include "GeoSourceInterface.H"


// -----------------------------------------------------------------------------
// A simple cylindrical coordinate system.
// Xi = r, Nu = theta, Zeta = z.
// -----------------------------------------------------------------------------
class CylindricalMap: public GeoSourceInterface
{
public:
    // Constructor
    CylindricalMap ();

    // Destructor
    virtual ~CylindricalMap ();


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
};


#endif //!__CYLINDRICALMAP_HPP__INCLUDED__

