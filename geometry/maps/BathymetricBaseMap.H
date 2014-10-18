/*******************************************************************************
 *    FILE: BathymetricBaseMap.H
 *    DESC: Abstract base class for all terrain-following coordinate maps.
 *    DATE: Sun 12 May 2013 12:08 PM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __BathymetricBaseMap_H__INCLUDED__
#define __BathymetricBaseMap_H__INCLUDED__

#include "GeoSourceInterface.H"


// -----------------------------------------------------------------------------
// Abstract base class for all terrain-following coordinate maps.
// -----------------------------------------------------------------------------
class BathymetricBaseMap: public GeoSourceInterface
{
public:
    // Constructor
    BathymetricBaseMap ();

    // Destructor
    virtual ~BathymetricBaseMap ();

    // It is assumed that no bathymetric map is uniform.
    virtual inline bool isUniform () const;

    // Fills a mapped box with Cartesian locations.
    virtual void fill_physCoor (FArrayBox&      a_dest,
                                const int       a_destComp,
                                const int       a_mu,
                                const RealVect& a_dXi) const;

    // Fills an FArrayBox with the Jacobian matrix elements d[x^mu] / d[xi^nu].
    virtual void fill_dxdXi (FArrayBox&      a_dest,
                             const int       a_destComp,
                             const int       a_mu,
                             const int       a_nu,
                             const RealVect& a_dXi,
                             const Real      a_scale = 1.0) const;

    // Fills an FArrayBox with J = det[Jacobian]
    virtual void fill_J (FArrayBox&      a_dest,
                         const int       a_destComp,
                         const RealVect& a_dXi,
                         const Real      a_scale = 1.0) const;

protected:
    // Fills a NodeFAB with the bathymetric data. a_dest must be flat in the
    // vertical. Upon return, each point in the horizontal (Xi,Eta) of a_dest
    // will contain the (positive) local elevation from the sea floor.
    // NOTE: This vertical distance is measured in a straight line perpendicular
    // the the surface. We are measuring this distance along the Cartesian
    // vertical coordinate line, not the mapped vertical coordinate line.
    virtual void fill_bathymetry (FArrayBox&       a_dest,
                                  const int        a_destComp,
                                  const FArrayBox& a_cartPos,
                                  const RealVect&  a_dXi) const = 0;

    RealVect m_L;
    RealVect m_lev0DXi;
};


// -----------------------------------------------------------------------------
// It is assumed that no bathymetric map is uniform.
// -----------------------------------------------------------------------------
inline bool BathymetricBaseMap::isUniform () const
{
    return false;
}


#endif //!__BathymetricBaseMap_H__INCLUDED__
