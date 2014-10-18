/*******************************************************************************
 *    FILE: LedgeMap.H
 *    DESC: A simple underwater ledge with a smooth, cubic transition region.
 *    DATE: Fri 07 Feb 2014 01:01:34 PM EST
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __LedgeMap_H__INCLUDED__
#define __LedgeMap_H__INCLUDED__

#include "BathymetricBaseMap.H"


class LedgeMap: public BathymetricBaseMap
{
public:
    // Constructor
    LedgeMap ();

    // Destructor
    virtual ~LedgeMap ();

    // Must return the name of the coordinate mapping
    virtual const char* getCoorMapName () const;

    // Must return whether or not this metric is diagonal
    virtual bool isDiagonal () const;

protected:
    // Fills a NodeFAB with the bathymetric data. a_dest must be flat in the
    // vertical. Upon return, each point in the horizontal (Xi,Eta) of a_dest
    // will contain the (positive) local depth.
    // NOTE: This vertical distance is measured in a straight line perpendicular
    // to the surface. We are measuring this distance along the Cartesian
    // vertical coordinate line, not the mapped vertical coordinate line.
    virtual void fill_bathymetry (FArrayBox&       a_dest,
                                  const int        a_destComp,
                                  const FArrayBox& a_cartPos,
                                  const RealVect&  a_dXi) const;

    int m_transitionOrder;
    Real m_hl, m_hr, m_xl, m_xr;
    Real m_coeff0, m_coeff1, m_coeff2, m_coeff3;
};


#endif //!__LedgeMap_H__INCLUDED__
