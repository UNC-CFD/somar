/*******************************************************************************
 *    FILE: BeamGeneratorMap.H
 *    DESC: Terrain-following coordinates for the beam generation problem.
 *    DATE: Tue 30 Jul 2013 08:18:05 PM PDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __BeamGeneratorMap_H__INCLUDED__
#define __BeamGeneratorMap_H__INCLUDED__

#include "BathymetricBaseMap.H"


class BeamGeneratorMap: public BathymetricBaseMap
{
public:
    // Constructor
    BeamGeneratorMap ();

    // Destructor
    virtual ~BeamGeneratorMap ();

    // Must return the name of the coordinate mapping
    virtual const char* getCoorMapName () const;

    // Must return whether or not this metric is diagonal
    virtual bool isDiagonal () const;

protected:
    // Fills a NodeFAB with the bathymetric data. a_dest must be flat in the
    // vertical. Upon return, each point in the horizontal (Xi,Eta) of a_dest
    // will contain the (positive) local depth.
    // NOTE: This vertical distance is measured in a straight line perpendicular
    // the the surface. We are measuring this distance along the Cartesian
    // vertical coordinate line, not the mapped vertical coordinate line.
    virtual void fill_bathymetry (FArrayBox&       a_dest,
                                  const int        a_destComp,
                                  const FArrayBox& a_cartPos,
                                  const RealVect&  a_dXi) const;

    Real m_alpha;
};


#endif //!__BeamGeneratorMap_H__INCLUDED__
