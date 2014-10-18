/*******************************************************************************
 *    FILE: MassBayMap.H
 *    DESC: Reads the Mass Bay topo from file.
 *    DATE: Sun 01 Jun 2014 08:12:02 PM CDT  (UTC-0400)
 *    MAIL: santilli@physics.unc.edu
 *
 *    NOTE: NOT FOR PUBLIC USE!!! This is a work in progress.
 *
 ******************************************************************************/

#ifndef MassBayMap_H__INCLUDED__
#define MassBayMap_H__INCLUDED__

#include "BathymetricBaseMap.H"
#include "NodeFArrayBox.H"
#include "LevelData.H"
#include "RefCountedPtr.H"


class MassBayMap: public BathymetricBaseMap
{
public:
    // Constructor
    MassBayMap ();

    // Destructor
    virtual ~MassBayMap ();

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

    // Reads the bathymetry from file and interpolates onto level 0 grids.
    // This also trims the region and smooths the edges.
    virtual void createBaseMap ();

    static bool                                 s_fileIsRead;
    static RefCountedPtr<LevelData<FArrayBox> > s_depthPtr;
};


#endif //!MassBayMap_H__INCLUDED__
