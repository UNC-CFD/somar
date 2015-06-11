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
#ifndef DEMMap_H__INCLUDED__
#define DEMMap_H__INCLUDED__

#include "BathymetricBaseMap.H"
#include "NodeFArrayBox.H"
#include "LevelData.H"
#include "RefCountedPtr.H"


// -----------------------------------------------------------------------------
// Reads the Mass Bay topo from file.
// NOT FOR PUBLIC USE!!! This is a work in progress.
// -----------------------------------------------------------------------------
class DEMMap: public BathymetricBaseMap
{
public:
    // Constructor
    DEMMap ();

    // Destructor
    virtual ~DEMMap ();

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
    static RefCountedPtr<FArrayBox>  s_depthPtr_level0;
    static RefCountedPtr<FArrayBox>  s_depthPtr_level1;
    static Vector<FArrayBox>  s_depthPtr; 
};


#endif //!MassBayMap_H__INCLUDED__
