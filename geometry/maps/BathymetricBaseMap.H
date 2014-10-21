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
