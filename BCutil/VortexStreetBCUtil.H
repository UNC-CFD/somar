/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Copyright (C) 2014 Edward Santilli & Alberto Scotti
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
#ifndef __VortexStreetBCUtil_H__INCLUDED__
#define __VortexStreetBCUtil_H__INCLUDED__

#include "PhysBCUtil.H"


// -----------------------------------------------------------------------------
// IBCs for the Von Karman vortex demonstration.
// -----------------------------------------------------------------------------
class VortexStreetBCUtil: public PhysBCUtil
{
public:

    // Default constructor
    VortexStreetBCUtil ();

    // Default destructor
    virtual ~VortexStreetBCUtil ();

    // Factory
    virtual PhysBCUtil* newPhysBCUtil () const;

    // Fills a FAB with the initial velocity.
    // Locations are in mapped space, but components are Cartesian.
    virtual void setVelIC (FArrayBox&           a_velFAB,
                           const int            a_velComp,
                           const LevelGeometry& a_levGeo,
                           const DataIndex&     a_di) const;

    // Sets the Cartesian-based target velocity for the sponge layer.
    // By default, this function throws an error.
    virtual void fillVelSpongeLayerTarget (FArrayBox&           a_target,
                                           const int            a_velComp,
                                           const int            a_spongeDir,
                                           const Side::LoHiSide a_spongeSide,
                                           const LevelGeometry& a_levGeo,
                                           const DataIndex&     a_di,
                                           const Real           a_time);

    // Sets the target values for the scalar sponge layer. If we are using a
    // background scalar, then this function set the perturbation to zero.
    // Otherwise, an error is thrown and this function will need to be overridden.
    virtual void fillScalarSpongeLayerTarget (FArrayBox&           a_target,
                                              const int            a_scalarComp,
                                              const int            a_spongeDir,
                                              const Side::LoHiSide a_spongeSide,
                                              const LevelGeometry& a_levGeo,
                                              const DataIndex&     a_di,
                                              const Real           a_time);

protected:
    // Sets physical BCs on velocities
    virtual BCMethodHolder basicVelFuncBC (int a_veldir, bool a_isViscous) const;

    static bool s_paramsRead;
    static RealVect s_inflowVel;
};


#endif //!__VortexStreetBCUtil_H__INCLUDED__
