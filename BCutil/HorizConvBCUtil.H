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
#ifndef __HorizConvBCUtil_H__INCLUDED__
#define __HorizConvBCUtil_H__INCLUDED__

#include "PhysBCUtil.H"


// -----------------------------------------------------------------------------
// Horizontal convection problem a la Scotti & White.
// -----------------------------------------------------------------------------
class HorizConvBCUtil: public PhysBCUtil
{
public:

    // Default constructor
    HorizConvBCUtil ();

    // Default destructor
    virtual ~HorizConvBCUtil ();

    // Factory
    virtual PhysBCUtil* newPhysBCUtil () const;

    // Fills a FAB with the initial scalars
    virtual void setScalarIC (FArrayBox&           a_scalarFAB,
                              const int            a_scalarComp,
                              const LevelGeometry& a_levGeo,
                              const DataIndex&     a_di) const;

    // Used in single-compoent velocity TGA solves
    virtual BCMethodHolder viscousSolveFuncBC (int a_dir) const;

    // used in scalar TGA solves
    virtual BCMethodHolder diffusiveSolveFuncBC () const;

    // Sets ghosts needed to calculate the diffusive term nu.L[scalar]
    virtual BCMethodHolder diffusiveSourceFuncBC () const;


protected:
    // Sets physical BCs on velocities
    virtual BCMethodHolder basicVelFuncBC (int a_veldir, bool a_isViscous) const;
};


#endif //!__HorizConvBCUtil_H__INCLUDED__

