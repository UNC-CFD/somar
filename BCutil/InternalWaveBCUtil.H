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
#ifndef __InternalWaveBCUtil_H__INCLUDED__
#define __InternalWaveBCUtil_H__INCLUDED__

#include "PhysBCUtil.H"


// -----------------------------------------------------------------------------
// Simple internal wave demo.
// -----------------------------------------------------------------------------
class InternalWaveBCUtil: public PhysBCUtil
{
public:

    // Default constructor
    InternalWaveBCUtil ();

    // Default destructor
    virtual ~InternalWaveBCUtil ();

    // Factory
    virtual PhysBCUtil* newPhysBCUtil () const;

    // This is in case the BC's have an effect on the timestep.
    // Pass in currently computed dt, along with the cfl and dx. If the effect
    // of the BCs requires a decreased timestep, then the newly reduced timestep
    // is returned.
    virtual void computeBoundaryDt (Real&                a_dt,
                                    const Real           a_cfl,
                                    const LevelGeometry& a_levGeo) const;

    // Fills a FAB with the initial scalars
    virtual void setScalarIC (FArrayBox&           a_scalarFAB,
                              const int            a_scalarComp,
                              const LevelGeometry& a_levGeo,
                              const DataIndex&     a_di) const;

    // Fills a FAB with the background scalar
    virtual void setBackgroundScalar (FArrayBox&           a_scalarFAB,
                                      const int            a_scalarComp,
                                      const LevelGeometry& a_levGeo,
                                      const DataIndex&     a_di,
                                      const Real           a_time) const;
protected:
    static void bscalBCValues (Real*           a_pos,
                               int*            a_dir,
                               Side::LoHiSide* a_side,
                               Real*           a_value,
                               Real            a_derivScale,
                               Real            a_time);
};



#endif //!__InternalWaveBCUtil_H__INCLUDED__
