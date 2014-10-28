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
#ifndef __HomogeneousCFInterp_H__INCLUDED__
#define __HomogeneousCFInterp_H__INCLUDED__

#include "CFRegion.H"
#include "IntVectSet.H"
#include "LevelData.H"
#include "RealVect.H"


// Utilities that interpolate ghosts assuming 0 on a coarser level.
// These functions are modifications of functions originally found in Chombo's
// AMRPoissonOp class.


// Interpolate ghosts at CF interface using zeros on coarser grids.
// This version acts on a set of IntVects.
void interpOnIVSHomo (LevelData<FArrayBox>& a_phif,
                      const DataIndex&      a_index,
                      const int             a_dir,
                      const Side::LoHiSide  a_side,
                      const IntVectSet&     a_interpIVS,
                      const Real            a_fineDxDir,
                      const Real            a_crseDxDir);

// Interpolate ghosts at CF interface using zeros on coarser grids.
// Only do one face.
void homogeneousCFInterp (LevelData<FArrayBox>& a_phif,
                          const DataIndex&      a_index,
                          const int             a_dir,
                          const Side::LoHiSide  a_side,
                          const Real            a_fineDxDir,
                          const Real            a_crseDxDir,
                          const CFRegion&       a_cfRegion);

// Interpolate ghosts at CF interface using zeros on coarser grids.
void homogeneousCFInterp (LevelData<FArrayBox>& a_phif,
                          const RealVect&       a_fineDx,
                          const RealVect&       a_crseDx,
                          const CFRegion&       a_cfRegion,
                          const IntVect&        a_applyDirs = IntVect::Unit);


#endif //!__HomogeneousCFInterp_H__INCLUDED__
