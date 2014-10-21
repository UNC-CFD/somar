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
#ifndef __ShapiroFilter_H__INCLUDED__
#define __ShapiroFilter_H__INCLUDED__

#include "LevelGeometry.H"
#include "LevelData.H"
#include "MappedQuadCFInterp.H"
#include "EllipticBCInterface.H"
#include "VelBCHolder.H"


// Shapiro filters are used by ROMS and mitGCM.
// For more info, see Ralph Shapiro's 1970 paper:
//     "Smoothing, filtering, and boundary effects",
//     DOI: 10.1029/RG008i002p00359.


// -----------------------------------------------------------------------------
// This uses a 3-point stencil on each pass. The parameters and number of passes
// are hard-coded into this funciton, but can be changed easily.
// -----------------------------------------------------------------------------
void ShapiroFilter1D (// Data stuff
                      LevelData<FArrayBox>&       a_data,
                      const Interval              a_interval,
                      // CF BC stuff
                      const LevelData<FArrayBox>* a_crseDataPtr,
                      const MappedQuadCFInterp&   a_cfInterp,
                      // Physical BC stuff
                      VelBCHolder                 a_bc,
                      const bool                  a_isHomogeneous,
                      const LevelGeometry&        a_levGeo,
                      const Real                  a_time,
                      // Filter stuff
                      const int                   a_filterDir);


// -----------------------------------------------------------------------------
// This uses a 5-point stencil on each pass. The parameters and number of passes
// are hard-coded into this funciton, but can be changed easily.
// Use this carefully. I found that the 9-point stencil that arises from using
// ShapiroFilter1D twice (it's commutative) produces better results. Shapiro
// talks about this in his paper.
// -----------------------------------------------------------------------------
void ShapiroFilter2D (// Data stuff
                      LevelData<FArrayBox>&       a_data,
                      const Interval              a_interval,
                      // CF BC stuff
                      const LevelData<FArrayBox>* a_crseDataPtr,
                      const MappedQuadCFInterp&   a_cfInterp,
                      // Physical BC stuff
                      VelBCHolder                 a_bc,
                      const bool                  a_isHomogeneous,
                      const LevelGeometry&        a_levGeo,
                      const Real                  a_time,
                      // Filter stuff
                      const int                   a_filterDir0,
                      const int                   a_filterDir1);


#endif //!__ShapiroFilter_H__INCLUDED__
