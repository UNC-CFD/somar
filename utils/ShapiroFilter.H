/*******************************************************************************
 *    FILE: ShapiroFilter.H
 *    DESC: Shapiro filters are used by ROMS and mitGCM.
 *          For more info, see Ralph Shapiro's 1970 paper:
 *          "Smoothing, filtering, and boundary effects",
 *          DOI: 10.1029/RG008i002p00359.
 *    DATE: Sun 21 Jul 2013 03:02:40 PM PDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __ShapiroFilter_H__INCLUDED__
#define __ShapiroFilter_H__INCLUDED__

#include "LevelGeometry.H"
#include "LevelData.H"
#include "MappedQuadCFInterp.H"
#include "EllipticBCInterface.H"
#include "VelBCHolder.H"


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
