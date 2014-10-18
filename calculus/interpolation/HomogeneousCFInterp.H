/*******************************************************************************
 *    FILE: HomogeneousCFInterp.H
 *    DESC: Assumes 0 on a coarser level and interpolates ghosts.
 *    DATE: Fri 22 Nov 2013 02:49:43 PM EST
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __FileName_H__INCLUDED__
#define __FileName_H__INCLUDED__

#include "CFRegion.H"
#include "IntVectSet.H"
#include "LevelData.H"
#include "RealVect.H"


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


#endif //!__FileName_H__INCLUDED__
