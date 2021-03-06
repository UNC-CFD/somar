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
#ifndef __EXTRAPOLATIONUTILSUTILS_H__INCLUDED__
#define __EXTRAPOLATIONUTILSUTILS_H__INCLUDED__

#include "LevelData.H"
class CFRegion;

// This file contains convenience functions that extrapolate ghosts...


// -----------------------------------------------------------------------------
// Linear extrapolation formula. Assumes all data points are equidistant.
// -----------------------------------------------------------------------------
inline Real linearExtrap (const Real& a_nearVal,
                          const Real& a_farVal)
{
    Real retval = 2.0*a_nearVal - a_farVal;
    return retval;
}


// -----------------------------------------------------------------------------
// Quadratic extrapolation formula. Assumes all data points are equidistant.
// -----------------------------------------------------------------------------
inline Real quadraticExtrap (const Real& a_nearVal,
                             const Real& a_midVal,
                             const Real& a_farVal)
{
    Real retval = 3.0*(a_nearVal - a_midVal) + a_farVal;
    return retval;
}


// -----------------------------------------------------------------------------
// Cubic extrapolation formula. Assumes all data points are equidistant.
// -----------------------------------------------------------------------------
inline Real cubicExtrap (const Real& a_val1,   // closest to extraped value
                         const Real& a_val2,
                         const Real& a_val3,
                         const Real& a_val4)   // farthest extraped value
{
    Real retval = 4.0*(a_val1 + a_val3) - 6.0*a_val2 - a_val4;
    return retval;
}


// -----------------------------------------------------------------------------
// Quartic extrapolation formula. Assumes all data points are equidistant.
// -----------------------------------------------------------------------------
inline Real quarticExtrap (const Real& a_val1, // closest to extraped value
                           const Real& a_val2,
                           const Real& a_val3,
                           const Real& a_val4,
                           const Real& a_val5) // farthest from extraped value
{
    Real retval = 5.0 * (a_val1 - a_val4) + 10.0 * (a_val3 - a_val2) + a_val5;
    return retval;
}


// -----------------------------------------------------------------------------
// Extrapolation boundary conditions for a side, specified component interval
// This is a special version that works out of place and copies an inner layer.
// For use in this file only.
// -----------------------------------------------------------------------------
void ExtrapolateFaceNoEV (FArrayBox&           a_dest,
                          const FArrayBox&     a_state,
                          const Box&           a_valid,
                          const int            a_dir,
                          const Side::LoHiSide a_side,
                          const int            a_order);


// -----------------------------------------------------------------------------
// Extrapolation BCs. Does not fill edges and vertices.
// -----------------------------------------------------------------------------
void ExtrapolateBCNoEV (FArrayBox&           a_state,
                        const Box&           a_valid,
                        const ProblemDomain& a_domain,
                        const int            a_order,
                        const IntVect&       a_doLoSides = IntVect::Unit,
                        const IntVect&       a_doHiSides = IntVect::Unit);


// -----------------------------------------------------------------------------
// Extrapolation boundary conditions for a side.
// This works out of place and also extrapolates onto faces a_numLayers deep.
// Edge and vertex ghosts will be extrapolated from ghosts adjacent to boundary.
// If a_dest != a_state, the inner layers will also be copied.
// -----------------------------------------------------------------------------
void ExtrapolateFaceAndCopy (FArrayBox&           a_dest,
                             const FArrayBox&     a_state,
                             const Box&           a_valid,
                             const int            a_dir,
                             const Side::LoHiSide a_side,
                             const int            a_order,
                             const int            a_numLayers = 1);

// -----------------------------------------------------------------------------
// Extrapolates edge and vertex ghosts at the fine side of the CF interface.
// This is used to finish the job that MappedQuadCFInterp started.
// For future reference, the coeff of val_k using an order n extrap is
//   (-1)^(k+1) * Binomial[n+1,k]
// where val_1 is the nearVal and val_(n+1) is the farVal.
// -----------------------------------------------------------------------------
void ExtrapolateCFEV (FArrayBox&           a_state,
                      const CFRegion&      a_region,
                      const DataIndex&     a_di,
                      const int            a_dir,
                      const Side::LoHiSide a_side,
                      const Interval&      a_interval,
                      const int            a_order,
                      const IntVect&       a_activeDirs = IntVect::Unit);

// -----------------------------------------------------------------------------
// LevelData version.
// -----------------------------------------------------------------------------
void ExtrapolateCFEV (LevelData<FArrayBox>& a_state,
                      const CFRegion&       a_region,
                      const int             a_order,
                      const IntVect&        a_activeDirs = IntVect::Unit);


// -----------------------------------------------------------------------------
// Debugging function.
// Extrapolate/exchange all ghosts. Useful for testing if ghosts need to be set.
// This function does no exchanges.
// -----------------------------------------------------------------------------
void extrapAllGhosts (LevelData<FArrayBox>& a_data,
                      const int             a_order = 2);


// -----------------------------------------------------------------------------
// Debugging function.
// Extrapolate/exchange all ghosts. Useful for testing if ghosts need to be set.
// This function does no exchanges.
// -----------------------------------------------------------------------------
void extrapAllGhosts (LevelData<FluxBox>& a_data,
                      const int           a_order = 2);


// -----------------------------------------------------------------------------
// Extrapolates ghosts at the fine side of the CF interface.
// This is a special version that works out of place.
// For use in MappedAMRPoissonOp.
// For future reference, the coeff of val_k using an order n extrap is
//   (-1)^(k+1) * Binomial[n+1,k]
// where val_1 is the nearVal and val_(n+1) is the farVal.
// -----------------------------------------------------------------------------
void ExtrapolateCF(FArrayBox&       a_dest,
                   const FArrayBox& a_state,
                   const CFRegion&  a_region,
                   const DataIndex& a_di,
                   int              a_dir,
                   Side::LoHiSide   a_side,
                   const Interval&  a_interval,
                   int              a_order);


// -----------------------------------------------------------------------------
// Extrapolates ghosts at the fine side of the CF interface.
// This works in place and does all dirs and sides.
// For use in MappedAMRPoissonOp.
// -----------------------------------------------------------------------------
void ExtrapolateAllCF(FArrayBox&       a_state,
                      const CFRegion&  a_region,
                      const DataIndex& a_di,
                      int              a_order);


#endif //!__EXTRAPOLATIONUTILSUTILS_H__INCLUDED__

