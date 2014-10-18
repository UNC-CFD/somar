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
#ifndef __GSRB_H__INCLUDED__
#define __GSRB_H__INCLUDED__

#include "RelaxationMethod.H"


// -----------------------------------------------------------------------------
// Abstract base class that performs tasks common to all GSRB methods.
// -----------------------------------------------------------------------------
class BaseGSRB: public RelaxationMethod
{
protected:
    virtual ~BaseGSRB () {;}

    // Applies one sweep of Red-Black Gauss-Seidel to a_phi on just one box.
    // This version uses the full stencil and does not take special care at the
    // box boundaries.
    //
    // NOTE: At the moment, m_activeDirs must be (1,1) or (1,0) in 2D and
    // (1,1,1) or (1,1,0) in 3D.
    virtual void fullStencilGSRB (FArrayBox&       a_phi,
                                  const FArrayBox& a_rhs,
                                  const Box&       a_region,
                                  const DataIndex& a_index,
                                  const int        a_whichPass) const;

    // Relaxes only on cells that abut boundaries and special care is taken when
    // choosing the stencils. If the a_doAllBoundaries flag is set, we iterate over
    // every boundary of every grid in the layout, not just at physical boundaries.
    //
    // NOTE: At the moment, m_activeDirs must be (1,1) or (1,0) in 2D and
    // (1,1,1) or (1,1,0) in 3D.
    virtual void boundaryGSRB (LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>& a_rhs,
                               const int                   a_whichPass,
                               const bool                  a_doAllBoundaries = false) const;
};


// -----------------------------------------------------------------------------
// The complete GSRB method with no shortcuts for speed.
// This performs two synchronous exchanges, one for red and one for black.
// -----------------------------------------------------------------------------
class LevelGSRB: public BaseGSRB
{
public:
    virtual ~LevelGSRB () {;}
    virtual void relax (LevelData<FArrayBox>&       a_phi,
                        const LevelData<FArrayBox>& a_rhs);
};


// -----------------------------------------------------------------------------
// GSRB that only performs one exchange instead of two, asynchronously.
// -----------------------------------------------------------------------------
class LooseGSRB: public BaseGSRB
{
public:
    virtual ~LooseGSRB () {;}
    virtual void relax (LevelData<FArrayBox>&       a_phi,
                        const LevelData<FArrayBox>& a_rhs);
};


// -----------------------------------------------------------------------------
// GSRB that uses a tridiagonal solver for vertical derivatives.
// WARINING: You better not have vertically decomposed grids!
// -----------------------------------------------------------------------------
class LineGSRB: public BaseGSRB
{
public:
    virtual ~LineGSRB () {;}
    virtual void relax (LevelData<FArrayBox>&       a_phi,
                        const LevelData<FArrayBox>& a_rhs);
};


#endif //!__LooseGSRB_H__INCLUDED__
