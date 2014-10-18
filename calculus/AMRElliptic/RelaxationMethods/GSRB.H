/*******************************************************************************
 *    FILE: GSRB.H
 *    DESC: Red-black Gauss-Seidel methods.
 *    DATE: Fri 22 Nov 2013 04:47:18 PM EST
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __GSRB_H__INCLUDED__
#define __GSRB_H__INCLUDED__

#include "RelaxationMethod.H"


// Abstract base class that performs tasks common to all GSRB methods.
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


// The complete GSRB method with no shortcuts for speed.
// This performs two synchronous exchanges, one for red and one for black.
class LevelGSRB: public BaseGSRB
{
public:
    virtual ~LevelGSRB () {;}
    virtual void relax (LevelData<FArrayBox>&       a_phi,
                        const LevelData<FArrayBox>& a_rhs);
};


// GSRB that only performs one exchange instead of two, asynchronously.
class LooseGSRB: public BaseGSRB
{
public:
    virtual ~LooseGSRB () {;}
    virtual void relax (LevelData<FArrayBox>&       a_phi,
                        const LevelData<FArrayBox>& a_rhs);
};


// GSRB that uses a tridiagonal solver for vertical derivatives.
// WARINING: You better not be vertically decomposing your grids!
class LineGSRB: public BaseGSRB
{
public:
    virtual ~LineGSRB () {;}
    virtual void relax (LevelData<FArrayBox>&       a_phi,
                        const LevelData<FArrayBox>& a_rhs);
};


#endif //!__LooseGSRB_H__INCLUDED__
