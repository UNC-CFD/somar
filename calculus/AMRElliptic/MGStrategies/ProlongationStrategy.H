/*******************************************************************************
 *    FILE: ProlongationStrategy.H
 *    DESC: The interface shared by all MG prolongation methods.
 *          This file also defines some common strategies.
 *    DATE: Wed 27 Nov 2013 03:13:32 PM EST
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __ProlongationStrategy_H__INCLUDED__
#define __ProlongationStrategy_H__INCLUDED__

#include "LevelData.H"
#include "FArrayBox.H"
#include "RefCountedPtr.H"


// -----------------------------------------------------------------------------
// The strategy interface.
// -----------------------------------------------------------------------------
class ProlongationStrategy
{
public:
    virtual ~ProlongationStrategy () {;}

    // This function is designed to have a signature identical to that of
    // an MGLevelOp. I am doing this because one day, this really should be
    // made into a policy that helps construct an MGLevelOp.
    //
    // From the comments of MGLevelOp::prolongIncrement:
    // This function corrects the fine solution based on a coarse correction.
    // a_phiThisLevel += I[2h->h](a_correctCoarse)
    virtual void prolongIncrement (LevelData<FArrayBox>&       a_phiThisLevel,
                                   const LevelData<FArrayBox>& a_correctCoarse) = 0;
};


// -----------------------------------------------------------------------------
// Prolongation strategy: constant interpolation.
// -----------------------------------------------------------------------------
class ConstInterpPS: public ProlongationStrategy
{
public:
    virtual ~ConstInterpPS () {;}

    // Adds coarse cell values directly to all overlying fine cells.
    virtual void prolongIncrement (LevelData<FArrayBox>&       a_phiThisLevel,
                                   const LevelData<FArrayBox>& a_correctCoarse);
};


// -----------------------------------------------------------------------------
// Prolongation strategy: constant interpolation + remove fine average.
// -----------------------------------------------------------------------------
class ZeroAvgConstInterpPS: public ProlongationStrategy
{
public:
    // Constructor -- If this is a horizontal solve, it is up to the caller to
    // set dxProduct accordingly.
    ZeroAvgConstInterpPS (const Real                                  a_dxProduct,
                          const RefCountedPtr<LevelData<FArrayBox> >& a_CCJinvPtr)
    : m_dxProduct(a_dxProduct),
      m_CCJinvPtr(a_CCJinvPtr)
    {;}

    // Destructor
    virtual ~ZeroAvgConstInterpPS () {;}

    // Adds coarse cell values directly to all overlying fine cells,
    // then removes the average from the fine result.
    virtual void prolongIncrement (LevelData<FArrayBox>&       a_phiThisLevel,
                                   const LevelData<FArrayBox>& a_correctCoarse);

protected:
    Real m_dxProduct;
    const RefCountedPtr<LevelData<FArrayBox> > m_CCJinvPtr;
};


#endif //!__ProlongationStrategy_H__INCLUDED__
