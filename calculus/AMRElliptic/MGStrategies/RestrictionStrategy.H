/*******************************************************************************
 *    FILE: RestrictionStrategy.H
 *    DESC: The interface shared by all MG restriction methods.
 *          This file also defines some common strategies.
 *    DATE: Wed 27 Nov 2013 05:07:58 PM EST
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __RestrictionStrategy_H__INCLUDED__
#define __RestrictionStrategy_H__INCLUDED__

#include "LinearSolver.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "RefCountedPtr.H"


// -----------------------------------------------------------------------------
// The strategy interface.
// -----------------------------------------------------------------------------
class RestrictionStrategy
{
public:
    virtual ~RestrictionStrategy () {;}

    // This will be used to restrict the fine-level residual. It is up to the
    // linear operator to compute the residual before calling this function.
    virtual void restrict (LevelData<FArrayBox>&       a_crseRes,
                           const LevelData<FArrayBox>& a_fineRes) = 0;
};


// -----------------------------------------------------------------------------
// Restriction strategy: full weighting.
// -----------------------------------------------------------------------------
class FullWeightingPS: public RestrictionStrategy
{
public:
    FullWeightingPS (const RefCountedPtr<LevelData<FArrayBox> >& a_JinvPtr);

    virtual ~FullWeightingPS () {;}

    // Set coarse cells to the average of the overlying fine cells.
    virtual void restrict (LevelData<FArrayBox>&       a_crseRes,
                           const LevelData<FArrayBox>& a_fineRes);
protected:
    RefCountedPtr<LevelData<FArrayBox> > m_JinvPtr;
};




#endif //!__RestrictionStrategy_H__INCLUDED__
