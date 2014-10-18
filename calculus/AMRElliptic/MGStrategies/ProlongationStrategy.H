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
#ifndef __ProlongationStrategy_H__INCLUDED__
#define __ProlongationStrategy_H__INCLUDED__

#include "LevelData.H"
#include "FArrayBox.H"
#include "RefCountedPtr.H"


// -----------------------------------------------------------------------------
// The interface shared by all MG prolongation methods.
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
