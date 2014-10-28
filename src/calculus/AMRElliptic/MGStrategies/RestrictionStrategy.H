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
#ifndef __RestrictionStrategy_H__INCLUDED__
#define __RestrictionStrategy_H__INCLUDED__

#include "LinearSolver.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "RefCountedPtr.H"


// -----------------------------------------------------------------------------
// The interface shared by all MG restriction methods.
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
