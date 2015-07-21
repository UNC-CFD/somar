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
#include "ProlongationStrategy.H"
#include "ProlongationStrategyF_F.H"


// This is temporary! CHF_CONST_FRA1_SHIFT has a bug.
#undef CHF_CONST_FRA1_SHIFT
#define CHF_CONST_FRA1_SHIFT( a, n, iv)      \
  a.dataPtr( n ),                            \
  D_DECL6( CHFPTR(a.loVect()[0] - iv[0]),    \
           CHFPTR(a.loVect()[1] - iv[1]),    \
           CHFPTR(a.loVect()[2] - iv[2]),    \
           CHFPTR(a.loVect()[3] - iv[3]),    \
           CHFPTR(a.loVect()[4] - iv[4]),    \
           CHFPTR(a.loVect()[5] - iv[5]) ),  \
  D_DECL6( CHFPTR(a.hiVect()[0] - iv[0]),    \
           CHFPTR(a.hiVect()[1] - iv[1]),    \
           CHFPTR(a.hiVect()[2] - iv[2]),    \
           CHFPTR(a.hiVect()[3] - iv[3]),    \
           CHFPTR(a.hiVect()[4] - iv[4]),    \
           CHFPTR(a.hiVect()[5] - iv[5]) )


// -----------------------------------------------------------------------------
// Adds coarse cell values directly to all overlying fine cells.
// -----------------------------------------------------------------------------
void ConstInterpPS::prolongIncrement (LevelData<FArrayBox>&       a_phiThisLevel,
                                      const LevelData<FArrayBox>& a_correctCoarse)
{
    CH_TIME("ConstInterpPS::prolongIncrement");

    // Gather grids, domains, refinement ratios...
    const DisjointBoxLayout& fineGrids = a_phiThisLevel.getBoxes();
    const DisjointBoxLayout& crseGrids = a_correctCoarse.getBoxes();
    CH_assert(fineGrids.compatible(crseGrids));

    const ProblemDomain& fineDomain = fineGrids.physDomain();
    const ProblemDomain& crseDomain = crseGrids.physDomain();

    const IntVect mgRefRatio = fineDomain.size() / crseDomain.size();
    CH_assert(mgRefRatio.product() > 1);

    DataIterator dit = fineGrids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        // Create references for convenience
        FArrayBox& fineFAB = a_phiThisLevel[dit];
        const FArrayBox& crseFAB = a_correctCoarse[dit];
        const Box& fineValid = fineGrids[dit];

        // To make things easier, we will offset the
        // coarse and fine data boxes to zero.
        const IntVect& fiv = fineValid.smallEnd();
        const IntVect civ = coarsen(fiv, mgRefRatio);

        // Correct the fine data
        FORT_CONSTINTERPPS (
            CHF_FRA_SHIFT(fineFAB, fiv),
            CHF_CONST_FRA_SHIFT(crseFAB, civ),
            CHF_BOX_SHIFT(fineValid, fiv),
            CHF_CONST_INTVECT(mgRefRatio));
    }
}


// -----------------------------------------------------------------------------
// Adds coarse cell values directly to all overlying fine cells,
// then removes the average from the fine result.
// -----------------------------------------------------------------------------
void ZeroAvgConstInterpPS::prolongIncrement (LevelData<FArrayBox>&       a_phiThisLevel,
                                             const LevelData<FArrayBox>& a_correctCoarse)
{
    CH_TIME("ZeroAvgConstInterpPS::prolongIncrement");

    // Gather grids, domains, refinement ratios...
    const DisjointBoxLayout& fineGrids = a_phiThisLevel.getBoxes();
    const DisjointBoxLayout& crseGrids = a_correctCoarse.getBoxes();
    CH_assert(fineGrids.compatible(crseGrids));

    const ProblemDomain& fineDomain = fineGrids.physDomain();
    const ProblemDomain& crseDomain = crseGrids.physDomain();

    const IntVect mgRefRatio = fineDomain.size() / crseDomain.size();
    CH_assert(mgRefRatio.product() > 1);

    // These will accumulate averaging data.
    Real localSum = 0.0;
    Real localVol = 0.0;
    CH_assert(!m_CCJinvPtr.isNull());
    CH_assert(m_dxProduct > 0.0);

    DataIterator dit = fineGrids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        // Create references for convenience
        FArrayBox& fineFAB = a_phiThisLevel[dit];
        const FArrayBox& crseFAB = a_correctCoarse[dit];
        const Box& fineValid = fineGrids[dit];
        const FArrayBox& JinvFAB = (*m_CCJinvPtr)[dit];

        // To make things easier, we will offset the
        // coarse and fine data boxes to zero.
        const IntVect& fiv = fineValid.smallEnd();
        const IntVect civ = coarsen(fiv, mgRefRatio);

        // Correct the fine data
        FORT_CONSTINTERPWITHAVGPS (
            CHF_FRA_SHIFT(fineFAB, fiv),
            CHF_CONST_FRA_SHIFT(crseFAB, civ),
            CHF_BOX_SHIFT(fineValid, fiv),
            CHF_CONST_INTVECT(mgRefRatio),
            CHF_CONST_FRA1_SHIFT(JinvFAB,0,fiv),
            CHF_CONST_REAL(m_dxProduct),
            CHF_REAL(localVol),
            CHF_REAL(localSum));
    }

    // Compute global sum (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalSum = 0.0;
    int result = MPI_Allreduce(&localSum, &globalSum, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in ZeroAvgConstInterpPS::prolongIncrement");
    }

    Real globalVol = 0.0;
    result = MPI_Allreduce(&localVol, &globalVol, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in ZeroAvgConstInterpPS::prolongIncrement");
    }

#else
    Real globalSum = localSum;
    Real globalVol = localVol;
#endif

    // Remove the average from phi.
    Real avgPhi = globalSum / globalVol;
    for (dit.reset(); dit.ok(); ++dit) {
        a_phiThisLevel[dit] -= avgPhi;
    }
}

