#include "ProlongationStrategy.H"
#include "ProlongationStrategyF_F.H"


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
        // Create references for convienience
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
        // Create references for convienience
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

