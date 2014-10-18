#include "RestrictionStrategy.H"
#include "MappedCoarseAverageF_F.H"


FullWeightingPS::FullWeightingPS (const RefCountedPtr<LevelData<FArrayBox> >& a_JinvPtr)
: m_JinvPtr(a_JinvPtr)
{;}


// -----------------------------------------------------------------------------
// Set coarse cells to the average of the overlying fine cells.
// -----------------------------------------------------------------------------
void FullWeightingPS::restrict (LevelData<FArrayBox>&       a_crseRes,
                                const LevelData<FArrayBox>& a_fineRes)
{
    CH_TIME("FullWeightingPS::restrictResidual");

    // Tripping an assert in C is easier to debug than an invalid memory access
    // in fortran!
    CH_assert(a_crseRes.nComp() == a_fineRes.nComp());

    // Gather grids, domains, refinement ratios...
    const DisjointBoxLayout& fineGrids = a_fineRes.getBoxes();
    const DisjointBoxLayout& crseGrids = a_crseRes.getBoxes();
    CH_assert(fineGrids.compatible(crseGrids));

    const ProblemDomain& fineDomain = fineGrids.physDomain();
    const ProblemDomain& crseDomain = crseGrids.physDomain();

    const IntVect mgRefRatio = fineDomain.size() / crseDomain.size();
    CH_assert(mgRefRatio.product() > 1);

    // Average residual to coarse level.
    const Box refBox(IntVect::Zero, mgRefRatio - IntVect::Unit);
    DataIterator dit = fineGrids.dataIterator();

#   define USE_J_WEIGHTING
#   ifdef USE_J_WEIGHTING
    {
        // Average with J weighting.
        CH_assert(a_fineRes .getBoxes().compatible(fineGrids));
        CH_assert(m_JinvPtr->getBoxes().compatible(fineGrids));

        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& crseFAB = a_crseRes[dit];
            const FArrayBox& fineFAB = a_fineRes[dit];
            const FArrayBox& fineJinvFAB = (*m_JinvPtr)[dit];
            const Box& crseValid = crseGrids[dit];

            FORT_MAPPEDAVERAGE2 (
                CHF_FRA(crseFAB),
                CHF_CONST_FRA(fineFAB),
                CHF_CONST_FRA1(fineJinvFAB,0),
                CHF_BOX(crseValid),
                CHF_CONST_INTVECT(mgRefRatio),
                CHF_BOX(refBox));
        }
    }
#   else
    {
        // Average without J weighting.
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& crseFAB = a_crseRes[dit];
            const FArrayBox& fineFAB = a_fineRes[dit];
            const Box& crseValid = crseGrids[dit];

            FORT_UNMAPPEDAVERAGE (
                CHF_FRA(crseFAB),
                CHF_CONST_FRA(fineFAB),
                CHF_BOX(crseValid),
                CHF_CONST_INTVECT(mgRefRatio),
                CHF_BOX(refBox));
        }
    }
#   endif
}
