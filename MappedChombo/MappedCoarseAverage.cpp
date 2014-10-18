#include "LevelData.H"
#include "BaseFab.H"
#include "REAL.H"
#include "FArrayBox.H"
#include "DataIterator.H"
#include "LayoutIterator.H"
#include "parstream.H"

#include "LevelGeometry.H"
#include "MappedCoarseAverageF_F.H"
#include "MappedCoarseAverage.H"
#include "AnisotropicRefinementTools.H"


// -----------------------------------------------------------------------------
// Default constructor -- leaves object unusable
// -----------------------------------------------------------------------------
MappedCoarseAverage::MappedCoarseAverage ()
: m_isDefined(false),
  m_isCopierDefined(false),
  m_coarsenedFineData(NULL)
{;}


// -----------------------------------------------------------------------------
// Full constructor
// -----------------------------------------------------------------------------
MappedCoarseAverage::MappedCoarseAverage (const DisjointBoxLayout& a_fineGrids,
                                          const int                a_numComps,
                                          const IntVect&           a_refRatio)
: m_isDefined(false),
  m_isCopierDefined(false),
  m_coarsenedFineData(NULL)
{
    this->define(a_fineGrids, a_numComps, a_refRatio);
}


// -----------------------------------------------------------------------------
// Full constructor -- caches a Copier.
// -----------------------------------------------------------------------------
MappedCoarseAverage::MappedCoarseAverage (const DisjointBoxLayout& a_fineGrids,
                                          const DisjointBoxLayout& a_crseGrids,
                                          const int                a_numComps,
                                          const IntVect&           a_refRatio,
                                          const IntVect&           a_ghostVect)
: m_isDefined(false),
  m_isCopierDefined(false),
  m_coarsenedFineData(NULL)
{
    this->define(a_fineGrids, a_crseGrids, a_numComps, a_refRatio, a_ghostVect);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
MappedCoarseAverage::~MappedCoarseAverage ()
{
    delete m_coarsenedFineData;
    m_coarsenedFineData = NULL;
}


// -----------------------------------------------------------------------------
// Full define constructor
// -----------------------------------------------------------------------------
void MappedCoarseAverage::define (const DisjointBoxLayout& a_fineGrids,
                                  const int                a_numComps,
                                  const IntVect&           a_refRatio)
{
    CH_TIME("MappedCoarseAverage::define_no_copier");
    CH_assert(0 < a_numComps);
    D_TERM(CH_assert(0 < a_refRatio[0]);,
           CH_assert(0 < a_refRatio[1]);,
           CH_assert(0 < a_refRatio[2]);)

    DisjointBoxLayout coarsenedFineGrids;
    coarsen(coarsenedFineGrids, a_fineGrids, a_refRatio);

    if (m_coarsenedFineData != NULL) {
        delete m_coarsenedFineData;
    }
    m_coarsenedFineData = new LevelData<FArrayBox>;
    m_coarsenedFineData->define(coarsenedFineGrids, a_numComps);

    m_refRatio = a_refRatio;
    m_isCopierDefined = false;
    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Full define constructor -- caches a Copier.
// -----------------------------------------------------------------------------
void MappedCoarseAverage::define (const DisjointBoxLayout& a_fineGrids,
                                  const DisjointBoxLayout& a_crseGrids,
                                  const int                a_numComps,
                                  const IntVect&           a_refRatio,
                                  const IntVect&           a_ghostVect)
{
    CH_TIME("MappedCoarseAverage::define");
    CH_assert(0 < a_numComps);
    D_TERM(CH_assert(0 < a_refRatio[0]);,
           CH_assert(0 < a_refRatio[1]);,
           CH_assert(0 < a_refRatio[2]);)

    DisjointBoxLayout coarsenedFineGrids;
    coarsen(coarsenedFineGrids, a_fineGrids, a_refRatio);

    if (m_coarsenedFineData != NULL) {
        delete m_coarsenedFineData;
    }
    m_coarsenedFineData = new LevelData<FArrayBox>;
    m_coarsenedFineData->define(coarsenedFineGrids, a_numComps);

    m_copier.define(coarsenedFineGrids, a_crseGrids); // Don't coarsen ghosts.

    m_refRatio = a_refRatio;
    m_isCopierDefined = true;
    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Is this object in a useable state?
// -----------------------------------------------------------------------------
bool MappedCoarseAverage::isDefined () const
{
    return m_isDefined;
}


// -----------------------------------------------------------------------------
// Do the averaging
// -----------------------------------------------------------------------------
void MappedCoarseAverage::averageToCoarse (LevelData<FArrayBox>&       a_crseData,
                                           const LevelData<FArrayBox>& a_fineData,
                                           const LevelGeometry* const  a_fineLevGeoPtr,
                                           const bool                  a_considerCellVol)
{
    this->computeAverages(a_crseData, a_fineData, ARITHMETIC,
                          a_fineLevGeoPtr, a_considerCellVol);
}


// -----------------------------------------------------------------------------
// Do the harmonic averaging
// -----------------------------------------------------------------------------
void MappedCoarseAverage::averageToCoarseHarmonic (LevelData<FArrayBox>&       a_crseData,
                                                   const LevelData<FArrayBox>& a_fineData,
                                                   const LevelGeometry* const  a_fineLevGeoPtr,
                                                   const bool                  a_considerCellVol)
{
    this->computeAverages(a_crseData, a_fineData, HARMONIC,
                          a_fineLevGeoPtr, a_considerCellVol);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void MappedCoarseAverage::computeAverages (LevelData<FArrayBox>&       a_crseData,
                                           const LevelData<FArrayBox>& a_fineData,
                                           const AverageType           a_averageType,
                                           const LevelGeometry* const  a_fineLevGeoPtr,
                                           const bool                  a_considerCellVol)
{
    CH_TIME("MappedCoarseAverage::computeAverages");
    CH_assert(m_isDefined);

    // Perform the averaging.
    DataIterator dit = a_fineData.boxLayout().dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
        BaseFab<Real>& coarsenedFineFAB = (*m_coarsenedFineData)[dit];
        const BaseFab<Real>& fineFAB = a_fineData[dit];

        this->averageGridData(coarsenedFineFAB,
                              fineFAB,
                              dit(),
                              a_averageType,
                              a_fineLevGeoPtr,
                              a_considerCellVol);
    }

    // Copy data to the user's holder.
    if (m_isCopierDefined) {
        // Use our cached copier to make things faster.
        m_coarsenedFineData->copyTo(m_coarsenedFineData->interval(),
                                    a_crseData,
                                    a_crseData.interval(),
                                    m_copier);
    } else {
        // We did not cache a copier.
        m_coarsenedFineData->copyTo(m_coarsenedFineData->interval(),
                                    a_crseData,
                                    a_crseData.interval());
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void MappedCoarseAverage::averageGridData (BaseFab<Real>&             a_coarse,
                                           const BaseFab<Real>&       a_fine,
                                           const DataIndex&           a_fineDI,
                                           const AverageType          a_averageType,
                                           const LevelGeometry* const a_fineLevGeoPtr,
                                           const bool                 a_considerCellVol) const
{
    CH_TIME("MappedCoarseAverage::averageGridData");
    CH_assert(m_isDefined);

    // Create the averaging regions.
    // refBox contains the offsets used to average each cell.
    const Box& destBox = a_coarse.box();
    Box refBox(IntVect::Zero, (m_refRatio-IntVect::Unit));

    // This is the region of fine data required to fill a_coarse.
    const IntVect srcRegionSmallEnd = destBox.smallEnd() * m_refRatio + refBox.smallEnd();
    const IntVect srcRegionBigEnd = destBox.bigEnd() * m_refRatio + refBox.bigEnd();
    const Box srcRegion(srcRegionSmallEnd, srcRegionBigEnd);

    // Make sure our source data is defined where we need it.
    CH_assert(a_fine.box().contains(srcRegion));

    if (a_considerCellVol) {
        // Gather metric data.
        CH_assert(a_fineLevGeoPtr != NULL);
        CH_assert(a_fineLevGeoPtr->isDefined());

        bool useCache = false;
        if (!a_fineLevGeoPtr->getCCJPtr().isNull()) {
            const LevelData<FArrayBox>& fineCCJ = a_fineLevGeoPtr->getCCJ();
            if (fineCCJ.getBoxes().check(a_fineDI)) {
                if (fineCCJ[a_fineDI].box().contains(srcRegion)) {
                    useCache = true;
                }
            }
        }

        if (useCache) {
            // Perform the weighted averaging using the cached J.
            const LevelData<FArrayBox>& fineCCJ = a_fineLevGeoPtr->getCCJ();

            if (a_averageType == ARITHMETIC) {
                FORT_MAPPEDAVERAGE(CHF_FRA(a_coarse),
                                   CHF_CONST_FRA(a_fine),
                                   CHF_CONST_FRA1(fineCCJ[a_fineDI],0),
                                   CHF_BOX(destBox),
                                   CHF_CONST_INTVECT(m_refRatio),
                                   CHF_BOX(refBox));

            } else if (a_averageType == HARMONIC) {
                MayDay::Error("MappedCoarseAverage::averageGridData -- Mapped harmonic averaging not complete");
            } else {
                MayDay::Error("MappedCoarseAverage::averageGridData -- bad averageType");
            }

        } else {
            // Fill a data holder with J.
            FArrayBox thisFineCCJ(srcRegion, 1);
            a_fineLevGeoPtr->fill_J(thisFineCCJ);

            // Perform the weighted averaging.
            if (a_averageType == ARITHMETIC) {
                FORT_MAPPEDAVERAGE(CHF_FRA(a_coarse),
                                   CHF_CONST_FRA(a_fine),
                                   CHF_CONST_FRA1(thisFineCCJ,0),
                                   CHF_BOX(destBox),
                                   CHF_CONST_INTVECT(m_refRatio),
                                   CHF_BOX(refBox));

            } else if (a_averageType == HARMONIC) {
                MayDay::Error("MappedCoarseAverage::averageGridData -- Mapped harmonic averaging not complete");
            } else {
                MayDay::Error("MappedCoarseAverage::averageGridData -- bad averageType");
            }
        }

    } else {
        // No need for metric data here. Just perform the unweighted averaging.
        if (a_averageType == ARITHMETIC) {
            FORT_UNMAPPEDAVERAGE(CHF_FRA(a_coarse),
                                 CHF_CONST_FRA(a_fine),
                                 CHF_BOX(destBox),
                                 CHF_CONST_INTVECT(m_refRatio),
                                 CHF_BOX(refBox));

        } else if (a_averageType == HARMONIC) {
            FORT_UNMAPPEDAVERAGEHARMONIC(CHF_FRA(a_coarse),
                                         CHF_CONST_FRA(a_fine),
                                         CHF_BOX(destBox),
                                         CHF_CONST_INTVECT(m_refRatio),
                                         CHF_BOX(refBox));

        } else {
            MayDay::Error("MappedCoarseAverage::averageGridData -- bad averageType");
        }
    }
}




// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
MappedCoarseAverageFace::MappedCoarseAverageFace ()
: m_isDefined(false),
  m_coarsenedFineData(NULL)
{;}


// -----------------------------------------------------------------------------
// Full constructor
// -----------------------------------------------------------------------------
MappedCoarseAverageFace::MappedCoarseAverageFace (const DisjointBoxLayout& a_fineGrids,
                                                  const int                a_numComps,
                                                  const IntVect&           a_refRatio)
: m_isDefined(false),
  m_coarsenedFineData(NULL)
{
    this->define(a_fineGrids, a_numComps, a_refRatio);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
MappedCoarseAverageFace::~MappedCoarseAverageFace ()
{
    delete m_coarsenedFineData;
    m_coarsenedFineData = NULL;
}


// -----------------------------------------------------------------------------
// Full define constructor
// -----------------------------------------------------------------------------
void MappedCoarseAverageFace::define (const DisjointBoxLayout& a_fineGrids,
                                      const int                a_numComps,
                                      const IntVect&           a_refRatio)
{
    CH_assert(0 < a_numComps);
    D_TERM(CH_assert(0 < a_refRatio[0]);,
           CH_assert(0 < a_refRatio[1]);,
           CH_assert(0 < a_refRatio[2]);)

    DisjointBoxLayout coarsened_fine_domain;
    coarsen(coarsened_fine_domain, a_fineGrids, a_refRatio);

    m_coarsenedFineData = new LevelData<FluxBox>;
    m_coarsenedFineData->define(coarsened_fine_domain, a_numComps);

    m_refRatio = a_refRatio;
    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Is this object in a useable state?
// -----------------------------------------------------------------------------
bool MappedCoarseAverageFace::isDefined () const
{
    return m_isDefined;
}


// -----------------------------------------------------------------------------
// Averages fine-level data to coarse level
// -----------------------------------------------------------------------------
void MappedCoarseAverageFace::averageToCoarse (LevelData<FluxBox>&        a_coarseData,
                                               const LevelData<FluxBox>&  a_fineData,
                                               const LevelGeometry* const a_levGeoPtr,
                                               const bool                 a_considerCellVol)
{
    // Put average into internal holder.
    this->computeAverages(a_coarseData, a_fineData, ARITHMETIC, a_levGeoPtr, a_considerCellVol);
}


// -----------------------------------------------------------------------------
// Averages fine-level data to coarse level using harmonic averaging
// -----------------------------------------------------------------------------
void MappedCoarseAverageFace::averageToCoarseHarmonic (LevelData<FluxBox>&        a_coarseData,
                                                       const LevelData<FluxBox>&  a_fineData,
                                                       const LevelGeometry* const a_levGeoPtr,
                                                       const bool                 a_considerCellVol)
{
    // Put average into internal holder.
    this->computeAverages(a_coarseData, a_fineData, HARMONIC, a_levGeoPtr, a_considerCellVol);
}


// -----------------------------------------------------------------------------
// Utility for averaging fine-level data to internal coarse representation
// of the fine grid. Sum of fine values is divided by
// a_refFactor^(CH_SPACEDIM-1). Normally this is the refinement ratio.
// -----------------------------------------------------------------------------
void MappedCoarseAverageFace::computeAverages (LevelData<FluxBox>&        a_coarseData,
                                               const LevelData<FluxBox>&  a_fineData,
                                               const AverageType          a_averageType,
                                               const LevelGeometry* const a_fineLevGeoPtr,
                                               const bool                 a_considerCellVol)
{
    CH_TIME("MappedCoarseAverageFace::computeAverages");
    CH_assert(isDefined());

    DataIterator dit = a_fineData.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        FluxBox& coarsenedFine = (*m_coarsenedFineData)[dit];
        const FluxBox& fine = a_fineData[dit];

        // coarsen from the entire fine grid onto the entire coarse grid
        this->averageGridData(coarsenedFine,
                              fine,
                              dit(),
                              a_averageType,
                              a_fineLevGeoPtr,
                              a_considerCellVol);
    }

    // If coarseData's DisjointBoxLayout is not a simple coarsenening of
    // the fine one, then it needs to have at least one ghost cell in
    // order to ensure that this copyTo is done correctly. In
    // particular, this is required in order to ensure that we handle
    // the case where the coarse-fine interface is coincident with a
    // coarse-coarse boundary. The other solution to this would be to
    // build a specialized Copier for LevelData<FluxBox>, but we're
    // hoping to avoid that for now...
    if ((a_coarseData.ghostVect() == IntVect::Zero) && !(a_coarseData.getBoxes().compatible(m_coarsenedFineData->getBoxes()))) {
        MayDay::Error("MappedCoarseAverageFace requires that coarse data which is not a coarsenening of the fine grids have at least one ghost cell");
    }

    // Copy average from internal holder to the user's holder.
    m_coarsenedFineData->copyTo(m_coarsenedFineData->interval(),
                                a_coarseData,
                                a_coarseData.interval());
}


// -----------------------------------------------------------------------------
// Averages single grid data from fine->crse.
// Sum of fine values is divided by a_refFactor^(CH_SPACEDIM-1).
// Normally this is the refinement ratio.
// -----------------------------------------------------------------------------
void MappedCoarseAverageFace::averageGridData (FluxBox&                   a_coarsenedFine,
                                               const FluxBox&             a_fine,
                                               const DataIndex&           a_fineDI,
                                               const AverageType          a_averageType,
                                               const LevelGeometry* const a_fineLevGeoPtr,
                                               const bool                 a_considerCellVol) const
{

    for (int dir = 0; dir < SpaceDim; ++dir) {
        FArrayBox& coarseFAB = a_coarsenedFine[dir];
        const FArrayBox& fineFAB = a_fine[dir];
        const Box& destBox = coarseFAB.box();
        const IntVect& boxType = destBox.type();

        // Set up refinement box remembering that we
        // don't want to index at all in dir direction --
        // instead, want to just march along face.
        IntVect hiVect(D_DECL(m_refRatio[0]-1,
                              m_refRatio[1]-1,
                              m_refRatio[2]-1));
        hiVect[dir] = 0;
        const Box refBox(IntVect::Zero, hiVect, boxType);

        // This is the region of fine data required to fill a_coarse.
        const IntVect srcRegionSmallEnd = destBox.smallEnd() * m_refRatio + refBox.smallEnd();
        const IntVect srcRegionBigEnd = destBox.bigEnd() * m_refRatio + refBox.bigEnd();
        const Box srcRegion(srcRegionSmallEnd, srcRegionBigEnd, boxType);

        // Make sure our source data is defined where we need it.
        CH_assert(fineFAB.box().contains(srcRegion));

        if (a_considerCellVol) {
            // Gather metric data.
            CH_assert(a_fineLevGeoPtr != NULL);
            CH_assert(a_fineLevGeoPtr->isDefined());
            FArrayBox thisFineFCJ(srcRegion, 1);
            a_fineLevGeoPtr->fill_J(thisFineFCJ);

            // Perform the weighted averaging.
            if (a_averageType == ARITHMETIC) {
                FORT_MAPPEDAVERAGEFACE(CHF_FRA(coarseFAB),
                                       CHF_CONST_FRA(fineFAB),
                                       CHF_CONST_FRA1(thisFineFCJ,0),
                                       CHF_BOX(destBox),
                                       CHF_CONST_INT(dir),
                                       CHF_CONST_INTVECT(m_refRatio),
                                       CHF_BOX(refBox));

            } else if (a_averageType == HARMONIC) {
                MayDay::Error("MappedCoarseAverageFace::averageGridData -- Mapped harmonic averaging not complete");
            } else {
                MayDay::Error("MappedCoarseAverageFace::averageGridData -- bad averageType");
            }
        } else {
            // No need for metric data here. Just perform the unweighted averaging.
            if (a_averageType == ARITHMETIC) {
                FORT_UNMAPPEDAVERAGEFACE(CHF_FRA(coarseFAB),
                                         CHF_CONST_FRA(fineFAB),
                                         CHF_BOX(destBox),
                                         CHF_CONST_INT(dir),
                                         CHF_CONST_INTVECT(m_refRatio),
                                         CHF_BOX(refBox));
            } else if (a_averageType == HARMONIC) {
                FORT_UNMAPPEDAVERAGEFACEHARMONIC(CHF_FRA(coarseFAB),
                                                 CHF_CONST_FRA(fineFAB),
                                                 CHF_BOX(destBox),
                                                 CHF_CONST_INT(dir),
                                                 CHF_CONST_INTVECT(m_refRatio),
                                                 CHF_BOX(refBox));
            } else {
                MayDay::Error("MappedCoarseAverageFace::averageGridData -- bad averageType");
            }
        }
    }
}

