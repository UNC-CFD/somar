#include "MappedFineInterp.H"
#include "InterpF_F.H"     // We can get away with most of Chombo's interp utils.
#include "MappedFineInterpF_F.H"
#include "MappedCoarseAverageF_F.H"
#include "AnisotropicRefinementTools.H"
#include "LevelGeometry.H"
#include "SetValLevel.H"
#include "Constants.H"


// -----------------------------------------------------------------------------
// Default constructor. User must subsequently call define().
// -----------------------------------------------------------------------------
MappedFineInterp::MappedFineInterp ()
: m_isDefined(false),
  m_fineLevGeoPtr(NULL),
  m_considerCellVol(false)
{;}


// -----------------------------------------------------------------------------
// Defining constructor. Constructs a valid object.
// Equivalent to default construction followed by define().
// -----------------------------------------------------------------------------
MappedFineInterp::MappedFineInterp (const DisjointBoxLayout&   a_fineGrids,
                                    const int                  a_numComps,
                                    const IntVect&             a_refRatio,
                                    const Box&                 a_fineDomainBox,
                                    const LevelGeometry* const a_fineLevGeoPtr,
                                    const bool                 a_considerCellVol)
: m_isDefined(false),
  m_fineLevGeoPtr(NULL),
  m_considerCellVol(false)
{
    ProblemDomain fineDomain(a_fineDomainBox);
    this->define(a_fineGrids, a_numComps, a_refRatio, fineDomain,
                 a_fineLevGeoPtr, a_considerCellVol);
}


// -----------------------------------------------------------------------------
// Defining constructor. Constructs a valid object.
// Equivalent to default construction followed by define().
// -----------------------------------------------------------------------------
MappedFineInterp::MappedFineInterp (const DisjointBoxLayout&   a_fineGrids,
                                    const int                  a_numComps,
                                    const IntVect&             a_refRatio,
                                    const ProblemDomain&       a_fineDomain,
                                    const LevelGeometry* const a_fineLevGeoPtr,
                                    const bool                 a_considerCellVol)
: m_isDefined(false),
  m_fineLevGeoPtr(NULL),
  m_considerCellVol(false)
{
    this->define(a_fineGrids, a_numComps, a_refRatio, a_fineDomain,
                 a_fineLevGeoPtr, a_considerCellVol);
}


// -----------------------------------------------------------------------------
// Destructor.
// -----------------------------------------------------------------------------
MappedFineInterp::~MappedFineInterp ()
{;}


// -----------------------------------------------------------------------------
// Defines this object. Existing information is overriden.
// -----------------------------------------------------------------------------
void MappedFineInterp::define (const DisjointBoxLayout&   a_fineGrids,
                               const int                  a_numComps,
                               const IntVect&             a_refRatio,
                               const Box&                 a_fineDomainBox,
                               const LevelGeometry* const a_fineLevGeoPtr,
                               const bool                 a_considerCellVol)
{
    ProblemDomain fineDomain(a_fineDomainBox);
    this->define(a_fineGrids, a_numComps, a_refRatio, fineDomain,
                 a_fineLevGeoPtr, a_considerCellVol);
}


// -----------------------------------------------------------------------------
// Defines this object. Existing information is overriden.
// -----------------------------------------------------------------------------
void MappedFineInterp::define (const DisjointBoxLayout&   a_fineGrids,
                               const int                  a_numComps,
                               const IntVect&             a_refRatio,
                               const ProblemDomain&       a_fineDomain,
                               const LevelGeometry* const a_fineLevGeoPtr,
                               const bool                 a_considerCellVol)
{
    CH_TIME("MappedFineInterp::define");
    CH_assert(a_fineGrids.checkPeriodic(a_fineDomain));
    CH_assert(!a_considerCellVol || a_fineLevGeoPtr);

    m_fineLevGeoPtr = a_fineLevGeoPtr;
    m_considerCellVol = a_considerCellVol;
    m_refRatio = a_refRatio;
    m_coarseProblemDomain = coarsen(a_fineDomain, a_refRatio);

    // Create the work array
    DisjointBoxLayout coarsenedFineGrids;
    coarsen(coarsenedFineGrids, a_fineGrids, a_refRatio);
    m_coarsenedFineData.define(coarsenedFineGrids, a_numComps, IntVect::Unit);

    // Redefine J if necessary.
    m_coarseJ.clear();

    if (a_considerCellVol) {
        // Make sure the fine LevGeo is well defined.
        CH_assert(a_fineLevGeoPtr != NULL);
        CH_assert(a_fineLevGeoPtr->getDomain() == a_fineDomain);

        // refBox contains the offsets used to average each cell.
        const Box refBox(IntVect::Zero, (a_refRatio-IntVect::Unit));

        // Define the coarse J holder.
        m_coarseJ.define(coarsenedFineGrids, 1, IntVect::Unit);
#ifndef NDEBUG
        setValLevel(m_coarseJ, quietNAN);
#endif

        // Fill the coarse J holder.
        DataIterator dit = m_coarseJ.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            // Reference the destination holder and region.
            FArrayBox& crseJFAB = m_coarseJ[dit];
            const Box& destBox = crseJFAB.box();

            // Create the source region and holder.
            Box srcBox = destBox;
            srcBox.refine(a_refRatio);
            srcBox.grow(1);
            FArrayBox fineJFAB(srcBox, 1);

            // Fill the source with the fine level J.
            a_fineLevGeoPtr->fill_J(fineJFAB);

            // Fill the dest with an average of the fine level J.
            FORT_UNMAPPEDAVERAGE(CHF_FRA(crseJFAB),
                                 CHF_CONST_FRA(fineJFAB),
                                 CHF_BOX(destBox),
                                 CHF_CONST_INTVECT(a_refRatio),
                                 CHF_BOX(refBox));
        }
    } // end if considering cell volumes

    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Returns true if this object was created with the defining
// constructor or if define() has been called.
// -----------------------------------------------------------------------------
bool MappedFineInterp::isDefined () const
{
    return m_isDefined;
}


// -----------------------------------------------------------------------------
// Replaces a_fineData with data interpolated from a_crseData. It
// is an error to call if not this->isDefined().  The domain of
// a_fineData should be the same as the fine domain specified in the
// most recent call to define().  It is expected that the coarse and
// fine level's domains are properly nested.  Both a_crseData and
// a_fineData should have the same number of components specified in
// the most recent call to define().
//
// a_fineData (modified): fine data.
// a_crseData (not modified): coarse data.
// a_averageFromDest: if true, first average data from a_fineData down
//                    to the resolution of a_crseData, then interp
//                    everything back up -- necessary when the coarse
//                    grids don't cover the fine grid (i.e when flattening
//                    an AMR hierarchy to a single resolution). Default is
//                    false.
// -----------------------------------------------------------------------------
void MappedFineInterp::interpToFine (LevelData<FArrayBox>&       a_fineData,
                                     const LevelData<FArrayBox>& a_crseData,
                                     const bool                  a_averageFromDest,
                                     const bool                  a_pwcAtBoxBoundaries)
{
    CH_TIME("MappedFineInterp::interpToFine");
    CH_assert(m_isDefined);

#ifndef NDEBUG
    setValLevel(m_coarsenedFineData, quietNAN);
#endif

    if (a_averageFromDest) {
        // Average down fine data -- this is a local operation
        DataIterator dit = a_fineData.dataIterator();
        for (dit.begin(); dit.ok(); ++dit) {
            // Create references for convenience.
            FArrayBox& crseFineFAB = m_coarsenedFineData[dit];
            const FArrayBox& fineFAB = a_fineData[dit];

            // The data destination box and averaging offsets.
            const Box& destBox = m_coarsenedFineData.getBoxes()[dit];
            Box refbox(IntVect::Zero, (m_refRatio-1)*IntVect::Unit);

            // Average appropriately.
            if (m_considerCellVol) {
                MayDay::Error("MappedFineInterp::interpToFine cannot yet consider cell vols and average from dest");
            } else {
                FORT_UNMAPPEDAVERAGE(CHF_FRA(crseFineFAB),
                                     CHF_CONST_FRA(fineFAB),
                                     CHF_BOX(destBox),
                                     CHF_CONST_INTVECT(m_refRatio),
                                     CHF_BOX(refbox));
            }
        }
    } // end if averaging from dest.

    // This should handle all the periodic BCs as well,
    // by filling in the ghost cells in an appropriate way
    a_crseData.copyTo(a_crseData.interval(),
                      m_coarsenedFineData,
                      m_coarsenedFineData.interval());

    ProblemDomain tmpCrseDomain;
    const BoxLayout& fineLayout = a_fineData.boxLayout();
    DataIterator dit = fineLayout.dataIterator();

    CH_assert(m_coarsenedFineData.getBoxes(). compatible(a_fineData.getBoxes()));

    for (dit.begin(); dit.ok(); ++dit) {
        // Create references for convenience
        FArrayBox& crseFineFAB = m_coarsenedFineData[dit];
        FArrayBox& fineFAB = a_fineData[dit];
        const Box& srcBox = m_coarsenedFineData.getBoxes()[dit];

        const ProblemDomain* crseDomPtr = &m_coarseProblemDomain;
        if (a_pwcAtBoxBoundaries) {
            tmpCrseDomain.define(srcBox);
            crseDomPtr = &tmpCrseDomain;
        }

        // Do we need to consider cell volumes?
        if (m_considerCellVol) {
            // Apply the volume weighting.
            for (int comp = 0; comp < crseFineFAB.nComp(); ++comp) {
                crseFineFAB.mult(m_coarseJ[dit], crseFineFAB.box(), 0, comp, 1);
            }


            // Interpolate to fine grid.
            this->interpGridData(fineFAB,
                                 crseFineFAB,
                                 srcBox,
                                 m_refRatio,
                                 *crseDomPtr);

            // Remove the volume weighting.
            // Using interpGridData ensures J is refined EXACTLY as a_crseData.
            FArrayBox fineJ(fineFAB.box(), 1);
            fineJ.setVal(quietNAN);
            this->interpGridData(fineJ,
                                 m_coarseJ[dit],
                                 srcBox,
                                 m_refRatio,
                                 *crseDomPtr);
            for (int comp = 0; comp < fineFAB.nComp(); ++comp) {
                fineFAB.divide(fineJ, fineFAB.box(), 0, comp, 1);
            }


        } else {
            // Don't consider cell volumes. Just perform the interpolation
            // on this grid box.
            this->interpGridData(fineFAB,
                                 crseFineFAB,
                                 srcBox,
                                 m_refRatio,
                                 *crseDomPtr);
        }
    }
}


// -----------------------------------------------------------------------------
// Just do piecewise-constant interpolation.
// -----------------------------------------------------------------------------
void MappedFineInterp::pwcinterpToFine (LevelData<FArrayBox>&       a_fineData,
                                        const LevelData<FArrayBox>& a_crseData,
                                        const bool                  a_averageFromDest)
{
    CH_TIME("MappedFineInterp::pwcinterpToFine");
    CH_assert(m_isDefined);

    if (a_averageFromDest) {
        // average down fine data -- this is a local operation
        DataIterator dit = a_fineData.dataIterator();
        for (dit.begin(); dit.ok(); ++dit) {
            FArrayBox& fineFab = a_fineData[dit];
            FArrayBox& crseFab = m_coarsenedFineData[dit];
            const Box& crseBox = m_coarsenedFineData.getBoxes()[dit];

            const Box refbox(IntVect::Zero, m_refRatio - IntVect::Unit);

            FORT_UNMAPPEDAVERAGE(CHF_FRA(crseFab),
                                 CHF_CONST_FRA(fineFab),
                                 CHF_BOX(crseBox),
                                 CHF_CONST_INTVECT(m_refRatio),
                                 CHF_BOX(refbox));
        }
    }

    // this should handle all the periodic BCs as well,
    // by filling in the ghost cells in an appropriate way
    a_crseData.copyTo(a_crseData.interval(),
                      m_coarsenedFineData,
                      m_coarsenedFineData.interval() );

    const BoxLayout fine_domain = a_fineData.boxLayout();
    DataIterator dit = fine_domain.dataIterator();

    for (dit.begin(); dit.ok(); ++dit) {
        const BaseFab<Real>& coarsened_fine = m_coarsenedFineData[dit()];
        const Box& coarsened_fine_box = m_coarsenedFineData.getBoxes()[dit()];
        BaseFab<Real>& fine = a_fineData[dit()];
        // interpGridData interpolates from an entire coarse grid onto an
        // entire fine grid.
        this->pwcinterpGridData(fine,
                                coarsened_fine,
                                coarsened_fine_box,
                                m_refRatio);
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void MappedFineInterp::pwcinterpGridData (BaseFab<Real>&       a_fine,
                                          const BaseFab<Real>& a_coarse,
                                          const Box&           a_coarsened_fine_box,
                                          const IntVect&       a_refRatio)
{
    CH_TIME("MappedFineInterp::pwcinterpGridData");
    // fill fine data with piecewise constant coarse data
    const Box& b = a_coarsened_fine_box;
    const Box refbox(IntVect::Zero, a_refRatio - IntVect::Unit);

    FORT_UNMAPPEDINTERPCONSTANT(CHF_FRA(a_fine),
                                CHF_CONST_FRA(a_coarse),
                                CHF_BOX(b),
                                CHF_CONST_INTVECT(a_refRatio),
                                CHF_BOX(refbox));
}


// -----------------------------------------------------------------------------
// Interpolate from fine grid to coarse grid.
// Prerequisite: coarsened.box contains coarsen(fine.box).
//
// Uses piecewise bilinear interpolation with multidimensional-limited
// slopes. See design document for details.
// -----------------------------------------------------------------------------
void MappedFineInterp::interpGridData (BaseFab<Real>&       a_fine,
                                       const BaseFab<Real>& a_coarse,
                                       const Box&           a_coarsened_fine_box,
                                       const IntVect&       a_refRatio,
                                       const ProblemDomain& a_coarseDomain)
{
    CH_TIME("MappedFineInterp::interpGridData");

    // fill fine data with piecewise constant coarse data
    const Box& b = a_coarsened_fine_box;
    const int num_comp = a_fine.nComp();
    const Box refbox(IntVect::Zero, (a_refRatio-IntVect::Unit));

    // This is used to effectively remove the non-contributing directions
    // from the interpolation. The non-contributors are the directions
    // in which the refRatio is 1.
    //
    // NOTE: To bring this back to the original, full dimensional interpolation,
    // just set refMask to IntVect::Unit;
    // const IntVect& refMask = IntVect::Unit;
    const IntVect refMask(D_DECL(
        (a_refRatio[0] > 1)? 1: 0,
        (a_refRatio[1] > 1)? 1: 0,
        (a_refRatio[2] > 1)? 1: 0));

    // Start with a piecewise constant interpolation (the prediction).
    FORT_UNMAPPEDINTERPCONSTANT(CHF_FRA(a_fine),
                                CHF_CONST_FRA(a_coarse),
                                CHF_BOX(b),
                                CHF_CONST_INTVECT(a_refRatio),
                                CHF_BOX(refbox));

    // hardwired to 3 due to lack of variable number of arguments in chfpp
    BaseFab<Real> slopes[3];
    for (int dir = 0; dir < 3; ++dir) {
        BaseFab<Real>& dir_slope = slopes[dir];
        dir_slope.resize(b, num_comp);
    }

    for (int dir = 0; dir < SpaceDim; ++dir) {
        BaseFab<Real>& dir_slope = slopes[dir];

        if (refMask[dir] > 0) {
            const Box interiorBox = grow(a_coarseDomain.domainBox(), -BASISV(dir));

            const Box bcenter = interiorBox & b;
            if (!bcenter.isEmpty()) {
                FORT_INTERPCENTRALSLOPE(CHF_FRA(dir_slope),
                                        CHF_CONST_FRA(a_coarse),
                                        CHF_BOX(bcenter),
                                        CHF_CONST_INT(dir));
            }
            const Box blo = b & adjCellLo(interiorBox, dir, 1);
            if (!blo.isEmpty()) {
                FORT_INTERPHISIDESLOPE(CHF_FRA(dir_slope),
                                       CHF_CONST_FRA(a_coarse),
                                       CHF_BOX(blo),
                                       CHF_CONST_INT(dir));
            }
            const Box bhi = b & adjCellHi(interiorBox, dir, 1);
            if (!bhi.isEmpty()) {
                FORT_INTERPLOSIDESLOPE(CHF_FRA(dir_slope),
                                       CHF_CONST_FRA(a_coarse),
                                       CHF_BOX(bhi),
                                       CHF_CONST_INT(dir));
            }
        } else {
            // Unrefined directions should not contribute. The entire
            // refinement process should reduce in dimensionality.
            dir_slope.setVal(0.0);
        }
    }

    // To do limits, we need to have a box which includes
    // the neighbors of a given point (to check for the
    // local maximum...
    Box neighborBox(-1*refMask, refMask);

    // GHM 7/12/01
    // interplimit iterates over box b_mod (was b), but cells within
    // 1 of the physical boundary never enter result (and this
    // wasted calculation may call upon uninitialized memory).
    // DFM 10/8/01
    // note that this turns off slope limiting for cells adjacent to the
    // boundary -- may want to revisit this in the future
    Box b_mod(b);
    b_mod.grow(refMask);
    b_mod = a_coarseDomain & b_mod;
    b_mod.grow(-refMask);

    // create a box grown big enough to remove periodic BCs from domain
    Box domBox = grow(b, 2*refMask);
    domBox = a_coarseDomain & domBox;

    FORT_INTERPLIMIT(CHF_FRA(slopes[0]),
                     CHF_FRA(slopes[1]),
                     CHF_FRA(slopes[2]),
                     CHF_CONST_FRA(a_coarse),
                     CHF_BOX(b_mod),
                     CHF_BOX(neighborBox),
                     CHF_BOX(domBox));

    // Update to a linear interpolation (the correction).
    for (int dir = 0; dir < SpaceDim; ++dir) {
        // Don't do unrefined directions since it should add nothing anyway.
        if (a_refRatio[dir] == 1) continue;

        BaseFab<Real>& dir_slope = slopes[dir];

        FORT_UNMAPPEDINTERPLINEAR(CHF_FRA(a_fine),
                                  CHF_CONST_FRA(dir_slope),
                                  CHF_BOX(b),
                                  CHF_CONST_INT(dir),
                                  CHF_CONST_INTVECT(a_refRatio),
                                  CHF_BOX(refbox));
    }
}



// -----------------------------------------------------------------------------
// interpolate from coarse grid to fine grid.  prerequisite:
// coarsened.box contains coarsen(fine.box).
//
// uses piecewise bilinear interpolation with multidimensional-limited
// slopes.  see design document for details.
// -----------------------------------------------------------------------------
void MappedFineInterpFace::interpGridData (BaseFab<Real>&       a_fine,
                                           const BaseFab<Real>& a_coarse,
                                           const Box&           a_coarsened_fine_box,
                                           const ProblemDomain& a_crseDomain,
                                           const IntVect&       a_ref_ratio,
                                           const int            a_faceDir)
{
    // Fill fine data with piecewise constant coarse data
    const Box& b = a_coarsened_fine_box;
    const int num_comp = a_fine.nComp();

    // define refbox for edges
    Box facerefbox(IntVect::Zero, a_ref_ratio - IntVect::Unit);
    facerefbox.surroundingNodes(a_faceDir);
    facerefbox.setBig(a_faceDir, 0);

    FORT_UNMAPPEDINTERPFACECONSTANT(CHF_FRA(a_fine),
                                    CHF_CONST_FRA(a_coarse),
                                    CHF_BOX(b),
                                    CHF_CONST_INTVECT(a_ref_ratio),
                                    CHF_BOX(facerefbox),
                                    CHF_CONST_INT(a_faceDir));

    //  Tuple<BaseFab<Real>, SpaceDim> slopes;
    //  for (int dir = 0; dir < SpaceDim; ++dir)
    // hardwired to 3 due to lack of variable number of arguments in chfpp
    // this is designed to accomodate the Chombo max-dimension of 3
#define MAXDIM 3
    BaseFab<Real> slopes[MAXDIM];

    // this is a trick to make domain face-centered while still
    // respecting periodic BC's; take domain, then grow
    // by one (to allow for centered differences), then intersect
    // with domain.  this will result in a box that is larger than
    // the computational domain in periodic directions, but is
    // the same as the periodic domain in non-periodic directions.
    // so when we later intersect boxes with this domain, it will
    // cut off cells in non-periodic directions while leaving
    // periodic directions intact.
    Box domainFaceBox(a_crseDomain.domainBox());
    // using a grow radius of two here, but probably could have used 1
    domainFaceBox.grow(2);
    domainFaceBox &= a_crseDomain;
    domainFaceBox.surroundingNodes(a_faceDir);

    for (int dir = 0; dir < SpaceDim; ++dir) {
        BaseFab<Real>& dir_slope = slopes[dir];
        dir_slope.resize(b, num_comp);
    }
    // define the extras over a unit box in order to avoid issues with null
    // pointers and undefined FAB's
    for (int dir = SpaceDim; dir < MAXDIM; dir++) {
        Box unitBox(IntVect::Zero, IntVect::Zero);
        BaseFab<Real>& dir_slope = slopes[dir];
        dir_slope.resize(unitBox, num_comp);
    }

    for (int dir = 0; dir < SpaceDim; ++dir) {
        // only do slopes if not doing normal direction
        if (dir != a_faceDir) {
            BaseFab<Real>& dir_slope = slopes[dir];

            const Box bcenter = b & grow(domainFaceBox, -BASISV(dir));
            if (!bcenter.isEmpty()) {
                FORT_INTERPCENTRALSLOPE(CHF_FRA(dir_slope),
                                        CHF_CONST_FRA(a_coarse),
                                        CHF_BOX(bcenter),
                                        CHF_CONST_INT(dir));
            }
            const Box blo = b & adjCellLo(grow(domainFaceBox, -BASISV(dir)), dir);
            if (!blo.isEmpty()) {
                FORT_INTERPHISIDESLOPE(CHF_FRA(dir_slope),
                                       CHF_CONST_FRA(a_coarse),
                                       CHF_BOX(blo),
                                       CHF_CONST_INT(dir));
            }
            const Box bhi = b & adjCellHi(grow(domainFaceBox,
                                               -BASISV(dir)), dir);
            if (!bhi.isEmpty()) {
                FORT_INTERPLOSIDESLOPE(CHF_FRA(dir_slope),
                                       CHF_CONST_FRA(a_coarse),
                                       CHF_BOX(bhi),
                                       CHF_CONST_INT(dir));
            }
        } else {
            // if this is the normal direction, just fill with 0's
            slopes[dir].setVal(0.0);
        }
    } // end loop over directions

    // need a box which covers all neighbors in the plane of the
    // coarse face. this should result in a box that's (-1,1)
    // in the plane of the face, and (0,0) normal to the face)
    Box neighborBox(-1 * IntVect::Unit,
                    IntVect::Unit);
    neighborBox.grow(a_faceDir, -1);

    // this will be a box which defines which cells in a_coarse
    // are valid (or, at least, which ones exist!)

    // first construct box of valid edges in physical domain
    Box validBox(a_crseDomain.domainBox());
    validBox.grow(2);
    validBox &= a_crseDomain;
    validBox.surroundingNodes(a_faceDir);
    // now intersect with existing cells in a_coarse
    validBox &= a_coarse.box();


    FORT_UNMAPPEDINTERPLIMITFACE(CHF_FRA(slopes[0]),
                                 CHF_FRA(slopes[1]),
                                 CHF_FRA(slopes[2]),
                                 CHF_CONST_FRA(a_coarse),
                                 CHF_BOX(b),
                                 CHF_BOX(neighborBox),
                                 CHF_BOX(validBox),
                                 CHF_CONST_INT(a_faceDir));

    // do (bi)linear interpolation on fine faces which overlie coarse
    // faces

    // need a cell-centered coarse box
    Box coarseCellBox = b;
    coarseCellBox.enclosedCells(a_faceDir);

    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (dir != a_faceDir) {
            BaseFab<Real>& dir_slope = slopes[dir];
            FORT_UNMAPPEDINTERPLINEARFACE(CHF_FRA(a_fine),
                                          CHF_CONST_FRA(dir_slope),
                                          CHF_BOX(b ),
                                          CHF_CONST_INT(dir),
                                          CHF_CONST_INTVECT(a_ref_ratio),
                                          CHF_BOX(facerefbox));
        }
    }

    // finally, do interiors

    // this box will contain the interior faces which do
    // not overlie a coarse face
    Box interiorRefBox(IntVect::Zero, a_ref_ratio - IntVect::Unit);
    interiorRefBox.surroundingNodes(a_faceDir);
    // remove outer faces from interior box...
    interiorRefBox.grow(a_faceDir, -1);

    FORT_UNMAPPEDINTERPLINEARINTERIORFACE(CHF_FRA(a_fine),
                                          CHF_BOX(coarseCellBox),
                                          CHF_CONST_INTVECT(a_ref_ratio),
                                          CHF_CONST_INT(a_faceDir),
                                          CHF_BOX(interiorRefBox));
}
