/*******************************************************************************
 *    FILE: MappedFineInterp.H
 *    DESC: Chombo's MappedFineInterp modified to handle anisotropic refinement and
 *          variable cell volumes.
 *    NOTE: Much of this code was taken from Chombo's FineInterp.
 *    DATE: Tue 04 Jun 2013 10:13:00 AM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __MappedFineInterp_H__INCLUDED__
#define __MappedFineInterp_H__INCLUDED__

#include "LevelData.H"
class DisjointBoxLayout;
class LevelGeometry;


// This class replaces data at a fine level of refinement with data
// interpolated from a coarser level of refinement.  Interpolation is
// piecewise bi(tri)linear, with van Leer slopes if there is room for
// the stencil, with lower-order slopes if there isn't.  See the
// design document.
class MappedFineInterp {
public:
    // Default constructor. User must subsequently call define().
    MappedFineInterp ();

    // Defining constructor. Constructs a valid object.
    // Equivalent to default construction followed by define().
    MappedFineInterp (const DisjointBoxLayout&    a_fineGrids,
                      const int                   a_numComps,
                      const IntVect&              a_refRatio,
                      const Box&                  a_fineDomainBox,
                      const LevelGeometry* const  a_fineLevGeoPtr = NULL,
                      const bool                  a_considerCellVol = false);

    // Defining constructor. Constructs a valid object.
    // Equivalent to default construction followed by define().
    MappedFineInterp (const DisjointBoxLayout&    a_fineGrids,
                      const int                   a_numComps,
                      const IntVect&              a_refRatio,
                      const ProblemDomain&        a_fineDomain,
                      const LevelGeometry* const  a_fineLevGeoPtr = NULL,
                      const bool                  a_considerCellVol = false);

    // Destructor.
    virtual ~MappedFineInterp ();

    // Defines this object. Existing information is overriden.
    virtual void define (const DisjointBoxLayout&    a_fineGrids,
                         const int                   a_numComps,
                         const IntVect&              a_refRatio,
                         const Box&                  a_fineDomainBox,
                         const LevelGeometry* const  a_fineLevGeoPtr = NULL,
                         const bool                  a_considerCellVol = false);

    // Defines this object. Existing information is overriden.
    virtual void define (const DisjointBoxLayout&    a_fineGrids,
                         const int                   a_numComps,
                         const IntVect&              a_refRatio,
                         const ProblemDomain&        a_fineDomain,
                         const LevelGeometry* const  a_fineLevGeoPtr = NULL,
                         const bool                  a_considerCellVol = false);

    // Returns true if this object was created with the defining
    // constructor or if define() has been called.
    virtual bool isDefined () const;

    // Replaces a_fine_data with data interpolated from a_coarse_data. It
    // is an error to call if not this->isDefined().  The domain of
    // a_fine_data should be the same as the fine domain specified in the
    // most recent call to define().  It is expected that the coarse and
    // fine level's domains are properly nested.  Both a_coarse_data and
    // a_fine_data should have the same number of components specified in
    // the most recent call to define().
    //
    // a_fine_data (modified): fine data.
    // a_coarse_data (not modified): coarse data.
    // a_averageFromDest: if true, first average data from a_fine_data down
    //                    to the resolution of a_coarse_data, then interp
    //                    everything back up -- necessary when the coarse
    //                    grids don't cover the fine grid (i.e when flattening
    //                    an AMR hierarchy to a single resolution). Default is
    //                    false.
    virtual void interpToFine (LevelData<FArrayBox>&       a_fineData,
                               const LevelData<FArrayBox>& a_coarseData,
                               const bool                  a_averageFromDest = false,
                               const bool                  a_pwcAtBoxBoundaries = false);

    // Just do piecewise-constant interpolation.
    virtual void pwcinterpToFine(LevelData<FArrayBox>&       a_fineData,
                                 const LevelData<FArrayBox>& a_coarseData,
                                 const bool                  a_averageFromDest = false);

    // Use these static utilities at your own risk!!!

    // Interpolate from fine grid to coarse grid.
    // Prerequisite: coarsened.box contains coarsen(fine.box).
    // Uses piecewise bilinear interpolation with multidimensional-limited
    // slopes. See design document for details.
    static void interpGridData (BaseFab<Real>&       a_fine,
                                const BaseFab<Real>& a_coarse,
                                const Box&           a_coarsenedFineBox,
                                const IntVect&       a_refRatio,
                                const ProblemDomain& a_coarseDomain);

    static void pwcinterpGridData(BaseFab<Real>&       a_fine,
                                  const BaseFab<Real>& a_coarse,
                                  const Box&           a_coarsenedFineBox,
                                  const IntVect&       a_refRatio);

protected:
    bool m_isDefined;

    IntVect m_refRatio;

    ProblemDomain m_coarseProblemDomain;

    // Work array for the coarse level data in a domain that is the
    // outline of the fine level domain on the coarse level
    LevelData<FArrayBox> m_coarsenedFineData;

    const LevelGeometry* m_fineLevGeoPtr;
    LevelData<FArrayBox> m_coarseJ;

    bool m_considerCellVol;
};



class MappedFineInterpFace {
public:
    // interpolate from coarse grid to fine grid.  prerequisite:
    // coarsened.box contains coarsen(fine.box).
    //
    // uses piecewise bilinear interpolation with multidimensional-limited
    // slopes.  see design document for details.
    static void interpGridData (BaseFab<Real>&       a_fine,
                                const BaseFab<Real>& a_coarse,
                                const Box&           a_coarsened_fine_box,
                                const ProblemDomain& a_crseDomain,
                                const IntVect&       a_ref_ratio,
                                const int            a_faceDir);
};


#endif //!__MappedFineInterp_H__INCLUDED__

