/*******************************************************************************
 *    FILE: MappedCoarseAverage.H
 *    DESC: Same as Chombo's CoarseAverage, but can account for varying cell
 *          areas and anisotropic refinement.
 *    NOTE: Much of this code was taken from Chombo's CoarseAverage.
 *    DATE: Thu 9 May 2013 10:02 AM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __MappedCoarseAverage_H__INCLUDED__
#define __MappedCoarseAverage_H__INCLUDED__


#include "REAL.H"
#include "BaseFab.H"
#include "FArrayBox.H"
#include "LevelData.H"
class DisjointBoxLayout;
class LevelGeometry;


// The averaging schemes
enum AverageType
{
    ARITHMETIC = 0,
    HARMONIC,
    NUM_AVERAGE_TYPES
};


// This class replaces data at a coarse level of refinement with an average of
// data at a finer level of refinement, in areas where fine data is present.
// Coarse level data is not modified where fine level data is not present.
class MappedCoarseAverage
{
public:
    // Constructor
    MappedCoarseAverage ();

    // Full constructor
    MappedCoarseAverage (const DisjointBoxLayout& a_fineGrids,
                         const int                a_numComps,
                         const IntVect&           a_refRatio);

    // Full constructor -- caches a Copier.
    MappedCoarseAverage (const DisjointBoxLayout& a_fineGrids,
                         const DisjointBoxLayout& a_crseGrids,
                         const int                a_numComps,
                         const IntVect&           a_refRatio,
                         const IntVect&           a_ghostVect = IntVect::Zero);

    // Destructor
    virtual ~MappedCoarseAverage ();

    // Full define constructor
    virtual void define (const DisjointBoxLayout& a_fineGrids,
                         const int                a_numComps,
                         const IntVect&           a_refRatio);

    // Full define constructor -- caches a Copier.
    virtual void define (const DisjointBoxLayout& a_fineGrids,
                         const DisjointBoxLayout& a_crseGrids,
                         const int                a_numComps,
                         const IntVect&           a_refRatio,
                         const IntVect&           a_ghostVect = IntVect::Zero);

    // Is this object in a usable state?
    virtual bool isDefined () const;

    // Do the averaging
    virtual void averageToCoarse (LevelData<FArrayBox>&       a_crseData,
                                  const LevelData<FArrayBox>& a_fineData,
                                  const LevelGeometry* const  a_fineLevGeoPtr = NULL,
                                  const bool                  a_considerCellVol = false);

    // Do the harmonic averaging
    virtual void averageToCoarseHarmonic (LevelData<FArrayBox>&       a_crseData,
                                          const LevelData<FArrayBox>& a_fineData,
                                          const LevelGeometry* const  a_fineLevGeoPtr = NULL,
                                          const bool                  a_considerCellVol = false);

protected:
    virtual void computeAverages (LevelData<FArrayBox>&       a_crseData,
                                  const LevelData<FArrayBox>& a_fineData,
                                  const AverageType           a_averageType,
                                  const LevelGeometry* const  a_fineLevGeoPtr = NULL,
                                  const bool                  a_considerCellVol = false);

    virtual void averageGridData (BaseFab<Real>&              a_coarse,
                                  const BaseFab<Real>&        a_fine,
                                  const DataIndex&            a_fineDI,
                                  const AverageType           a_averageType,
                                  const LevelGeometry* const  a_fineLevGeoPtr = NULL,
                                  const bool                  a_considerCellVol = false) const;

    // Did you call the full define constructor?
    bool m_isDefined;

    // Has a copier been defined to transfer data to coarse-grid layout?
    bool m_isCopierDefined;

    // Cached copier to handle transfer to coarse-grid layout.
    Copier m_copier;

    // The refinement ratio
    IntVect m_refRatio;

    // Work array for the coarsening of the fine data, of the same "shape"
    // as the fine data.
    LevelData<FArrayBox>* m_coarsenedFineData;
};



// This class replaces face-centered data on a coarse level of refinement
// with the average of the finer-level data which overlays the face.
// This class is similar to CoarseAverage.
class MappedCoarseAverageFace {
public:
    // Default constructor
    MappedCoarseAverageFace ();

    // Full constructor
    MappedCoarseAverageFace (const DisjointBoxLayout& a_fineGrids,
                             const int                a_numComps,
                             const IntVect&           a_refRatio);

    // Destructor
    ~MappedCoarseAverageFace ();

    // Full define constructor
    virtual void define (const DisjointBoxLayout& a_fineGrids,
                         const int                a_numComps,
                         const IntVect&           a_refRatio);

    // Is this object in a usable state?
    virtual bool isDefined () const;

    // Averages fine-level data to coarse level.
    virtual void averageToCoarse (LevelData<FluxBox>&        a_coarse_data,
                                  const LevelData<FluxBox>&  a_fine_data,
                                  const LevelGeometry* const a_levGeoPtr = NULL,
                                  const bool                 a_considerCellVol = false);

    // Averages fine-level data to coarse level using harmonic averaging.
    virtual void averageToCoarseHarmonic (LevelData<FluxBox>&        a_coarse_data,
                                          const LevelData<FluxBox>&  a_fine_data,
                                          const LevelGeometry* const a_levGeoPtr = NULL,
                                          const bool                 a_considerCellVol = false);

protected:
    // Utility for averaging fine-level data to internal coarse representation
    // of the fine grid. Sum of fine values is divided by
    // a_refFactor^(CH_SPACEDIM-1). Normally this is the refinement ratio.
    virtual void computeAverages (LevelData<FluxBox>&        a_coarseData,
                                  const LevelData<FluxBox>&  a_fineData,
                                  const AverageType          a_averageType,
                                  const LevelGeometry* const a_fineLevGeoPtr = NULL,
                                  const bool                 a_considerCellVol = false);

    // Averages single grid data from fine->crse.
    // Sum of fine values is divided by a_refFactor^(CH_SPACEDIM-1).
    // Normally this is the refinement ratio.
    virtual void averageGridData (FluxBox&                   a_coarsenedFine,
                                  const FluxBox&             a_fine,
                                  const DataIndex&           a_fineDI,
                                  const AverageType          a_averageType,
                                  const LevelGeometry* const a_fineLevGeoPtr = NULL,
                                  const bool                 a_considerCellVol = false) const;

    // Did you call the full define constructor?
    bool m_isDefined;

    // Was the average function called to set up the internal representation
    // of the coarse data?
    bool m_isAveraged;

    // The refinement ratio
    IntVect m_refRatio;

    // Work array for coarsening of fine data, same "shape" as fine data
    LevelData<FluxBox>* m_coarsenedFineData;
};



#endif //!__MappedCoarseAverage_H__INCLUDED__
