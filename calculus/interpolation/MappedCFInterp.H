/*******************************************************************************
 *    FILE: MappedCFInterp.H
 *    DESC: Interpolates ALL fine ghosts using only coarse data.
 *    DATE: Tue 22 Apr 2014 05:16:39 PM EDT  (UTC-0400)
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef MappedCFInterp_H__INCLUDED__
#define MappedCFInterp_H__INCLUDED__

#include "LevelGeometry.H"
#include "CFRegion.H"


class MappedCFInterp
{
public:
    // Default constructor
    MappedCFInterp ();

    // Full constructor
    MappedCFInterp (const LevelGeometry& a_levGeo);

    // Full constructor
    MappedCFInterp (const DisjointBoxLayout& a_grids,
                    const DisjointBoxLayout* a_crseGridsPtr,
                    const RealVect&          a_dx,
                    const IntVect&           a_refRatio);

    // Destructor
    virtual ~MappedCFInterp ();

    // Full define
    void define (const LevelGeometry& a_levGeo);

    // Full define
    void define (const DisjointBoxLayout& a_grids,
                 const DisjointBoxLayout* a_crseGridsPtr,
                 const RealVect&          a_dx,
                 const IntVect&           a_refRatio);

    // Frees memory, makes object unusable
    void undefine ();

    // Is this object ready for use?
    virtual inline bool isDefined () const;

    // Coarse / fine interpolation in time and space.
    virtual void coarseFineInterp (LevelData<FArrayBox>&       a_phif,
                                   const LevelData<FArrayBox>& a_oldPhic,
                                   const LevelData<FArrayBox>& a_newPhic,
                                   const Real                  a_timeInterpCoeff,
                                   const int                   a_order,
                                   const bool                  a_useLinearLimiting = false) const;

    // Coarse / fine interpolation in space only.
    virtual void coarseFineInterp (LevelData<FArrayBox>&       a_phif,
                                   const LevelData<FArrayBox>& a_phic,
                                   const int                   a_order,
                                   const bool                  a_useLinearLimiting = false) const;

protected:

    // Van-Leer minmod limiter
    static inline Real minmod (const Real a, const Real b, const Real c);

    bool m_isDefined;
    const LevelGeometry* m_levGeoPtr;
    const LevelGeometry* m_crseLevGeoPtr;
    bool m_deleteLevGeoPtrs;

    CFRegion          m_cfregion;
    DisjointBoxLayout m_crseGrids;

    ProblemDomain     m_domain;
    DisjointBoxLayout m_grids;
    IntVect           m_refRatio;

    IntVect           m_crseGhosts;
    Copier            m_crseCopier;
};


// -----------------------------------------------------------------------------
// Is this object ready for use?
// -----------------------------------------------------------------------------
bool MappedCFInterp::isDefined () const
{
    return m_isDefined;
}


// -----------------------------------------------------------------------------
// Van-Leer minmod limiter
// -----------------------------------------------------------------------------
Real MappedCFInterp::minmod (const Real a, const Real b, const Real c)
{
    if (a > 0.0 && b > 0.0 & c > 0.0) {
        return Min(a,Min(b,c));
    } else if (a < 0.0 && b < 0.0 & c < 0.0) {
        return Max(a,Max(b,c));
    }
    return 0.0;
}


#endif //!MappedCFInterp_H__INCLUDED__
