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
#ifndef MappedCFInterp_H__INCLUDED__
#define MappedCFInterp_H__INCLUDED__

#include "LevelGeometry.H"
#include "CFRegion.H"


// -----------------------------------------------------------------------------
// Interpolates ALL fine ghosts using only coarse data. This is a modification
// of Chombo's QuadCFInterp class.
// -----------------------------------------------------------------------------
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
    if ((a > 0.0) && (b > 0.0) & (c > 0.0)) {
        return Min(a,Min(b,c));
    } else if ((a < 0.0) && (b < 0.0) & (c < 0.0)) {
        return Max(a,Max(b,c));
    }
    return 0.0;
}


#endif //!MappedCFInterp_H__INCLUDED__
