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
#include "MappedCFInterp.H"
#include "MappedFineInterpF_F.H"
#include "AnisotropicRefinementTools.H"
#include "BoxIterator.H"
#include "timeInterp.H"
#include "MappedQuadCFInterp.H"
#include "Printing.H"
#include "Constants.H"
#include "MiscUtils.H"
#include "CornerCopier.H"
#include "MappedCoarseAverage.H"


// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
MappedCFInterp::MappedCFInterp ()
: m_isDefined(false),
  m_levGeoPtr(NULL),
  m_crseLevGeoPtr(NULL),
  m_deleteLevGeoPtrs(false)
{;}


// -----------------------------------------------------------------------------
// Full constructor
// -----------------------------------------------------------------------------
MappedCFInterp::MappedCFInterp (const LevelGeometry& a_levGeo)
: m_isDefined(false),
  m_levGeoPtr(NULL),
  m_crseLevGeoPtr(NULL),
  m_deleteLevGeoPtrs(false)
{
    define(a_levGeo);
}


// -----------------------------------------------------------------------------
// Full constructor
// -----------------------------------------------------------------------------
MappedCFInterp::MappedCFInterp (const DisjointBoxLayout& a_grids,
                                const DisjointBoxLayout* a_crseGridsPtr,
                                const RealVect&          a_dx,
                                const IntVect&           a_refRatio)
: m_isDefined(false),
  m_levGeoPtr(NULL),
  m_crseLevGeoPtr(NULL),
  m_deleteLevGeoPtrs(false)
{
    define(a_grids, a_crseGridsPtr, a_dx, a_refRatio);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
MappedCFInterp::~MappedCFInterp ()
{
    undefine();
}


// -----------------------------------------------------------------------------
// Full define
// -----------------------------------------------------------------------------
void MappedCFInterp::define (const LevelGeometry& a_levGeo)
{
    if (m_isDefined) undefine();

    m_levGeoPtr = &a_levGeo;
    m_crseLevGeoPtr = a_levGeo.getCoarserPtr();
    CH_assert(m_crseLevGeoPtr != NULL);
    m_deleteLevGeoPtrs = false;

    m_domain = m_levGeoPtr->getDomain();
    m_grids = m_levGeoPtr->getBoxes();
    m_refRatio = m_levGeoPtr->getCrseRefRatio();

    m_cfregion.define(m_grids, m_domain);
    coarsen(m_crseGrids, m_grids, m_refRatio);

    m_crseGhosts = IntVect(D_DECL(2,2,2));
    m_crseCopier.define(m_crseLevGeoPtr->getBoxes(),
                        m_crseGrids,
                        m_crseGrids.physDomain(),
                        m_crseGhosts);

    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Full define
// -----------------------------------------------------------------------------
void MappedCFInterp::define (const DisjointBoxLayout& a_grids,
                             const DisjointBoxLayout* a_crseGridsPtr,
                             const RealVect&          a_dx,
                             const IntVect&           a_refRatio)
{
    if (m_isDefined) undefine();

    m_levGeoPtr = new LevelGeometry(a_dx);
    CH_assert(m_levGeoPtr != NULL);
    m_crseLevGeoPtr = new LevelGeometry(a_dx*RealVect(a_refRatio));
    CH_assert(m_crseLevGeoPtr != NULL);
    m_deleteLevGeoPtrs = true;

    m_domain = a_grids.physDomain();
    m_grids = a_grids;
    m_refRatio = a_refRatio;

    if (a_crseGridsPtr == NULL) MayDay::Error("Just checking");
    m_cfregion.define(m_grids, m_domain);
    coarsen(m_crseGrids, m_grids, m_refRatio);

    m_crseGhosts = IntVect(D_DECL(2,2,2));
    m_crseCopier.define(*a_crseGridsPtr,
                        m_crseGrids,
                        m_crseGrids.physDomain(),
                        m_crseGhosts);

    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Frees memory, makes object unusable
// -----------------------------------------------------------------------------
void MappedCFInterp::undefine ()
{
    if (!m_isDefined) return;

    if (m_deleteLevGeoPtrs) {
        delete m_levGeoPtr;
        delete m_crseLevGeoPtr;
    }
    m_levGeoPtr = NULL;
    m_crseLevGeoPtr = NULL;
    m_deleteLevGeoPtrs = false;

    m_cfregion = CFRegion();
    m_crseGrids = DisjointBoxLayout();

    m_domain = ProblemDomain();
    m_grids = DisjointBoxLayout();
    m_refRatio = IntVect::Zero;

    m_crseGhosts = IntVect(D_DECL(-1,-1,-1));
    m_crseCopier.clear();

    m_isDefined = false;
}


// -----------------------------------------------------------------------------
// Coarse / fine interpolation in time and space.
// -----------------------------------------------------------------------------
void MappedCFInterp::coarseFineInterp (LevelData<FArrayBox>&       a_phif,
                                       const LevelData<FArrayBox>& a_oldPhic,
                                       const LevelData<FArrayBox>& a_newPhic,
                                       const Real                  a_timeInterpCoeff,
                                       const int                   a_order,
                                       const bool                  a_useLinearLimiting) const
{
    CH_assert(m_isDefined);
    CH_assert(a_phif.nComp() == a_oldPhic.nComp());
    CH_assert(a_phif.nComp() == a_newPhic.nComp());
    CH_assert(0.0 <= a_timeInterpCoeff && a_timeInterpCoeff <= 1.0);

    if (a_timeInterpCoeff == 0.0) {
        // Just use old coarse data
        coarseFineInterp(a_phif, a_oldPhic, a_order, a_useLinearLimiting);

    } else if (a_timeInterpCoeff == 1.0) {
        // Just use new coarse data
        coarseFineInterp(a_phif, a_newPhic, a_order, a_useLinearLimiting);

    } else {
        // We need to interpolate
        LevelData<FArrayBox> phic(a_oldPhic.getBoxes(), a_oldPhic.nComp(), IntVect::Unit);
        timeInterp(phic, a_timeInterpCoeff, a_oldPhic, 0.0, a_newPhic, 1.0, phic.interval());
        coarseFineInterp(a_phif, phic, a_order, a_useLinearLimiting);
    }
}


// -----------------------------------------------------------------------------
// Coarse / fine interpolation operator.
// -----------------------------------------------------------------------------
void MappedCFInterp::coarseFineInterp (LevelData<FArrayBox>&       a_phif,
                                       const LevelData<FArrayBox>& a_phic,
                                       const int                   a_order,
                                       const bool                  a_useLinearLimiting) const
{
    CH_assert(m_isDefined);
    CH_assert(a_phif.nComp() == a_phic.nComp());
    CH_assert(a_phif.ghostVect() >= IntVect::Unit); // Changed from a_phic
    CH_assert(0 <= a_order && a_order <= 2);

    // Send coarse data to grids that are compatible with the finer level.
    const int ncomp = a_phif.nComp();
    LevelData<FArrayBox> crseData(m_crseGrids, ncomp, m_crseGhosts);
    a_phic.copyTo(crseData, m_crseCopier);

    MappedCoarseAverage avgObj;
    avgObj.define(m_grids, m_crseGrids, ncomp, m_refRatio, IntVect::Zero);
    avgObj.averageToCoarse(crseData, a_phif, m_levGeoPtr, true);

    // TODO: Is this needed?
    Copier exCopier(m_crseGrids, m_crseGrids, m_crseGrids.physDomain(), m_crseGhosts, true);
    crseData.exchange(exCopier);

    // Loop over grids, directions, and sides.
    DataIterator dit = m_grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& fineFAB = a_phif[dit];
        const FArrayBox& crseFAB = crseData[dit];

        for (int dir = 0; dir < SpaceDim; ++dir) {
            for (int iside = 0; iside < 2; ++iside) {
                // We are now at a specific grid boundary.
                // Interpolate!

                // Find the ghost region
                CFRegion& cfRef = (CFRegion&)m_cfregion;
                Box ghostBox = (iside == 0)?
                               cfRef.loCFIVS(dit(), dir).minBox():
                               cfRef.hiCFIVS(dit(), dir).minBox();
                if (ghostBox.isEmpty()) continue;

                // Include corner cells
                IntVect growVect = a_phif.ghostVect();
                growVect[dir] = 0;
                ghostBox.grow(growVect);

                // This is where the underlying coarse points live.
                Box crseGhostBox = coarsen(ghostBox, m_refRatio);
                CH_assert(!crseGhostBox.isEmpty());

                // Gather physical cell locations
                FArrayBox fineX(ghostBox, SpaceDim);
                m_levGeoPtr->fill_physCoor(fineX);

                FArrayBox crseX(grow(crseGhostBox,1), SpaceDim);
                m_crseLevGeoPtr->fill_physCoor(crseX);

                // Constant interpolation
                for (int n = 0; n < ncomp; ++n) {
                    BoxIterator bit(ghostBox);
                    for (bit.reset(); bit.ok(); ++bit) {
                        const IntVect& fineIV = bit();
                        const IntVect  crseIV = coarsen(fineIV, m_refRatio);
                        fineFAB(fineIV,n) = crseFAB(crseIV,n);
                    } // end loop over box (bit)
                } // end loop over comps

                // Is this all we needed?
                if (a_order < 1) continue;

                // Compute first derivatives
                FArrayBox slopes[CH_SPACEDIM];
                for (int slopeDir = 0; slopeDir < SpaceDim; ++slopeDir) {
                    slopes[slopeDir].define(crseGhostBox, ncomp);
                    const IntVect e = BASISV(slopeDir);

                    for (int n = 0; n < ncomp; ++n) {
                        BoxIterator bit(crseGhostBox);
                        for (bit.reset(); bit.ok(); ++bit) {
                            const IntVect& cc = bit();

                            Real dphi = crseFAB(cc+e,n) - crseFAB(cc-e,n);
                            Real dx = crseX(cc+e,slopeDir) - crseX(cc-e,slopeDir);
                            Real dphidx = dphi/dx;

                            // Linear limiting
                            if (a_useLinearLimiting) {
                                Real dphil = crseFAB(cc,n) - crseFAB(cc-e,n);
                                Real dxl = crseX(cc,slopeDir) - crseX(cc-e,slopeDir);
                                Real dphidxl = dphil/dxl;

                                Real dphir = crseFAB(cc+e,n) - crseFAB(cc,n);
                                Real dxr = crseX(cc+e,slopeDir) - crseX(cc,slopeDir);
                                Real dphidxr = dphir/dxr;

                                Real theta = 1.0; // (1 = most dissipative, 2 = least dissipative)
                                dphidx = minmod(theta*dphidxl, dphidx, theta*dphidxr);
                            }

                            slopes[slopeDir](cc,n) = dphidx;

                        } // end loop over box (bit)
                    } // end loop over comps (n)
                } // end loop over slope dirs (slopeDir)

                // Upgrade to linear interpolation
                for (int slopeDir = 0; slopeDir < SpaceDim; ++slopeDir) {
                    for (int n = 0; n < ncomp; ++n) {
                        BoxIterator bit(ghostBox);
                        for (bit.reset(); bit.ok(); ++bit) {
                            const IntVect& fineIV = bit();
                            const IntVect  crseIV = coarsen(fineIV, m_refRatio);

                            Real dphidx = slopes[slopeDir](crseIV,n);
                            Real dx = fineX(fineIV,slopeDir) - crseX(crseIV,slopeDir);
                            fineFAB(fineIV,n) += dphidx*dx;
                        } // end loop over box (bit)
                    } // end loop over comps (n)
                } // end loop over slope dirs (slopeDir)

                // Is this all we needed?
                if (a_order < 2) continue;

                // Compute second derivatives
                FArrayBox curves[CH_SPACEDIM][CH_SPACEDIM];
                for (int adir = 0; adir < SpaceDim; ++adir) {
                    for (int bdir = 0; bdir < SpaceDim; ++bdir) {
                        if (adir == bdir) {
                            // Reuse the slope holder
                            curves[adir][adir].define(Interval(0,ncomp-1), slopes[adir]);
                        } else {
                            curves[adir][bdir].define(crseGhostBox, ncomp);
                        }
                        const IntVect ea = BASISV(adir);
                        const IntVect eb = BASISV(bdir);

                        for (int n = 0; n < ncomp; ++n) {
                            BoxIterator bit(crseGhostBox);
                            for (bit.reset(); bit.ok(); ++bit) {
                                const IntVect& cc = bit();

                                Real dSlopeAdxB;
                                if (adir == bdir) {
                                    Real dxA = 0.5 * (crseX(cc+ea,adir) - crseX(cc-ea,adir));
                                    dSlopeAdxB = crseFAB(cc+ea,n) - 2.0*crseFAB(cc,n) + crseFAB(cc-ea,n);
                                    dSlopeAdxB /= (dxA*dxA);
                                } else {
                                    // Technically, dxA should be evaluated at cc+/-eb.
                                    Real dxA = crseX(cc+ea,adir) - crseX(cc-ea,adir);
                                    Real dxB = crseX(cc+eb,bdir) - crseX(cc-eb,bdir);
                                    dSlopeAdxB = (crseFAB(cc+ea+eb,n) - crseFAB(cc-ea+eb,n))
                                               - (crseFAB(cc+ea-eb,n) - crseFAB(cc-ea-eb,n));
                                    dSlopeAdxB /= (dxA*dxB);
                                }

                                curves[adir][bdir](cc,n) = dSlopeAdxB;
                            } // end loop over box (bit)
                        } // end loop over comps (n)
                    } // (bdir)
                } // (adir)

                // Upgrade to quadratic interpolation
                for (int adir = 0; adir < SpaceDim; ++adir) {
                    for (int bdir = 0; bdir < SpaceDim; ++bdir) {
                        for (int n = 0; n < ncomp; ++n) {
                            BoxIterator bit(ghostBox);
                            for (bit.reset(); bit.ok(); ++bit) {
                                const IntVect& fineIV = bit();
                                const IntVect  crseIV = coarsen(fineIV, m_refRatio);

                                Real curveAB = curves[adir][bdir](crseIV,n);
                                Real dxA = fineX(fineIV,adir) - crseX(crseIV,adir);
                                Real dxB = fineX(fineIV,bdir) - crseX(crseIV,bdir);
                                fineFAB(fineIV,n) += curveAB*dxA*dxB / 2.0;
                            } // end loop over box (bit)
                        } // end loop over comps (n)
                    } // (bdir)
                } // (adir)

            } // end loop over sides (iside)
        } // end loop over dirs (dir)
    } // end loop over fine grids (dit)

    Copier excp(m_grids, m_grids, m_domain, a_phif.ghostVect(), true);
    a_phif.exchange(excp);
}
