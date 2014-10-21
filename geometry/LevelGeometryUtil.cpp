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
#include "LevelGeometry.H"
#include "LevelGeometryF_F.H"


// -----------------------------------------------------------------------------
// Return a pointer to the coarser LevelGeometry
// -----------------------------------------------------------------------------
const LevelGeometry* LevelGeometry::getCoarserPtr () const
{
    CH_assert(this->isDefined());
    return m_coarserPtr;
}


// -----------------------------------------------------------------------------
// Return a pointer to the finer LevelGeometry
// -----------------------------------------------------------------------------
const LevelGeometry* LevelGeometry::getFinerPtr () const
{
    CH_assert(this->isDefined());
    return m_finerPtr;
}


// -----------------------------------------------------------------------------
// Returns all levgGeos over the AMR hierarchy
// -----------------------------------------------------------------------------
Vector<const LevelGeometry*> LevelGeometry::getAMRLevGeos () const
{
    Vector<const LevelGeometry*> vLevGeo(0);

    // Find the coarsest levGeo
    const LevelGeometry* levGeoPtr = this;
    while (levGeoPtr->getCoarserPtr() != NULL) {
        levGeoPtr = levGeoPtr->getCoarserPtr();
    }

    // Work up the levels, collecting pointers
    while (levGeoPtr != NULL) {
        vLevGeo.push_back(levGeoPtr);
        levGeoPtr = levGeoPtr->getFinerPtr();
    }

    return vLevGeo;
}


// -----------------------------------------------------------------------------
// Returns all levgGeos over the AMR hierarchy (private, non-const version)
// -----------------------------------------------------------------------------
Vector<LevelGeometry*> LevelGeometry::getAMRLevGeos ()
{
    Vector<LevelGeometry*> vLevGeo(0);

    // Find the coarsest levGeo
    LevelGeometry* levGeoPtr = this;
    while (levGeoPtr->getCoarserPtr() != NULL) {
        levGeoPtr = const_cast<LevelGeometry*>(levGeoPtr->getCoarserPtr());
    }

    // Work up the levels, collecting pointers
    while (levGeoPtr != NULL) {
        vLevGeo.push_back(levGeoPtr);
        levGeoPtr = const_cast<LevelGeometry*>(levGeoPtr->getFinerPtr());
    }

    return vLevGeo;
}


// -----------------------------------------------------------------------------
// Returns all fineRefRatios over the AMR hierarchy
// -----------------------------------------------------------------------------
Vector<IntVect> LevelGeometry::getAMRRefRatios () const
{
    Vector<IntVect> vRefRatios(0);

    // Find the coarsest levGeo
    const LevelGeometry* levGeoPtr = this;
    while (levGeoPtr->getCoarserPtr() != NULL) {
        levGeoPtr = levGeoPtr->getCoarserPtr();
    }

    // Work up the levels, collecting crseRefRatios
    while (levGeoPtr != NULL) {
        vRefRatios.push_back(levGeoPtr->getFineRefRatio());
        levGeoPtr = levGeoPtr->getFinerPtr();
    }

    return vRefRatios;
}


// -----------------------------------------------------------------------------
// Returns all DisjointBoxLayouts over the AMR hierarchy
// -----------------------------------------------------------------------------
Vector<DisjointBoxLayout> LevelGeometry::getAMRGrids () const
{
    Vector<DisjointBoxLayout> vGrids(0);

    // Find the coarsest levGeo
    const LevelGeometry* levGeoPtr = this;
    while (levGeoPtr->getCoarserPtr() != NULL) {
        levGeoPtr = levGeoPtr->getCoarserPtr();
    }

    // Work up the levels, collecting crseRefRatios
    while (levGeoPtr != NULL) {
        vGrids.push_back(levGeoPtr->getBoxes());
        levGeoPtr = levGeoPtr->getFinerPtr();
    }

    return vGrids;
}


// -----------------------------------------------------------------------------
// Sends a vector from the scaled, mapped basis to an unscaled, Cartesian basis.
// uCart^{a} = [dx^{a} / dXi^{b}] * uMapped^{b}
// -----------------------------------------------------------------------------
void LevelGeometry::sendToCartesianBasis (LevelData<FArrayBox>& a_CCvect,
                                          bool                  a_doGhosts) const
{
    const DisjointBoxLayout& grids = a_CCvect.getBoxes();

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        FArrayBox& vectFAB = a_CCvect[dit];
        const Box& region = (a_doGhosts ? vectFAB.box() : grids[dit]);

        FArrayBox dxdXiFAB(region, SpaceDim*SpaceDim);
        this->fill_dxdXi(dxdXiFAB);

        FORT_CONTRACTMATRIXVECTORCC(
            CHF_FRA(vectFAB),
            CHF_CONST_FRA(dxdXiFAB),
            CHF_BOX(region));
    }
}


// -----------------------------------------------------------------------------
// Sends a vector from the scaled, mapped basis to an unscaled, Cartesian basis.
// uCart^{a} = [dx^{a} / dXi^{b}] * uMapped^{b}
// -----------------------------------------------------------------------------
void LevelGeometry::sendToCartesianBasis (LevelData<FluxBox>& a_FCvect,
                                          bool                a_doGhosts) const
{
    const DisjointBoxLayout& grids = a_FCvect.getBoxes();
    const int ncomp = a_FCvect.nComp();

    CH_assert(ncomp == SpaceDim);

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        // Convert to mapped basis
        for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
            FArrayBox& vectFAB = a_FCvect[dit][dir];
            const Box& region = (a_doGhosts ? vectFAB.box() : surroundingNodes(grids[dit], dir));

            // Fill a holder with the inverse Jacobian matrix
            FArrayBox dxdXiFAB(region, SpaceDim*SpaceDim);
            this->fill_dxdXi(dxdXiFAB);

            FORT_CONTRACTMATRIXVECTORCC(
                CHF_FRA(vectFAB),
                CHF_CONST_FRA(dxdXiFAB),
                CHF_BOX(region));
        } // end loop over FC dirs (dir)
    } // end loop over grids (dit)
}


// -----------------------------------------------------------------------------
// Sends a vector from the unscaled, Cartesian basis to an scaled, mapped basis.
// uMapped^{a} = [dXi^{a} / dx^{b}] * uCart^{b}
// -----------------------------------------------------------------------------
void LevelGeometry::sendToMappedBasis (LevelData<FArrayBox>& a_CCvect,
                                       bool                  a_doGhosts) const
{
    CH_assert(a_CCvect.nComp() == SpaceDim);

    const DisjointBoxLayout& grids = a_CCvect.getBoxes();

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        FArrayBox& vectFAB = a_CCvect[dit];
        const Box& region = (a_doGhosts ? vectFAB.box() : grids[dit]);

        FArrayBox dXidxFAB(region, SpaceDim*SpaceDim);
        this->fill_dXidx(dXidxFAB);

        FORT_CONTRACTMATRIXVECTORCC(
            CHF_FRA(vectFAB),
            CHF_CONST_FRA(dXidxFAB),
            CHF_BOX(region));
    }
}


// -----------------------------------------------------------------------------
// Sends a vector from the unscaled, Cartesian basis to an scaled, mapped basis.
// uMapped^{a} = [dXi^{a} / dx^{b}] * uCart^{b}
// -----------------------------------------------------------------------------
void LevelGeometry::sendToMappedBasis (LevelData<FluxBox>& a_FCvect,
                                       bool                a_doGhosts) const
{
    const DisjointBoxLayout& grids = a_FCvect.getBoxes();
    const int ncomp = a_FCvect.nComp();

    CH_assert(ncomp == SpaceDim);

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        // Convert to mapped basis
        for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
            FArrayBox& vectFAB = a_FCvect[dit][dir];
            const Box& region = (a_doGhosts ? vectFAB.box() : surroundingNodes(grids[dit], dir));

            // Fill a holder with the inverse Jacobian matrix
            FArrayBox dXidxFAB(region, SpaceDim*SpaceDim);
            this->fill_dXidx(dXidxFAB);

            FORT_CONTRACTMATRIXVECTORCC(
                CHF_FRA(vectFAB),
                CHF_CONST_FRA(dXidxFAB),
                CHF_BOX(region));
        } // end loop over FC dirs (dir)
    } // end loop over grids (dit)
}


// -----------------------------------------------------------------------------
// Multiplies every element by J (Entire AMR hierarchy version)
// -----------------------------------------------------------------------------
void LevelGeometry::multByJ(Vector<LevelData<FArrayBox>*>& a_data) const
{

    // Begin at the bottom level
    const LevelGeometry* levGeoPtr = this;
    while (levGeoPtr->m_coarserPtr != NULL) {
        levGeoPtr = levGeoPtr->getCoarserPtr();
    }

    // Iterate over data vector elements
    const int numLevels = a_data.size();
    for (int ilev = 0; ilev < numLevels; ++ilev) {
        // Do the multiplication if this level has data
        if (a_data[ilev] != NULL) {
            CH_assert(a_data[ilev]->getBoxes().physDomain().domainBox()
                      == levGeoPtr->getDomain().domainBox());
            levGeoPtr->multByJ(*a_data[ilev]);
        }

        // Move on to the next finer level
        levGeoPtr = levGeoPtr->getFinerPtr();
    }
}


// -----------------------------------------------------------------------------
// Multiplies every element by J (Entire level version)
// -----------------------------------------------------------------------------
void LevelGeometry::multByJ(LevelData<FArrayBox>& a_data) const
{
    CH_assert(a_data.getBoxes().physDomain().domainBox() == this->getDomain().domainBox());

    // Just loop over this level's grids and call the single grid version.
    DataIterator dit = a_data.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        this->multByJ(a_data[dit], dit());
    }
}


// -----------------------------------------------------------------------------
// Multiplies every element by J (Single grid version)
// -----------------------------------------------------------------------------
void LevelGeometry::multByJ(FArrayBox& a_data, const DataIndex& a_di, int a_comp) const
{

    CH_assert(m_grids.check(a_di));
    CH_assert(-1 <= a_comp && a_comp < a_data.nComp());

    // If a_comp is -1, then do all comps.
    int startcomp, numcomp;
    if (a_comp == -1) {
        startcomp = 0;
        numcomp = a_data.nComp();
    } else {
        startcomp = a_comp;
        numcomp = 1;
    }

    const Box& region = a_data.box();
    if (region.type() == IntVect::Zero) {
        // We can use the Jinv stored in the cache.
        const FArrayBox& JFAB = this->getCCJ()[a_di];
        CH_assert(JFAB.box().contains(region));

        // Scale each component requested
        for (int n = startcomp; n < startcomp + numcomp; ++n) {
            a_data.mult(JFAB, 0, n, 1);
        }
    } else {
        // The cache does not support this centering, so we need
        // to fill a custom holder in this region with J.
        FArrayBox JFAB(region, 1);
        this->fill_J(JFAB);

        // Scale each component requested
        for (int n = startcomp; n < startcomp + numcomp; ++n) {
            a_data.mult(JFAB, 0, n, 1);
        }
    }
}


// -----------------------------------------------------------------------------
// Multiplies every element by J (Entire AMR hierarchy version)
// -----------------------------------------------------------------------------
void LevelGeometry::multByJ(Vector<LevelData<FluxBox>*>& a_data) const
{

    // Begin at the bottom level
    const LevelGeometry* levGeoPtr = this;
    while (levGeoPtr->m_coarserPtr != NULL) {
        levGeoPtr = levGeoPtr->getCoarserPtr();
    }

    // Iterate over data vector elements
    const int numLevels = a_data.size();
    for (int ilev = 0; ilev < numLevels; ++ilev) {
        // Do the multiplication if this level has data
        if (a_data[ilev] != NULL) {
            CH_assert(a_data[ilev]->getBoxes().physDomain().domainBox() == levGeoPtr->getDomain().domainBox());
            levGeoPtr->multByJ(*a_data[ilev]);
        }

        // Move on to the next finer level
        levGeoPtr = levGeoPtr->getFinerPtr();
    }
}


// -----------------------------------------------------------------------------
// Multiplies every element by J (Entire level version)
// -----------------------------------------------------------------------------
void LevelGeometry::multByJ(LevelData<FluxBox>& a_data) const
{
    CH_assert(a_data.getBoxes().physDomain().domainBox() == this->getDomain().domainBox());

    // Just loop over this level's grids and call the single grid version.
    DataIterator dit = a_data.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        this->multByJ(a_data[dit], dit());
    }
}


// -----------------------------------------------------------------------------
// Multiplies every element by J (Single grid version)
// -----------------------------------------------------------------------------
void LevelGeometry::multByJ(FluxBox& a_data, const DataIndex& a_di, int a_comp) const
{
    CH_TIME("LevelGeometry::multByJ");

    CH_assert(-1 <= a_comp && a_comp < a_data.nComp());

    // If a_comp is -1, then do all comps.
    int startcomp, numcomp;
    if (a_comp == -1) {
        startcomp = 0;
        numcomp = a_data.nComp();
    } else {
        startcomp = a_comp;
        numcomp = 1;
    }

    // Loop over FABs in the FluxBox
    for (int FCdir = 0; FCdir < SpaceDim; ++FCdir) {
        FArrayBox& dataFAB = a_data[FCdir];

        // Fill a holder in this region with J.
        FArrayBox JFAB(dataFAB.box(), 1);
        this->fill_J(JFAB);

        // Scale each component requested
        for (int n = startcomp; n < startcomp + numcomp; ++n) {
            dataFAB.mult(JFAB, 0, n, 1);
        }
    }
}


// -----------------------------------------------------------------------------
// Divides every element by J (Entire AMR hierarchy version)
// -----------------------------------------------------------------------------
void LevelGeometry::divByJ(Vector<LevelData<FArrayBox>*>& a_data) const
{

    // Begin at the bottom level
    const LevelGeometry* levGeoPtr = this;
    while (levGeoPtr->m_coarserPtr != NULL) {
        levGeoPtr = levGeoPtr->getCoarserPtr();
    }

    // Iterate over data vector elements
    const int numLevels = a_data.size();
    for (int ilev = 0; ilev < numLevels; ++ilev) {
        // Do the multiplication if this level has data
        if (a_data[ilev] != NULL) {
            CH_assert(a_data[ilev]->getBoxes().physDomain().domainBox() == levGeoPtr->getDomain().domainBox());
            levGeoPtr->divByJ(*a_data[ilev]);
        }

        // Move on to the next finer level
        levGeoPtr = levGeoPtr->getFinerPtr();
    }
}


// -----------------------------------------------------------------------------
// Divides every element by J (Entire level version)
// -----------------------------------------------------------------------------
void LevelGeometry::divByJ(LevelData<FArrayBox>& a_data) const
{
    CH_assert(a_data.getBoxes().physDomain().domainBox() == this->getDomain().domainBox());

    // Just loop over this level's grids and call the single grid version.
    DataIterator dit = a_data.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        this->divByJ(a_data[dit], dit());
    }
}


// -----------------------------------------------------------------------------
// Divides every element by J (Single grid version)
// -----------------------------------------------------------------------------
void LevelGeometry::divByJ(FArrayBox& a_data, const DataIndex& a_di, int a_comp) const
{

    CH_assert(m_grids.check(a_di));
    CH_assert(-1 <= a_comp && a_comp < a_data.nComp());

    // If a_comp is -1, then do all comps.
    int startcomp, numcomp;
    if (a_comp == -1) {
        startcomp = 0;
        numcomp = a_data.nComp();
    } else {
        startcomp = a_comp;
        numcomp = 1;
    }

    const Box& region = a_data.box();
    if (region.type() == IntVect::Zero) {
        // We can use the Jinv stored in the cache.
        const FArrayBox& JinvFAB = this->getCCJinv()[a_di];
        CH_assert(JinvFAB.box().contains(region));

        // Scale each component requested
        for (int n = startcomp; n < startcomp + numcomp; ++n) {
            a_data.mult(JinvFAB, 0, n, 1);
        }
    } else {
        // The cache does not support this centering, so we need
        // to fill a custom holder in this region with Jinv.
        FArrayBox JinvFAB(region, 1);
        this->fill_Jinv(JinvFAB);

        // Scale each component requested
        for (int n = startcomp; n < startcomp + numcomp; ++n) {
            a_data.mult(JinvFAB, 0, n, 1);
        }
    }
}


// -----------------------------------------------------------------------------
// Divides every element by J (Entire AMR hierarchy version)
// -----------------------------------------------------------------------------
void LevelGeometry::divByJ(Vector<LevelData<FluxBox>*>& a_data) const
{

    // Begin at the bottom level
    const LevelGeometry* levGeoPtr = this;
    while (levGeoPtr->m_coarserPtr != NULL) {
        levGeoPtr = levGeoPtr->getCoarserPtr();
    }

    // Iterate over data vector elements
    const int numLevels = a_data.size();
    for (int ilev = 0; ilev < numLevels; ++ilev) {
        // Do the multiplication if this level has data
        if (a_data[ilev] != NULL) {
            CH_assert(a_data[ilev]->getBoxes().physDomain().domainBox() == levGeoPtr->getDomain().domainBox());
            levGeoPtr->divByJ(*a_data[ilev]);
        }

        // Move on to the next finer level
        levGeoPtr = levGeoPtr->getFinerPtr();
    }
}


// -----------------------------------------------------------------------------
// Divides every element by J (Entire level version)
// -----------------------------------------------------------------------------
void LevelGeometry::divByJ(LevelData<FluxBox>& a_data) const
{
    CH_assert(a_data.getBoxes().physDomain().domainBox() == this->getDomain().domainBox());

    // Just loop over this level's grids and call the single grid version.
    DataIterator dit = a_data.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        this->divByJ(a_data[dit], dit());
    }
}


// -----------------------------------------------------------------------------
// Divides every element by J (Single grid version)
// -----------------------------------------------------------------------------
void LevelGeometry::divByJ(FluxBox& a_data, const DataIndex& a_di, int a_comp) const
{
    CH_TIME("LevelGeometry::divByJ");

    CH_assert(-1 <= a_comp && a_comp < a_data.nComp());

    // If a_comp is -1, then do all comps.
    int startcomp, numcomp;
    if (a_comp == -1) {
        startcomp = 0;
        numcomp = a_data.nComp();
    } else {
        startcomp = a_comp;
        numcomp = 1;
    }

    // Loop over FABs in the FluxBox
    for (int FCdir = 0; FCdir < SpaceDim; ++FCdir) {
        FArrayBox& dataFAB = a_data[FCdir];

        // Fill a holder in this region with Jinv.
        FArrayBox JinvFAB(dataFAB.box(), 1);
        this->fill_Jinv(JinvFAB);

        // Scale each component requested
        for (int n = startcomp; n < startcomp + numcomp; ++n) {
            dataFAB.mult(JinvFAB, 0, n, 1);
        }
    }
}


// -----------------------------------------------------------------------------
// Contracts a CC contra-vector with gdn, making it covariant.
// (Single-grid version)
// -----------------------------------------------------------------------------
void LevelGeometry::makeCovariant (FArrayBox&       a_coVect,
                                   const FArrayBox& a_contraVect,
                                   const Box&       a_region,
                                   const DataIndex& a_di) const
{
    CH_TIME("LevelGeometry::makeCovariant");

    CH_assert(a_coVect.nComp() == SpaceDim);
    CH_assert(a_contraVect.nComp() == SpaceDim);
    CH_assert(a_coVect.box().contains(a_region));
    CH_assert(a_contraVect.box().contains(a_region));

    const FArrayBox& gdnFAB = this->getCCgdn()[a_di];

    FORT_CONTRACTWITHMETRIC(
        CHF_FRA(a_coVect),
        CHF_CONST_FRA(a_contraVect),
        CHF_BOX(a_region),
        CHF_CONST_FRA(gdnFAB));
}


// -----------------------------------------------------------------------------
// Contracts a CC covariant one-form with gup, making it a contravariant vector.
// (static, single-grid version)
// -----------------------------------------------------------------------------
void LevelGeometry::makeContravariant (FArrayBox&       a_contraVect,
                                       const FArrayBox& a_coVect,
                                       const Box&       a_region,
                                       const DataIndex& a_di,
                                       const RealVect&  a_dXi)
{
    CH_TIME("LevelGeometry::makeContravariant (static version)");

    CH_assert(a_coVect.nComp() == SpaceDim);
    CH_assert(a_contraVect.nComp() == SpaceDim);
    CH_assert(a_coVect.box().contains(a_region));
    CH_assert(a_contraVect.box().contains(a_region));

    // const FArrayBox& gup = this->getCCgup()[a_di];
    FArrayBox gupFAB(a_region, SpaceDim*(SpaceDim+1)/2);
    fill_gup(gupFAB, a_dXi);

    FORT_CONTRACTWITHMETRIC(
        CHF_FRA(a_contraVect),
        CHF_CONST_FRA(a_coVect),
        CHF_BOX(a_region),
        CHF_CONST_FRA(gupFAB));
}


// -----------------------------------------------------------------------------
// Contracts a CC covariant one-form with gup, making it a contravariant vector.
// (Single-grid version)
// -----------------------------------------------------------------------------
void LevelGeometry::makeContravariant (FArrayBox&       a_contraVect,
                                       const FArrayBox& a_coVect,
                                       const Box&       a_region,
                                       const DataIndex& a_di) const
{
    CH_TIME("LevelGeometry::makeContravariant");

    CH_assert(a_coVect.nComp() == SpaceDim);
    CH_assert(a_contraVect.nComp() == SpaceDim);
    CH_assert(a_coVect.box().contains(a_region));
    CH_assert(a_contraVect.box().contains(a_region));

    // const FArrayBox& gup = this->getCCgup()[a_di];
    FArrayBox gupFAB(a_region, SpaceDim*(SpaceDim+1)/2);
    this->fill_gup(gupFAB);

    FORT_CONTRACTWITHMETRIC(
        CHF_FRA(a_contraVect),
        CHF_CONST_FRA(a_coVect),
        CHF_BOX(a_region),
        CHF_CONST_FRA(gupFAB));
}


// -----------------------------------------------------------------------------
// Compute the magnitude of a vector
// (Single-grid version)
// -----------------------------------------------------------------------------
void LevelGeometry::contractVectors (FArrayBox&       a_res,
                                     const FArrayBox& a_vec1,
                                     const FArrayBox& a_vec2,
                                     const DataIndex& a_di)
{
    CH_TIME("LevelGeometry::contractVectors");

    // Sanity checks
    CH_assert(a_res.nComp() == 1);
    CH_assert(a_vec1.nComp() == SpaceDim);
    CH_assert(a_vec2.nComp() == SpaceDim);

    Box region = a_res.box();
    CH_assert(a_vec1.box().contains(region));
    CH_assert(a_vec2.box().contains(region));

    const FArrayBox& gdn = this->getCCgdn()[a_di];
    CH_assert(gdn.box().contains(region));

    // Compute
    FORT_CONTRACTVECTORS(
        CHF_FRA1(a_res,0),
        CHF_CONST_FRA(a_vec1),
        CHF_CONST_FRA(a_vec2),
        CHF_BOX(region),
        CHF_CONST_FRA(gdn));
}


#ifndef NDEBUG
// -----------------------------------------------------------------------------
// Dumps debugging info to pout()
// -----------------------------------------------------------------------------
void LevelGeometry::dump () const
{
    pout() << "LevelGeometry dump:\n"
           << "\tm_finerPtr   = " << (void*)m_finerPtr << "\n"
           << "\tthis         = " << (void*)this << "\n"
           << "\tm_coarserPtr = " << (void*)m_coarserPtr << "\n"
           << "\tCoordinate map = " << getCoorMapName() << "\n"
           << "\tis uniform  = " << (isUniform()? "true": "false") << "\n"
           << "\tis diagonal = " << (isDiagonal()? "true": "false") << "\n"
           << "\tm_dXi            = " << m_dXi << "\n"
           << "\tm_fineRefRatio   = " << m_fineRefRatio << "\n"
           << "\tm_crseRefRatio   = " << m_crseRefRatio << "\n"
           << "\tgrids.physDomain() = " << m_grids.physDomain() << "\n"
           << "\tgrids = " << m_grids << std::flush;
}
#endif

