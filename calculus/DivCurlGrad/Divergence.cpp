/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Copyright (C) 2014 Edward Santilli & Alberto Scotti
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
#include "Divergence.H"
#include "DivCurlGradF_F.H"
#include "LevelGeometry.H"
#include "EllipticBCInterface.H"
#include "MappedQuadCFInterp.H"
#include "MappedLevelFluxRegister.H"
#include "CellToEdge.H"
#include "AnisotropicRefinementTools.H"
#include "Debug.H"


// With this defined, the CC divergence functions will use a smaller stencil
// that do not require corner cells.
// #define USE_SIMPLE_STENCIL


// Single-level divergence functions...

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Divergence::levelDivergenceMAC (LevelData<FArrayBox>&      a_div,
                                     const LevelData<FluxBox>&  a_uEdge,
                                     const LevelGeometry&       a_levGeo,
                                     const Real                 a_time,
                                     const BC_type*             a_fluxBC)
{
    CH_TIME("Divergence::levelDivergenceMAC");

    // Sanity check
    const RealVect& dx = a_levGeo.getDx();
    CH_assert(dx.product() > 0.0);

    const DisjointBoxLayout& grids = a_div.getBoxes();
    DataIterator dit = a_div.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        // Gather data holders
        FArrayBox& divFAB = a_div[dit];
        FluxBox& fluxFB = (FluxBox&)a_uEdge[dit];

        // Gather CC region to perform calculation
        // NOTE: Changed from grids[dit] & divFAB.box() on Nov 4, 2013 to fix bkgdSrc bug.
        const Box& region = divFAB.box();
        CH_assert(fluxFB.box().contains(region));

        // Gather reference to 1/J
        CH_assert(a_levGeo.getCCJinv().getBoxes().check(dit()));
        const FArrayBox& JinvFAB = a_levGeo.getCCJinv()[dit];
        CH_assert(JinvFAB.box().contains(region));

        // Set boundary fluxes if needed
        // WARNING: This only fills ghosts that abut grids[dit].
        if (a_fluxBC != NULL) {
            const ProblemDomain& domain = a_levGeo.getDomain();
            const FluxBox& JgupFB = a_levGeo.getFCJgup()[dit];

            for (int dir = 0; dir < SpaceDim; ++dir) {
                // NOTE: Changed from valid = region on Nov 4, 2013 to fix bkgdSrc bug.
                Box valid = grids[dit] & divFAB.box();
                valid.surroundingNodes(dir);
                Interval interv = Interval(0, fluxFB[dir].nComp() - 1);

                (*a_fluxBC)[dir].setGhosts(fluxFB[dir],    // state
                                           NULL,           // extrap
                                           valid,          // valid
                                           domain,         // domain
                                           dx,             // dx
                                           dit(),          // dataIndex
                                           &JgupFB,        // &JgupFB
                                           false,          // isHomogeneous
                                           a_time,         // time
                                           interv);        // interval
            }
        }

        // Compute divergence
#if CH_SPACEDIM == 2
        FORT_MAPPEDFLUXDIVERGENCE2D(
            CHF_FRA(divFAB),
            CHF_CONST_FRA(fluxFB[0]),
            CHF_CONST_FRA(fluxFB[1]),
            CHF_CONST_FRA1(JinvFAB,0),
            CHF_BOX(region),
            CHF_CONST_REALVECT(dx));

#elif CH_SPACEDIM == 3
        FORT_MAPPEDFLUXDIVERGENCE3D(
            CHF_FRA(divFAB),
            CHF_CONST_FRA(fluxFB[0]),
            CHF_CONST_FRA(fluxFB[1]),
            CHF_CONST_FRA(fluxFB[2]),
            CHF_CONST_FRA1(JinvFAB,0),
            CHF_BOX(region),
            CHF_CONST_REALVECT(dx));

#else
#   error Bad CH_SPACEDIM
#endif
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Divergence::simpleDivergenceMAC(FArrayBox&            a_div,
                                     const FluxBox&        a_uEdge,
                                     const DataIndex&      a_di,
                                     const LevelGeometry&  a_levGeo,
                                     const Real            a_time,
                                     const BC_type*        a_fluxBC)
{
    CH_TIME("Divergence::simpleDivergenceMAC");

    // Sanity checks
    const RealVect& dx = a_levGeo.getDx();
    CH_assert(dx.product() > 0.0);

    const DisjointBoxLayout& grids = a_levGeo.getCCJinv().getBoxes();
    CH_assert(grids.check(a_di));

    // Gather CC region to perform calculation
    const Box& region = grids[a_di] & a_div.box();
    CH_assert(a_uEdge.box().contains(region));

    // Gather reference to 1/J
    const FArrayBox& JinvFAB = a_levGeo.getCCJinv()[a_di];

    // Set boundary fluxes if needed
    if (a_fluxBC != NULL) {
        FluxBox& fluxFB = (FluxBox&)a_uEdge;
        const ProblemDomain& domain = a_levGeo.getDomain();
        const FluxBox& JgupFB = a_levGeo.getFCJgup()[a_di];

        for (int dir = 0; dir < SpaceDim; ++dir) {
            Box valid = surroundingNodes(region, dir);
            Interval interv = Interval(0, a_uEdge[dir].nComp() - 1);

            (*a_fluxBC)[dir].setGhosts(fluxFB[dir],    // state
                                       NULL,           // extrap
                                       valid,          // valid
                                       domain,         // domain
                                       dx,             // dx
                                       a_di,           // dataIndex
                                       &JgupFB,        // &JgupFB
                                       false,          // isHomogeneous
                                       a_time,         // time
                                       interv);        // interval
        }
    }

    // Compute divergence
#if CH_SPACEDIM == 2
    FORT_MAPPEDFLUXDIVERGENCE2D(
        CHF_FRA(a_div),
        CHF_CONST_FRA(a_uEdge[0]),
        CHF_CONST_FRA(a_uEdge[1]),
        CHF_CONST_FRA1(JinvFAB,0),
        CHF_BOX(region),
        CHF_CONST_REALVECT(dx));

#elif CH_SPACEDIM == 3
    FORT_MAPPEDFLUXDIVERGENCE3D(
        CHF_FRA(a_div),
        CHF_CONST_FRA(a_uEdge[0]),
        CHF_CONST_FRA(a_uEdge[1]),
        CHF_CONST_FRA(a_uEdge[2]),
        CHF_CONST_FRA1(JinvFAB,0),
        CHF_BOX(region),
        CHF_CONST_REALVECT(dx));

#else
#   error Bad CH_SPACEDIM
#endif
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Divergence::levelDivergenceCC (LevelData<FArrayBox>&       a_div,
                                    LevelData<FArrayBox>&       a_u,
                                    const LevelData<FArrayBox>* a_uCrsePtr,
                                    const bool                  a_quadInterp,
                                    const LevelGeometry&        a_levGeo,
                                    const Real                  a_time,
                                    const BC_type*              a_fluxBC)
{
    CH_TIME("Divergence::levelDivergenceCC_1");

    const RealVect& dx = a_levGeo.getDx();
    CH_assert(dx.product() > 0.0);

    const IntVect& nRefCrse = a_levGeo.getCrseRefRatio();
    const ProblemDomain& domain = a_levGeo.getDomain();

    if (a_uCrsePtr != NULL) {
        if (a_quadInterp) {

            // Need to compute coarse-fine BCs
            const DisjointBoxLayout& boxes = a_div.getBoxes();
            const DisjointBoxLayout& coarseBoxes = a_uCrsePtr->getBoxes();
            int nComp = a_u.nComp();

            MappedQuadCFInterp cfInterp(boxes,
                                        &coarseBoxes,
                                        dx,
                                        nRefCrse,
                                        nComp,
                                        domain);

            levelDivergenceCC(a_div,
                              a_u,
                              a_uCrsePtr,
                              a_quadInterp,
                              cfInterp,
                              a_levGeo,
                              a_time,
                              a_fluxBC);

        } else {
            // non-quadratic CF interp bc's not implemented yet
            CH_assert(a_quadInterp);
        }

    } else {

        // For single-level, won't need C/F interpolation anyway...
        MappedQuadCFInterp cfInterpBogus;
        levelDivergenceCC(a_div,
                          a_u,
                          a_uCrsePtr,
                          a_quadInterp,
                          cfInterpBogus,
                          a_levGeo,
                          a_time,
                          a_fluxBC);
    }
}


#ifdef USE_SIMPLE_STENCIL
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Divergence::levelDivergenceCC (LevelData<FArrayBox>&       a_div,
                                    LevelData<FArrayBox>&       a_u,
                                    const LevelData<FArrayBox>* a_uCrsePtr,
                                    const bool                  a_quadInterp,
                                    MappedQuadCFInterp&         a_cfInterp,
                                    const LevelGeometry&        a_levGeo,
                                    const Real                  a_time,
                                    const BC_type*              a_fluxBC)
{
    CH_TIME("Divergence::levelDivergenceCC_2 (simple stencil");

    // for now, hardwire to simple single-component case
    CH_assert(a_div.nComp() == 1);
    CH_assert(a_u.nComp() == SpaceDim);

    // Compute coarse-fine BC data
    if (a_uCrsePtr != NULL) {
        if (a_quadInterp) {
            CH_assert(a_cfInterp.isDefined());
            a_cfInterp.coarseFineInterp(a_u, *a_uCrsePtr);
        } else {
            // non-quadInterp not implemented yet
            CH_assert(a_quadInterp);
        }
    } // end if coarser level exists

    // Gather geometric info
    const RealVect& dx = a_levGeo.getDx();
    CH_assert(dx.product() > 0.0);
    const DisjointBoxLayout& grids = a_div.getBoxes();
    DataIterator dit = a_div.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        // Gather data holders
        FArrayBox& divFAB = a_div[dit];
        const FArrayBox& fluxFAB = a_u[dit];

        // Gather CC region to perform calculation
        // NOTE: Changed from grids[dit] & divFAB.box() on Nov 4, 2013 to fix bkgdSrc bug.
        const Box& region = divFAB.box();
        CH_assert(fluxFAB.box().contains(region));

        // Gather reference to 1/J
        CH_assert(a_levGeo.getCCJinv().getBoxes().check(dit()));
        const FArrayBox& JinvFAB = a_levGeo.getCCJinv()[dit];
        CH_assert(JinvFAB.box().contains(region));


        // Set boundary fluxes if needed
        if (a_fluxBC != NULL) {
            FArrayBox& castFluxFAB = (FArrayBox&)fluxFAB;
            const ProblemDomain& domain = a_levGeo.getDomain();
            const FluxBox& JgupFB = a_levGeo.getFCJgup()[dit];

            for (int dir = 0; dir < SpaceDim; ++dir) {
                Box valid = surroundingNodes(region, dir);
                // Interval interv = Interval(0, castFluxFAB.nComp() - 1);
                Interval interv = Interval(dir, dir);

                (*a_fluxBC)[dir].setGhosts(castFluxFAB,    // state
                                           NULL,           // extrap
                                           valid,          // valid
                                           domain,         // domain
                                           dx,             // dx
                                           dit(),          // dataIndex
                                           &JgupFB,        // &JgupFB
                                           false,          // isHomogeneous
                                           a_time,         // time
                                           interv);        // interval
            }
        }


        // Compute divergence
        const RealVect pdScale = RealVect::Unit;

        FORT_MAPPEDCCDIVSCALE(
            CHF_FRA1(divFAB,0),
            CHF_CONST_FRA(fluxFAB),
            CHF_CONST_FRA1(JinvFAB,0),
            CHF_BOX(region),
            CHF_CONST_REALVECT(dx),
            CHF_CONST_REALVECT(pdScale));
    }
}

#else
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Divergence::levelDivergenceCC(LevelData<FArrayBox>&       a_div,
                                   LevelData<FArrayBox>&       a_u,
                                   const LevelData<FArrayBox>* a_uCrsePtr,
                                   const bool                  a_quadInterp,
                                   MappedQuadCFInterp&         a_cfInterp,
                                   const LevelGeometry&        a_levGeo,
                                   const Real                  a_time,
                                   const BC_type*              a_fluxBC)
{
    CH_TIME("Divergence::levelDivergenceCC_2");

    // for now, hardwire to simple single-component case
    CH_assert (a_div.nComp() == 1);
    CH_assert (a_u.nComp() == SpaceDim);

    // Compute coarse-fine BC data
    if (a_uCrsePtr != NULL) {
        if (a_quadInterp) {
            CH_assert(a_cfInterp.isDefined());
            a_cfInterp.coarseFineInterp(a_u, *a_uCrsePtr);
        } else {
            // non-quadInterp not implemented yet
            CH_assert(a_quadInterp);
        }
    } // end if coarser level exists

    // Average CC->FC
    const DisjointBoxLayout& boxes = a_u.getBoxes();
    LevelData<FluxBox> uEdge(boxes, 1);
    CellToEdge(a_u, uEdge);

    // Take FC divergence
    levelDivergenceMAC(a_div,
                       uEdge,
                       a_levGeo,
                       a_time,
                       a_fluxBC);
}
#endif //USE_SIMPLE_STENCIL


// Composite divergence functions...

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Divergence::compDivergenceMAC(LevelData<FArrayBox>&     a_div,
                                   LevelData<FluxBox>&       a_uEdge,
                                   const LevelData<FluxBox>* a_uEdgeFinePtr,
                                   const LevelGeometry&      a_levGeo,
                                   const Real                a_time,
                                   const BC_type*            a_fluxBC)
{
    CH_TIME("Divergence::compDivergenceMAC_1");

    if (a_uEdgeFinePtr != NULL) {
        // define a MappedLevelFluxRegister to do coarse-fine mismatch accounting
        const DisjointBoxLayout& dblCrse = a_div.getBoxes();
        const DisjointBoxLayout& dblFine = a_uEdgeFinePtr->getBoxes();
        int ncomp = 1;

        ProblemDomain fineDomain = a_levGeo.getDomain();
        const IntVect& nRefine = a_levGeo.getFineRefRatio();
        fineDomain.refine(nRefine);

        MappedLevelFluxRegister FR(dblFine, dblCrse, fineDomain, nRefine, ncomp);

        // Then compute the composite divergence
        compDivergenceMAC(a_div,
                          a_uEdge,
                          a_uEdgeFinePtr,
                          &FR,
                          a_levGeo,
                          a_time,
                          a_fluxBC);
    } else {
        // Do single-level version
        levelDivergenceMAC(a_div,
                           a_uEdge,
                           a_levGeo,
                           a_time,
                           a_fluxBC);
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Divergence::compDivergenceMAC(LevelData<FArrayBox>&     a_div,
                                   LevelData<FluxBox>&       a_uEdge,
                                   const LevelData<FluxBox>* a_uEdgeFinePtr,
                                   MappedLevelFluxRegister*  a_fluxRegPtr,
                                   const LevelGeometry&      a_levGeo,
                                   const Real                a_time,
                                   const BC_type*            a_fluxBC)
{
    CH_TIME("Divergence::compDivergenceMAC_2");

    // For now, hardwire to simple single-component case
    CH_assert(a_div.nComp() == 1);
    int comp = 0;

    // First do simple single-level divergence, then fix up if necessary
    levelDivergenceMAC(a_div, a_uEdge, a_levGeo, a_time, a_fluxBC);

    // Now adjust for effect of finer level (if applicable)
    if (a_uEdgeFinePtr != NULL) {
        const RealVect& dx = a_levGeo.getDx();

        // Initialize flux register
        CH_assert(a_fluxRegPtr != NULL);
        MappedLevelFluxRegister& FR = *a_fluxRegPtr;
        FR.setToZero();

        // Do coarse side of FR
        DataIterator dit = a_div.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            FluxBox& thisFlux = a_uEdge[dit()];
            const Interval compInterval(comp, comp);

            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real scale = 1.0 / dx[dir];
                FR.incrementCoarse(thisFlux[dir], scale, dit(),
                                   compInterval, compInterval, dir);
            }
        }

        // Do fine side of FR
        DataIterator ditFine = a_uEdgeFinePtr->dataIterator();
        for (ditFine.reset(); ditFine.ok(); ++ditFine) {
            const FluxBox& thisFineFlux = (*a_uEdgeFinePtr)[ditFine()];
            Real scale = 1.0;
            Interval srcComps(comp,comp);

            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real scale = 1.0 / dx[dir];
                FR.incrementFine(thisFineFlux[dir], scale,
                                 ditFine(), srcComps,
                                 srcComps, dir, Side::Lo);

                FR.incrementFine(thisFineFlux[dir], scale,
                                 ditFine(), srcComps,
                                 srcComps, dir, Side::Hi);
            }
        }

        // Reflux
        FR.reflux(a_div, a_levGeo);

    } // end correction for finer level
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Divergence::compDivergenceCC(LevelData<FArrayBox>&       a_div,
                                  LevelData<FArrayBox>&       a_u,
                                  const LevelData<FArrayBox>* a_uCrsePtr,
                                  LevelData<FArrayBox>*       a_uFinePtr,
                                  const bool                  a_quadInterp,
                                  const LevelGeometry&        a_levGeo,
                                  const Real                  a_time,
                                  const BC_type*              a_fluxBC)
{
    CH_TIME("Divergence::compDivergenceCC_1");

    // for now, hardwire to single-component
    CH_assert(a_div.nComp() == 1);
    int nComp = 1;

    const RealVect& dx = a_levGeo.getDx();
    const IntVect& nRefCrse = a_levGeo.getCrseRefRatio();
    const IntVect& nRefFine = a_levGeo.getFineRefRatio();
    const ProblemDomain& domain = a_levGeo.getDomain();
    const DisjointBoxLayout& thisLevelBoxes = a_u.getBoxes();
    MappedQuadCFInterp cfInterpCrse;
    MappedQuadCFInterp cfInterpFine;
    MappedLevelFluxRegister FR;

    // Define coarse-level CF-BC object
    if (a_uCrsePtr != NULL) {
        D_TERM(CH_assert(nRefCrse[0] >= 1);,
               CH_assert(nRefCrse[1] >= 1);,
               CH_assert(nRefCrse[2] >= 1);)
        const DisjointBoxLayout& crseLevelBoxes = a_uCrsePtr->getBoxes();

        cfInterpCrse.define(thisLevelBoxes,
                            &crseLevelBoxes,
                            dx,
                            nRefCrse,
                            SpaceDim * nComp,
                            domain);
    }

    // Define fine-level CF-BC object and flux register
    if (a_uFinePtr != NULL) {
        D_TERM(CH_assert(nRefFine[0] >= 1);,
               CH_assert(nRefFine[1] >= 1);,
               CH_assert(nRefFine[2] >= 1);)
        const RealVect dxFine = dx / RealVect(nRefFine);

        const DisjointBoxLayout& fineLevelBoxes = a_uFinePtr->getBoxes();
        ProblemDomain fineDomain(domain);
        refine(fineDomain, domain, nRefFine);

        cfInterpFine.define(fineLevelBoxes,
                            &thisLevelBoxes,
                            dxFine,
                            nRefFine,
                            SpaceDim * nComp,
                            fineDomain);

        FR.define(fineLevelBoxes, thisLevelBoxes, fineDomain, nRefFine, nComp);
    }

    // Call the composite divergence operator
    compDivergenceCC(a_div,
                     a_u,
                     a_uCrsePtr,
                     a_uFinePtr,
                     a_quadInterp,
                     &FR,
                     cfInterpCrse,
                     cfInterpFine,
                     a_levGeo,
                     a_time,
                     a_fluxBC);
}


#ifdef USE_SIMPLE_STENCIL
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Divergence::compDivergenceCC(LevelData<FArrayBox>&       a_div,
                                  LevelData<FArrayBox>&       a_u,
                                  const LevelData<FArrayBox>* a_uCrsePtr,
                                  LevelData<FArrayBox>*       a_uFinePtr,
                                  const bool                  a_quadInterp,
                                  MappedLevelFluxRegister*    a_fluxRegFinePtr,
                                  MappedQuadCFInterp&         a_cfInterpCrse,
                                  MappedQuadCFInterp&         a_cfInterpFine,
                                  const LevelGeometry&        a_levGeo,
                                  const Real                  a_time,
                                  const BC_type*              a_fluxBC)
{
    CH_TIME("Divergence::compDivergenceCC_2");

    // for now, hardwire to simplest single-component case
    CH_assert(a_div.nComp() == 1);
    CH_assert(a_u.nComp() == SpaceDim);
    int comp = 0;

    // First, compute the no-fine-level CC divergence.
    levelDivergenceCC(a_div,
                      a_u,
                      a_uCrsePtr,
                      a_quadInterp,
                      a_cfInterpCrse,
                      a_levGeo,
                      a_time,
                      a_fluxBC);

    // Then, if a finer level exists, fix up at C/F interface.
    if (a_uFinePtr != NULL) {
        const RealVect& dx = a_levGeo.getDx();

        // Initialize flux register
        MappedLevelFluxRegister& FR = *a_fluxRegFinePtr;
        FR.setToZero();

        // Do CF-BC's on finer level data
        if (a_quadInterp) {
            CH_assert(a_cfInterpFine.isDefined());
            a_cfInterpFine.coarseFineInterp(*a_uFinePtr, a_u);
        } else {
            // non-quadInterp BC's not implemented at this point
            CH_assert(a_quadInterp);
        }

        // Do coarse side of FR
        const DisjointBoxLayout& grids = a_div.getBoxes();
        DataIterator dit = a_div.dataIterator();

        for (dit.reset(); dit.ok(); ++dit) {
            const Interval comps(comp,comp);

            for (int dir = 0; dir < SpaceDim; ++dir) {
                FArrayBox fluxFAB(surroundingNodes(grids[dit], dir), 1);
                CellToEdge(a_u[dit], fluxFAB, dir);

                const Real scale = 1.0 / dx[dir];
                FR.incrementCoarse(fluxFAB, scale, dit(),
                                   comps, comps, dir);
            }
        }

        // Do fine side of FR
        // Only averaging cell-> edge for the appropriate regions.
        const DisjointBoxLayout& fineGrids = a_uFinePtr->getBoxes();
        DataIterator ditFine = a_uFinePtr->dataIterator();
        FArrayBox cellData;
        FArrayBox edgeData;
        const Interval comps(comp,comp);

        for (ditFine.reset(); ditFine.ok(); ++ditFine) {
            FArrayBox& fineCCFab = (*a_uFinePtr)[ditFine()];
            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real scale = 1.0 / dx[dir];
                SideIterator sit;
                for (sit.begin(); sit.ok(); sit.next()) {
                    Side::LoHiSide hiorlo = sit();
                    const Box& fineValid = fineGrids[ditFine];
                    Box ccEdgeBox;
                    Box edgeBox;

                    if (hiorlo == Side::Lo) {
                        ccEdgeBox = adjCellLo(fineValid, dir, 2);
                        ccEdgeBox.shift(dir, 1);
                        edgeBox = bdryLo(fineValid, dir, 1);
                    } else {
                        ccEdgeBox = adjCellHi(fineValid, dir, 2);
                        ccEdgeBox.shift(dir, -1);
                        edgeBox = bdryHi(fineValid, dir, 1);
                    }

                    CH_assert(!ccEdgeBox.isEmpty());
                    CH_assert(!edgeBox.isEmpty());
                    CH_assert(fineCCFab.contains(ccEdgeBox));

                    cellData.resize(ccEdgeBox, 1);
                    edgeData.resize(edgeBox, 1);

                    // now copy cell-centered fine-level data into cellData
                    cellData.copy(fineCCFab, ccEdgeBox, dir, ccEdgeBox, 0, 1);

                    // now need to average cellData->edgeData
                    CellToEdge(cellData, 0, edgeData, 0, dir);

                    // now increment flux register
                    FR.incrementFine(edgeData, scale, ditFine(),
                                     comps, comps, dir, hiorlo);

                } // end iteration over high-lo
            } // end iteration over directions
        } // end iteration over fine boxes

        // Reflux
        FR.reflux(a_div, a_levGeo);

    } // end correction for finer levels
}

#else
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Divergence::compDivergenceCC(LevelData<FArrayBox>&       a_div,
                                  LevelData<FArrayBox>&       a_u,
                                  const LevelData<FArrayBox>* a_uCrsePtr,
                                  LevelData<FArrayBox>*       a_uFinePtr,
                                  const bool                  a_quadInterp,
                                  MappedLevelFluxRegister*    a_fluxRegFinePtr,
                                  MappedQuadCFInterp&         a_cfInterpCrse,
                                  MappedQuadCFInterp&         a_cfInterpFine,
                                  const LevelGeometry&        a_levGeo,
                                  const Real                  a_time,
                                  const BC_type*              a_fluxBC)
{
    CH_TIME("Divergence::compDivergenceCC_2");

    // for now, hardwire to simplest single-component case
    CH_assert(a_div.nComp() == 1);
    CH_assert(a_u.nComp() == SpaceDim);
    int comp = 0;

    // First do coarse-level BC's
    if (a_uCrsePtr != NULL) {
        if (a_quadInterp) {
            CH_assert(a_cfInterpCrse.isDefined());
            a_cfInterpCrse.coarseFineInterp(a_u, *a_uCrsePtr);
        } else {
            // non-quadInterp BC's not yet implemented
            CH_assert(a_quadInterp);
        }
    }

    // now average cells->edges
    const DisjointBoxLayout& boxes = a_u.getBoxes();
    LevelData<FluxBox> uEdge(boxes, 1);
    CellToEdge(a_u, uEdge);

    // Compute no-fine-level divergence
    levelDivergenceMAC(a_div,
                       uEdge,
                       a_levGeo,
                       a_time,
                       a_fluxBC);

    // if a fine level exists, fix up at C/F interface:
    if (a_uFinePtr != NULL) {
        const RealVect& dx = a_levGeo.getDx();

        // Initialize flux register
        MappedLevelFluxRegister& FR = *a_fluxRegFinePtr;
        FR.setToZero();

        // Do CF-BC's on finer level data
        if (a_quadInterp) {
            CH_assert(a_cfInterpFine.isDefined());
            a_cfInterpFine.coarseFineInterp(*a_uFinePtr, a_u);
        } else {
            // non-quadInterp BC's not implemented at this point
            CH_assert(a_quadInterp);
        }

        // Do coarse side of FR
        DataIterator dit = a_div.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            FluxBox& thisFlux = uEdge[dit()];
            const Interval comps(comp,comp);

            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real scale = 1.0 / dx[dir];
                FR.incrementCoarse(thisFlux[dir], scale, dit(),
                                   comps, comps, dir);
            }
        }

        // Do fine side of FR
        // Only averaging cell-> edge for the appropriate regions.
        const DisjointBoxLayout& fineGrids = a_uFinePtr->getBoxes();
        DataIterator ditFine = a_uFinePtr->dataIterator();
        FArrayBox cellData;
        FArrayBox edgeData;
        const Interval comps(comp,comp);

        for (ditFine.reset(); ditFine.ok(); ++ditFine) {
            FArrayBox& fineCCFab = (*a_uFinePtr)[ditFine()];
            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real scale = 1.0 / dx[dir];
                SideIterator sit;
                for (sit.begin(); sit.ok(); sit.next()) {
                    Side::LoHiSide hiorlo = sit();
                    const Box& fineValid = fineGrids[ditFine];
                    Box ccEdgeBox;
                    Box edgeBox;

                    if (hiorlo == Side::Lo) {
                        ccEdgeBox = adjCellLo(fineValid, dir, 2);
                        ccEdgeBox.shift(dir, 1);
                        edgeBox = bdryLo(fineValid, dir, 1);
                    } else {
                        ccEdgeBox = adjCellHi(fineValid, dir, 2);
                        ccEdgeBox.shift(dir, -1);
                        edgeBox = bdryHi(fineValid, dir, 1);
                    }

                    CH_assert(!ccEdgeBox.isEmpty());
                    CH_assert(!edgeBox.isEmpty());
                    CH_assert(fineCCFab.contains(ccEdgeBox));

                    cellData.resize(ccEdgeBox, 1);
                    edgeData.resize(edgeBox, 1);

                    // now copy cell-centered fine-level data into cellData
                    cellData.copy(fineCCFab, ccEdgeBox, dir, ccEdgeBox, 0, 1);

                    // now need to average cellData->edgeData
                    CellToEdge(cellData, 0, edgeData, 0, dir);

                    // now increment flux register
                    FR.incrementFine(edgeData, scale, ditFine(),
                                     comps, comps, dir, hiorlo);

                } // end iteration over high-lo
            } // end iteration over directions
        } // end iteration over fine boxes

        // Reflux
        FR.reflux(a_div, a_levGeo);

    } // end correction for finer levels
}
#endif //USE_SIMPLE_STENCIL

