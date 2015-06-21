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
#include "Gradient.H"
#include "DivCurlGradF_F.H"
#include "Mask.H"
#include "LevelGeometry.H"
#include "ExtrapolationUtils.H"
#include "EllipticBCInterface.H"
#include "MappedQuadCFInterp.H"
#include "IntVectSet.H"
#include "EdgeToCell.H"
#include "AnisotropicRefinementTools.H"
#include "Debug.H"
#include "Constants.H"


// You should leave these undefined. They exist for debugging purposes.
//
// Old method: Handle all corner ghosts the same.
// New method: Handle domain corners specially and treat all exchange
//             ghosts the same.
// #define USE_OLD_EXTRAP
//
// With this defined, the CC gradient functions will use a smaller stencil
// that do not require corner cells.
// #define USE_SIMPLE_STENCIL


// Single-level gradient fuunctions ============================================
// -----------------------------------------------------------------------------
void Gradient::levelGradientMAC (LevelData<FluxBox>&         a_edgeGrad,
                                 LevelData<FArrayBox>&       a_phi,
                                 const LevelData<FArrayBox>* a_phiCrsePtr,
                                 const LevelGeometry&        a_levGeo,
                                 const Real                  a_time,
                                 const BC_type*              a_fluxBC)
{
    CH_TIME("Gradient::levelGradientMAC_1");

    // Create the QuadFCInterp object
    MappedQuadCFInterp cfInterpCrse;

    if (a_phiCrsePtr != NULL) {
        cfInterpCrse.define(a_phi.getBoxes(),
                            &(a_phiCrsePtr->getBoxes()),
                            a_levGeo.getDx(),
                            a_levGeo.getCrseRefRatio(),
                            a_phi.nComp(),
                            a_levGeo.getDomain());
    }

    // Call the full gradient function
    levelGradientMAC(a_edgeGrad,
                     a_phi,
                     a_phiCrsePtr,
                     cfInterpCrse,
                     a_levGeo,
                     a_time,
                     a_fluxBC);
}


// -----------------------------------------------------------------------------
void Gradient::levelGradientMAC (LevelData<FluxBox>&           a_edgeGrad,
                                 LevelData<FArrayBox>&         a_phi,
                                 const LevelData<FArrayBox>*   a_phiCrsePtr,
                                 const MappedQuadCFInterp&     a_cfInterpCrse,
                                 const LevelGeometry&          a_levGeo,
                                 const Real                    a_time,
                                 const BC_type*                a_fluxBC)
{
    CH_TIME("Gradient::levelGradientMAC_2");

    CH_assert(a_edgeGrad.nComp() >= a_phi.nComp());
    CH_assert(a_edgeGrad.getBoxes() == a_phi.getBoxes());

    const ProblemDomain& domain = a_levGeo.getDomain();
    const DisjointBoxLayout grids = a_phi.getBoxes();
    DataIterator dit = a_phi.dataIterator();

    // Do coarse-fine BC's
    if (a_phiCrsePtr != NULL) {
        CH_assert(a_phi.nComp() == a_phiCrsePtr->nComp());
        CH_assert(a_cfInterpCrse.isDefined());

        a_cfInterpCrse.coarseFineInterp(a_phi, *a_phiCrsePtr);
    }

    Box validDomain = domain.domainBox();
#   ifndef USE_OLD_EXTRAP
    {
        // Without this copier, I've been having problems when
        // telling MeshRefine to vertically span the domain.
        Copier excp(grids, grids, domain, a_phi.ghostVect(), true);
        a_phi.exchange(excp);

        for (int dir = 0; dir < SpaceDim; ++dir) {
            if (domain.isPeriodic(dir)) {
                validDomain.grow(dir, 1);
            }
        }
    }
#   endif

    if (a_edgeGrad.nComp() == a_phi.nComp()) {
        // loop over boxes and compute gradient
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& thisPhi = a_phi[dit];
            FluxBox& thisEdgeGrad = a_edgeGrad[dit];

            Box validPhi = thisPhi.box();
#           ifdef USE_OLD_EXTRAP
                validPhi &= grids[dit];
#           else
                validPhi &= validDomain;
#           endif

            // loop over directions
            for (int dir = 0; dir < SpaceDim; ++dir) {
                FArrayBox& edgeGradDirFab = thisEdgeGrad[dir];

                // only do this in interior of grid
                Box edgeBox(grids[dit]);
                edgeBox.surroundingNodes(dir);

                // Only do one component
                int gradComp = 0;
                int phiComp = 0;
                int numComp = thisPhi.nComp();

                // for one-component gradient, only do normal gradients
                int edgeDir = dir;

                singleBoxMacGrad(edgeGradDirFab, thisPhi,
                                 gradComp, phiComp, numComp,
                                 edgeBox,
                                 validPhi,
                                 dir, edgeDir,
                                 dit(), a_levGeo,
                                 a_time, a_fluxBC);
            } // end loop over dir
        } // end loop over grids
    } // end if only doing normal directions
    else {
        // multicomponent gradPhi means that we also need to do
        // transverse directions.
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& thisPhi = a_phi[dit];
            FluxBox& thisEdgeGrad = a_edgeGrad[dit];

            Box validPhi = thisPhi.box();
#           ifdef USE_OLD_EXTRAP
                validPhi &= grids[dit];
#           else
                validPhi &= validDomain;
#           endif

            // loop over edges
            for (int edgeDir = 0; edgeDir < SpaceDim; ++edgeDir) {
                FArrayBox& thisEdgeGradDirFab = thisEdgeGrad[edgeDir];
                const Box& edgeBox = thisEdgeGradDirFab.box();

                // loop over component directions in edgeGrad
                // (this will be direction of gradient)
                for (int dir = 0; dir < thisEdgeGrad.nComp(); ++dir) {
                    int gradComp = dir;
                    int phiComp = 0;
                    int numComp = 1;
                    singleBoxMacGrad(thisEdgeGradDirFab, thisPhi,
                                     dir, phiComp, numComp,
                                     edgeBox,
                                     validPhi,
                                     dir, edgeDir,
                                     dit(), a_levGeo,
                                     a_time, a_fluxBC);
                }
            } // end loop over edgeDir
        } // end loop over grids
    } // end if we're also computing transverse direcctions
}


// -----------------------------------------------------------------------------
void Gradient::levelGradientMAC (LevelData<FluxBox>&         a_edgeGrad,
                                 const LevelData<FArrayBox>& a_phi,
                                 const LevelGeometry&        a_levGeo,
                                 const Real                  a_time,
                                 const BC_type*              a_fluxBC)

{
    CH_TIME("Gradient::levelGradientMAC_3");

    CH_assert(a_edgeGrad.nComp() >= a_phi.nComp());

    // Gather grid info
    const DisjointBoxLayout& grids = a_levGeo.getBoxes();
    DataIterator dit = grids.dataIterator();
    const ProblemDomain& domain = grids.physDomain();
    Box validDomain = domain.domainBox();

#   ifndef USE_OLD_EXTRAP
        ((LevelData<FArrayBox>&)a_phi).exchange(); // TODO: Is this needed?

        for (int dir = 0; dir < SpaceDim; ++dir) {
            if (domain.isPeriodic(dir)) {
                validDomain.grow(dir, 1);
            }
        }
#   endif

    if (a_edgeGrad.nComp() == a_phi.nComp()) {
        // We are only computing the normal dir

        for (dit.reset(); dit.ok(); ++dit) {
            const FArrayBox& phiFAB = a_phi[dit];

            Box validPhi = phiFAB.box();
#           ifdef USE_OLD_EXTRAP
                validPhi &= grids[dit];
#           else
                validPhi &= validDomain;
#           endif

            // loop over directions
            for (int dir = 0; dir < SpaceDim; ++dir) {
                FArrayBox& gradFAB = a_edgeGrad[dit][dir];
                const Box& edgeBox = gradFAB.box();

                // for one-component gradient, only do normal gradients
                int edgeDir = dir;
                int gradComp = 0;
                int phiComp = 0;
                int numComp = phiFAB.nComp();

                singleBoxMacGrad(gradFAB,
                                 phiFAB,
                                 gradComp, phiComp, numComp,
                                 edgeBox,
                                 validPhi,
                                 dir, edgeDir,
                                 dit(), a_levGeo,
                                 a_time, a_fluxBC);

            } // end loop over dir
        } // end loop over grids

    } else if (a_edgeGrad.nComp() == SpaceDim * a_phi.nComp()) {

        // multicomponent gradPhi means that we also need to do
        // transverse directions.

        CH_assert(a_phi.nComp() == 1); // Just to see if this is ever false.

        DataIterator dit = a_phi.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            const FArrayBox& phiFAB = a_phi[dit];

            Box validPhi = phiFAB.box();
#           ifdef USE_OLD_EXTRAP
                validPhi &= grids[dit];
#           else
                validPhi &= validDomain;
#           endif

            // Loop over FC directions
            for (int edgeDir = 0; edgeDir < SpaceDim; ++edgeDir) {
                FArrayBox& gradFAB = a_edgeGrad[dit][edgeDir];
                const Box& edgeBox = gradFAB.box();

                // Loop over gradient directions
                for (int gradDir = 0; gradDir < SpaceDim; ++gradDir) {

                    // Loop over the components of phi
                    for (int phiComp = 0; phiComp < a_phi.nComp(); ++phiComp) {
                        int gradComp = phiComp * SpaceDim + gradDir;
                        int numComp = 1;

                        singleBoxMacGrad(gradFAB, phiFAB,
                                         gradComp, phiComp, numComp,
                                         edgeBox,
                                         validPhi,
                                         gradDir, edgeDir,
                                         dit(), a_levGeo,
                                         a_time, a_fluxBC);

                    } // end loop over components
                } // end loop over gradient directions
            } // end loop over edge orientation
        } // end loop over grids
    } else {
        // bad number of components in either phi or gradPhi
        MayDay::Error("levelGradientMAC: bad number of components!");
    }
}


// -----------------------------------------------------------------------------
void Gradient::levelGradientCC (LevelData<FArrayBox>&       a_grad,
                                LevelData<FArrayBox>&       a_phi,
                                const LevelData<FArrayBox>* a_phiCrsePtr,
                                const LevelGeometry&        a_levGeo,
                                const Real                  a_time,
                                const BC_type*              a_fluxBC)
{
    CH_TIME("Gradient::levelGradientCC_1");

    // Create MappedQuadCFInterp object
    MappedQuadCFInterp cfInterpCrse;

    if (a_phiCrsePtr != NULL) {
        cfInterpCrse.define(a_phi.getBoxes(),
                            &(a_phiCrsePtr->getBoxes()),
                            a_levGeo.getDx(),
                            a_levGeo.getCrseRefRatio(),
                            a_phi.nComp(),
                            a_levGeo.getDomain());
    }

    // Call full gradient function
    levelGradientCC(a_grad,
                    a_phi,
                    a_phiCrsePtr,
                    cfInterpCrse,
                    a_levGeo,
                    a_time,
                    a_fluxBC);
}


#ifdef USE_SIMPLE_STENCIL
// -----------------------------------------------------------------------------
void Gradient::levelGradientCC (LevelData<FArrayBox>&         a_grad,
                                LevelData<FArrayBox>&         a_phi,
                                const LevelData<FArrayBox>*   a_phiCrsePtr,
                                MappedQuadCFInterp&           a_cfInterpCrse,
                                const LevelGeometry&          a_levGeo,
                                const Real                    a_time,
                                const BC_type*                a_fluxBC)
{
    CH_TIME("Gradient::levelGradientCC_2 (simple stencil)");

    // Collect geometric info
    const RealVect& dx = a_levGeo.getDx();
    const ProblemDomain& domain = a_levGeo.getDomain();
    const DisjointBoxLayout& grids = a_levGeo.getBoxes();
    DataIterator dit = grids.dataIterator();
    int ncomp = a_phi.nComp();

    // Sanity checks
    CH_assert(a_grad.nComp() == ncomp * SpaceDim);
    CH_assert(a_grad.getBoxes() == grids);
    CH_assert(a_phi .getBoxes() == grids);

    // Set CFBCs
    if (a_phiCrsePtr != NULL) {
        const LevelData<FArrayBox>& crseData = *a_phiCrsePtr;
        CH_assert(a_phi.nComp() == crseData.nComp());
        CH_assert(a_cfInterpCrse.isDefined());

        a_cfInterpCrse.coarseFineInterp(a_phi, crseData);
    }

    // Exchange
    a_phi.exchange();

    // Loop over grids and compute gradient
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& gradFAB = a_grad[dit];
        const FArrayBox& phiFAB = a_phi[dit];
        const FluxBox& JgupFB = a_levGeo.getFCJgup()[dit];
        const Box& valid = grids[dit];
        const RealVect& pdScale = RealVect::Unit;

        // Set physical BCs.
        if (a_fluxBC != NULL) {
            FArrayBox& phiFABRef = (FArrayBox&)phiFAB;

#ifdef USE_OLD_EXTRAP
            // Extrapolate ghosts for non-diagonal derivatives
            const int extrapOrder = 2;
            FArrayBox extrapFAB(phiFABRef.box(), phiFABRef.nComp());
            extrapFAB.copy(phiFABRef);

            Box thisValid = valid & a_phiFab.box();
            for (int fdir = 0; fdir < SpaceDim; ++fdir) {
                SideIterator fsit;
                for (fsit.reset(); fsit.ok(); ++fsit) {
                    ExtrapolateFaceAndCopy(extrapFAB, extrapFAB, thisValid,
                                           fdir, fsit(), extrapOrder);
                }
                thisValid.grow(fdir, 1);
            }
#else
            // Extrapolate ghosts for transverse derivatives
            FArrayBox extrapFAB(phiFABRef.box(), phiFABRef.nComp());
            if (!LevelGeometry::isDiagonal()) {
                const int extrapOrder = 2;
                extrapFAB.copy(phiFABRef);

                Box thisValid = valid;
                for (int fdir = 0; fdir < SpaceDim; ++fdir) {
                    SideIterator fsit;
                    for (fsit.reset(); fsit.ok(); ++fsit) {
                        ExtrapolateFaceAndCopy(extrapFAB, extrapFAB, thisValid,
                                               fdir, fsit(), extrapOrder);
                    }
                    thisValid.grow(fdir, 1);
                }
            }
#endif

            a_fluxBC->setGhosts(phiFABRef,              // stateFAB
                                &extrapFAB,             // &extrapFAB
                                valid,                  // valid box
                                domain,                 // ProblemDomain
                                dx,                     // derivScale
                                dit(),                  // DataIndex
                                &JgupFB,                // &JgupFB
                                false,                  // isHomogeneous
                                a_time,                 // time
                                phiFABRef.interval());  // interval
        }

        for (int dir = 0; dir < SpaceDim; ++dir) {
            FORT_MAPPEDCCGRADSCALE (
                CHF_FRA1(gradFAB,dir),
                CHF_CONST_FRA1(phiFAB,0),
                CHF_CONST_FRA1(phiFAB,0),
                CHF_CONST_FRA(JgupFB[dir]),
                CHF_BOX(valid),
                CHF_CONST_REALVECT(dx),
                CHF_CONST_REALVECT(pdScale),
                CHF_CONST_INT(dir));
        } // end loop over gradient directions (dir)

        // TODO: Fix at physical boundary using fluxBCs?

    } // end loop over grids (dit)
}

#else
// -----------------------------------------------------------------------------
void Gradient::levelGradientCC (LevelData<FArrayBox>&         a_grad,
                                LevelData<FArrayBox>&         a_phi,
                                const LevelData<FArrayBox>*   a_phiCrsePtr,
                                MappedQuadCFInterp&           a_cfInterpCrse,
                                const LevelGeometry&          a_levGeo,
                                const Real                    a_time,
                                const BC_type*                a_fluxBC)
{
    CH_TIME("Gradient::levelGradientCC_2");

    // First, compute FC gradient
    const DisjointBoxLayout& grids = a_grad.getBoxes();
    int ncompgrad = a_grad.nComp() / SpaceDim;
    LevelData<FluxBox> edgeGrad(grids, ncompgrad);

    levelGradientMAC(edgeGrad,
                     a_phi,
                     a_phiCrsePtr,
                     a_cfInterpCrse,
                     a_levGeo,
                     a_time,
                     a_fluxBC);

    // Then average edges->cells
    EdgeToCell(edgeGrad, a_grad);
}
#endif //USE_SIMPLE_STENCIL


// -----------------------------------------------------------------------------
void Gradient::levelGradientCC (LevelData<FArrayBox>&       a_grad,
                                const LevelData<FArrayBox>& a_phi,
                                const LevelGeometry&        a_levGeo,
                                const Real                  a_time,
                                const BC_type*              a_fluxBC)
{
    levelGradientCC(a_grad,
                    (LevelData<FArrayBox>&)a_phi,
                    NULL, //a_phiCrsePtr,
                    a_levGeo,
                    a_time,
                    a_fluxBC);
}


// Composite gradient functions ================================================

// -----------------------------------------------------------------------------
void Gradient::compGradientMAC (LevelData<FluxBox>&         a_edgeGrad,
                                LevelData<FArrayBox>&       a_phi,
                                const LevelData<FArrayBox>* a_phiCrsePtr,
                                const LevelData<FArrayBox>* a_phiFinePtr,
                                const LevelGeometry&        a_levGeo,
                                const Real                  a_time,
                                const BC_type*              a_fluxBC)
{
    CH_TIME("Gradient::compGradientMAC_1");

    levelGradientMAC(a_edgeGrad,
                     a_phi,
                     a_phiCrsePtr,
                     a_levGeo,
                     a_time,
                     a_fluxBC);
}


// -----------------------------------------------------------------------------
void Gradient::compGradientMAC (LevelData<FluxBox>&           a_edgeGrad,
                                LevelData<FArrayBox>&         a_phi,
                                const LevelData<FArrayBox>*   a_phiCrsePtr,
                                const LevelData<FArrayBox>*   a_phiFinePtr,
                                MappedQuadCFInterp&           a_cfInterpCrse,
                                const LevelGeometry&          a_levGeo,
                                const Real                    a_time,
                                const BC_type*                a_fluxBC)
{
    CH_TIME("Gradient::compGradientMAC_2");

    levelGradientMAC(a_edgeGrad,
                     a_phi,
                     a_phiCrsePtr,
                     a_cfInterpCrse,
                     a_levGeo,
                     a_time,
                     a_fluxBC);
}


// -----------------------------------------------------------------------------
void Gradient::compGradientCC (LevelData<FArrayBox>&       a_grad,
                               LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>* a_phiCrsePtr,
                               const LevelData<FArrayBox>* a_phiFinePtr,
                               const LevelGeometry&        a_levGeo,
                               const Real                  a_time,
                               const BC_type*              a_fluxBC)
{
    CH_TIME("Gradient::compGradientCC_1");

    // Create QuadCFInterp object
    MappedQuadCFInterp cfInterpCrse;

    if (a_phiCrsePtr != NULL) {
        cfInterpCrse.define(a_phi.getBoxes(),
                            &(a_phiCrsePtr->getBoxes()),
                            a_levGeo.getDx(),
                            a_levGeo.getCrseRefRatio(),
                            a_phi.nComp(),
                            a_levGeo.getDomain());
    }

    // Call full gradient function
    compGradientCC(a_grad,
                   a_phi,
                   a_phiCrsePtr,
                   a_phiFinePtr,
                   cfInterpCrse,
                   a_levGeo,
                   a_time,
                   a_fluxBC);
}


#ifdef USE_SIMPLE_STENCIL
// -----------------------------------------------------------------------------
void Gradient::compGradientCC (LevelData<FArrayBox>&         a_grad,
                               LevelData<FArrayBox>&         a_phi,
                               const LevelData<FArrayBox>*   a_phiCrsePtr,
                               const LevelData<FArrayBox>*   a_phiFinePtr,
                               MappedQuadCFInterp&           a_cfInterpCrse,
                               const LevelGeometry&          a_levGeo,
                               const Real                    a_time,
                               const BC_type*                a_fluxBC)
{
    CH_TIME("Gradient::compGradientCC_2 (simple stencil");

    // First, compute level gradient.
    levelGradientCC(a_grad,
                    a_phi,
                    a_phiCrsePtr,
                    a_cfInterpCrse,
                    a_levGeo,
                    a_time,
                    a_fluxBC);

    // Then, fix up the interface with the finer level...
    const DisjointBoxLayout& grids = a_levGeo.getBoxes();
    const ProblemDomain& domain = a_levGeo.getDomain();
    const IntVect& nRefFine = a_levGeo.getFineRefRatio();

    D_TERM(CH_assert(nRefFine[0] >= 1 || a_phiFinePtr == NULL);,
           CH_assert(nRefFine[1] >= 1 || a_phiFinePtr == NULL);,
           CH_assert(nRefFine[2] >= 1 || a_phiFinePtr == NULL);)

    if (a_phiFinePtr != NULL) {
        // do one-sided differencing at interfaces w/ finer levels
        // loop over this level grids
        const DisjointBoxLayout& fineGrids = a_phiFinePtr->getBoxes();

        // Do exchange to ensure that all edges are filled with appropriate data
        Copier exCopier;
        // exCopier.define(grids, grids, domain, a_grad.ghostVect(), true);
        exCopier.exchangeDefine(grids, a_grad.ghostVect());
        exCopier.trimEdges(grids, a_grad.ghostVect());
        a_grad.exchangeBegin(exCopier);

        // will need a mask for this one
        IntVect maskGrow(IntVect::Unit);
        maskGrow *= 2;
        LevelData<BaseFab<int> > masks(grids, 1, maskGrow);
        Mask thisMask;
        thisMask.buildMasks(masks, domain, grids, &fineGrids, nRefFine);

        // will need coarsened fine grids
        CH_assert(fineGrids.isClosed());
        BoxLayout coarsenedFineGrids;
        coarsen(coarsenedFineGrids, fineGrids, nRefFine);

        // End the exchange and perform the grad calculation on each grid.
        a_grad.exchangeEnd();

        DataIterator dit = a_grad.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& gradFAB = a_grad[dit];
            //const Box& thisGradBox = thisEdgeGrad.box();
            const Box& thisGradBox = a_grad.getBoxes()[dit];
            BaseFab<int>& thisMask = masks[dit];

            // now loop over (coarsened) fine boxes
            LayoutIterator litFine = coarsenedFineGrids.layoutIterator();

            for (litFine.reset(); litFine.ok(); ++litFine) {
                Box overlapBox(thisGradBox);
                // grow fine-grid box by one to make sure we catch case
                // where coarse-fine and coarse-coarse interfaces meet.
                const Box& crseFineBox = coarsenedFineGrids.get(litFine());
                Box testBox(crseFineBox);
                testBox.grow(1);
                overlapBox &= testBox;

                if (!overlapBox.isEmpty()) {
                    // fine grid overlays coarse grid, so we need to modify grad.

                    // loop over directions
                    for (int dir = 0; dir<SpaceDim; dir++) {
                        // figure out which edge faces we need to correct
                        Box loEdgeBox = adjCellLo(crseFineBox,dir,1);
                        Box hiEdgeBox = adjCellHi(crseFineBox,dir,1);

                        loEdgeBox &= thisGradBox;
                        hiEdgeBox &= thisGradBox;

                        loEdgeBox.shiftHalf(dir, 1);
                        hiEdgeBox.shiftHalf(dir,-1);

                        // check to see if we want to do correction
                        // in this direction, both for low and hi
                        int do_lo = 1;
                        int do_hi = 1;
                        // don't do correction if fine-grid and crse-grid
                        // edges coincide
                        if (overlapBox.smallEnd(dir)<=thisGradBox.smallEnd(dir))
                            do_lo = 0;
                        if (overlapBox.bigEnd(dir) >= thisGradBox.bigEnd(dir))
                            do_hi = 0;

                        FORT_SIMPLECRSEONESIDEGRAD(
                            CHF_FRA1(gradFAB,dir),
                            CHF_FIA1(thisMask,0),
                            CHF_BOX(loEdgeBox),
                            CHF_BOX(hiEdgeBox),
                            CHF_INT(dir),
                            CHF_INT(do_lo),
                            CHF_INT(do_hi));

                    } // end loop over directions
                } // end if fine grid overlays this grid
            } // end loop over fine grids
        } // end loop over this level's grids
    }  // end if a finer level exists
}

#else
// -----------------------------------------------------------------------------
void Gradient::compGradientCC (LevelData<FArrayBox>&         a_grad,
                               LevelData<FArrayBox>&         a_phi,
                               const LevelData<FArrayBox>*   a_phiCrsePtr,
                               const LevelData<FArrayBox>*   a_phiFinePtr,
                               MappedQuadCFInterp&           a_cfInterpCrse,
                               const LevelGeometry&          a_levGeo,
                               const Real                    a_time,
                               const BC_type*                a_fluxBC)
{
    CH_TIME("Gradient::compGradientCC_2");

    const ProblemDomain& domain = a_levGeo.getDomain();
    const IntVect& nRefFine = a_levGeo.getFineRefRatio();

    D_TERM(CH_assert(nRefFine[0] >= 1 || a_phiFinePtr == NULL);,
           CH_assert(nRefFine[1] >= 1 || a_phiFinePtr == NULL);,
           CH_assert(nRefFine[2] >= 1 || a_phiFinePtr == NULL);)

    // first compute edge-centered gradient
    const DisjointBoxLayout& grids = a_phi.getBoxes();
    int ncompGrad = a_grad.nComp() / SpaceDim;

    // nGhost is necessary to ensure that all edges for
    // extrapolation of face-centered gradient at
    // coarse-fine interpolation are present on the same
    // grid
    int nGhost = 1;
    LevelData<FluxBox> edgeGrad(grids, ncompGrad, nGhost*IntVect::Unit);

    compGradientMAC(edgeGrad,
                    a_phi,
                    a_phiCrsePtr,
                    a_phiFinePtr,
                    a_cfInterpCrse,
                    a_levGeo,
                    a_time,
                    a_fluxBC);

    if (a_phiFinePtr != NULL) {
        // do one-sided differencing at interfaces w/ finer levels
        // loop over this level grids
        const DisjointBoxLayout& fineGrids = a_phiFinePtr->getBoxes();

        // Do exchange to ensure that all edges are filled with appropriate data
        Copier exCopier;
        exCopier.define(grids, grids, domain, edgeGrad.ghostVect(), true);
        TODO();
        // exCopier.exchangeDefine(grids, edgeGrad.ghostVect());
        // exCopier.trimEdges(grids, edgeGrad.ghostVect());
        edgeGrad.exchangeBegin(exCopier);

        // will need a mask for this one
        IntVect maskGrow(IntVect::Unit);
        maskGrow *= 2;
        LevelData<BaseFab<int> > masks(grids, 1, maskGrow);
        Mask thisMask;
        thisMask.buildMasks(masks, domain, grids, &fineGrids, nRefFine);

        // will need coarsened fine grids
        CH_assert(fineGrids.isClosed());
        BoxLayout coarsenedFineGrids;
        coarsen(coarsenedFineGrids, fineGrids, nRefFine);

        // End the exchange and perform the grad calculation on each grid.
        edgeGrad.exchangeEnd();
        DataIterator dit = edgeGrad.dataIterator();

        for (dit.reset(); dit.ok(); ++dit) {
            FluxBox& thisEdgeGrad = edgeGrad[dit];
            //const Box& thisGradBox = thisEdgeGrad.box();
            const Box& thisGradBox = edgeGrad.getBoxes()[dit];
            BaseFab<int>& thisMask = masks[dit];

            // now loop over (coarsened) fine boxes
            LayoutIterator litFine = coarsenedFineGrids.layoutIterator();

            for (litFine.reset(); litFine.ok(); ++litFine) {
                Box overlapBox(thisGradBox);
                // grow fine-grid box by one to make sure we catch case
                // where coarse-fine and coarse-coarse interfaces meet.
                const Box& crseFineBox = coarsenedFineGrids.get(litFine());
                Box testBox(crseFineBox);
                testBox.grow(1);
                overlapBox &= testBox;

                if (!overlapBox.isEmpty()) {
                    // fine grid overlays coarse grid, so we need to modify grad.

                    // loop over directions
                    for (int dir = 0; dir<SpaceDim; dir++) {
                        FArrayBox& thisGradDir = thisEdgeGrad[dir];

                        // figure out which edge faces we need to correct
                        Box loEdgeBox = adjCellLo(crseFineBox,dir,1);
                        Box hiEdgeBox = adjCellHi(crseFineBox,dir,1);

                        loEdgeBox &= thisGradBox;
                        hiEdgeBox &= thisGradBox;

                        loEdgeBox.shiftHalf(dir, 1);
                        hiEdgeBox.shiftHalf(dir,-1);

                        // check to see if we want to do correction
                        // in this direction, both for low and hi
                        int do_lo = 1;
                        int do_hi = 1;
                        // don't do correction if fine-grid and crse-grid
                        // edges coincide
                        if (overlapBox.smallEnd(dir)<=thisGradBox.smallEnd(dir))
                            do_lo = 0;
                        if (overlapBox.bigEnd(dir) >= thisGradBox.bigEnd(dir))
                            do_hi = 0;

                        FORT_CRSEONESIDEGRAD(CHF_FRA1(thisGradDir,0),
                                             CHF_FIA1(thisMask,0),
                                             CHF_BOX(loEdgeBox),
                                             CHF_BOX(hiEdgeBox),
                                             CHF_INT(dir),
                                             CHF_INT(do_lo),
                                             CHF_INT(do_hi));
                    } // end loop over directions
                } // end if fine grid overlays this grid
            } // end loop over fine grids
        } // end loop over this level's grids
    }  // end if a finer level exists

    // now average to cells
    EdgeToCell(edgeGrad, a_grad);
}
#endif //USE_SIMPLE_STENCIL


// -----------------------------------------------------------------------------
void Gradient::compGradientCC (LevelData<FArrayBox>&       a_grad,
                               const LevelData<FArrayBox>& a_phi,
                               const LevelData<FArrayBox>* a_phiFinePtr,
                               const LevelGeometry&        a_levGeo,
                               const Real                  a_time,
                               const BC_type*              a_fluxBC)
{
    MappedQuadCFInterp dummyInterp;
    compGradientCC(a_grad,
                   (LevelData<FArrayBox>&)a_phi,
                   NULL, //a_phiCrsePtr,
                   a_phiFinePtr,
                   dummyInterp,
                   a_levGeo,
                   a_time,
                   a_fluxBC);
}


// Additional utilities ========================================================

// -----------------------------------------------------------------------------
// Use this to compute (u.Del)U...that is, CC.grad(FC)
// This function makes no assumptions on the scaling of a_CC or a_FC.
// -----------------------------------------------------------------------------
void Gradient::levelCCDotGradFC (LevelData<FArrayBox>&       a_CCDotGradFC,
                                 const LevelData<FArrayBox>& a_CC,
                                 const LevelData<FluxBox>&   a_FC,
                                 const LevelGeometry&        a_levGeo)
{
    CH_TIME("Gradient::levelCCDotGradFC");

    // Collect grids and dx
    const RealVect& dx = a_levGeo.getDx();
    const DisjointBoxLayout& grids = a_levGeo.getBoxes();

    // Sanity checks
    CH_assert(a_CCDotGradFC.nComp() == SpaceDim);
    CH_assert(a_CC.nComp() == SpaceDim);
    CH_assert(a_FC.nComp() == SpaceDim); // SpaceDim*SpaceDim total comps

    CH_assert(a_CCDotGradFC.getBoxes().compatible(grids));
    CH_assert(a_CC.getBoxes().compatible(grids));
    CH_assert(a_FC.getBoxes().compatible(grids));
    CH_assert(grids == a_CCDotGradFC.getBoxes());
    CH_assert(grids == a_CC.getBoxes());
    CH_assert(grids == a_FC.getBoxes());

    // Loop over grids and calculate
    for (DataIterator dit(grids); dit.ok(); ++dit) {
        FArrayBox& uDelU = a_CCDotGradFC[dit];
        const FArrayBox& CC = a_CC[dit];
        const FluxBox& FC = a_FC[dit];
        const Box& region = grids[dit];

#if CH_SPACEDIM == 2
        for (int dir0 = 0; dir0 < SpaceDim; ++dir0) {
            const int dir1 = (dir0 + 1) % SpaceDim;

            const FArrayBox& FC0 = FC[dir0];
            const FArrayBox& FC1 = FC[dir1];

            FORT_CCGRADFC2D(CHF_FRA1(uDelU,dir0),
                            CHF_CONST_FRA(CC),
                            CHF_CONST_FRA1(FC0,dir0),
                            CHF_CONST_FRA1(FC1,dir0),
                            CHF_BOX(region),
                            CHF_CONST_REALVECT(dx),
                            CHF_CONST_INT(dir0));
        } // end loop over dir0
#elif CH_SPACEDIM == 3

        for (int dir0 = 0; dir0 < SpaceDim; ++dir0) {
            const int dir1 = (dir0 + 1) % SpaceDim;
            const int dir2 = (dir0 + 2) % SpaceDim;

            const FArrayBox& FC0 = FC[dir0];
            const FArrayBox& FC1 = FC[dir1];
            const FArrayBox& FC2 = FC[dir2];

            FORT_CCGRADFC3D(CHF_FRA1(uDelU,dir0),
                            CHF_CONST_FRA(CC),
                            CHF_CONST_FRA1(FC0,dir0),
                            CHF_CONST_FRA1(FC1,dir0),
                            CHF_CONST_FRA1(FC2,dir0),
                            CHF_BOX(region),
                            CHF_CONST_REALVECT(dx),
                            CHF_CONST_INT(dir0));
        } // end loop over dir0
#else
#   error Bad CH_SPACEDIM
#endif
    } // end loop over grids (dit)
}


// -----------------------------------------------------------------------------
// Utility function to do the actual computin' (to reduce code duplication)
// -----------------------------------------------------------------------------
void Gradient::singleBoxMacGrad (FArrayBox&            a_gradFab,
                                 const FArrayBox&      a_phiFab,
                                 int                   a_gradComp,
                                 int                   a_phiComp,
                                 int                   a_numComp,
                                 const Box&            a_edgeBox,
                                 const Box&            a_validPhi,
                                 int                   a_gradDir,
                                 int                   a_edgeDir,
                                 const DataIndex&      a_di,
                                 const LevelGeometry&  a_levGeo,
                                 const Real            a_time,
                                 const BC_type*        a_fluxBC)
{
    CH_TIME("Gradient::singleBoxMacGrad");

#ifndef NDEBUG
   // Set gradient to bogus value
   a_gradFab.setVal(quietNAN, a_gradFab.box(), a_gradComp, a_numComp);
#endif

    const RealVect& dx = a_levGeo.getDx();

#ifdef USE_OLD_EXTRAP
    // Extrapolate ghosts for non-diagonal derivatives
    const int extrapOrder = 2;
    FArrayBox extrapFAB(a_phiFab.box(), a_phiFab.nComp());
    extrapFAB.copy(a_phiFab);

    Box valid = a_levGeo.getBoxes()[a_di] & a_phiFab.box();
    for (int fdir = 0; fdir < SpaceDim; ++fdir) {
        SideIterator fsit;
        for (fsit.reset(); fsit.ok(); ++fsit) {
            ExtrapolateFaceAndCopy(extrapFAB, extrapFAB, valid,
                                   fdir, fsit(), extrapOrder);
        }
        valid.grow(fdir, 1);
    }
#else
    // Extrapolate ghosts for transverse derivatives
    FArrayBox extrapFAB(a_phiFab.box(), a_phiFab.nComp());
    if (!a_levGeo.isDiagonal() || a_edgeDir != a_gradDir) {
        const int extrapOrder = 2;
        extrapFAB.copy(a_phiFab);

        Box valid = a_validPhi;
        for (int fdir = 0; fdir < SpaceDim; ++fdir) {
            SideIterator fsit;
            for (fsit.reset(); fsit.ok(); ++fsit) {
                ExtrapolateFaceAndCopy(extrapFAB, extrapFAB, valid,
                                       fdir, fsit(), extrapOrder);
            }
            valid.grow(fdir, 1);
        }
    }
#endif

    // Fill all ghosts of a_phi, if requested to do so.
    // This is more efficient than filling ghosts before calling grad
    // because it avoids an extra extrapolation.
    if (a_fluxBC != NULL) {
        const DisjointBoxLayout& grids = a_levGeo.getBoxes();
        const ProblemDomain& domain = a_levGeo.getDomain();
        const FluxBox& JgupFB = a_levGeo.getFCJgup()[a_di];
        FArrayBox& phiFABRef = (FArrayBox&)a_phiFab;

        a_fluxBC->setGhosts(phiFABRef,              // stateFAB
                            &extrapFAB,             // &extrapFAB
                            grids[a_di],            // valid box
                            domain,                 // ProblemDomain
                            dx,                     // derivScale
                            a_di,                   // DataIndex
                            &JgupFB,                // &JgupFB
                            false,                  // isHomogeneous
                            a_time,                 // time
                            phiFABRef.interval());  // interval
    }

    // We need a SpaceDim-component FAB to be filled with Jgup^{a_gradDir,*}
    // at every point a_gradFab is to be calculated.
    //FArrayBox JgupFAB(a_gradFab.box(), SpaceDim);
    //a_levGeo.fill_Jgup(JgupFAB, a_gradDir);
    FArrayBox JgupFAB;
    {
        const FArrayBox& JgupCachedFAB = a_levGeo.getFCJgup()[a_di][a_gradDir];
        const Box& JgupCacheBox = JgupCachedFAB.box();

        // Can we use the cached metric?
        if (JgupCacheBox.type() == a_gradFab.box().type()) {
            if (JgupCacheBox.contains(a_gradFab.box())) {
                // Use the cache
                JgupFAB.define(JgupCachedFAB.interval(), (FArrayBox&)JgupCachedFAB);
            } else {
                // The region of interest lies outside of the cached region.
                // Use an expensive fill function. (This never hits.)
                JgupFAB.define(a_gradFab.box(), SpaceDim);
                a_levGeo.fill_Jgup(JgupFAB, a_gradDir);
            }
        } else {
            // This hits many times.
            JgupFAB.define(a_gradFab.box(), SpaceDim);
            a_levGeo.fill_Jgup(JgupFAB, a_gradDir);
        }
    }

    // Loop over gradFab components
    int phiComp = a_phiComp;
    for (int comp = a_gradComp; comp < a_gradComp + a_numComp; ++comp) {

        // Compute gradient on this grid
        if (!a_levGeo.isDiagonal()) {
            // Use full non-orthogonal version
            FORT_MAPPEDMACGRAD(CHF_FRA1(a_gradFab, comp),
                               CHF_CONST_FRA1(a_phiFab, phiComp),
                               CHF_CONST_FRA1(extrapFAB, phiComp),
                               CHF_CONST_FRA(JgupFAB),
                               CHF_BOX(a_edgeBox),
                               CHF_CONST_REALVECT(dx),
                               CHF_INT(a_gradDir),
                               CHF_INT(a_edgeDir));
        } else {
            // Use simpler orthogonal version
            const Real dxDir = dx[a_gradDir];
            FORT_MAPPEDMACGRADORTHO(CHF_FRA1(a_gradFab, comp),
                                    CHF_CONST_FRA1(a_phiFab, phiComp),
                                    CHF_CONST_FRA1(extrapFAB, phiComp),
                                    CHF_CONST_FRA(JgupFAB),
                                    CHF_BOX(a_edgeBox),
                                    CHF_CONST_REAL(dxDir),
                                    CHF_INT(a_gradDir),
                                    CHF_INT(a_edgeDir));
        }

        // Set boundary fluxes on this component, if needed
        if (a_fluxBC != NULL) {
            const Box& valid = a_levGeo.getBoxes()[a_di];
            const ProblemDomain& domain = a_levGeo.getDomain();
            Interval interv = Interval(a_phiComp, a_phiComp + a_numComp - 1);

            a_fluxBC->setFluxes(a_gradFab,      // state
                                NULL,           // extrap
                                valid,          // valid
                                domain,         // domain
                                dx,             // dx
                                a_di,           // dataIndex
                                NULL,           // &JgupFB
                                a_gradDir,      // edge dir
                                false,          // isHomogeneous
                                a_time,         // time
                                interv);        // interval (comps of state)
        } // end if applying BCs

        phiComp++;
    } // end loop over gradFab components
}

