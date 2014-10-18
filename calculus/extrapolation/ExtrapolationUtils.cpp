#include "ExtrapolationUtils.H"
#include "ExtrapolationUtilsF_F.H"
#include "FluxBox.H"
#include "CFRegion.H"
#include "Debug.H"


// -----------------------------------------------------------------------------
// Extrapolation boundary conditions for a side.
// -----------------------------------------------------------------------------
void ExtrapolateFaceNoEV (FArrayBox&           a_dest,
                          const FArrayBox&     a_state,
                          const Box&           a_valid,
                          const int            a_dir,
                          const Side::LoHiSide a_side,
                          const int            a_order)
{
    CH_TIME("ExtrapolateFaceNoEV");

    // Sanity checks
    CH_assert(0 <= a_order && a_order <= 4);
    CH_assert(a_dest.box().type() == a_valid.type());
    CH_assert(a_dest.nComp() == a_state.nComp());
    CH_assert(a_dest.box().type() == a_state.box().type());

    int isign = sign(a_side);

    // Find ghost region
    Box toRegion = adjCellBox(a_valid, a_dir, a_side, 1);
    if (a_valid.type(a_dir) == IndexType::NODE) {
        toRegion.shiftHalf(a_dir, -isign);
    }
    toRegion &= a_state.box();
    if (toRegion.isEmpty()) return;

    // Do extrapolation
    FORT_EXTRAPOLATEFACENOEV(CHF_FRA(a_dest),
                             CHF_FRA(a_state),
                             CHF_BOX(toRegion),
                             CHF_CONST_INT(a_dir),
                             CHF_CONST_INT(isign),
                             CHF_CONST_INT(a_order));
}


// -----------------------------------------------------------------------------
// Extrapolation BCs. Does not fill edges and vertices.
// -----------------------------------------------------------------------------
void ExtrapolateBCNoEV (FArrayBox&           a_state,
                        const Box&           a_valid,
                        const ProblemDomain& a_domain,
                        const int            a_order,
                        const IntVect&       a_doLoSides,
                        const IntVect&       a_doHiSides)
{
    CH_TIME("ExtrapolateBCNoEV");

    Box thisDom = a_domain.domainBox();
    thisDom.convert(a_valid.type());

    for (int fdir = 0; fdir < SpaceDim; ++fdir) {
        if (a_domain.isPeriodic(fdir)) continue;

        if (a_doLoSides[fdir] != 0) {
            const Side::LoHiSide fside = Side::Lo;
            if (thisDom.sideEnd(fside)[fdir] == a_valid.sideEnd(fside)[fdir]) {
                ExtrapolateFaceNoEV(a_state, a_state, a_valid, fdir, fside, a_order);
            }
        }
        if (a_doHiSides[fdir] != 0) {
            const Side::LoHiSide fside = Side::Hi;
            if (thisDom.sideEnd(fside)[fdir] == a_valid.sideEnd(fside)[fdir]) {
                ExtrapolateFaceNoEV(a_state, a_state, a_valid, fdir, fside, a_order);
            }
        }
    } // end loop over FC dirs (fdir)
}


// -----------------------------------------------------------------------------
// Extrapolation boundary conditions for a side.
// This works out of place and also extrapolates onto faces a_numLayers deep.
// Edge and vertex ghosts will be extrapolated from ghosts adjacent to boundary.
// If a_dest != a_state, the inner layers will also be copied.
// -----------------------------------------------------------------------------
void ExtrapolateFaceAndCopy (FArrayBox&           a_dest,
                             const FArrayBox&     a_state,
                             const Box&           a_valid,
                             const int            a_dir,
                             const Side::LoHiSide a_side,
                             const int            a_order,
                             const int            a_numLayers)
{
    CH_TIME("ExtrapolateFaceAndCopy");

    const int ncomps = a_state.nComp();
    CH_assert(a_dest.nComp() == ncomps);

    // Extrapolate face
    ExtrapolateFaceNoEV(a_dest, a_state, a_valid, a_dir, a_side, a_order);

    Box ghostBox = adjCellBox(a_valid, a_dir, a_side, 1);
    ghostBox &= a_state.box();

    Box nearBox = ghostBox;
    nearBox.shift(a_dir, -sign(a_side));
    nearBox.growDir(a_dir, flip(a_side), a_numLayers-1);

    // Copy the inner layer, no ghosts
    if (&a_dest != &a_state) {
        a_dest.copy(a_state, nearBox);
    }

    // Extrapolate edges and vertices around face and inner layer
    for (int edir = 0; edir < SpaceDim; ++edir) {
        if (edir == a_dir) continue;

        SideIterator esit;
        for (esit.reset(); esit.ok(); ++esit) {
            Side::LoHiSide eside = esit();

            ExtrapolateFaceNoEV(a_dest, a_dest, ghostBox, edir, eside, a_order);
            ExtrapolateFaceNoEV(a_dest, a_dest, nearBox, edir, eside, a_order);
        } // end loop over sides (esit)

        ghostBox.grow(BASISV(edir));
        ghostBox &= a_state.box();

        nearBox.grow(BASISV(edir));
        nearBox &= a_state.box();
    } // end loop over dirs (edir)
}


// -----------------------------------------------------------------------------
// Extrapolates edge and vertex ghosts at the fine side of the CF interface.
// This is used to finish the job that MappedQuadCFInterp started.
// For future reference, the coeff of val_k using an order n extrap is
//   (-1)^(k+1) * Binomial[n+1,k]
// where val_1 is the nearVal and val_(n+1) is the farVal.
// -----------------------------------------------------------------------------
void ExtrapolateCFEV (FArrayBox&           a_state,
                      const CFRegion&      a_region,
                      const DataIndex&     a_di,
                      const int            a_dir,
                      const Side::LoHiSide a_side,
                      const Interval&      a_interval,
                      const int            a_order,
                      const IntVect&       a_activeDirs)
{
    CH_TIME("ExtrapolateCFEV");
    CH_assert(0 <= a_order && a_order <= 4);

    const int isign = sign(a_side);

    const CFIVS* cfivs = NULL;
    CFRegion& cfregionRef = const_cast<CFRegion&>(a_region);

    if (a_side == Side::Lo) {
        cfivs = &(cfregionRef.loCFIVS(a_di, a_dir));
    } else {
        cfivs = &(cfregionRef.hiCFIVS(a_di, a_dir));
    }

    TODO();
    // if (!cfivs->isPacked()) {
    //     const Box& stateBox = a_state.box();
    //     Box faceBox = cfivs->minBox();
    //     faceBox &= stateBox;
    //     if (faceBox.isEmpty()) return;

    //     const Box& edgeBox = faceBox;
    //     const int vdir = 1 - a_dir;
    //     IntVect v = BASISV(vdir);

    //     IntVect cc = edgeBox.sideEnd(Side::Lo);
    //     a_state(cc-v) = -1000.0;
    //     cc = edgeBox.sideEnd(Side::Hi);
    //     a_state(cc+v) = 1000.0;
    //     return;
    // }

    const Box& stateBox = a_state.box();

    Box faceBox = cfivs->minBox();
    // Box faceBox = cfivs->packedBox();
    faceBox &= stateBox;
    if (faceBox.isEmpty()) return;

#if CH_SPACEDIM == 2
    const Box& edgeBox = faceBox;
    const int vdir = 1 - a_dir;

    if (a_activeDirs[vdir] == 0) return;

    IntVect v = BASISV(vdir);
    IntVect cc = edgeBox.sideEnd(Side::Lo);

    if (a_order == 0) {
        a_state(cc-v) = a_state(cc);

        cc = edgeBox.sideEnd(Side::Hi);
        a_state(cc+v) = a_state(cc);

    } else if (a_order == 1) {
        a_state(cc-v) = linearExtrap(a_state(cc),
                                     a_state(cc+v));

        cc = edgeBox.sideEnd(Side::Hi);
        a_state(cc+v) = linearExtrap(a_state(cc),
                                     a_state(cc-v));

    } else if (a_order == 2) {
        a_state(cc-v) = quadraticExtrap(a_state(cc),
                                        a_state(cc+v),
                                        a_state(cc+2*v));

        cc = edgeBox.sideEnd(Side::Hi);
        a_state(cc+v) = quadraticExtrap(a_state(cc),
                                        a_state(cc-v),
                                        a_state(cc-2*v));

    } else if (a_order == 3) {
        a_state(cc-v) = cubicExtrap(a_state(cc),
                                    a_state(cc+v),
                                    a_state(cc+2*v),
                                    a_state(cc+3*v));

        cc = edgeBox.sideEnd(Side::Hi);
        a_state(cc+v) = cubicExtrap(a_state(cc),
                                    a_state(cc-v),
                                    a_state(cc-2*v),
                                    a_state(cc-3*v));

    } else if (a_order == 4) {
        a_state(cc-v) = quarticExtrap(a_state(cc),
                                      a_state(cc+v),
                                      a_state(cc+2*v),
                                      a_state(cc+3*v),
                                      a_state(cc+4*v));

        cc = edgeBox.sideEnd(Side::Hi);
        a_state(cc+v) = quarticExtrap(a_state(cc),
                                      a_state(cc-v),
                                      a_state(cc-2*v),
                                      a_state(cc-3*v),
                                      a_state(cc-4*v));

    } else {
        MayDay::Error("ExtrapolateCFEV: Bad order");
    }

#else
    for (int edir = 0; edir < SpaceDim; ++edir) {
        if (edir == a_dir) continue;
        if (a_activeDirs[edir] == 0) continue;

        for (SideIterator esit; esit.ok(); ++esit) {
            const Side::LoHiSide eside = esit();
            const int esign = sign(eside);

            Box edgeBox = adjCellBox(faceBox, edir, eside, 1);
            edgeBox &= stateBox;
            if (edgeBox.isEmpty()) continue;

            FORT_EXTRAPOLATEFACENOEV (
                CHF_FRA(a_state),
                CHF_FRA(a_state),
                CHF_BOX(edgeBox),
                CHF_CONST_INT(edir),
                CHF_CONST_INT(esign),
                CHF_CONST_INT(a_order));

            const int vdir = SpaceDim - edir - a_dir;
            if (a_activeDirs[vdir] == 0) continue;

            IntVect v = BASISV(vdir);
            IntVect cc = edgeBox.sideEnd(Side::Lo);

            if (a_order == 0) {
                a_state(cc-v) = a_state(cc);

                cc = edgeBox.sideEnd(Side::Hi);
                a_state(cc+v) = a_state(cc);

            } else if (a_order == 1) {
                a_state(cc-v) = linearExtrap(a_state(cc),
                                             a_state(cc+v));

                cc = edgeBox.sideEnd(Side::Hi);
                a_state(cc+v) = linearExtrap(a_state(cc),
                                             a_state(cc-v));

            } else if (a_order == 2) {
                a_state(cc-v) = quadraticExtrap(a_state(cc),
                                                a_state(cc+v),
                                                a_state(cc+2*v));

                cc = edgeBox.sideEnd(Side::Hi);
                a_state(cc+v) = quadraticExtrap(a_state(cc),
                                                a_state(cc-v),
                                                a_state(cc-2*v));

            } else if (a_order == 3) {
                a_state(cc-v) = cubicExtrap(a_state(cc),
                                            a_state(cc+v),
                                            a_state(cc+2*v),
                                            a_state(cc+3*v));

                cc = edgeBox.sideEnd(Side::Hi);
                a_state(cc+v) = cubicExtrap(a_state(cc),
                                            a_state(cc-v),
                                            a_state(cc-2*v),
                                            a_state(cc-3*v));

            } else if (a_order == 4) {
                a_state(cc-v) = quarticExtrap(a_state(cc),
                                              a_state(cc+v),
                                              a_state(cc+2*v),
                                              a_state(cc+3*v),
                                              a_state(cc+4*v));

                cc = edgeBox.sideEnd(Side::Hi);
                a_state(cc+v) = quarticExtrap(a_state(cc),
                                              a_state(cc-v),
                                              a_state(cc-2*v),
                                              a_state(cc-3*v),
                                              a_state(cc-4*v));

            } else {
                MayDay::Error("ExtrapolateCFEV: Bad order");
            }
        } // end loop over edge sides (esit)
    } // end loop over edge directions (edir)
#endif
}


// -----------------------------------------------------------------------------
// LevelData version.
// -----------------------------------------------------------------------------
void ExtrapolateCFEV (LevelData<FArrayBox>& a_phi,
                      const CFRegion&       a_region,
                      const int             a_order,
                      const IntVect&        a_activeDirs)
{
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (a_activeDirs[dir] == 0) continue;

        DataIterator dit = a_phi.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& phiFAB = a_phi[dit];
            const Interval& interv = phiFAB.interval();

            ExtrapolateCFEV(phiFAB, a_region, dit(), dir, Side::Lo, interv, a_order, a_activeDirs);
            ExtrapolateCFEV(phiFAB, a_region, dit(), dir, Side::Hi, interv, a_order, a_activeDirs);
        }
    }
}


// -----------------------------------------------------------------------------
// Debugging function.
// Extrapolate/exchange all ghosts. Useful for testing if ghosts need to be set.
// This function does no exchanges.
// -----------------------------------------------------------------------------
void extrapAllGhosts (LevelData<FArrayBox>& a_data,
                      const int             a_order)
{
    const IntVect ghostVect = a_data.ghostVect();
    const DisjointBoxLayout& grids = a_data.getBoxes();
    DataIterator dit = a_data.dataIterator();

#if CH_SPACEDIM == 1
    const int numghosts = ghostVect[0];
#else
#if CH_SPACEDIM == 2
    const int numghosts = Min(ghostVect[0],ghostVect[1]);
#else
    const int numghosts = Min(Min(ghostVect[0],ghostVect[1]),ghostVect[2]);
#endif
#endif

    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& dataFAB = a_data[dit];
        Box valid = grids[dit];

        for (int ghost = 0; ghost < numghosts; ++ghost) {
            for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
                SideIterator sit;
                for (sit.reset(); sit.ok(); ++sit) {
                    Side::LoHiSide iside = sit();
                    ExtrapolateFaceNoEV(dataFAB, dataFAB, valid, dir, iside, a_order);
                }
                valid.grow(dir, 1);
            }
        }
    }
}


// -----------------------------------------------------------------------------
// Debugging function.
// Extrapolate/exchange all ghosts. Useful for testing if ghosts need to be set.
// This function does no exchanges.
// -----------------------------------------------------------------------------
void extrapAllGhosts (LevelData<FluxBox>& a_data,
                      const int           a_order)
{
    const IntVect ghostVect = a_data.ghostVect();
    const DisjointBoxLayout& grids = a_data.getBoxes();
    DataIterator dit = a_data.dataIterator();

#if CH_SPACEDIM == 1
    const int numghosts = ghostVect[0];
#else
#if CH_SPACEDIM == 2
    const int numghosts = Min(ghostVect[0],ghostVect[1]);
#else
    const int numghosts = Min(Min(ghostVect[0],ghostVect[1]),ghostVect[2]);
#endif
#endif

    for (dit.reset(); dit.ok(); ++dit) {
        for (int fdir = 0; fdir < CH_SPACEDIM; ++fdir) {
            FArrayBox& dataFAB = a_data[dit][fdir];

            Box valid = grids[dit];
            valid.surroundingNodes(fdir);

            for (int ghost = 0; ghost < numghosts; ++ghost) {
                for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
                    SideIterator sit;
                    for (sit.reset(); sit.ok(); ++sit) {
                        Side::LoHiSide iside = sit();
                        int isign = sign(iside);

                        Box toRegion = bdryBox(valid, dir, iside, 1);
                        if (dir == fdir) {
                            toRegion.shift(dir, isign);
                        } else {
                            toRegion.shiftHalf(dir, isign);
                        }
                        CH_assert(toRegion.sameType(dataFAB.box()));
                        toRegion &= dataFAB.box();
                        CH_assert(!toRegion.isEmpty());

#ifndef NDEBUG
                        Box stencilBox = toRegion;
                        stencilBox.growDir(dir, flip(iside), a_order + 1);
                        CH_assert(dataFAB.box().contains(stencilBox));
#endif

                        FORT_EXTRAPOLATEFACENOEV(CHF_FRA(dataFAB),
                                                 CHF_FRA(dataFAB),
                                                 CHF_BOX(toRegion),
                                                 CHF_CONST_INT(dir),
                                                 CHF_CONST_INT(isign),
                                                 CHF_CONST_INT(a_order));
                    }
                    valid.grow(dir, 1);
                }
            }
        }
    }
}


// -----------------------------------------------------------------------------
// Extrapolates ghosts at the fine side of the CF interface.
// This is a special version that works out of place.
// For use in MappedAMRPoissonOp.
// For future reference, the coeff of val_k using an order n extrap is
//   (-1)^(k+1) * Binomial[n+1,k]
// where val_1 is the nearVal and val_(n+1) is the farVal.
// -----------------------------------------------------------------------------
void ExtrapolateCF (FArrayBox&       a_dest,
                    const FArrayBox& a_state,
                    const CFRegion&  a_region,
                    const DataIndex& a_di,
                    int              a_dir,
                    Side::LoHiSide   a_side,
                    const Interval&  a_interval,
                    int              a_order)
{
    CH_TIME("ExtrapolateCF 2");
    CH_assert(0 <= a_order && a_order <= 4);

    int isign = sign(a_side);

    const CFIVS* cfivs = NULL;
    CFRegion& cfregionRef = const_cast<CFRegion&>(a_region);

    if (a_side == Side::Lo) {
        cfivs = &(cfregionRef.loCFIVS(a_di, a_dir));
    } else {
        cfivs = &(cfregionRef.hiCFIVS(a_di, a_dir));
    }


    if (cfivs->isPacked()) {
        Box toRegion = cfivs->packedBox();
        if (!toRegion.isEmpty()) {
            FORT_EXTRAPOLATEFACENOEV(CHF_FRA(a_dest),
                                     CHF_FRA(a_state),
                                     CHF_BOX(toRegion),
                                     CHF_CONST_INT(a_dir),
                                     CHF_CONST_INT(isign),
                                     CHF_CONST_INT(a_order));
        }
    } else {
        if (a_order == 1) {
            const IntVectSet& region = cfivs->getIVS();
            for (IVSIterator ivsit(region); ivsit.ok(); ++ivsit) {
                const IntVect& ivTo = ivsit();
                const IntVect ivNear = ivTo -   isign*BASISV(a_dir);
                const IntVect ivFar  = ivTo - 2*isign*BASISV(a_dir);

                for (int icomp = a_interval.begin(); icomp <= a_interval.end(); ++icomp) {
                    const Real nearVal = a_state(ivNear, icomp);
                    const Real farVal  = a_state(ivFar,  icomp);

                    a_dest(ivTo, icomp) = linearExtrap(nearVal, farVal);
                }
            }
        } else if (a_order == 2) {
            const IntVectSet& region = cfivs->getIVS();
            for (IVSIterator ivsit(region); ivsit.ok(); ++ivsit) {
                const IntVect& ivTo = ivsit();
                const IntVect ivNear = ivTo -   isign*BASISV(a_dir);
                const IntVect ivMid  = ivTo - 2*isign*BASISV(a_dir);
                const IntVect ivFar  = ivTo - 3*isign*BASISV(a_dir);

                for (int icomp = a_interval.begin(); icomp <= a_interval.end(); ++icomp) {
                    const Real nearVal = a_state(ivNear, icomp);
                    const Real midVal  = a_state(ivMid,  icomp);
                    const Real farVal  = a_state(ivFar,  icomp);

                    a_dest(ivTo, icomp) = quadraticExtrap(nearVal, midVal, farVal);
                }
            }
        } else if (a_order == 3) {
            const IntVectSet& region = cfivs->getIVS();
            for (IVSIterator ivsit(region); ivsit.ok(); ++ivsit) {
                const IntVect& ivTo = ivsit();
                const IntVect ivNear      = ivTo -   isign*BASISV(a_dir);
                const IntVect ivMid       = ivTo - 2*isign*BASISV(a_dir);
                const IntVect ivFar       = ivTo - 3*isign*BASISV(a_dir);
                const IntVect ivReallyFar = ivTo - 4*isign*BASISV(a_dir);

                for (int icomp = a_interval.begin(); icomp <= a_interval.end(); ++icomp) {
                    const Real nearVal      = a_state(ivNear,      icomp);
                    const Real midVal       = a_state(ivMid,       icomp);
                    const Real farVal       = a_state(ivFar,       icomp);
                    const Real reallyFarVal = a_state(ivReallyFar, icomp);

                    a_dest(ivTo, icomp) = cubicExtrap(nearVal, midVal, farVal, reallyFarVal);
                }
            }
        } else if (a_order == 4) {
            const IntVectSet& region = cfivs->getIVS();
            for (IVSIterator ivsit(region); ivsit.ok(); ++ivsit) {
                const IntVect& ivTo = ivsit();
                const IntVect iv1 = ivTo -   isign*BASISV(a_dir);
                const IntVect iv2 = ivTo - 2*isign*BASISV(a_dir);
                const IntVect iv3 = ivTo - 3*isign*BASISV(a_dir);
                const IntVect iv4 = ivTo - 4*isign*BASISV(a_dir);
                const IntVect iv5 = ivTo - 5*isign*BASISV(a_dir);

                for (int icomp = a_interval.begin(); icomp <= a_interval.end(); ++icomp) {
                    const Real val1 = a_state(iv1, icomp);
                    const Real val2 = a_state(iv2, icomp);
                    const Real val3 = a_state(iv3, icomp);
                    const Real val4 = a_state(iv4, icomp);
                    const Real val5 = a_state(iv5, icomp);

                    a_dest(ivTo, icomp) = quarticExtrap(val1, val2, val3, val4, val5);
                }
            }
        } else {
            MayDay::Error("bogus order argument");
        }
    }
}


// -----------------------------------------------------------------------------
// Extrapolates ghosts at the fine side of the CF interface.
// This works in place and does all dirs and sides.
// For use in MappedAMRPoissonOp.
// -----------------------------------------------------------------------------
void ExtrapolateAllCF(FArrayBox&       a_state,
                      const CFRegion&  a_region,
                      const DataIndex& a_di,
                      int              a_order)
{
    CH_TIME("ExtrapolateAllCF");

    Interval interv = Interval(0, a_state.nComp()-1);
    for (int idir = 0; idir < SpaceDim; ++idir) {
        for (SideIterator sit; sit.ok(); ++sit) {
            ExtrapolateCF(a_state, a_state, a_region, a_di,
                          idir, sit(), interv, a_order);
        }
    }
}
