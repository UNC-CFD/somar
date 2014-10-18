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
#include "EllipticBCUtils.H"
#include "EllipticBCUtilsF_F.H"
#include "ExtrapolationUtils.H"
#include "LevelGeometry.H"
#include "BoxIterator.H"
#include "Debug.H"



// ***************************** Single-side BCs *******************************

// -----------------------------------------------------------------------------
// Sets Dirichlet BCs on a side
// For a CC state, this fills ghosts.
// For a FC state, this sets faces at the physical boundary.
// a_state and a_valid must have the same centering.
// -----------------------------------------------------------------------------
void setSideDiriBC (FArrayBox&           a_state,
                    const Box&           a_valid,
                    const ProblemDomain& a_domain,
                    const Real           a_value,
                    const int            a_dir,
                    const Side::LoHiSide a_side,
                    const bool           a_homogeneous,
                    const int            a_order)
{
    // Sanity checks
    CH_assert(0 <= a_dir);
    CH_assert(a_dir < SpaceDim);

    // We shouldn't touch periodic dirs.
    if (a_domain.isPeriodic(a_dir)) return;

    // Set up some basic parameters.
    const Real bcval = a_homogeneous? 0.0: a_value;
    const int isign = sign(a_side);

    // Ensure valid and state boxes are of the same type
    // and that valid is a subset of the state box.
    const IntVect boxType = a_state.box().type();
    CH_assert(a_valid.type() == boxType);
    CH_assert(a_state.box().contains(a_valid));

    // Create a domain box with the same centering as state.
    Box thisDom = a_domain.domainBox();
    thisDom.convert(boxType);

#ifndef NDEBUG
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (!a_domain.isPeriodic(dir)) {
            // Check that the valid box is within the domain extents.
            CH_assert(thisDom.smallEnd(dir) <= a_valid.smallEnd(dir));
            CH_assert(a_valid.bigEnd(dir) <= thisDom.bigEnd(dir));
        } else {
            // Do nothing in periodic dirs. We are already assured
            // that the valid box is contained within the state box.
            // That's good enough since all points, even the ghosts in
            // this dir, are in the valid domain.
        }
    }
#endif

    // Compute the physical boundary faces contained within
    // the region of interest.
    Box destBox = bdryBox(thisDom, a_dir, a_side, 1)
                & surroundingNodes(a_valid, a_dir);

    // If destBox is empty, then we are not at a physical boundary
    // and have nothing to do.
    if (destBox.isEmpty()) return;

    if (boxType[a_dir] == IndexType::NODE) {
        // a_state is FC. Just set the values directly...

        CH_assert(a_state.box().contains(destBox));
        a_state.setVal(bcval, destBox, 0, a_state.nComp());
    } else {
        // a_state is CC. Set ghost values...

        // Shift destBox to overlap with the ghost cells of interest.
        destBox.shiftHalf(a_dir, isign);

        // It is possible this state box has no ghosts.
        // If so, just leave.
        destBox &= a_state.box();
        if (destBox.isEmpty()) return;

        // Fill the ghosts.
        FORT_ELLIPTICCONSTDIRIBCGHOST(CHF_FRA(a_state),
                                      CHF_BOX(destBox),
                                      CHF_CONST_REAL(bcval),
                                      CHF_CONST_INT(a_dir),
                                      CHF_CONST_INT(isign),
                                      CHF_CONST_INT(a_order));
    }
}


// -----------------------------------------------------------------------------
// Sets the outward pointing normal covariant derivative of a_state on a side.
// For a CC state, this fills ghosts.
// For a FC state, this throws an error.
// a_state and a_valid must have the same centering.
// -----------------------------------------------------------------------------
void setSideNeumBC (FArrayBox&           a_state,
                    const Box&           a_valid,
                    const ProblemDomain& a_domain,
                    const Real           a_value,
                    const int            a_dir,
                    const Side::LoHiSide a_side,
                    const bool           a_homogeneous,
                    const FArrayBox&     a_Jgupi,
                    const FArrayBox*     a_extrapPtr,
                    const RealVect&      a_dx,
                    const DataIndex&     a_index)
{
    CH_TIME("setSideNeumBC");

    // Sanity checks
    CH_assert(a_state.box().type() == IntVect::Zero);
    CH_assert(a_valid.type() == IntVect::Zero);
    CH_assert(a_state.nComp() == 1); // For now, this is all that is needed.

    // If the domain is periodic in this direction, scram.
    if (a_domain.isPeriodic(a_dir)) return;

    // Find the boundary locations. If we are not at a boundary, scram.
    const Box& domBox = a_domain.domainBox();
    const Box faceBox = bdryBox(a_valid, a_dir, a_side, 1)
                      & bdryBox(domBox, a_dir, a_side, 1);
    if (faceBox.isEmpty()) return;

    // Find the ghost cells
    const int isign = sign(a_side);
    Box ghostBox = faceBox;
    ghostBox.shiftHalf(a_dir, isign);
    ghostBox &= a_state.box();
    if (ghostBox.isEmpty()) return;

    if (LevelGeometry::isDiagonal()) {
        // This metric is diagonal. Just set the ghosts.
        const Real dxDir = a_dx[a_dir];

        FORT_ELLIPTICCONSTNEUMBCGHOSTORTHO (
            CHF_FRA(a_state),
            CHF_CONST_FRA1(a_Jgupi, a_dir),
            CHF_BOX(ghostBox),
            CHF_CONST_REAL(a_value),
            CHF_CONST_INT(a_dir),
            CHF_CONST_INT(isign),
            CHF_CONST_REAL(dxDir));

    } else {
        // Grab the boundary normals.
        const FArrayBox& nhat = a_Jgupi;

        // Allocate space for extrapolated values...
        // Used to compute non-diagonal derivatives at boundaries.
        // This is identical to using one-sided derivative stencils.
        FArrayBox extrap;
        FArrayBox* thisExtrapPtr = NULL;
        if (a_extrapPtr == NULL) {
            thisExtrapPtr = new FArrayBox(a_state.box(), a_state.nComp());
            extrap.define(thisExtrapPtr->interval(), *thisExtrapPtr);
        } else {
            Interval interv(0, a_state.nComp()-1);
            thisExtrapPtr = const_cast<FArrayBox*>(a_extrapPtr);
            extrap.define(interv, *thisExtrapPtr);
        }

        // Extrapolate the boundary values for cross derivatives
        const int extrapOrder = 2;
        ExtrapolateFaceAndCopy(extrap, a_state, a_valid, a_dir, a_side, extrapOrder);

        // Set BCs
        FORT_ELLIPTICCONSTNEUMBCGHOST (
            CHF_FRA(a_state),
            CHF_CONST_FRA(extrap),
            CHF_CONST_FRA(nhat),
            CHF_BOX(ghostBox),
            CHF_CONST_REAL(a_value),
            CHF_CONST_INT(a_dir),
            CHF_CONST_INT(isign),
            CHF_CONST_REALVECT(a_dx));

        // Delete extrap space if we created it
        if (a_extrapPtr == NULL) {
            delete thisExtrapPtr;
        }
    }
}


// -----------------------------------------------------------------------------
// Sets extrapolation BCs on a side
// For a CC state, this fills ghosts.
// For a FC state, this sets faces at the physical boundary.
// a_state and a_valid must have the same centering.
// This function does NOT check for periodicity!
// -----------------------------------------------------------------------------
void setSideExtrapBC (FArrayBox&           a_state,
                      const Box&           a_valid,
                      const ProblemDomain& a_domain,
                      const int            a_dir,
                      const Side::LoHiSide a_side,
                      const int            a_order)
{
    // Set up some basic parameters.
    const int isign = sign(a_side);

    // Ensure valid and state boxes are of the same type
    // and that valid is a subset of the state box.
    const IntVect boxType = a_state.box().type();
    CH_assert(a_valid.type() == boxType);
    CH_assert(a_state.box().contains(a_valid));

    // Create a domain box with the same centering as state.
    Box thisDom = a_domain.domainBox();
    thisDom.convert(boxType);

#ifndef NDEBUG
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (!a_domain.isPeriodic(dir)) {
            // Check that the valid box is within the domain extents.
            CH_assert(thisDom.smallEnd(dir) <= a_valid.smallEnd(dir));
            CH_assert(a_valid.bigEnd(dir) <= thisDom.bigEnd(dir));
        } else {
            // Do nothing in periodic dirs. We are already assured
            // that the valid box is contained within the state box.
            // That's good enough since all points, even the ghosts in
            // this dir, are in the valid domain.
        }
    }
#endif

    // Compute the physical boundary faces contained within
    // the region of interest.
    Box destBox = bdryBox(thisDom, a_dir, a_side, 1)
                & surroundingNodes(a_valid, a_dir);

    // If destBox is empty, then we are not at a physical boundary
    // and have nothing to do.
    if (destBox.isEmpty()) return;

    if (boxType[a_dir] == IndexType::NODE) {
        // a_state is FC. Just set the values directly...

        // Since we are setting faces at the physical boundary and not in a
        // ghost layer, destBox had better be a subset of the state box!
        CH_assert(a_state.box().contains(destBox));

        // Fill the ghosts.
        FORT_ELLIPTICEXTRAPBCGHOST(
            CHF_FRA(a_state),
            CHF_BOX(destBox),
            CHF_CONST_INT(a_dir),
            CHF_CONST_INT(isign),
            CHF_CONST_INT(a_order));

    } else {
        // a_state is CC. Set ghost values...

        // Shift destBox to overlap with the ghost cells of interest.
        destBox.shiftHalf(a_dir, isign);

        // It is possible this state box has no ghosts.
        // If so, just leave.
        destBox &= a_state.box();
        if (destBox.isEmpty()) return;

        // Fill the ghosts.
        FORT_ELLIPTICEXTRAPBCGHOST(
            CHF_FRA(a_state),
            CHF_BOX(destBox),
            CHF_CONST_INT(a_dir),
            CHF_CONST_INT(isign),
            CHF_CONST_INT(a_order));
    }
}


// -----------------------------------------------------------------------------
// Sets extrapolation BCs on a side.
// This version extrapolates CC data to faces and does not work in place.
// a_dest must be FC in a_dir.
// a_src and a_valid must be CC a_dir.
// a_dest, a_src, and a_valid must have the same centering in other directions.
// This function does NOT check for periodicity!
// -----------------------------------------------------------------------------
void setSideExtrapBC (FArrayBox&           a_dest,
                      const FArrayBox&     a_src,
                      const Box&           a_valid,
                      const ProblemDomain& a_domain,
                      const int            a_dir,
                      const Side::LoHiSide a_side,
                      const int            a_order)
{
    // Set up some basic parameters.
    const int isign = sign(a_side);

    // In this function, it is assumed that a_dest is FC in a_dir and a_src
    // is CC in a_dir.
    CH_assert(a_dest.box().type()[a_dir] == 1);
    CH_assert(a_src.box().type()[a_dir] == 0);

    // Ensure valid and src boxes are of the same type
    // and that valid is a subset of the src box.
    CH_assert(a_src.box().type() == a_valid.type());
    CH_assert(a_src.box().contains(a_valid));

    // Create a domain box with the same centering as dest.
    const IntVect destBoxType = a_dest.box().type();
    Box thisDom = a_domain.domainBox();
    thisDom.convert(destBoxType);

    // Create a valid box with the same centering as dest.
    Box thisValid = a_valid;
    thisValid.convert(destBoxType);

    // The valid box had better be in the physical domain!
    CH_assert(thisDom.contains(thisValid));

    // Compute the physical boundary faces contained within
    // the region of interest.
    Box destBox = bdryBox(thisDom, a_dir, a_side, 1) & thisValid;

    // If destBox is empty, then we are not at a physical boundary
    // and have nothing to do.
    if (destBox.isEmpty()) return;

    // Make sure we have enough room to perform the calculation.
#ifndef NDEBUG
    {
        CH_assert(a_dest.box().contains(destBox));
        Box srcBox = destBox;
        srcBox.shiftHalf(a_dir, -isign);
        srcBox.growDir(a_dir, flip(a_side), a_order);
        CH_assert(a_src.box().contains(srcBox));
    }
#endif

    // Fill the ghosts.
    FORT_ELLIPTICEXTRAPBCCELLTOFACE(
        CHF_FRA(a_dest),
        CHF_CONST_FRA(a_src),
        CHF_BOX(destBox),
        CHF_CONST_INT(a_dir),
        CHF_CONST_INT(isign),
        CHF_CONST_INT(a_order));
}


// *************************** Complete BC methods *****************************

// -----------------------------------------------------------------------------
// A simple override of BCGhostClass to set constant Diri BCs.
// For CC dirs, fills ghosts. For FC dirs, sets bdry data directly.
// -----------------------------------------------------------------------------
void EllipticConstDiriBCGhostClass::operator() (FArrayBox&           a_state,
                                                const FArrayBox*     a_extrapPtr,   // Just a dummy
                                                const Box&           a_valid,
                                                const ProblemDomain& a_domain,
                                                const RealVect&      a_dx,          // Just a dummy
                                                const DataIndex&     a_index,       // Just a dummy
                                                const FluxBox*       a_JgupPtr,     // Just a dummy
                                                bool                 a_homogeneous,
                                                Real                 a_time,        // Just a dummy
                                                const Interval&      a_interval) const
{
    CH_TIME("EllipticConstDiriBCGhostClass::operator()");

    // This is the only thing the setSide* functions care about.
    CH_assert(a_state.box().type() == a_valid.type());

    // Alias the components we want to alter
    const Interval interv = (a_interval == Interval())? a_state.interval(): a_interval;
    CH_assert(a_state.nComp() > interv.end());
    FArrayBox stateAlias(interv, a_state);

    // Set BCs in each non-periodic direction.
    for (int idir = 0; idir < CH_SPACEDIM; ++idir) {
        if (a_domain.isPeriodic(idir)) continue;

        // For this direction, set BCs on each requested side.
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide iside = sit();
            if (m_doSide[iside][idir] == 0) continue;

            // Set the BCs
            setSideDiriBC(stateAlias,
                          a_valid,
                          a_domain,
                          m_BCValue[iside][idir],
                          idir,
                          iside,
                          a_homogeneous,
                          1); // order
        }
    }
}


// -----------------------------------------------------------------------------
// A simple override of BCGhostClass that fills ghosts using constant Neum BCs.
// NOTE: It only makes sense to use this on CC data!
// -----------------------------------------------------------------------------
void EllipticConstNeumBCGhostClass::operator() (FArrayBox&           a_state,
                                                const FArrayBox*     a_extrapPtr,
                                                const Box&           a_valid,
                                                const ProblemDomain& a_domain,
                                                const RealVect&      a_dx,
                                                const DataIndex&     a_index,       // Just a dummy
                                                const FluxBox*       a_JgupPtr,
                                                bool                 a_homogeneous,
                                                Real                 a_time,        // Just a dummy
                                                const Interval&      a_interval) const
{
    CH_TIME("EllipticConstNeumBCGhostClass::operator()");

    // We will definitely need the metric to set Neumann BCs.
    CH_assert(a_JgupPtr != NULL);

    // This is the only thing the setSide* functions care about.
    CH_assert(a_state.box().type() == a_valid.type());

    // Alias the components we want to alter
    const Interval interv = (a_interval == Interval())? a_state.interval(): a_interval;
    CH_assert(a_state.nComp() > interv.end());
    FArrayBox stateAlias(interv, a_state);

    // Set BCs in each non-periodic direction.
    for (int idir = 0; idir < SpaceDim; ++idir) {
        if (a_domain.isPeriodic(idir)) continue;

        // For this direction, set BCs on each requested side.
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            Side::LoHiSide iside = sit();
            if (m_doSide[iside][idir] == 0) continue;

            // Set the BCs
            setSideNeumBC(stateAlias,
                          a_valid,
                          a_domain,
                          m_BCValue[iside][idir],
                          idir,
                          iside,
                          a_homogeneous,
                          (*a_JgupPtr)[idir],
                          a_extrapPtr,
                          a_dx,
                          a_index);
        }
    }
}


// -----------------------------------------------------------------------------
// A simple override of BCFluxClass that sets fluxes using constant Neum BCs.
// This is used to set BCs after taking the gradient of a field.
//
// NOTE: This function was not really made to apply zero-normal BCs to a
// velocity. For that, consider using a zero EllipticConstDiriBCGhostClass.
//
// NOTE: This can only set face values, not ghosts!
// -----------------------------------------------------------------------------
void EllipticConstNeumBCFluxClass::operator() (FArrayBox&           a_state,
                                               const FArrayBox*     a_extrapPtr,    // Just a dummy
                                               const Box&           a_valid,
                                               const ProblemDomain& a_domain,
                                               const RealVect&      a_dx,           // Just a dummy
                                               const DataIndex&     a_index,        // Just a dummy
                                               const FluxBox*       a_Jgup,         // Just a dummy
                                               int                  a_dir,
                                               bool                 a_homogeneous,
                                               Real                 a_time,         // Just a dummy
                                               const Interval&      a_interval) const
{
    CH_TIME("EllipticConstNeumBCFluxClass::operator()");

    // For the *BCFluxClass stuff, we require the state box to be FC and
    // the valid box to be CC.
    CH_assert(a_state.box().type() == BASISV(a_dir));
    CH_assert(a_valid.type() == IntVect::Zero);

    // We can think of this function as setting Dirichlet BCs on a gradient.
    // So, we can adjust the centering of the valid box and use the Diri code.
    Box FCValid = a_valid;
    FCValid.surroundingNodes(a_dir);

    // Alias the components we want to alter
    const Interval interv = (a_interval == Interval())? a_state.interval(): a_interval;
    CH_assert(a_state.nComp() > interv.end());
    FArrayBox stateAlias(interv, a_state);

    // Set BCs in each non-periodic direction.
    for (int idir = 0; idir < CH_SPACEDIM; ++idir) {
        if (a_domain.isPeriodic(idir)) continue;

        // For this direction, set BCs on each requested side.
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide iside = sit();
            if (m_doSide[iside][idir] == 0) continue;

            // Set the BCs
            setSideDiriBC(stateAlias,
                          FCValid,
                          a_domain,
                          m_BCValue[iside][idir],
                          idir,
                          iside,
                          a_homogeneous,
                          1); // order (this is just a dummy value)

        }
    }
}


// -----------------------------------------------------------------------------
// A simple override of BCGhostClass that fills ghosts using general Diri BCs.
// -----------------------------------------------------------------------------
void EllipticDiriBCGhostClass::operator() (FArrayBox&           a_state,
                                           const FArrayBox*     a_extrapPtr,    // Just a dummy
                                           const Box&           a_valid,
                                           const ProblemDomain& a_domain,
                                           const RealVect&      a_dx,
                                           const DataIndex&     a_index,        // Just a dummy
                                           const FluxBox*       a_JgupPtr,      // Just a dummy
                                           bool                 a_homogeneous,
                                           Real                 a_time,
                                           const Interval&      a_interval) const
{
    CH_TIME("EllipticDiriBCGhostClass::operator()");

    // Sanity checks
    const IntVect boxType = a_state.box().type();
    CH_assert(boxType == a_valid.type());

    Box thisDom = a_domain.domainBox();
    thisDom.convert(boxType);

    const Interval interv = (a_interval == Interval())? a_state.interval(): a_interval;
    const int scomp = interv.begin();
    const int ncomp = interv.size();
    CH_assert(a_state.nComp() > interv.end());

    for (int idir = 0; idir < SpaceDim; ++idir) {
        if (a_domain.isPeriodic(idir)) continue;

        const IntVect e = BASISV(idir);
        Box extDom = thisDom;
        extDom.grow(4*(IntVect::Unit-BASISV(idir)));

        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            Side::LoHiSide iside = sit();
            Real isign = sign(iside);

            if (m_doSide[iside][idir] == 0) continue;

            // Find the bdry locations
            Box faceBox = bdryBox(a_valid, idir, iside, 1)
                        & bdryBox(extDom, idir, iside, 1);
            if (faceBox.isEmpty()) continue;

            if (boxType[idir] == 0) {
                // CC state: Fill the ghost cell

                // Calculate the FC offsets
                RealVect offset = 0.5 * (RealVect::Unit - RealVect(boxType));
                offset[idir] = 0.0;

                // Find the source data locations
                Box srcBox = faceBox;
                srcBox.shiftHalf(idir, -isign);
                srcBox &= a_state.box();
                CH_assert(!srcBox.isEmpty());

                // Find the ghost locations
                Box ghostBox = faceBox;
                ghostBox.shiftHalf(idir, isign);
                ghostBox &= a_state.box();
                if (ghostBox.isEmpty()) continue;

                // Set the ghosts (1st order)
                a_state.copy(a_state, srcBox, scomp, ghostBox, scomp, ncomp);
                a_state.negate(ghostBox, scomp, ncomp);

                if (!a_homogeneous) {
                    Real* bcval = new Real[ncomp];
                    Real pos[CH_SPACEDIM];
                    BoxIterator bit(ghostBox);

                    for (bit.reset(); bit.ok(); ++bit) {
                        const IntVect& ghostCC = bit();
                        const IntVect fc = ghostCC + (1-iside)*e;

                        D_TERM(pos[0] = (Real(fc[0]) + offset[0]) * a_dx[0];,
                               pos[1] = (Real(fc[1]) + offset[1]) * a_dx[1];,
                               pos[2] = (Real(fc[2]) + offset[2]) * a_dx[2];)

                        if (m_funcPtr != NULL) {
                            m_funcPtr(pos, &idir, &iside, bcval, a_dx, a_time);
                        } else {
                            (*m_objPtr)(pos, &idir, &iside, bcval, a_dx, a_time);
                        }

                        for (int n = scomp; n < ncomp; ++n) {
                            // 1st-order ghost extrapolation
                            a_state(ghostCC,n) += 2.0 * bcval[n];

                            // 2nd-order ghost extrapolation
                            // IntVect ncc = ghostCC - iside * e;
                            // IntVect fcc = ghostCC - 2.0 * iside * e;
                            // a_state(ghostCC,n) = (8.0*bcval[n] - 6.0*a_state(ncc,n) + a_state(fcc,n)) / 3.0;
                        }
                    }

                    delete[] bcval;
                } // end if not homogeneous

            } else {
                // FC state: Set the BCs directly

                // Get BC val
                if (a_homogeneous) {
                    // Set the state BCs
                    faceBox &= a_state.box();
                    a_state.setVal(0.0, faceBox, scomp, ncomp);

                } else {
                    Real* bcval = new Real[ncomp];
                    Real pos[CH_SPACEDIM];

                    RealVect offset = 0.5 * (RealVect::Unit - RealVect(boxType));
                    offset[idir] = 0.0;

                    faceBox &= a_state.box();
                    BoxIterator bit(faceBox);

                    for (bit.reset(); bit.ok(); ++bit) {
                        const IntVect& fc = bit();

                        D_TERM(pos[0] = (Real(fc[0]) + offset[0]) * a_dx[0];,
                               pos[1] = (Real(fc[1]) + offset[1]) * a_dx[1];,
                               pos[2] = (Real(fc[2]) + offset[2]) * a_dx[2];)

                        if (m_funcPtr != NULL) {
                            m_funcPtr(pos, &idir, &iside, bcval, a_dx, a_time);
                        } else {
                            (*m_objPtr)(pos, &idir, &iside, bcval, a_dx, a_time);
                        }

                        for (int n = scomp; n < ncomp; ++n) {
                            a_state(fc,n) = bcval[n];
                        }
                    }

                    delete[] bcval;
                } // end if not homogeneous
            } // end if CC or FC
        } // end loop over sides
    } // end loop over dirs
}


// -----------------------------------------------------------------------------
// A simple override of BCGhostClass that fills ghosts using general Neum BCs.
// -----------------------------------------------------------------------------
void EllipticNeumBCGhostClass::operator() (FArrayBox&           a_state,
                                           const FArrayBox*     a_extrapPtr,
                                           const Box&           a_valid,
                                           const ProblemDomain& a_domain,
                                           const RealVect&      a_dx,
                                           const DataIndex&     a_index,
                                           const FluxBox*       a_JgupPtr,
                                           bool                 a_homogeneous,
                                           Real                 a_time,
                                           const Interval&      a_interval) const
{
    CH_TIME("EllipticNeumBCGhostClass::operator()");

    // Sanity checks
    CH_assert(a_state.box().type() == IntVect::Zero);
    CH_assert(a_valid.type() == IntVect::Zero);
    CH_assert(a_JgupPtr != NULL);

    // Alias only the components we want to alter
    const Interval interv = (a_interval == Interval())? a_state.interval(): a_interval;
    const int scomp = interv.begin();
    const int ncomp = interv.size();
    CH_assert(a_state.nComp() > interv.end());
    FArrayBox stateAlias(interv, a_state);

    // Allocate space for extrapolated values.
    // Used to compute non-diagonal derivatives at boundaries.
    // This is identical to using one-sided derivative stencils.
    FArrayBox extrapAlias;
    FArrayBox* thisExtrapPtr = NULL;
    if (a_extrapPtr == NULL) {
        thisExtrapPtr = new FArrayBox(a_state.box(), ncomp);
        extrapAlias.define(thisExtrapPtr->interval(), *thisExtrapPtr);
    } else {
        thisExtrapPtr = const_cast<FArrayBox*>(a_extrapPtr);
        extrapAlias.define(interv, *thisExtrapPtr);
    }

    // Calculate off-diagonal derivatives to 2nd order
    const int extrapOrder = 2;

    Real* bcval = new Real[a_state.nComp()];
    Real pos[CH_SPACEDIM];
    if (a_homogeneous) {
        for (int icomp = 0; icomp < a_state.nComp(); ++icomp) {
            bcval[icomp] = 0.0;
        }
    }

    for (int idir = 0; idir < SpaceDim; ++idir) {
        if (a_domain.isPeriodic(idir)) continue;

        const FArrayBox& Jgupi = (*a_JgupPtr)[idir];
        D_TERM(
            ;,
            const int jdir = (idir + 1) % SpaceDim;,
            const int kdir = (idir + 2) % SpaceDim;
        )
        D_TERM(
            ;,
            const IntVect j = BASISV(jdir);,
            const IntVect k = BASISV(kdir);
        )
        CH_assert(idir != jdir);
#if CH_SPACEDIM == 3
        CH_assert(idir != kdir);
        CH_assert(jdir != kdir);
#endif

        // Calculate the FC offsets
        RealVect offset = 0.5 * RealVect::Unit;
        offset[idir] -= 0.5;

        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            Side::LoHiSide iside = sit();
            int isign = sign(iside);

            if (m_doSide[iside][idir] == 0) continue;

            // Find the bdry locations
            const Box faceBox = bdryBox(a_valid, idir, iside, 1)
                              & bdryBox(a_domain.domainBox(), idir, iside, 1);
            if (faceBox.isEmpty()) continue;

            // Find the source data locations
            Box srcBox = faceBox;
            srcBox.shiftHalf(idir, -isign);
            CH_assert(a_state.box().contains(srcBox));

            // Find the ghost locations
            Box ghostBox = srcBox;
            ghostBox.shift(idir, isign);
            CH_assert(a_state.box().contains(ghostBox));

            // Extrapolate the boundary values for cross derivatives
            ExtrapolateFaceAndCopy(extrapAlias, stateAlias, a_valid, idir, iside, extrapOrder);

            // Gather BC values
            FArrayBox bcValFAB;
            if (m_bdryDataPtr != NULL) {
                const FArrayBox& fieldFAB = m_bdryDataPtr->getData(a_index, idir, iside);

                CH_assert(!fieldFAB.box().isEmpty());   // If this trips, we have no values at this boundary.
                CH_assert(fieldFAB.box() == faceBox);
                CH_assert(ncomp == 1); // For now.

                bcValFAB.define(Interval(0,0), (FArrayBox&)fieldFAB);
            }

            for (BoxIterator bit(ghostBox); bit.ok(); ++bit) {
                const IntVect& ccGhost = bit();
                const IntVect ccNear = ccGhost - isign*BASISV(idir);
                const IntVect fcBdry = ((iside == Side::Hi)? ccGhost: (ccGhost + BASISV(idir)));

                if (!a_homogeneous) {
                    D_TERM(pos[0] = (Real(fcBdry[0]) + offset[0]) * a_dx[0];,
                           pos[1] = (Real(fcBdry[1]) + offset[1]) * a_dx[1];,
                           pos[2] = (Real(fcBdry[2]) + offset[2]) * a_dx[2];)

                    if (m_funcPtr != NULL) {
                        m_funcPtr(pos, &idir, &iside, bcval, a_dx, a_time);
                    } else if (m_bdryDataPtr != NULL) {
                        bcval[0] = bcValFAB(fcBdry);
                    } else {
                        (*m_objPtr)(pos, &idir, &iside, bcval, a_dx, a_time);
                    }
                }

                for (int icomp = 0; icomp < ncomp; icomp++) {
                    Real cross = 0.25 * (D_TERM(
                        0.0,

                        + (+ extrapAlias(ccGhost + j, icomp)
                           - extrapAlias(ccGhost - j, icomp)
                           + extrapAlias(ccNear + j, icomp)
                           - extrapAlias(ccNear - j, icomp)) * Jgupi(fcBdry,jdir) / a_dx[jdir],

                        + (+ extrapAlias(ccGhost + k, icomp)
                           - extrapAlias(ccGhost - k, icomp)
                           + extrapAlias(ccNear + k, icomp)
                           - extrapAlias(ccNear - k, icomp)) * Jgupi(fcBdry,kdir) / a_dx[kdir]
                    ));

                    Real bcvalComp = bcval[icomp + scomp];

                    stateAlias(ccGhost, icomp) = stateAlias(ccNear, icomp)
                                               + Real(isign) * a_dx[idir] * (bcvalComp - cross) / Jgupi(fcBdry,idir);

                } // end loop over comps
            } // end loop over ghost box
        } // end loop over sides
    } // end loop over dirs

    delete[] bcval;
    // Delete extrap space if we created it
    if (a_extrapPtr == NULL) {
        delete thisExtrapPtr;
    }
}


// -----------------------------------------------------------------------------
// A simple override of BCFluxClass that sets fluxes using general Neum BCs.
// -----------------------------------------------------------------------------
void EllipticNeumBCFluxClass::operator() (FArrayBox&           a_state,
                                          const FArrayBox*     a_extrapPtr,     // Just a dummy
                                          const Box&           a_valid,
                                          const ProblemDomain& a_domain,
                                          const RealVect&      a_dx,
                                          const DataIndex&     a_index,         // Just a dummy
                                          const FluxBox*       a_JgupPtr,       // Just a dummy
                                          int                  a_dir,
                                          bool                 a_homogeneous,
                                          Real                 a_time,
                                          const Interval&      a_interval) const
{
    CH_TIME("EllipticNeumBCFluxClass::operator()");

    if (a_domain.isPeriodic(a_dir)) return;

    // Sanity checks
    CH_assert(a_state.box().type() == BASISV(a_dir));
    CH_assert(a_valid.type() == IntVect::Zero);

    const Interval interv = (a_interval == Interval())? a_state.interval(): a_interval;
    const int scomp = interv.begin();
    const int ncomp = interv.size();
    CH_assert(a_state.nComp() > interv.end());

    // Calculate the FC offsets
    RealVect offset = 0.5 * RealVect::Unit;
    offset[a_dir] -= 0.5;

    SideIterator sit;
    for (sit.reset(); sit.ok(); ++sit) {
        Side::LoHiSide iside = sit();
        Real isign = sign(iside);

        if (m_doSide[iside][a_dir] == 0) continue;

        // Find the bdry locations
        const Box faceBox = bdryBox(a_valid, a_dir, iside, 1)
                          & bdryBox(a_domain.domainBox(), a_dir, iside, 1)
                          & a_state.box();

        if (!faceBox.isEmpty()) {
            // Set the fluxes
            if (a_homogeneous) {
                a_state.setVal(0.0, faceBox, scomp, ncomp);
            } else {
                if (m_fluxPtr != NULL) {
                    const FArrayBox& fieldFAB = (*m_fluxPtr)[a_index][a_dir];
                    a_state.copy(fieldFAB, faceBox, 0, faceBox, scomp, ncomp);

                } else if (m_bdryDataPtr != NULL) {
                    const FArrayBox& fieldFAB = m_bdryDataPtr->getData(a_index, a_dir, iside);
                    const Box fieldBox = fieldFAB.box();

                    CH_assert(!fieldBox.isEmpty());                 // If this trips, we have no values at this boundary.
                    CH_assert(fieldBox.size() == faceBox.size());   // If this trips, the src and dest grids may differ.
                    CH_assert(fieldBox.type() == BASISV(a_dir));    // If this trips, BoundaryData isn't returning a FC box.

                    a_state.copy(fieldFAB, fieldBox, 0, faceBox, scomp, ncomp);

                } else {
                    Real* bcval = new Real[a_state.nComp()];
                    Real pos[CH_SPACEDIM];

                    BoxIterator bit(faceBox);
                    for (bit.reset(); bit.ok(); ++bit) {
                        const IntVect& fc = bit();
                        D_TERM(pos[0] = (Real(fc[0]) + offset[0]) * a_dx[0];,
                               pos[1] = (Real(fc[1]) + offset[1]) * a_dx[1];,
                               pos[2] = (Real(fc[2]) + offset[2]) * a_dx[2];)

                        if (m_funcPtr != NULL) {
                            m_funcPtr(pos, &a_dir, &iside, bcval, a_dx, a_time);
                        } else {
                            (*m_objPtr)(pos, &a_dir, &iside, bcval, a_dx, a_time);
                        }

                        for (int n = scomp; n < scomp + ncomp; ++n) {
                            a_state(fc,n) = bcval[n];
                        } // end loop over comps
                    } // end loop over faceBox

                    delete[] bcval;
                }
            } // end if homogeneous
        } // end if faceBox not empty
    } // end loop over sides
}


// -----------------------------------------------------------------------------
// A simple override of BCFluxClass that sets fluxes to zero if a_homogeneous.
// Does nothing if !a_homogeneous.
//
// NOTE: This function was designed to be used in an iterative solver to
// eliminate errors introduced by taking the gradient over extrapolated ghosts.
// -----------------------------------------------------------------------------
void EllipticDoNothingUnlessHomogNeumBCFluxClass::operator() (FArrayBox&           a_state,
                                                              const FArrayBox*     a_extrapPtr,     // Just a dummy
                                                              const Box&           a_valid,
                                                              const ProblemDomain& a_domain,
                                                              const RealVect&      a_dx,            // Just a dummy
                                                              const DataIndex&     a_index,         // Just a dummy
                                                              const FluxBox*       a_JgupPtr,       // Just a dummy
                                                              int                  a_dir,
                                                              bool                 a_homogeneous,   // Just a dummy
                                                              Real                 a_time,          // Just a dummy
                                                              const Interval&      a_interval) const
{
    CH_TIME("EllipticDoNothingUnlessHomogNeumBCFluxClass::operator()");

    // Do we have anything to do?
    if (!a_homogeneous || a_domain.isPeriodic(a_dir)) return;

    // For the *BCFluxClass stuff, we require the state box to be FC and
    // the valid box to be CC.
    CH_assert(a_state.box().type() == BASISV(a_dir));
    CH_assert(a_valid.type() == IntVect::Zero);

    // We can think of this function as setting Dirichlet BCs on a gradient.
    // So, we can adjust the centering of the valid box and use the Diri code.
    Box FCValid = a_valid;
    FCValid.surroundingNodes(a_dir);

    // Alias the components we want to alter
    const Interval interv = (a_interval == Interval())? a_state.interval(): a_interval;
    CH_assert(a_state.nComp() > interv.end());
    FArrayBox stateAlias(interv, a_state);

    // For this direction, set BCs on each requested side.
    SideIterator sit;
    for (sit.reset(); sit.ok(); ++sit) {
        const Side::LoHiSide iside = sit();
        if (m_doSide[iside][a_dir] == 0) continue;

        // Set the BCs
        setSideDiriBC(stateAlias,
                      FCValid,
                      a_domain,
                      0.0,  // The bcval
                      a_dir,
                      iside,
                      a_homogeneous,
                      1); // order (this is just a dummy value)
    }
}


// -----------------------------------------------------------------------------
// A simple override of BCGhostClass that fills ghosts using extrapolation BCs.
// -----------------------------------------------------------------------------
void EllipticExtrapBCGhostClass::operator() (FArrayBox&           a_state,
                                             const FArrayBox*     a_extrapPtr,      // Just a dummy
                                             const Box&           a_valid,
                                             const ProblemDomain& a_domain,
                                             const RealVect&      a_dx,             // Just a dummy
                                             const DataIndex&     a_index,          // Just a dummy
                                             const FluxBox*       a_JgupPtr,        // Just a dummy
                                             bool                 a_homogeneous,    // Just a dummy
                                             Real                 a_time,           // Just a dummy
                                             const Interval&      a_interval) const
{
    CH_TIME("EllipticExtrapBCGhostClass::operator()");

    // Alias the state components we want to alter
    const Interval interv = (a_interval == Interval())? a_state.interval(): a_interval;
    const int scomp = interv.begin();
    const int ncomp = interv.size();
    CH_assert(a_state.nComp() > interv.end());
    FArrayBox stateAlias(interv, a_state);

    for (int idir = 0; idir < SpaceDim; ++idir) {
        if (a_domain.isPeriodic(idir)) continue;

        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide iside = sit();
            if (m_doSide[iside][idir] == 0) continue;

            setSideExtrapBC(stateAlias,
                            a_valid,
                            a_domain,
                            idir,
                            iside,
                            m_order);
        }
    }
}


// -----------------------------------------------------------------------------
// A simple override of BCFluxClass that fills ghosts using extrapolation BCs.
// a_state is FC.
// a_valid is CC.
// -----------------------------------------------------------------------------
void EllipticExtrapBCFluxClass::operator() (FArrayBox&           a_state,
                                            const FArrayBox*     a_extrapPtr,  // Can be CC!
                                            const Box&           a_valid,
                                            const ProblemDomain& a_domain,
                                            const RealVect&      a_dx,         // Just a dummy
                                            const DataIndex&     a_index,      // Just a dummy
                                            const FluxBox*       a_Jgup,       // Just a dummy
                                            int                  a_dir,
                                            bool                 a_homogeneous,
                                            Real                 a_time,
                                            const Interval&      a_interval) const
{
    CH_TIME("EllipticExtrapBCFluxClass::operator()");

    // For the *BCFluxClass stuff, we require the state box to be FC and
    // the valid box to be CC.
    CH_assert(a_state.box().type() == BASISV(a_dir));
    CH_assert(a_valid.type() == IntVect::Zero);

    // Alias the components we want to alter
    const Interval interv = (a_interval == Interval())? a_state.interval(): a_interval;
    CH_assert(a_state.nComp() > interv.end());
    FArrayBox stateAlias(interv, a_state);

    if (m_useExtrapPtrAsSource) {
        // Use extrapPtr...

        // If this dir is periodic, just leave.
        if (a_domain.isPeriodic(a_dir)) return;

        CH_assert(a_extrapPtr != NULL);
        const FArrayBox& srcFAB = *a_extrapPtr;

        // Make sure src is CC
        CH_assert(srcFAB.box().type() == IntVect::Zero);

        // For this direction, set BCs on each requested side.
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide iside = sit();
            if (m_doSide[iside][a_dir] == 0) continue;

            // Set the BCs
            setSideExtrapBC(stateAlias,
                            srcFAB,
                            a_valid,
                            a_domain,
                            a_dir,
                            iside,
                            m_order);
        }

    } else {
        // Don't use extrapPtr...

        // Adjust the centering of the valid box to conform to extrap function.
        Box thisValid = a_valid;
        thisValid.surroundingNodes(a_dir);

        // Set BCs in each non-periodic direction.
        for (int idir = 0; idir < CH_SPACEDIM; ++idir) {
            if (a_domain.isPeriodic(idir)) continue;

            // For this direction, set BCs on each requested side.
            SideIterator sit;
            for (sit.reset(); sit.ok(); ++sit) {
                const Side::LoHiSide iside = sit();
                if (m_doSide[iside][idir] == 0) continue;

                // Set the BCs
                setSideExtrapBC(stateAlias,
                                thisValid,
                                a_domain,
                                idir,
                                iside,
                                m_order);
            }
        }
    }
}


// -----------------------------------------------------------------------------
// A simple override of BCGhostClass that fills ghosts using zero Diri in normal
// m_normDir and extrapolation BCs in tangential dirs.
// -----------------------------------------------------------------------------
void EllipticInviscidBCGhostClass::operator() (FArrayBox&           a_state,
                                               const FArrayBox*     a_extrapPtr,   // Just a dummy
                                               const Box&           a_valid,
                                               const ProblemDomain& a_domain,
                                               const RealVect&      a_dx,          // Just a dummy
                                               const DataIndex&     a_index,       // Just a dummy
                                               const FluxBox*       a_JgupPtr,     // Just a dummy
                                               bool                 a_homogeneous,
                                               Real                 a_time,        // Just a dummy
                                               const Interval&      a_interval) const
{
    CH_TIME("EllipticInviscidBCGhostClass::operator()");

    // This is the only thing the setSide* functions care about.
    CH_assert(a_state.box().type() == a_valid.type());

    // Alias the components we want to alter
    const Interval interv = (a_interval == Interval())? a_state.interval(): a_interval;
    CH_assert(a_state.nComp() > interv.end());
    FArrayBox stateAlias(interv, a_state);

    // Set BCs in each non-periodic direction.
    for (int idir = 0; idir < CH_SPACEDIM; ++idir) {
        if (a_domain.isPeriodic(idir)) continue;

        if (idir == m_normDir) {
            // No flow through the solid wall
            SideIterator sit;
            for (sit.reset(); sit.ok(); ++sit) {
                const Side::LoHiSide iside = sit();
                if (m_doSide[iside][idir] == 0) continue;

                setSideDiriBC(stateAlias,
                              a_valid,
                              a_domain,
                              0.0, // bcval
                              idir,
                              iside,
                              a_homogeneous,
                              1); // order

            }
        } else {
            // Free-slip
            SideIterator sit;
            for (sit.reset(); sit.ok(); ++sit) {
                const Side::LoHiSide iside = sit();
                if (m_doSide[iside][idir] == 0) continue;

                setSideExtrapBC(stateAlias,
                                a_valid,
                                a_domain,
                                idir,
                                iside,
                                1); // extrap order
            }
        } // end if this is the normal or tangential dir
    } // end loop over directions (idir)
}


// -----------------------------------------------------------------------------
// A simple override of BCGhostClass that fills ghosts using zero Diri in normal
// m_normDir and extrapolation BCs in tangential dirs.
// -----------------------------------------------------------------------------
void BasicVelocityBCGhostClass::operator() (FArrayBox&           a_state,
                                            const FArrayBox*     a_extrapPtr,   // Just a dummy
                                            const Box&           a_valid,
                                            const ProblemDomain& a_domain,
                                            const RealVect&      a_dx,          // Just a dummy
                                            const DataIndex&     a_index,       // Just a dummy
                                            const FluxBox*       a_JgupPtr,     // Just a dummy
                                            bool                 a_homogeneous,
                                            Real                 a_time,
                                            const Interval&      a_interval) const
{
    CH_TIME("BasicVelocityBCGhostClass::operator()");

    // Get the inflow velocity
    Real           inflowVel   = m_inflowVel;
    Side::LoHiSide inflowSide  = m_inflowSide;
    Side::LoHiSide outflowSide = m_outflowSide;

    if (m_inflowVelFuncPtr != NULL) {
        inflowVel = m_inflowVelFuncPtr(a_time) * m_inflowVelScale;
        if (0.0 < inflowVel) {
            inflowSide = Side::Lo;
            outflowSide = Side::Hi;
        } else {
            inflowSide = Side::Hi;
            outflowSide = Side::Lo;
        }
    }

    // This is the only thing the setSide* functions care about.
    CH_assert(a_state.box().type() == a_valid.type());

    // Alias the components we want to alter
    const Interval interv = (a_interval == Interval())? a_state.interval(): a_interval;
    CH_assert(a_state.nComp() > interv.end());
    FArrayBox stateAlias(interv, a_state);

    // Set BCs in each non-periodic direction.
    for (int idir = 0; idir < CH_SPACEDIM; ++idir) {
        if (a_domain.isPeriodic(idir)) continue;

        // For this direction, set BCs on each requested side.
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide iside = sit();
            if (m_doSide[iside][idir] == 0) continue;

            if (idir == m_inflowDir && iside == inflowSide) {
                // Inflow
                if (m_velComp == m_inflowDir) {
                    // inflow val Diri
                    setSideDiriBC(stateAlias,
                                  a_valid,
                                  a_domain,
                                  inflowVel,
                                  idir,
                                  iside,
                                  a_homogeneous,
                                  1); // order
                } else {
                    // zero Diri
                    setSideDiriBC(stateAlias,
                                  a_valid,
                                  a_domain,
                                  0.0,
                                  idir,
                                  iside,
                                  a_homogeneous,
                                  1); // order
                }
            } else if (idir == m_outflowDir && iside == outflowSide) {
                // Outflow
                // TODO: Should be zero Neumann
                setSideExtrapBC(stateAlias,
                                a_valid,
                                a_domain,
                                idir,
                                iside,
                                0); // extrap order
            } else {
                // Solid wall
                if (m_isViscous) {
                    if (m_velComp == m_inflowDir) {
                        // inflow val Diri
                        setSideDiriBC(stateAlias,
                                      a_valid,
                                      a_domain,
                                      inflowVel,
                                      idir,
                                      iside,
                                      a_homogeneous,
                                      1); // order
                    } else {
                        // zero Diri
                        setSideDiriBC(stateAlias,
                                      a_valid,
                                      a_domain,
                                      0.0,
                                      idir,
                                      iside,
                                      a_homogeneous,
                                      1); // order
                    }
                } else {
                    // inviscid
                    if (idir == m_velComp) {
                        // No flow through wall - zero Diri
                        setSideDiriBC(stateAlias,
                                      a_valid,
                                      a_domain,
                                      0.0,
                                      idir,
                                      iside,
                                      a_homogeneous,
                                      1); // order
                    } else {
                        setSideExtrapBC(stateAlias,
                                        a_valid,
                                        a_domain,
                                        idir,
                                        iside,
                                        1); // extrap order
                    }
                }
            } // end if inflow, outflow, or solid wall
        } // end loop over sides (sit)
    } // end loop over directions (idir)
}

