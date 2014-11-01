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
#include "AdvectUtil.H"
#include "AdvectUtilF_F.H"
#include "LevelGeometry.H"
#include "Constants.H"
#include "MappedGodunovUtilitiesF_F.H"
#include "PeriodicLoHiCenter.H"
#include "Debug.H"

// #define nanCheck(x) checkForValidNAN(x)
#define nanCheck(x)


// -----------------------------------------------------------------------------
// Default constructor - leaves object in a useless state.
// -----------------------------------------------------------------------------
MappedAdvectionUtil::MappedAdvectionUtil ()
: m_levGeoPtr(NULL)
{;}


// -----------------------------------------------------------------------------
// Full constructor - just calls the define function
// -----------------------------------------------------------------------------
MappedAdvectionUtil::MappedAdvectionUtil (const LevelGeometry* a_levGeoPtr,
                                          const int            a_normalPredOrder,
                                          const bool           a_useFourthOrderSlopes,
                                          const bool           a_useLimiting,
                                          const bool           a_useHighOrderLimiter,
                                          const bool           a_useUpwinding)
: m_levGeoPtr(NULL)
{

    this->define(a_levGeoPtr,
                 a_normalPredOrder,
                 a_useFourthOrderSlopes,
                 a_useLimiting,
                 a_useHighOrderLimiter,
                 a_useUpwinding);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
MappedAdvectionUtil::~MappedAdvectionUtil ()
{
    m_levGeoPtr = NULL;
}


// -----------------------------------------------------------------------------
// define -- Essentially, the full constructor
//
// a_normalPredOrder is CTU_NORMAL_PRED, PLM_NORMAL_PRED, or PPM_NORMAL_PRED.
// If a_useLimiting is true, then a traditional van Leer slope limiter is used.
// If a_useHighOrderLimiter is true, then an extremum-preserving van Leer slope
//  limiter and special version of the PPM face interpolator is used.
// -----------------------------------------------------------------------------
void MappedAdvectionUtil::define (const LevelGeometry* a_levGeoPtr,
                                  const int            a_normalPredOrder,
                                  const bool           a_useFourthOrderSlopes,
                                  const bool           a_useLimiting,
                                  const bool           a_useHighOrderLimiter,
                                  const bool           a_useUpwinding)
{

    CH_assert(a_levGeoPtr != NULL);
    m_levGeoPtr = a_levGeoPtr;

    CH_assert(0 <= a_normalPredOrder);
    CH_assert(a_normalPredOrder < _NUM_NORMAL_PRED);
    m_normalPredOrder = a_normalPredOrder;

    m_useFourthOrderSlopes = a_useFourthOrderSlopes;
    m_useLimiting = a_useLimiting;
    m_useHighOrderLimiter = a_useHighOrderLimiter;
    m_useUpwinding = a_useUpwinding;

    const ProblemDomain& domain = a_levGeoPtr->getDomain();
    const RealVect& dx = a_levGeoPtr->getDx();

    // GodunovUtilities asked for a scalar dx, but never used it. Lucky me.
    m_util.define(a_levGeoPtr);
    m_util.highOrderLimiter(a_useHighOrderLimiter);
}


#define USE_OLD_METHOD
#ifdef USE_OLD_METHOD
// -----------------------------------------------------------------------------
// predictScalar

// Performs the characteristic tracing of a_Wold and predicts the time-centered,
// FC a_Whalf. This function does not perform exchanges. If a_returnFlux is
// false, the user must set final BCs.

// a_Whalf is FC and time-centered.
// a_Wold is CC.
// a_sourceTermPtr is CC and can be NULL.
// a_oldVel is CC.
// a_advVel is FC.
// a_oldTime and a_fluxBCPtr only need to be given if a_returnFlux is true.
// -----------------------------------------------------------------------------
void MappedAdvectionUtil::predictScalar (LevelData<FluxBox>&             a_Whalf,
                                         const LevelData<FArrayBox>&     a_Wold,
                                         const LevelData<FArrayBox>*     a_sourceTermPtr,
                                         const LevelData<FArrayBox>&     a_oldVel,
                                         const LevelData<FluxBox>&       a_advVel,
                                         const Real                      a_dt,
                                         const LevelGeometry&            a_levGeo,
                                         Tuple<BCMethodHolder,SpaceDim>  a_BCValues,
                                         Tuple<BCMethodHolder,SpaceDim>  a_BCSlopes,
                                         const Real                      a_oldTime,
                                         const bool                      a_returnFlux)
{
    CH_TIME("MappedAdvectionUtil::predictScalar");

    // Sanity checks
    CH_assert(m_levGeoPtr != NULL);

    CH_assert(a_Whalf .nComp() == 1);
    CH_assert(a_Wold  .nComp() == 1);
    CH_assert(a_oldVel.nComp() == SpaceDim);
    CH_assert(a_advVel.nComp() == 1);
    if (a_sourceTermPtr != NULL) {
        CH_assert(a_sourceTermPtr->nComp() == 1);
    }

    // Gather some basic data
    const ProblemDomain& domain = m_levGeoPtr->getDomain();
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    DataIterator dit = grids.dataIterator();

    const RealVect& dx = m_levGeoPtr->getDx();
    const RealVect dtondx = a_dt / dx;
    const Real halfTime = a_oldTime + 0.5 * a_dt;
    const IntVect tracingGhosts(D_DECL(ADVECT_GROW, ADVECT_GROW, ADVECT_GROW));
    CH_assert(a_Wold.ghostVect() == tracingGhosts);

    // Domain checks
    CH_assert(a_Whalf.getBoxes().compatible(grids));
    CH_assert(a_Wold .getBoxes().compatible(grids));
    CH_assert(a_oldVel .getBoxes().compatible(grids));
    CH_assert(a_advVel .getBoxes().compatible(grids));
    if (a_sourceTermPtr != NULL) {
        CH_assert(a_sourceTermPtr->getBoxes().compatible(grids));
    }

    // Perform the calculation over one grid at a time.
    for (dit.reset(); dit.ok(); ++dit) {
        // Create references for convenience.
        FluxBox& WhalfFB = a_Whalf[dit];
        const FArrayBox& WoldFAB = a_Wold[dit];
        const FArrayBox& oldVelFAB = a_oldVel[dit];
        const FluxBox& advVelFB = a_advVel[dit];
        const FluxBox* JgupPtr = &(m_levGeoPtr->getFCJgup()[dit]);

        // The current box of valid cells
        Box curBox = grids[dit];

        // Boxes for face-centered state - used for the riemann() and
        // artificialViscosity() calls
        Box faceBox[SpaceDim];

        // Boxes for face-centered fluxes
        Box fluxBox[SpaceDim];

        // Boxes for cell-centered state - used for the updatePrim() calls
        Box ccBox[SpaceDim];

        for (int dir1 = 0; dir1 < SpaceDim; ++dir1) {
            // Usually, we employ the Riemann solver at valid faces that are not
            // at a physical boundary, then apply the BCs. However, if the BCs
            // are periodic, then the ghost cells have valid data and can be
            // used by the Riemann solver at the physical boundary.
            Box domGrow = domain.domainBox();
            if (domain.isPeriodic(dir1)) domGrow.grow(dir1, 1);

            // faceBox[dir1] is face-centered in direction "dir1",
            // is valid "curBox" grown by 1 in all directions except "dir1",
            // and stays one cell away, in "dir1", from the domain boundary.
            faceBox[dir1] = curBox; // valid cell-centered box
            faceBox[dir1].grow(1);
            faceBox[dir1] &= domGrow; //domain;
            faceBox[dir1].grow(dir1, -1);
            faceBox[dir1].surroundingNodes(dir1);

            // ccBox[dir1] is cell-centered,
            // is valid "curBox" grown by 1 in all directions except "dir1",
            // but only within the domain, "m_domain".
            ccBox[dir1] = curBox;
            ccBox[dir1].grow(1);
            ccBox[dir1].grow(dir1, -1);
            ccBox[dir1] &= domain;

            // fluxBox[dir1] is face-centered in direction "dir1",
            // consisting of all those faces of cells of "ccBox[dir1]".
            fluxBox[dir1] = ccBox[dir1];
            fluxBox[dir1].surroundingNodes(dir1);

            // The difference between faceBox[dir1] and fluxBox[dir1] is:
            // if curBox abuts the boundary of the domain in direction "dir1",
            // then fluxBox[dir1] contains faces along that boundary,
            // but faceBox[dir1] does not.
        }

        // cell-centered box:  but restrict it to within domain
        Box WBox = WoldFAB.box();
        WBox &= domain;

        // slopeBox is the cell-centered box where slopes will be needed:
        // it is one larger than the final update box "curBox" of valid cells,
        // but within domain.
        // On slopeBox we define the FABs flattening, WMinus, and WPlus.
        Box slopeBox = curBox;
        slopeBox.grow(1);
        slopeBox &= domain;

        // Intermediate, extrapolated primitive variables
        FArrayBox WMinus[SpaceDim];
        FArrayBox WPlus [SpaceDim];

        // Initial fluxes
        FArrayBox WHalf1[SpaceDim];

        // Compute initial fluxes
        for (int dir1 = 0; dir1 < SpaceDim; ++dir1) {
            // Size the intermediate, extrapolated primitive variables
            WMinus[dir1].resize(slopeBox, 1); // cell-centered
            WPlus [dir1].resize(slopeBox, 1); // cell-centered

            // We can just use a_Whalf's memory since we won't use it until
            // we no longer need WHalf1.
            WHalf1[dir1].define(Interval(0,0), WhalfFB[dir1]);

            // Compute predictor step to obtain extrapolated primitive variables
            switch (m_normalPredOrder) {
                case CTU_NORMAL_PRED:
                    CTUNormalPred(WMinus[dir1],
                                  WPlus[dir1],
                                  a_dt,
                                  WoldFAB,
                                  oldVelFAB,
                                  dir1,
                                  slopeBox,
                                  dit());
                    break;
                case PLM_NORMAL_PRED:
                    PLMNormalPred(WMinus[dir1],
                                  WPlus[dir1],
                                  a_dt,
                                  WoldFAB,
                                  oldVelFAB,
                                  dir1,
                                  slopeBox,
                                  dit());
                    break;
                case PPM_NORMAL_PRED:
                    PPMNormalPred(WMinus[dir1],
                                  WPlus[dir1],
                                  a_dt,
                                  WoldFAB,
                                  oldVelFAB,
                                  dir1,
                                  slopeBox,
                                  dit());
                    break;
                default:
                    MayDay::Error("MappedAdvectionUtil::traceScalar: Normal predictor order must be 0 (CTU), 1 (PLM), or 2 (PPM)");
            }

            // If the source term is valid add it to the primitive quantities
            if (a_sourceTermPtr != NULL) {
                const FArrayBox& srcFAB = (*a_sourceTermPtr)[dit];
                const Real scale = 0.5 * a_dt;

                WMinus[dir1].plus(srcFAB, scale);
                WPlus [dir1].plus(srcFAB, scale);
            }

            // Solve the Riemann problem
            WHalf1[dir1].resize(fluxBox[dir1],1);
            RiemannSolver(WHalf1[dir1],
                          WPlus[dir1],
                          WMinus[dir1],
                          advVelFB,
                          dir1,
                          faceBox[dir1]);

            if (!domain.isPeriodic(dir1)) {
                Box valid = surroundingNodes(domain.domainBox(), dir1);
                valid.grow(ADVECT_GROW);
                valid.grow(dir1, -ADVECT_GROW);
                CH_assert(WHalf1[dir1].box().type() == valid.type());
                valid &= WHalf1[dir1].box();
                CH_assert(!valid.isEmpty());

                a_BCValues[dir1].setGhosts(WHalf1[dir1],
                                           NULL,
                                           valid,
                                           domain,
                                           dx,
                                           dit(),
                                           JgupPtr,
                                           false,
                                           halfTime);
            }

            // WhalfFB[dir1].copy(WHalf1[dir1]);
        } // end loop over FC dir (dir1)

#if (CH_SPACEDIM == 3)
        // In 3D, compute some additional intermediate fluxes
        //
        // NOTE:  The diagonal entries of this array of fluxes are not
        // used and will not be defined.
        FArrayBox WHalf2[SpaceDim][SpaceDim];

        // Compute the intermediate, corrected fluxes in each direction
        for (int dir1 = 0; dir1 < SpaceDim; dir1++) {
            // Correct fluxes using fluxes from a different direction
            for (int dir2 = 0; dir2 < SpaceDim; dir2++) {
                if (dir2 == dir1) continue;

                // Temporary primitive variables
                FArrayBox WTempMinus(WMinus[dir1].box(),1);
                FArrayBox WTempPlus (WPlus [dir1].box(),1);
                FArrayBox AdWdx(WPlus[dir1].box(),1);

                // preventing uninitialized memory reads which cause
                // FPE's on some machines
                // AdWdx.setVal(666.666);

                // Copy data for in place modification
                WTempMinus.copy(WMinus[dir1]);
                WTempPlus .copy(WPlus [dir1]);

                // Update the current, extrapolated primitive variable using a flux
                // in a different direction
                quasilinearUpdate(AdWdx,
                                  WHalf1[dir2],
                                  oldVelFAB,
                                  -(1.0/3.0) * dtondx[dir2],
                                  dir2,
                                  ccBox[dir2]);

                WTempMinus += AdWdx;
                WTempPlus  += AdWdx;

                // Solve the Riemann problem.
                WHalf2[dir1][dir2].resize(fluxBox[dir1],1);
                RiemannSolver(WHalf2[dir1][dir2],
                              WTempPlus,
                              WTempMinus,
                              advVelFB,
                              dir1,
                              faceBox[dir1]);

                if (!domain.isPeriodic(dir1)) {
                    Box valid = surroundingNodes(domain.domainBox(), dir1);
                    valid.grow(ADVECT_GROW);
                    valid.grow(dir1, -ADVECT_GROW);
                    CH_assert(WHalf2[dir1][dir2].box().type() == valid.type());
                    valid &= WHalf2[dir1][dir2].box();
                    CH_assert(!valid.isEmpty());

                    a_BCValues[dir1].setGhosts(WHalf2[dir1][dir2],  // stateFAB
                                               NULL,                // extrapFABPtr
                                               valid,               // valid
                                               domain,              // domain
                                               dx,                  // dx
                                               dit(),               // DataIndex
                                               JgupPtr,             // JgupFBPtr
                                               false,               // isHomogeneous
                                               halfTime);           // time
                }
            } // end loop over transverse dir (dir2)
        } // end loop over FC dir (dir1)
#endif

        // faceBox and fluxBox are now a bit smaller for the final corrections
        for (int dir1 = 0; dir1 < SpaceDim; ++dir1) {
            Box domGrow = domain.domainBox();
            if (domain.isPeriodic(dir1)) domGrow.grow(dir1, 1);

            faceBox[dir1] = curBox;
            faceBox[dir1].grow(dir1,1);
            faceBox[dir1] &= domGrow;
            faceBox[dir1].grow(dir1,-1);
            faceBox[dir1].surroundingNodes(dir1);

            fluxBox[dir1] = curBox;
            fluxBox[dir1].surroundingNodes(dir1);
        }

        // Do the final corrections to the fluxes
        for (int dir1 = 0; dir1 < SpaceDim; ++dir1) {
            // Correct the flux using fluxes in the remaining direction(s)
            for (int dir2 = 0; dir2 < SpaceDim; ++dir2) {
                if (dir2 == dir1) continue;

#if (CH_SPACEDIM == 2)
                // In 2D, the current primitive state is updated by a flux in
                // the other direction
                FArrayBox AdWdx(WPlus[dir1].box(), 1);
                quasilinearUpdate(AdWdx,
                                  WHalf1[dir2],
                                  oldVelFAB,
                                  -(1.0/2.0) * dtondx[dir2],
                                  dir2,
                                  ccBox[dir2]);

                WMinus[dir1] += AdWdx;
                WPlus[dir1] += AdWdx;

#elif (CH_SPACEDIM == 3)
                // In 3D, find a direction different from the two above
                int dir3 = 3 - dir1 - dir2;

                // Update the conservative state using both corrected fluxes in
                // the other two directions
                FArrayBox AdWdx(WPlus[dir1].box(),1);

                quasilinearUpdate(AdWdx,
                                  WHalf2[dir2][dir3],
                                  oldVelFAB,
                                  -(1.0/2.0) * dtondx[dir2],
                                  dir2,
                                  ccBox[dir2]);

                WMinus[dir1].plus(AdWdx,0,0,1);
                WPlus[dir1].plus(AdWdx,0,0,1);
#else
                // Only 2D and 3D should be possible
                MayDay::Error("traceScalar: CH_SPACEDIM not 2 or 3!");
#endif
            } // end loop over transverse dir (dir2)

            // Solve the Riemann problem to obtain time-centered face values.
            RiemannSolver(WhalfFB[dir1],
                          WPlus[dir1],
                          WMinus[dir1],
                          advVelFB,
                          dir1,
                          faceBox[dir1]);

            if (!domain.isPeriodic(dir1)) {
                Box valid = surroundingNodes(domain.domainBox(), dir1);
                // valid.grow(ADVECT_GROW);
                // valid.grow(dir1, -ADVECT_GROW);
                CH_assert(WhalfFB[dir1].box().type() == valid.type());
                valid &= WhalfFB[dir1].box();
                CH_assert(!valid.isEmpty());

                a_BCValues[dir1].setGhosts(WhalfFB[dir1],   // stateFAB
                                           NULL,            // extrapFABPtr
                                           valid,           // valid
                                           domain,          // domain
                                           dx,              // dx
                                           dit(),           // DataIndex
                                           JgupPtr,         // JgupFBPtr
                                           false,           // isHomogeneous
                                           halfTime);       // time
            }

            // Set boundary fluxes if requested
            if (a_returnFlux) {
                // NOTE: If this is ever needed, you'll have to decide if it's really
                // FCgup that you want due to the mult-by-advVel move.
                // const FluxBox& JgupFB = a_levGeo.getFCJgup()[dit()];

                if (!domain.isPeriodic(dir1)) {
                    Box valid = WhalfFB[dir1].box();
                    valid.enclosedCells();
                    valid &= domain.domainBox();

                    a_BCValues[dir1].setFluxes(WhalfFB[dir1],   // stateFAB
                                               &WoldFAB,        // extrapFABPtr
                                               valid,           // valid
                                               domain,          // domain
                                               dx,              // dx
                                               dit(),           // DataIndex
                                               NULL,            // JgupFBPtr
                                               dir1,            // dir
                                               false,           // isHomogeneous
                                               halfTime);       // time
                } // end if domain is not periodic in dir1

                // Multiply by Uadv (J * advecting velocity)
                // NOTE: This used to happen before the setFluxes call.
                WhalfFB[dir1].mult(advVelFB[dir1], 0, 0, 1);

            } // end if returning fluxes
        } // end loop over FC dir (dir1)
    } // end loop over grids (dit)
}

#else

void MappedAdvectionUtil::predictScalar (LevelData<FluxBox>&             a_Whalf,
                                         const LevelData<FArrayBox>&     a_Wold,
                                         const LevelData<FArrayBox>*     a_sourceTermPtr,
                                         const LevelData<FArrayBox>&     a_oldVel,
                                         const LevelData<FluxBox>&       a_advVel,
                                         const Real                      a_dt,
                                         const LevelGeometry&            a_levGeo,
                                         Tuple<BCMethodHolder,SpaceDim>  a_BCValues,
                                         Tuple<BCMethodHolder,SpaceDim>  a_BCSlopes,
                                         const Real                      a_oldTime,
                                         const bool                      a_returnFlux)
{
    CH_TIME("MappedAdvectionUtil::predictScalar");

    // Sanity checks
    CH_assert(m_levGeoPtr != NULL);

    CH_assert(a_Whalf .nComp() == 1);
    CH_assert(a_Wold  .nComp() == 1);
    CH_assert(a_oldVel.nComp() == SpaceDim);
    CH_assert(a_advVel.nComp() == 1);
    if (a_sourceTermPtr != NULL) {
        CH_assert(a_sourceTermPtr->nComp() == 1);
    }

    // Gather some basic data
    const ProblemDomain& domain = m_levGeoPtr->getDomain();
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    DataIterator dit = grids.dataIterator();

    const RealVect& dx = m_levGeoPtr->getDx();
    const RealVect dtondx = a_dt / dx;
    const Real halfTime = a_oldTime + 0.5 * a_dt;
    const IntVect tracingGhosts(D_DECL(ADVECT_GROW, ADVECT_GROW, ADVECT_GROW));
    CH_assert(a_Wold.ghostVect() == tracingGhosts);

    // Domain checks
    CH_assert(a_Whalf.getBoxes().compatible(grids));
    CH_assert(a_Wold .getBoxes().compatible(grids));
    CH_assert(a_oldVel .getBoxes().compatible(grids));
    CH_assert(a_advVel .getBoxes().compatible(grids));
    if (a_sourceTermPtr != NULL) {
        CH_assert(a_sourceTermPtr->getBoxes().compatible(grids));
    }

    // Make sure we were given usable data.
    nanCheck(a_Wold);
    nanCheck(a_oldVel);
    nanCheck(a_advVel);

    // This is a domain box that includes ghosts in periodic directions.
    Box perDomBox = domain.domainBox();
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (domain.isPeriodic(dir)) {
            perDomBox.grow(dir, ADVECT_GROW);
        }
    }

    // Perform the calculation over one grid at a time.
    for (dit.reset(); dit.ok(); ++dit) {
        // Create references for convenience.
        FluxBox& WhalfFB = a_Whalf[dit];
        const FArrayBox& WoldFAB = a_Wold[dit];
        const FArrayBox& oldVelFAB = a_oldVel[dit];
        const FluxBox& advVelFB = a_advVel[dit];
        const FluxBox* JgupPtr = &(m_levGeoPtr->getFCJgup()[dit]);

        // The current box of valid cells
        Box curBox = grids[dit];

        // Boxes for face-centered state - used for the riemann() and
        // artificialViscosity() calls
        Box faceBox[SpaceDim];

        // Boxes for face-centered fluxes
        Box fluxBox[SpaceDim];

        // Boxes for cell-centered state - used for the updatePrim() calls
        Box ccBox[SpaceDim];

        for (int dir1 = 0; dir1 < SpaceDim; ++dir1) {

            // faceBox[dir1] is face-centered in direction "dir1",
            // is valid "curBox" grown by 1 in all directions except "dir1",
            // and stays one cell away, in "dir1", from the domain boundary.
            faceBox[dir1] = curBox; // valid cell-centered box
            faceBox[dir1].grow(1);
            faceBox[dir1] &= perDomBox;
            faceBox[dir1].grow(dir1, -1);
            faceBox[dir1].surroundingNodes(dir1);

            // ccBox[dir1] is cell-centered,
            // is valid "curBox" grown by 1 in all directions except "dir1",
            // but only within the domain, "m_domain".
            ccBox[dir1] = curBox;
            ccBox[dir1].grow(1);
            ccBox[dir1].grow(dir1, -1);
            ccBox[dir1] &= perDomBox;

            // fluxBox[dir1] is face-centered in direction "dir1",
            // consisting of all those faces of cells of "ccBox[dir1]".
            fluxBox[dir1] = ccBox[dir1];
            fluxBox[dir1].surroundingNodes(dir1);

            // The difference between faceBox[dir1] and fluxBox[dir1] is:
            // if curBox abuts the boundary of the domain in direction "dir1",
            // then fluxBox[dir1] contains faces along that boundary,
            // but faceBox[dir1] does not.
        }

        // cell-centered box:  but restrict it to within domain
        Box WBox = WoldFAB.box();
        WBox &= perDomBox;

        // slopeBox is the cell-centered box where slopes will be needed:
        // it is one larger than the final update box "curBox" of valid cells,
        // but within domain.
        // On slopeBox we define the FABs flattening, WMinus, and WPlus.
        Box slopeBox = curBox;
        slopeBox.grow(1);
        slopeBox &= perDomBox;

        // Intermediate, extrapolated primitive variables
        FArrayBox WMinus[SpaceDim];
        FArrayBox WPlus [SpaceDim];

        // Initial fluxes
        FArrayBox WHalf1[SpaceDim];

        // Compute initial fluxes
        for (int dir1 = 0; dir1 < SpaceDim; ++dir1) {
            // Size the intermediate, extrapolated primitive variables
            WMinus[dir1].resize(slopeBox, 1); // cell-centered
            WPlus [dir1].resize(slopeBox, 1); // cell-centered

            // We can just use a_Whalf's memory since we won't use it until
            // we no longer need WHalf1.
            WHalf1[dir1].define(Interval(0,0), WhalfFB[dir1]);

            // Compute predictor step to obtain extrapolated primitive variables
            switch (m_normalPredOrder) {
                case CTU_NORMAL_PRED:
                    CTUNormalPred(WMinus[dir1],
                                  WPlus[dir1],
                                  a_dt,
                                  WoldFAB,
                                  oldVelFAB,
                                  dir1,
                                  slopeBox,
                                  dit());
                    break;
                case PLM_NORMAL_PRED:
                    PLMNormalPred(WMinus[dir1],
                                  WPlus[dir1],
                                  a_dt,
                                  WoldFAB,
                                  oldVelFAB,
                                  dir1,
                                  slopeBox,
                                  dit());
                    break;
                case PPM_NORMAL_PRED:
                    PPMNormalPred(WMinus[dir1],
                                  WPlus[dir1],
                                  a_dt,
                                  WoldFAB,
                                  oldVelFAB,
                                  dir1,
                                  slopeBox,
                                  dit());
                    break;
                default:
                    MayDay::Error("MappedAdvectionUtil::traceScalar: Normal predictor order must be 0 (CTU), 1 (PLM), or 2 (PPM)");
            }

            // If the source term is valid add it to the primitive quantities
            if (a_sourceTermPtr != NULL) {
                const FArrayBox& srcFAB = (*a_sourceTermPtr)[dit];
                const Real scale = 0.5 * a_dt;

                WMinus[dir1].plus(srcFAB, scale);
                WPlus [dir1].plus(srcFAB, scale);
            }

            // Solve the Riemann problem
            WHalf1[dir1].resize(fluxBox[dir1],1);
            RiemannSolver(WHalf1[dir1],
                          WPlus[dir1],
                          WMinus[dir1],
                          advVelFB,
                          dir1,
                          faceBox[dir1]);

            if (!domain.isPeriodic(dir1)) {
                TODO();
                Box valid = surroundingNodes(domain.domainBox(), dir1);
                valid.grow(ADVECT_GROW);
                valid.grow(dir1, -ADVECT_GROW);
                CH_assert(WHalf1[dir1].box().type() == valid.type());
                valid &= WHalf1[dir1].box();
                CH_assert(!valid.isEmpty());

                // Box valid = faceBox[dir1];
                // valid &= surroundingNodes(perDomBox, dir1);
                // valid &= WHalf1[dir1].box();
                // CH_assert(!valid.isEmpty());

                a_BCValues[dir1].setGhosts(WHalf1[dir1],
                                           NULL,
                                           valid,
                                           domain,
                                           dx,
                                           dit(),
                                           JgupPtr,
                                           false,
                                           halfTime);
            }
        } // end loop over FC dir (dir1)

#if (CH_SPACEDIM == 3)
        // In 3D, compute some additional intermediate fluxes
        //
        // NOTE:  The diagonal entries of this array of fluxes are not
        // used and will not be defined.
        FArrayBox WHalf2[SpaceDim][SpaceDim];

        // Compute the intermediate, corrected fluxes in each direction
        for (int dir1 = 0; dir1 < SpaceDim; dir1++) {
            // Correct fluxes using fluxes from a different direction
            for (int dir2 = 0; dir2 < SpaceDim; dir2++) {
                if (dir2 == dir1) continue;

                // Temporary primitive variables
                FArrayBox WTempMinus(WMinus[dir1].box(),1);
                FArrayBox WTempPlus (WPlus [dir1].box(),1);
                FArrayBox AdWdx(WPlus[dir1].box(),1);

                // preventing uninitialized memory reads which cause
                // FPE's on some machines
                // AdWdx.setVal(666.666);

                // Copy data for in place modification
                WTempMinus.copy(WMinus[dir1]);
                WTempPlus .copy(WPlus [dir1]);

                // Update the current, extrapolated primitive variable using a flux
                // in a different direction
                quasilinearUpdate(AdWdx,
                                  WHalf1[dir2],
                                  oldVelFAB,
                                  -(1.0/3.0) * dtondx[dir2],
                                  dir2,
                                  ccBox[dir2]);

                WTempMinus += AdWdx;
                WTempPlus  += AdWdx;

                // Solve the Riemann problem.
                WHalf2[dir1][dir2].resize(fluxBox[dir1],1);
                RiemannSolver(WHalf2[dir1][dir2],
                              WTempPlus,
                              WTempMinus,
                              advVelFB,
                              dir1,
                              faceBox[dir1]);

                if (!domain.isPeriodic(dir1)) {
                    TODO();
                    Box valid = surroundingNodes(domain.domainBox(), dir1);
                    valid.grow(ADVECT_GROW);
                    valid.grow(dir1, -ADVECT_GROW);
                    CH_assert(WHalf2[dir1][dir2].box().type() == valid.type());
                    valid &= WHalf2[dir1][dir2].box();
                    CH_assert(!valid.isEmpty());

                    a_BCValues[dir1].setGhosts(WHalf2[dir1][dir2],  // stateFAB
                                               NULL,                // extrapFABPtr
                                               valid,               // valid
                                               domain,              // domain
                                               dx,                  // dx
                                               dit(),               // DataIndex
                                               JgupPtr,             // JgupFBPtr
                                               false,               // isHomogeneous
                                               halfTime);           // time
                }
            } // end loop over transverse dir (dir2)
        } // end loop over FC dir (dir1)
#endif

        // faceBox and fluxBox are now a bit smaller for the final corrections
        for (int dir1 = 0; dir1 < SpaceDim; ++dir1) {
            faceBox[dir1] = curBox;
            faceBox[dir1].grow(dir1,1);
            faceBox[dir1] &= perDomBox;
            faceBox[dir1].grow(dir1,-1);
            faceBox[dir1].surroundingNodes(dir1);

            fluxBox[dir1] = curBox;
            fluxBox[dir1].surroundingNodes(dir1);
        }

        // Do the final corrections to the fluxes
        for (int dir1 = 0; dir1 < SpaceDim; ++dir1) {
            // Correct the flux using fluxes in the remaining direction(s)
            for (int dir2 = 0; dir2 < SpaceDim; ++dir2) {
                if (dir2 == dir1) continue;

#if (CH_SPACEDIM == 2)
                // In 2D, the current primitive state is updated by a flux in
                // the other direction
                CH_assert(WPlus[dir1].box() == WMinus[dir1].box());
                FArrayBox AdWdx(WPlus[dir1].box(), 1);
                quasilinearUpdate(AdWdx,
                                  WHalf1[dir2],
                                  oldVelFAB,
                                  -(1.0/2.0) * dtondx[dir2],
                                  dir2,
                                  ccBox[dir2]);

                WMinus[dir1] += AdWdx;
                WPlus[dir1] += AdWdx;

#elif (CH_SPACEDIM == 3)
                // In 3D, find a direction different from the two above
                int dir3 = 3 - dir1 - dir2;

                // Update the conservative state using both corrected fluxes in
                // the other two directions
                FArrayBox AdWdx(WPlus[dir1].box(),1);

                quasilinearUpdate(AdWdx,
                                  WHalf2[dir2][dir3],
                                  oldVelFAB,
                                  -(1.0/2.0) * dtondx[dir2],
                                  dir2,
                                  ccBox[dir2]);

                WMinus[dir1].plus(AdWdx,0,0,1);
                WPlus[dir1].plus(AdWdx,0,0,1);
#else
                // Only 2D and 3D should be possible
                MayDay::Error("traceScalar: CH_SPACEDIM not 2 or 3!");
#endif
            } // end loop over transverse dir (dir2)

            // Solve the Riemann problem to obtain time-centered face values.
            RiemannSolver(WhalfFB[dir1],
                          WPlus[dir1],
                          WMinus[dir1],
                          advVelFB,
                          dir1,
                          faceBox[dir1]);

            if (!domain.isPeriodic(dir1)) {
                TODO();
                Box valid = surroundingNodes(domain.domainBox(), dir1);
                // valid.grow(ADVECT_GROW);
                // valid.grow(dir1, -ADVECT_GROW);
                CH_assert(WhalfFB[dir1].box().type() == valid.type());
                valid &= WhalfFB[dir1].box();
                CH_assert(!valid.isEmpty());

                // Box valid = faceBox[dir1];
                // valid &= surroundingNodes(perDomBox, dir1);
                // valid &= WhalfFB[dir1].box();
                // CH_assert(!valid.isEmpty());

                a_BCValues[dir1].setGhosts(WhalfFB[dir1],   // stateFAB
                                           NULL,            // extrapFABPtr
                                           valid,           // valid
                                           domain,          // domain
                                           dx,              // dx
                                           dit(),           // DataIndex
                                           JgupPtr,         // JgupFBPtr
                                           false,           // isHomogeneous
                                           halfTime);       // time
            }

            // Set boundary fluxes if requested
            if (a_returnFlux) {
                // NOTE: If this is ever needed, you'll have to decide if it's really
                // FCgup that you want due to the mult-by-advVel move.
                // const FluxBox& JgupFB = a_levGeo.getFCJgup()[dit()];

                if (!domain.isPeriodic(dir1)) {
                    TODO();
                    Box valid = WhalfFB[dir1].box();
                    valid.enclosedCells();
                    valid &= domain.domainBox();

                    a_BCValues[dir1].setFluxes(WhalfFB[dir1],   // stateFAB
                                               &WoldFAB,        // extrapFABPtr
                                               valid,           // valid
                                               domain,          // domain
                                               dx,              // dx
                                               dit(),           // DataIndex
                                               NULL,            // JgupFBPtr
                                               dir1,            // dir
                                               false,           // isHomogeneous
                                               halfTime);       // time
                } // end if domain is not periodic in dir1

                // Multiply by Uadv (J * advecting velocity)
                // NOTE: This used to happen before the setFluxes call.
                WhalfFB[dir1].mult(advVelFB[dir1], 0, 0, 1);

            } // end if returning fluxes
        } // end loop over FC dir (dir1)
    } // end loop over grids (dit)

    // How did we do?
    nanCheck(a_Whalf);
}
#endif


// -----------------------------------------------------------------------------
// RiemannSolver
//
// a_WGdnv is the FC solution to the Riemann problem.
// a_WLeft is the state predicted from the left. It is CC.
// a_WRight is the state predicted from the right. It is CC.
// a_advVel is the advecting velocity. It is FC.
// -----------------------------------------------------------------------------
void MappedAdvectionUtil::RiemannSolver (FArrayBox&       a_WGdnv,
                                         const FArrayBox& a_WLeft,
                                         const FArrayBox& a_WRight,
                                         const FluxBox&   a_advVel,
                                         const int        a_dir,
                                         const Box&       a_box)
{

    // Sanity checks
    CH_assert(a_WGdnv.box().contains(a_box));
    CH_assert(a_WGdnv .nComp() == 1);
    CH_assert(a_WLeft .nComp() == 1);
    CH_assert(a_WRight.nComp() == 1);

    // Cast away "const" inputs so their boxes can be shifted left or right
    // 1/2 cell and then back again (no net change is made!)
    FArrayBox& shiftWLeft  = (FArrayBox&)a_WLeft;
    FArrayBox& shiftWRight = (FArrayBox&)a_WRight;

    // Solution to the Riemann problem

    // Shift the left and right primitive variable boxes 1/2 cell so they are
    // face centered
    shiftWLeft .shiftHalf(a_dir, 1);
    shiftWRight.shiftHalf(a_dir,-1);

    CH_assert(shiftWLeft .box().contains(a_box));
    CH_assert(shiftWRight.box().contains(a_box));

    const FArrayBox& advVelDir = a_advVel[a_dir];
    CH_assert(advVelDir.box().contains(a_box));

    // Riemann solver computes WGdnv all edges that are not on the physical
    // boundary.
    if (m_useUpwinding) {
        FORT_RIEMANNSOLVER (
            CHF_FRA(a_WGdnv),
            CHF_CONST_FRA(shiftWLeft),
            CHF_CONST_FRA(shiftWRight),
            CHF_CONST_FRA1(advVelDir,0),
            CHF_CONST_INT(a_dir),
            CHF_BOX(a_box));
    } else {
        FORT_AVGSTATES (
            CHF_FRA(a_WGdnv),
            CHF_CONST_FRA(shiftWLeft),
            CHF_CONST_FRA(shiftWRight),
            CHF_CONST_INT(a_dir),
            CHF_BOX(a_box));
    }

    // Shift the left and right primitive variable boxes back to their original
    // position
    shiftWLeft .shiftHalf(a_dir,-1);
    shiftWRight.shiftHalf(a_dir, 1);
}


// -----------------------------------------------------------------------------
// quasilinearUpdate
//
// Computes A * dW/dx where A = dF/dW
// a_AdWdx is CC
// a_WHalf is FC in a_dir
// -----------------------------------------------------------------------------
void MappedAdvectionUtil::quasilinearUpdate (FArrayBox&       a_AdWdx,
                                             const FArrayBox& a_WHalf,
                                             const FArrayBox& a_cellVel,
                                             const Real&      a_scale,
                                             const int&       a_dir,
                                             const Box&       a_box)
{

    // Sanity checks
    CH_assert(a_AdWdx  .box().contains(a_box));
    CH_assert(a_cellVel.box().contains(a_box));
    // CH_assert(a_cellVel.nComp() == SpaceDim);

    FArrayBox cellVelAlias;
    if (a_cellVel.nComp() == 1) {
        cellVelAlias.define(Interval(0,0), (FArrayBox&)a_cellVel);
    } else {
        cellVelAlias.define(Interval(a_dir,a_dir), (FArrayBox&)a_cellVel);
    }

    FORT_QUASILINEARUPDATEF(CHF_FRA(a_AdWdx),
                            CHF_CONST_FRA(a_WHalf),
                            CHF_CONST_FRA1(cellVelAlias,0),
                            CHF_CONST_REAL(a_scale),
                            CHF_CONST_INT(a_dir),
                            CHF_BOX(a_box));
}


// -----------------------------------------------------------------------------
// CTUNormalPred
//
// All inputs and outputs are CC.
// -----------------------------------------------------------------------------
void MappedAdvectionUtil::CTUNormalPred (FArrayBox&       a_WMinus,
                                         FArrayBox&       a_WPlus,
                                         const Real&      a_dt,
                                         const FArrayBox& a_W,
                                         const FArrayBox& a_cellVel,
                                         const int&       a_dir,
                                         const Box&       a_box,
                                         const DataIndex  a_di)
{

    // for CTU, increments are 0 -- straight copy from cells to faces
    a_WMinus.copy(a_W);
    a_WPlus.copy(a_W);
}


// -----------------------------------------------------------------------------
// PLMNormalPred
//
// All inputs and outputs are CC.
// a_oldVel needs SpaceDim comps or an error will be thrown, but only the a_dir
// comp will be used.
// -----------------------------------------------------------------------------
void MappedAdvectionUtil::PLMNormalPred (FArrayBox&       a_WMinus,
                                         FArrayBox&       a_WPlus,
                                         const Real&      a_dt,
                                         const FArrayBox& a_W,
                                         const FArrayBox& a_oldVel,
                                         const int&       a_dir,
                                         const Box&       a_box,
                                         const DataIndex  a_di)
{
    // Sanity checks
    CH_assert(m_levGeoPtr != NULL);
    CH_assert(a_oldVel.nComp() == SpaceDim);

    // Gather needed data
    const ProblemDomain& domain = m_levGeoPtr->getDomain();
    const RealVect& dx = m_levGeoPtr->getDx();

    // This is a domain box that includes ghosts in periodic directions.
    Box perDomBox = domain.domainBox();
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (domain.isPeriodic(dir)) {
            perDomBox.grow(dir, ADVECT_GROW);
        }
    }

    // This will hold 2nd or 4th order slopes
    FArrayBox dW(a_box,1);

    if (m_useFourthOrderSlopes) {
        // 2nd order slopes need to be computed over a larger box to accommodate
        // the 4th order slope computation
        Box boxVL = a_box;
        boxVL.grow(a_dir,1);
        boxVL &= perDomBox; //domain;

        // Compute 2nd order (van Leer) slopes
        FArrayBox dWvL(boxVL, 1);
        m_util.vanLeerSlopes(dWvL, a_W, 1, m_useLimiting, a_dir, boxVL);

        // TODO: setBdrySlopes

        // Compute 4th order slopes, without limiting.
        m_util.fourthOrderSlopes(dW, a_W, dWvL, 1, a_dir, a_box);
    } else {
        // Compute 2nd order (van Leer) slopes
        m_util.vanLeerSlopes(dW, a_W, 1, m_useLimiting, a_dir, a_box);

        // TODO: setBdrySlopes
    }

    // To save on storage, we use the input values as temporaries for the delta's
    a_WMinus.setVal(0.0);
    a_WPlus .setVal(0.0);

    if (m_useFourthOrderSlopes) {
        // Compute one-sided differences as inputs for limiting.
        m_util.oneSidedDifferences(a_WMinus, a_WPlus, a_W, a_dir, a_box);
    }

    // Limiting is already done for 2nd order slopes, so don't do it again
    if (m_useLimiting || m_useFourthOrderSlopes) {
        m_util.slopeLimiter(dW, a_WMinus, a_WPlus, 1, a_box);
    }

    {
        // Alias the a_dir velocity component.
        const Interval velInt(a_dir, a_dir);
        const FArrayBox oldVelDir(velInt, (FArrayBox&)a_oldVel);

        // To the normal prediction
        m_util.PLMNormalPred(a_WMinus,
                             a_WPlus,
                             dW,
                             oldVelDir,
                             a_dt / dx[a_dir],
                             a_box,
                             a_di);
    }

    // Compute the state from the increments
    a_WMinus += a_W;
    a_WPlus  += a_W;
}


// -----------------------------------------------------------------------------
// PPMNormalPred
//
// All inputs and outputs are CC.
// a_oldVel needs SpaceDim comps or an error will be thrown, but only the a_dir
// comp will be used.
// -----------------------------------------------------------------------------
void MappedAdvectionUtil::PPMNormalPred (FArrayBox&       a_WMinus,
                                         FArrayBox&       a_WPlus,
                                         const Real&      a_dt,
                                         const FArrayBox& a_W,
                                         const FArrayBox& a_oldVel,
                                         const int&       a_dir,
                                         const Box&       a_box,
                                         const DataIndex  a_di)
{
    // Sanity checks
    CH_assert(m_levGeoPtr != NULL);

    const int numSlopes = 1;
    CH_assert(a_WMinus.nComp() == numSlopes);
    CH_assert(a_WPlus .nComp() == numSlopes);
    CH_assert(a_W     .nComp() == numSlopes);
    CH_assert(a_oldVel.nComp() == SpaceDim);

    // Gather needed data
    const RealVect& dx = m_levGeoPtr->getDx();

    Box faceBox = a_box;
    // added by petermc, 22 Sep 2008:
    // for 4th order, need extra faces in all the directions
    if (m_useHighOrderLimiter) faceBox.grow(1);
    faceBox.surroundingNodes(a_dir);
    FArrayBox WFace(faceBox, numSlopes);

    // Return WFace on face-centered faceBox.
    this->PPMFaceValues(WFace,
                        a_W,
                        a_dir,
                        faceBox,
                        a_di);

    // To save on storage, we use the input values as temporaries for the deltas
    a_WMinus.setVal(0.0);
    a_WPlus .setVal(0.0);

    a_WMinus -= a_W;
    a_WPlus  -= a_W;

    WFace.shiftHalf(a_dir,1);
    a_WMinus += WFace;

    WFace.shift(a_dir,-1);
    a_WPlus  += WFace;

    if (m_useLimiting) {
        m_util.PPMLimiter(a_WMinus,
                          a_WPlus,
                          a_W,
                          numSlopes,
                          a_dir,
                          a_box);
    }

    {
        // Alias the a_dir velocity component.
        const Interval velInt(a_dir, a_dir);
        const FArrayBox oldVelDir(velInt, (FArrayBox&)a_oldVel);

        // To the normal prediction in characteristic variables
        m_util.PPMNormalPred(a_WMinus,
                             a_WPlus,
                             oldVelDir,
                             a_dt / dx[a_dir],
                             numSlopes,
                             a_box,
                             a_di);
    }

    // Compute the state from the increments
    a_WMinus += a_W;
    a_WPlus  += a_W;
}


// -----------------------------------------------------------------------------
// PPMFaceValues
//
// Given the cell average a_W, compute fourth-order accurate FC values WFace on
// a_box by differentiating the indefinite integral. Limiting is performed in a
// separate pass.
//
// a_WFace will be FC and a_W is CC.
// a_box is the FC box on which a_WFace is computed.
// -----------------------------------------------------------------------------
void MappedAdvectionUtil::PPMFaceValues (FArrayBox&       a_WFace,
                                         const FArrayBox& a_W,
                                         const int        a_dir,
                                         const Box&       a_box,
                                         const DataIndex  a_di)
{

    // Sanity checks
    CH_assert(m_levGeoPtr != NULL);

    const ProblemDomain& domain = m_levGeoPtr->getDomain();
    const int numSlopes = 1;
    CH_assert(a_WFace.nComp() == numSlopes);
    CH_assert(a_W    .nComp() == numSlopes);
    CH_assert(a_WFace.box().contains(a_box));

    // Gather CC J
    FArrayBox ccJ(a_W.box(), 1);
    m_levGeoPtr->fill_J(ccJ);

    // Compute face-averaged J using the exact same averaging weights
    // we will use on the state variable. This ensures conservation.
    FArrayBox edgeJ(a_box, 1);
    if (m_useHighOrderLimiter) {
        // petermc, 7 Jan 2010: changed from a_box to grown a_box.
        Box box1cells = grow(a_box, BASISV(a_dir));
        box1cells.enclosedCells();

        Box loFaces, nextLoFaces;
        Box hiFaces, nextHiFaces;
        Box centerFaces, innerCenterFaces, entireFaces;
        int hasLoFaces, hasHiFaces;
        PeriodicLoHiCenterFace4(loFaces, nextLoFaces, hasLoFaces,
                                hiFaces, nextHiFaces, hasHiFaces,
                                centerFaces, innerCenterFaces, entireFaces,
                                box1cells, domain, a_dir);

        // For i-e/2 in innerCenterFaces, set a_WFace[i-e/2] from a_W[i-2e:i+e].
        // If hasLoFaces:
        // For i-e/2 in loFaces, set a_WFace[i-e/2] from a_W[i:i+3e].
        //           in nextLoFaces,                from a_W[i-e:i+2e].
        // If hasHiFaces:
        // For i-e/2 in hiFaces, set a_WFace[i-e/2] from a_W[i-4e:i-e].
        //           in nextHiFaces,                from a_W[i-3e:i].
        FORT_FOURTHINTERPFACES(CHF_FRA(edgeJ),
                               CHF_CONST_FRA(ccJ),
                               CHF_CONST_INT(numSlopes),
                               CHF_CONST_INT(a_dir),
                               CHF_BOX(loFaces),
                               CHF_BOX(nextLoFaces),
                               CHF_CONST_INT(hasLoFaces),
                               CHF_BOX(hiFaces),
                               CHF_BOX(nextHiFaces),
                               CHF_CONST_INT(hasHiFaces),
                               CHF_BOX(innerCenterFaces));

        // WAS if (a_a_useLimiting) call face limiter;
        // this removed by petermc, 7 Oct 2010

        // dummy statement in order to get around gdb bug
        int dummy_unused = 0; dummy_unused = 0;

    } else { // !m_highOrderLimiter :  this is the old method

        // A box one larger (in direction "a_dir") than the final result box
        // petermc, 14 Aug 2009:  first grow(), then enclosedCells(),
        // rather than reverse order, so you don't end up with empty box1cells.
        Box box1cells(a_box);
        int ghostbox1 = 1;
        box1cells.grow(a_dir, ghostbox1);
        box1cells.enclosedCells();

        FArrayBox dW(box1cells, numSlopes);
        m_util.vanLeerSlopes(dW, ccJ, numSlopes, m_useLimiting, a_dir, box1cells);

        // TODO: setBdrySlopes

        Box loBox,hiBox,centerBox,entireBox;
        int hasLo,hasHi;

        PeriodicLoHiCenterFace(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                               box1cells, domain, a_dir);

        // a_Wface[i-e/2] = (a_W[i-e] + dW[i-e]/3)/2 + (a_W[i] - dW[i]/3)/2
        FORT_PPMFACEVALUESF(CHF_FRA(edgeJ),
                            CHF_CONST_FRA(ccJ),
                            CHF_CONST_FRA(dW),
                            CHF_CONST_INT(numSlopes),
                            CHF_CONST_INT(a_dir),
                            CHF_BOX(loBox),
                            CHF_CONST_INT(hasLo),
                            CHF_BOX(hiBox),
                            CHF_CONST_INT(hasHi),
                            CHF_BOX(centerBox));
    }

    // Weight each cc state element with ccJ.
    FArrayBox JW(a_W.box(), 1);
    JW.copy(a_W);
    JW *= ccJ;

    // Interpolate state to faces.
    if (m_useHighOrderLimiter) {
        // petermc, 7 Jan 2010: changed from a_box to grown a_box.
        Box box1cells = grow(a_box, BASISV(a_dir));
        box1cells.enclosedCells();

        Box loFaces, nextLoFaces;
        Box hiFaces, nextHiFaces;
        Box centerFaces, innerCenterFaces, entireFaces;
        int hasLoFaces, hasHiFaces;
        PeriodicLoHiCenterFace4(loFaces, nextLoFaces, hasLoFaces,
                                hiFaces, nextHiFaces, hasHiFaces,
                                centerFaces, innerCenterFaces, entireFaces,
                                box1cells, domain, a_dir);

        // For i-e/2 in innerCenterFaces, set a_WFace[i-e/2] from a_W[i-2e:i+e].
        // If hasLoFaces:
        // For i-e/2 in loFaces, set a_WFace[i-e/2] from a_W[i:i+3e].
        //           in nextLoFaces,                from a_W[i-e:i+2e].
        // If hasHiFaces:
        // For i-e/2 in hiFaces, set a_WFace[i-e/2] from a_W[i-4e:i-e].
        //           in nextHiFaces,                from a_W[i-3e:i].
        FORT_FOURTHINTERPFACES(CHF_FRA(a_WFace),
                               CHF_CONST_FRA(JW),
                               CHF_CONST_INT(numSlopes),
                               CHF_CONST_INT(a_dir),
                               CHF_BOX(loFaces),
                               CHF_BOX(nextLoFaces),
                               CHF_CONST_INT(hasLoFaces),
                               CHF_BOX(hiFaces),
                               CHF_BOX(nextHiFaces),
                               CHF_CONST_INT(hasHiFaces),
                               CHF_BOX(innerCenterFaces));

        // WAS if (a_a_useLimiting) call face limiter;
        // this removed by petermc, 7 Oct 2010

        // dummy statement in order to get around gdb bug
        int dummy_unused = 0; dummy_unused = 0;

    } else { // !m_highOrderLimiter :  this is the old method

        // A box one larger (in direction "a_dir") than the final result box
        // petermc, 14 Aug 2009:  first grow(), then enclosedCells(),
        // rather than reverse order, so you don't end up with empty box1cells.
        Box box1cells(a_box);
        int ghostbox1 = 1;
        box1cells.grow(a_dir, ghostbox1);
        box1cells.enclosedCells();

        FArrayBox dW(box1cells, numSlopes);
        m_util.vanLeerSlopes(dW, JW, numSlopes, m_useLimiting, a_dir, box1cells);

        // TODO: setBdrySlopes

        Box loBox,hiBox,centerBox,entireBox;
        int hasLo,hasHi;

        PeriodicLoHiCenterFace(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                               box1cells, domain, a_dir);

        // a_Wface[i-e/2] = (a_W[i-e] + dW[i-e]/3)/2 + (a_W[i] - dW[i]/3)/2
        FORT_PPMFACEVALUESF(CHF_FRA(a_WFace),
                            CHF_CONST_FRA(JW),
                            CHF_CONST_FRA(dW),
                            CHF_CONST_INT(numSlopes),
                            CHF_CONST_INT(a_dir),
                            CHF_BOX(loBox),
                            CHF_CONST_INT(hasLo),
                            CHF_BOX(hiBox),
                            CHF_CONST_INT(hasHi),
                            CHF_BOX(centerBox));
    }

    // Remove J scaling, retrieving the face-averaged state.
    a_WFace.divide(edgeJ, a_box, 0, 0, 1);
}

