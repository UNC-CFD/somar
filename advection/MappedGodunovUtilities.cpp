#include "MappedGodunovUtilities.H"
#include "MappedGodunovUtilitiesF_F.H"
#include "PeriodicLoHiCenter.H"
#include "LoHiCenter.H"
#include "Debug.H"


// -----------------------------------------------------------------------------
// Constructor.
// Flag everything as not defined or set
// -----------------------------------------------------------------------------
MappedGodunovUtilities::MappedGodunovUtilities ()
: m_levGeoPtr(NULL)
{
    m_highOrderLimiter = false;
    m_isDefined = false;
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
MappedGodunovUtilities::~MappedGodunovUtilities ()
{
    m_levGeoPtr = NULL;
}


// -----------------------------------------------------------------------------
// Define this object and the boundary condition object
// -----------------------------------------------------------------------------
void MappedGodunovUtilities::define (const LevelGeometry* a_levGeoPtr)
{
    CH_assert(a_levGeoPtr != NULL);
    m_levGeoPtr = a_levGeoPtr;

    // Store the domain and grid spacing
    m_domain    = a_levGeoPtr->getDomain();
    m_dx        = a_levGeoPtr->getDx();
    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Compute componentwise van Leer slopes.
//  Given cell averages W, compute van Leer slopes dW on a
//  component-by-component basis.
// -----------------------------------------------------------------------------
void MappedGodunovUtilities::vanLeerSlopes (FArrayBox&       a_dW,
                                            const FArrayBox& a_W,
                                            const int&       a_numSlopes,
                                            const bool&      a_useLimiting,
                                            const int&       a_dir,
                                            const Box&       a_box)
{
    // A box one larger (in direction "a_dir") than the final result box
    // 2 Sep 2008:  For vanLeerSlopesExtPreserving, expand by 2?  17 Sep 2008
    Box box1 = a_box;
    int ghostbox1 = 1; // FIX, 19 sep 2008 (m_highOrderLimiter) ? 2 : 1;
    // int ghostbox1 = (m_highOrderLimiter) ? 2 : 1;
    box1.grow(a_dir, ghostbox1);

    // Compute where centered differences can be used and where one-sided
    // differences need to be used.
    Box loBox,hiBox,centerBox,entireBox;
    int hasLo,hasHi;

    PeriodicLoHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                       box1,m_domain,a_dir);

    // Compute 2nd order slopes - including one-sided differences
    FArrayBox dWMinus(entireBox,a_numSlopes);
    FArrayBox dWPlus (entireBox,a_numSlopes);

    // We calculate a_dW and dWMinus and dWPlus on centerBox,
    // and we need a_W on centerBox grown by 1 in direction a_dir.
    slopes(a_dW, dWMinus, dWPlus, a_W, a_numSlopes, a_dir,
           loBox, hiBox, centerBox, entireBox, hasLo, hasHi);

    // Apply the slope limiter if requested
    if (a_useLimiting) {
        // Apply slopeLimiter only on centerBox; elsewhere, a_dW is unchanged.

        // 2 Sep 2008:  replace slopeLimiter with slopeLimiterExtPreserving

        if (m_highOrderLimiter) {
            // Modifies a_dW on centerBox,
            // and needs dWMinus on centerBox including shift down by 1 in direction a_dir,
            // and needs dWPlus on centerBox including shift up by 1 in direction a_dir.
            slopeLimiterExtPreserving(a_dW, dWMinus, dWPlus, a_numSlopes, centerBox, a_dir);
        } else {
            slopeLimiter(a_dW, dWMinus, dWPlus, a_numSlopes, centerBox);
        }
    }
}


// -----------------------------------------------------------------------------
// Compute fourth-order slopes.
//  Given cell averages W and van Leer slopes dWvL, compute fourth-order
//  slopes dW4. Limiting is performed in a separate pass.
// -----------------------------------------------------------------------------
void MappedGodunovUtilities::fourthOrderSlopes (FArrayBox&       a_dW,
                                                const FArrayBox& a_W,
                                                const FArrayBox& a_dWvL,
                                                const int&       a_numSlopes,
                                                const int&       a_dir,
                                                const Box&       a_box)
{
    // Number of slopes to compute
    int numSlope = a_numSlopes;

    CH_assert(a_dW.nComp() == numSlope);
    CH_assert(a_W.nComp() >= numSlope);

    // A box one larger (in direction "a_dir") than the final result box
    Box box1 = a_box;
    box1.grow(a_dir,1);

    // Compute where centered differences can be used and where one sided
    // differences need to be used.
    Box loBox,hiBox,centerBox,entireBox;
    int hasLo,hasHi;

    PeriodicLoHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                       box1,m_domain,a_dir);

    CH_assert(a_dW.box().contains(entireBox));

    FORT_FOURTHSLOPEDIFFSF (
        CHF_FRA(a_dW),
        CHF_CONST_FRA(a_W),
        CHF_CONST_FRA(a_dWvL),
        CHF_CONST_INT(numSlope),
        CHF_CONST_INT(a_dir),
        CHF_BOX(loBox),
        CHF_CONST_INT(hasLo),
        CHF_BOX(hiBox),
        CHF_CONST_INT(hasHi),
        CHF_BOX(centerBox));
}


// -----------------------------------------------------------------------------
// Compute slopes (dW- and dW+) using one sided differences
// -----------------------------------------------------------------------------
void MappedGodunovUtilities::oneSidedDifferences (FArrayBox&       a_dWMinus,
                                                  FArrayBox&       a_dWPlus,
                                                  const FArrayBox& a_W,
                                                  const int&       a_dir,
                                                  const Box&       a_box)
{
    const int numSlopes = a_dWMinus.nComp();

    // A box one larger (in direction "a_dir") than the final result box
    Box box1 = a_box;
    box1.grow(a_dir,1);

    // Compute where centered differences can be used and where one sided
    // differences need to be used.
    Box loBox,hiBox,centerBox,entireBox;
    int hasLo,hasHi;

    PeriodicLoHiCenter(loBox,hasLo,hiBox,hasHi,centerBox,entireBox,
                       box1,m_domain,a_dir);

    // Compute 2nd order slopes - including one sided differences
    FArrayBox deltaWC(entireBox, numSlopes);

    slopes(deltaWC, a_dWMinus, a_dWPlus, a_W, numSlopes, a_dir,
           loBox, hiBox, centerBox, entireBox, hasLo, hasHi);
}


// -----------------------------------------------------------------------------
// Compute slopes (dW (center), dW-, and dW+)
//  a_dwCent, a_dwMinus, a_dWPlus, and a_dW all live on
//  cell-centered a_entireBox.

//  For i in a_centerBox:
//  a_dwCent[i] = (a_W[i+1] - a_W[i-1])/2;
//  a_dWMinus[i] = a_W[i] - a_W[i-1];
//  a_dWPlus[i] = a_W[i+1] - a_W[i].

//  For i in a_loBox, set only a_dWPlus[i] = a_W[i+1] - a_W[i]
//  and copy it to a_dWCent[i] and a_dWMinus[i].

//  For i in a_hiBox, set only a_dWMinus[i] = a_W[i] - a_W[i-1]
//  and copy it to a_dWCent[i] and a_dWPlus[i].
// -----------------------------------------------------------------------------
void MappedGodunovUtilities::slopes (FArrayBox&       a_dWCent,
                                     FArrayBox&       a_dWMinus,
                                     FArrayBox&       a_dWPlus,
                                     const FArrayBox& a_W,
                                     const int&       a_numSlopes,
                                     const int&       a_dir,
                                     const Box&       a_loBox,
                                     const Box&       a_hiBox,
                                     const Box&       a_centerBox,
                                     const Box&       a_entireBox,
                                     const int&       a_hasLo,
                                     const int&       a_hasHi)
{
    CH_assert(a_dWCent .nComp() == a_numSlopes);
    CH_assert(a_dWMinus.nComp() == a_numSlopes);
    CH_assert(a_dWPlus .nComp() == a_numSlopes);
    CH_assert(a_W.nComp() >= a_numSlopes);

    CH_assert(a_dWCent .box().contains(a_entireBox));
    CH_assert(a_dWMinus.box().contains(a_entireBox));
    CH_assert(a_dWPlus .box().contains(a_entireBox));
    CH_assert(a_W.box().contains( a_entireBox));

    CH_assert((a_dir >= 0 )&& (a_dir < SpaceDim));

    FORT_SECONDSLOPEDIFFSF (
        CHF_FRA(a_dWCent),
        CHF_FRA(a_dWMinus),
        CHF_FRA(a_dWPlus),
        CHF_CONST_FRA(a_W),
        CHF_CONST_INT(a_numSlopes),
        CHF_CONST_INT(a_dir),
        CHF_BOX(a_loBox),
        CHF_CONST_INT(a_hasLo),
        CHF_BOX(a_hiBox),
        CHF_CONST_INT(a_hasHi),
        CHF_BOX(a_centerBox));
}


// -----------------------------------------------------------------------------
// Apply a van Leer limiter directly to the slopes.
//  On input, dW contains the centered, unlimited slopes, and
//  dW(Minus,Plus) contain the one-sided slopes from the minus, plus sides.
//  On output, dW contains the limited slopes.
//  slopes dW4. Limiting is performed in a separate pass.
// -----------------------------------------------------------------------------
void MappedGodunovUtilities::slopeLimiter(FArrayBox&       a_dW,
                                          const FArrayBox& a_dWLeft,
                                          const FArrayBox& a_dWRigh,
                                          const int&       a_numSlopes,
                                          const Box&       a_box)
{
    CH_assert(m_isDefined);
    CH_assert(a_dW.nComp()     == a_numSlopes);
    CH_assert(a_dWLeft.nComp() == a_numSlopes);
    CH_assert(a_dWRigh.nComp() == a_numSlopes);
    CH_assert(a_dW.box().contains(a_box));
    CH_assert(a_dWLeft.box().contains(a_box));
    CH_assert(a_dWRigh.box().contains(a_box));

    FORT_VANLEERLIMITERF (
        CHF_FRA(a_dW),
        CHF_CONST_FRA(a_dWLeft),
        CHF_CONST_FRA(a_dWRigh),
        CHF_CONST_INT(a_numSlopes),
        CHF_BOX(a_box));
}


// -----------------------------------------------------------------------------
// Apply an extremum-preserving van Leer limiter directly to the slopes.
//  On input, dW contains the centered, unlimited slopes, and
//  dW(Minus,Plus) contain the one-sided slopes from the minus, plus sides.
//  On output, dW contains the limited slopes.
//  slopes dW4. Limiting is performed in a separate pass.
// -----------------------------------------------------------------------------
void MappedGodunovUtilities::slopeLimiterExtPreserving (FArrayBox&       a_dW,
                                                        const FArrayBox& a_dWLeft,
                                                        const FArrayBox& a_dWRigh,
                                                        const int&       a_numSlopes,
                                                        const Box&       a_box,
                                                        const int&       a_dir)
{
    CH_assert(m_isDefined);
    CH_assert(a_dW.nComp()     == a_numSlopes);
    CH_assert(a_dWLeft.nComp() == a_numSlopes);
    CH_assert(a_dWRigh.nComp() == a_numSlopes);
    CH_assert(a_dW.box().contains(a_box));
    CH_assert(a_dWLeft.box().contains(a_box));
    CH_assert(a_dWRigh.box().contains(a_box));

    // Compute where centered differences can be used and where one-sided
    // differences need to be used.
    Box loBox, hiBox, centerBox, entireBox;
    int hasLo, hasHi;

    PeriodicLoHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                       a_box, m_domain, a_dir);

    if (hasLo) {
        FORT_VANLEERLIMITERF (
            CHF_FRA(a_dW),
            CHF_CONST_FRA(a_dWLeft),
            CHF_CONST_FRA(a_dWRigh),
            CHF_CONST_INT(a_numSlopes),
            CHF_BOX(loBox));
    }
    if (hasHi) {
        FORT_VANLEERLIMITERF (
            CHF_FRA(a_dW),
            CHF_CONST_FRA(a_dWLeft),
            CHF_CONST_FRA(a_dWRigh),
            CHF_CONST_INT(a_numSlopes),
            CHF_BOX(hiBox));
    }
    if (!centerBox.isEmpty()) {
        // Modifies a_dW on centerBox,
        // and needs a_dWLeft on centerBox including shift down by 1 in direction a_dir,
        // and needs a_dWRigh on centerBox including shift up by 1 in direction a_dir.
        FORT_EXTPRESERVINGVANLEERLIMITERF (
            CHF_FRA(a_dW),
            CHF_CONST_FRA(a_dWLeft),
            CHF_CONST_FRA(a_dWRigh),
            CHF_CONST_INT(a_numSlopes),
            CHF_CONST_INT(a_dir),
            CHF_BOX(centerBox));
    }
    // dummy statement in order to get around gdb bug
    int dummy_unused = 0; dummy_unused = 0;
}


// -----------------------------------------------------------------------------
// Piecewise linear normal predictor. Computes increments in the
// characteristic amplitudes.
// -----------------------------------------------------------------------------
void MappedGodunovUtilities::PLMNormalPred(FArrayBox&       a_dWCharLo,
                                           FArrayBox&       a_dWCharHi,
                                           const FArrayBox& a_dWChar,
                                           const FArrayBox& a_Lambda,
                                           const Real&      a_dtbydx,
                                           const Box&       a_box,
                                           const DataIndex  a_di)
{
    const int numPrim = a_dWChar.nComp();

    CH_assert(a_dWCharLo.nComp() == numPrim);
    CH_assert(a_dWCharHi.nComp() == numPrim);
    CH_assert(a_Lambda.nComp() == numPrim);
    CH_assert(a_dWCharLo.box().contains(a_box));
    CH_assert(a_dWCharHi.box().contains(a_box));
    CH_assert(a_dWChar.box().contains(a_box));
    CH_assert(a_Lambda.box().contains(a_box));

    FORT_PLMNORMALPREDF (
        CHF_FRA(a_dWCharLo),
        CHF_FRA(a_dWCharHi),
        CHF_CONST_FRA(a_dWChar),
        CHF_CONST_FRA(a_Lambda),
        CHF_CONST_REAL(a_dtbydx),
        CHF_CONST_INT(numPrim),
        CHF_BOX(a_box));
}


// -----------------------------------------------------------------------------
// PPM Limiter.
//  a_dWMinus and a_dWPlus are the differences between the face values
//  on the minus and plus sides of cells and the average in the cell.
//  That is,
//  a_dWMinus[i] = WFace[i - e/2] - a_W[i]
//  a_dWPlus[i] = WFace[i + e/2] - a_W[i]
//  where e is the unit vector in dimension a_dir.
//  The PPM limiter is applied to these values to obtain a monotone
//  interpolant in the cell.
//  The function returns the limited a_dWMinus and a_dWPlus on a_box.
//  petermc, 4 Sep 2008:  included a_W in argument list

//  If m_highOrderLimiter,
//  then need a_dWMinus and a_dWPlus on a_box,
//  and need a_W on on a_box grown by 3 in dimension a_dir.
//  Returns limited a_dWMinus and a_dWPlus on a_box.
// -----------------------------------------------------------------------------
void MappedGodunovUtilities::PPMLimiter (FArrayBox& a_dWMinus,
                                         FArrayBox& a_dWPlus,
                                         const FArrayBox& a_W,
                                         const int& a_numSlopes,
                                         const int& a_dir,
                                         const Box& a_box)
{
    // Called by PatchGodunov::PPMNormalPred,
    // which is called by PatchGodunov::computeWHalf.
    // a_dWMinus[i] = WFace[i - e/2] - a_W[i]
    // a_dWPlus[i] = WFace[i + e/2] - a_W[i]
    // where e is unit vector in dimension a_dir.

    if (m_highOrderLimiter) {
        // this option added by petermc, 5 Sep 2008
        // Will need to recalculate some D^2's.

        // We calculate d2Wfcf on a_box,
        // and we need a_dWMinus on a_box,
        // and we need a_dWPlus on a_box.
        // d2Wfcf[i] = 6 * (a_dWMinus[i] + a_dWPlus[i])
        //           = 6 * (thisFaceWDir[i-e/2] - a_cellW[i] +
        //                  thisFaceWDir[i+e/2] - a_cellW[i])
        FArrayBox d2Wfcf(a_box, a_numSlopes);
        d2Wfcf.copy(a_dWMinus);
        d2Wfcf.plus(a_dWPlus, 0, 0, a_numSlopes);
        d2Wfcf *= 6.;

        // petermc, 21 Sep 2010:
        // In order to get a_dWMinus and a_dWPlus on a_box,
        // we need d2W on a_box grown by 3 in a_dir directions.
        Box box3 = a_box;
        box3.grow(a_dir, 3);

        Box loBox, hiBox, centerBox, entireBox;
        int hasLo, hasHi;
        PeriodicLoHiCenter(loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                           box3, m_domain, a_dir);

        FArrayBox d2W(entireBox, a_numSlopes);
        // On centerBox, use 3-point stencil for d2W;
        // on loBox and hiBox, copy result from neighbor.
        // In the case of centerBox == entireBox == box1,
        // where box1 is a_box grown by 1 in dimension a_dir:
        // We calculate d2W on a_box grown by 2 (was 1) in dimension a_dir,
        // and we need a_W on a_box grown by 3 (was 2) in dimension a_dir.

        // petermc, 21 Sep 2010, changed layer of d2W from 2 to 1,
        // and of a_W from 3 to 2.
        FORT_GETSECONDDIFF (
            CHF_FRA(d2W),
            CHF_CONST_FRA(a_W),
            CHF_CONST_INT(a_numSlopes),
            CHF_CONST_INT(a_dir),
            CHF_BOX(loBox),
            CHF_CONST_INT(hasLo),
            CHF_BOX(hiBox),
            CHF_CONST_INT(hasHi),
            CHF_BOX(centerBox));

        Box box1 = a_box;
        box1.grow(a_dir, 1);
        Box nextLoBox, nextHiBox, innerCenterBox;
        PeriodicLoHiCenter5(loBox, nextLoBox, hasLo,
                            hiBox, nextHiBox, hasHi,
                            centerBox, innerCenterBox, entireBox,
                            box1, m_domain, a_dir);
        CH_assert(entireBox == a_box);

        Real limitC = 1.25;
        Real eps = 1.0e-12;
        Real c3 = 0.1;
        FORT_CHECKCUBICLIMITERF (
            CHF_FRA(a_dWMinus), // <W>_(i-1/2) - <W>_i, on a_box
            CHF_FRA(a_dWPlus),  // <W>_(i+1/2) - <W>_i, on a_box
            CHF_CONST_FRA(a_W),
            CHF_CONST_FRA(d2W),
            CHF_CONST_FRA(d2Wfcf),
            CHF_CONST_INT(a_numSlopes),
            CHF_CONST_INT(a_dir),
            CHF_BOX(loBox),
            CHF_BOX(nextLoBox),
            CHF_CONST_INT(hasLo),
            CHF_BOX(hiBox),
            CHF_BOX(nextHiBox),
            CHF_CONST_INT(hasHi),
            CHF_BOX(innerCenterBox),
            CHF_CONST_REAL(limitC),
            CHF_CONST_REAL(c3),
            CHF_CONST_REAL(eps));

        // dummy statement in order to get around gdb bug
        int dummy_unused = 0; dummy_unused = 0;

    } else {
        FORT_PPMLIMITERF (
            CHF_FRA(a_dWMinus), // <W>_(i-1/2) - <W>_i, on a_box
            CHF_FRA(a_dWPlus),  // <W>_(i+1/2) - <W>_i, on a_box
            CHF_CONST_INT(a_numSlopes),
            CHF_BOX(a_box));
    }
}


// -----------------------------------------------------------------------------
// Piecewise Parabolic Method normal predictor.
//  On input, dW(Minus,Plus), contain the characteristic
//  expansions of the differences between the (minus, plus) face values
//  and the cell average. On output, dW(Minus,Plus) contain the
//  characteristic amplitudes of the corrections required to compute
//  the normal predictor.
// -----------------------------------------------------------------------------
void MappedGodunovUtilities::PPMNormalPred (FArrayBox&       a_dWMinus,
                                            FArrayBox&       a_dWPlus,
                                            const FArrayBox& a_Lambda,
                                            const Real&      a_dtbydx,
                                            const int&       a_numSlopes,
                                            const Box&       a_box,
                                            const DataIndex  a_di)
{
    FORT_PPMNORMALPREDF (
        CHF_FRA(a_dWMinus),
        CHF_FRA(a_dWPlus),
        CHF_CONST_FRA(a_Lambda),
        CHF_CONST_REAL(a_dtbydx),
        CHF_CONST_INT(a_numSlopes),
        CHF_BOX(a_box));
}


// -----------------------------------------------------------------------------
// Set whether to use high-order limiter.
// -----------------------------------------------------------------------------
void MappedGodunovUtilities::highOrderLimiter (bool a_highOrderLimiter)
{
    m_highOrderLimiter = a_highOrderLimiter;
}

