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
#ifndef __MappedGodunovUtilities_H__INCLUDED__
#define __MappedGodunovUtilities_H__INCLUDED__

#include "LevelGeometry.H"
#include "ProblemDomain.H"
#include "FArrayBox.H"
#include "RealVect.H"
#include "FluxBox.H"


// -----------------------------------------------------------------------------
// Utility class for higher-order Godunov methods: slopes, parabolic
// interpolants, limiters. Contains no physics-dependent methods, but
// one of the member functions (PPMFaceValues()) may require a pointer
// to a GodunovPhysics class in order to set boundary slopes.
// Basically, this is a rewrite of Chombo's GodunovUtilities.
// -----------------------------------------------------------------------------
class MappedGodunovUtilities
{
public:
    // Default constructor
    MappedGodunovUtilities ();

    // Destructor
    virtual ~MappedGodunovUtilities ();

    // Full define constructor
    virtual void define (const LevelGeometry* a_levGeoPtr);

    // Slopes...

    // Compute componentwise van Leer slopes.
    //  Given cell averages W, compute van Leer slopes dW on a
    //  component-by-component basis.
    virtual void vanLeerSlopes (FArrayBox&       a_dW,
                                const FArrayBox& a_W,
                                const int&       a_numSlopes,
                                const bool&      a_useLimiting,
                                const int&       a_dir,
                                const Box&       a_box);

    // Compute fourth-order slopes.
    //  Given cell averages W and van Leer slopes dWvL, compute fourth-order
    //  slopes dW4. Limiting is performed in a separate pass.
    virtual void fourthOrderSlopes (FArrayBox&       a_dW4,
                                    const FArrayBox& a_W,
                                    const FArrayBox& a_dWvL,
                                    const int&       a_numSlopes,
                                    const int&       a_dir,
                                    const Box&       a_box);

    // Compute slopes (dW- and dW+) using one sided differences
    virtual void oneSidedDifferences (FArrayBox&       a_dWMinus,
                                      FArrayBox&       a_dWPlus,
                                      const FArrayBox& a_W,
                                      const int&       a_dir,
                                      const Box&       a_box);

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
    virtual void slopes (FArrayBox&       a_dWCent,
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
                         const int&       a_hasHi);


    // Parabolic interpolants...

    // PLM normal predictor.
    // Compute the increments in the characteristic amplitudes.
    virtual void PLMNormalPred (FArrayBox&       a_dWCharMinus,
                                FArrayBox&       a_dWCharPlus,
                                const FArrayBox& a_dWChar,
                                const FArrayBox& a_Lambda,
                                const Real&      a_dtbydx,
                                const Box&       a_box,
                                const DataIndex  a_di);

    // PPM normal predictor.
    //  On input, dW(Minus,Plus), contain the characteristic
    //  expansions of the differences between the (minus, plus) face values
    //  and the cell average. On output, dW(Minus,Plus) contain the
    //  characteristic amplitudes of the corrections required to compute
    //  the normal predictor.
    virtual void PPMNormalPred (FArrayBox&       a_dWMinus,
                                FArrayBox&       a_dWPlus,
                                const FArrayBox& a_Lambda,
                                const Real&      a_dtbydx,
                                const int&       a_numSlopes,
                                const Box&       a_box,
                                const DataIndex  a_di);


    // Limiters...

    // van Leer slope limiter.
    //  On input, dW contains the centered, unlimited slopes, and
    //  dW(Minus,Plus) contain the one-sided slopes from the minus, plus sides.
    //  On output, dW contains the limited slopes.
    //  slopes dW4. Limiting is performed in a separate pass.
    virtual void slopeLimiter (FArrayBox&       a_dW,
                               const FArrayBox& a_dWMinus,
                               const FArrayBox& a_dWPlus,
                               const int&       a_numSlopes,
                               const Box&       a_box);

    // extremum-preserving van Leer slope limiter.
    //  On input, dW contains the centered, unlimited slopes, and
    //  dW(Minus,Plus) contain the one-sided slopes from the minus, plus sides.
    //  On output, dW contains the limited slopes.
    //  slopes dW4. Limiting is performed in a separate pass.
    virtual void slopeLimiterExtPreserving (FArrayBox&       a_dW,
                                            const FArrayBox& a_dWMinus,
                                            const FArrayBox& a_dWPlus,
                                            const int&       a_numSlopes,
                                            const Box&       a_box,
                                            const int&       a_dir);

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
    virtual void PPMLimiter (FArrayBox& a_dWMinus,
                             FArrayBox& a_dWPlus,
                             const FArrayBox& a_W,
                             const int& a_numSlopes,
                             const int& a_dir,
                             const Box& a_box);

    // Set whether to use high-order limiter.
    virtual void highOrderLimiter (bool a_highOrderLimiter);

protected:
    // Problem domain and grid spacing
    ProblemDomain        m_domain;
    RealVect             m_dx;
    const LevelGeometry* m_levGeoPtr;

    // Has this object been defined
    bool m_isDefined;

    // Use a high-order limiter?  (default false)
    bool m_highOrderLimiter;

private:
    // We may allow copying and assignment later.
    // Disallowed for all the usual reasons
    void operator= (const MappedGodunovUtilities&);
    MappedGodunovUtilities (const MappedGodunovUtilities&);
};


#endif //!__MappedGodunovUtilities_H__INCLUDED__
