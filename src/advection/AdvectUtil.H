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
#ifndef AdvectUtil_H__INCLUDED__
#define AdvectUtil_H__INCLUDED__

#include "LevelGeometry.H"
#include "PhysBCUtil.H"
#include "Copier.H"
#include "MappedGodunovUtilities.H"
class LevelGeometry;


// This needs to be at least 3 for PLM and 4 for PPM
#define ADVECT_GROW 4


// -----------------------------------------------------------------------------
// Calculates the upwinded fluxes for the advection terms using the
// CTU, PLM, or PPM methods.
// -----------------------------------------------------------------------------
class MappedAdvectionUtil
{
public:
    enum {
        CTU_NORMAL_PRED = 0,
        PLM_NORMAL_PRED,
        PPM_NORMAL_PRED,
        _NUM_NORMAL_PRED
    };

    // Default constructor - leaves object in a useless state.
    MappedAdvectionUtil ();

    // Full constructor - just calls the define function
    MappedAdvectionUtil (const LevelGeometry* a_levGeoPtr,
                         const int            a_normalPredOrder      = PLM_NORMAL_PRED,
                         const bool           a_useFourthOrderSlopes = true,
                         const bool           a_useLimiting          = false,
                         const bool           a_useHighOrderLimiter  = false,  // Only used in PPM method
                         const bool           a_useUpwinding         = true);

    // Destructor
    virtual ~MappedAdvectionUtil ();

    // Essentially, this is the full constructor.
    // a_normalPredOrder is CTU_NORMAL_PRED, PLM_NORMAL_PRED, or PPM_NORMAL_PRED.
    // If a_useLimiting is true, then a traditional van Leer slope limiter is used.
    // If a_useHighOrderLimiter is true, then an extremum-preserving van Leer slope
    //  limiter and special version of the PPM face interpolator is used.
    virtual void define (const LevelGeometry* a_levGeoPtr,
                         const int            a_normalPredOrder      = PLM_NORMAL_PRED,
                         const bool           a_useFourthOrderSlopes = true,
                         const bool           a_useLimiting          = false,
                         const bool           a_useHighOrderLimiter  = false,
                         const bool           a_useUpwinding         = true);

    // Performs the characteristic tracing of a_Wold and predicts the time-centered,
    // FC a_Whalf. This function does not perform exchanges. If a_returnFlux is
    // false, the user must set final BCs.
    //
    // a_Whalf is FC and time-centered.
    // a_Wold is CC.
    // a_sourceTermPtr is CC and can be NULL.
    // a_oldVel is CC.
    // a_advVel is FC.
    // a_oldTime and a_fluxBCPtr only need to be given if a_returnFlux is true.
    virtual void predictScalar (LevelData<FluxBox>&             a_Whalf,
                                const LevelData<FArrayBox>&     a_Wold,
                                const LevelData<FArrayBox>*     a_sourceTermPtr,
                                const LevelData<FArrayBox>&     a_oldVel,
                                const LevelData<FluxBox>&       a_advVel,
                                const Real                      a_dt,
                                const LevelGeometry&            a_levGeo,
                                Tuple<BCMethodHolder,SpaceDim>  a_BCValues,
                                Tuple<BCMethodHolder,SpaceDim>  a_BCSlopes,
                                const Real                      a_oldTime,
                                const bool                      a_returnFlux = false);

protected:
    // NOTE: This returns the scalar, not the flux!
    // a_WGdnv is the FC solution to the Riemann problem.
    // a_WLeft is the state predicted from the left. It is CC.
    // a_WRight is the state predicted from the right. It is CC.
    // a_advVel is the advecting velocity. It is FC.
    virtual void RiemannSolver (FArrayBox&       a_WGdnv,
                                const FArrayBox& a_WLeft,
                                const FArrayBox& a_WRight,
                                const FluxBox&   a_advVel,
                                const int        a_dir,
                                const Box&       a_box);

    // Computes A * dW/dx where A = dF/dW
    // a_AdWdx is CC
    // a_WHalf is FC in a_dir
    virtual void quasilinearUpdate (FArrayBox&       a_AdWdx,
                                    const FArrayBox& a_WHalf,
                                    const FArrayBox& a_cellVel,
                                    const Real&      a_scale,
                                    const int&       a_dir,
                                    const Box&       a_box);

    // All inputs and outputs are CC.
    virtual void CTUNormalPred (FArrayBox&       a_WMinus,
                                FArrayBox&       a_WPlus,
                                const Real&      a_dt,
                                const FArrayBox& a_W,
                                const FArrayBox& a_cellVel,
                                const int&       a_dir,
                                const Box&       a_box,
                                const DataIndex  a_di);

    // All inputs and outputs are CC.
    // a_oldVel needs SpaceDim comps or an error will be thrown, but only the
    // a_dir comp will be used.
    virtual void PLMNormalPred (FArrayBox&       a_WMinus,
                                FArrayBox&       a_WPlus,
                                const Real&      a_dt,
                                const FArrayBox& a_W,
                                const FArrayBox& a_oldVel,
                                const int&       a_dir,
                                const Box&       a_box,
                                const DataIndex  a_di);

    // All inputs and outputs are CC.
    // a_oldVel needs SpaceDim comps or an error will be thrown, but only the
    // a_dir comp will be used.
    virtual void PPMNormalPred (FArrayBox&       a_WMinus,
                                FArrayBox&       a_WPlus,
                                const Real&      a_dt,
                                const FArrayBox& a_W,
                                const FArrayBox& a_oldVel,
                                const int&       a_dir,
                                const Box&       a_box,
                                const DataIndex  a_di);

    // Given the cell average a_W, compute fourth-order accurate FC values WFace on
    // a_box by differentiating the indefinite integral. Limiting is performed in a
    // separate pass.
    //
    // a_WFace will be FC and a_W is CC.
    // a_box is the FC box on which a_WFace is computed.
    virtual void PPMFaceValues (FArrayBox&       a_WFace,
                                const FArrayBox& a_W,
                                const int        a_dir,
                                const Box&       a_box,
                                const DataIndex  a_di);

    // Member variables
    const LevelGeometry*   m_levGeoPtr;
    MappedGodunovUtilities m_util;

    int  m_normalPredOrder;
    bool m_useFourthOrderSlopes;
    bool m_useLimiting;
    bool m_useHighOrderLimiter;
    bool m_useUpwinding;
};


#endif //!AdvectUtil_H__INCLUDED__

