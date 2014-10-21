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
#ifndef __DIVERGENCE_H__INCLUDED__
#define __DIVERGENCE_H__INCLUDED__

#include "LevelData.H"
#include "Tuple.H"
class BCMethodHolder;
class LevelGeometry;
class MappedQuadCFInterp;
class MappedLevelFluxRegister;

extern const Real BOGUS_TIME;


// -----------------------------------------------------------------------------
// Class to encapsulate Divergence functions
// -----------------------------------------------------------------------------
class Divergence
{
public:
    // The divergence functions will call setGhosts() on a_u before the
    // calculation is performed. It would be more efficient, however, for
    // the ghosts to be set before calling the div function in order to
    // cache the extrapolated values. This is why the default a_fluxBC
    // is set to NULL.
    typedef Tuple<BCMethodHolder,SpaceDim> BC_type;


    // Single-level divergence functions ---------------------------------------

    /** computes cell-centered level-operator divergence of
        edge-centered vector field uEdge; assumes all coarse-fine BC's
        have already been set. */
    static void levelDivergenceMAC (LevelData<FArrayBox>&       a_div,
                                    const LevelData<FluxBox>&   a_uEdge,
                                    const LevelGeometry&        a_levGeo,
                                    const Real                  a_time = BOGUS_TIME,
                                    const BC_type*              a_fluxBC = NULL);


    /// really basic single-fab MAC divergence
    static void simpleDivergenceMAC (FArrayBox&             a_div,
                                     const FluxBox&         a_uEdge,
                                     const DataIndex&       a_di,
                                     const LevelGeometry&   a_levGeo,
                                     const Real             a_time = BOGUS_TIME,
                                     const BC_type*         a_fluxBC = NULL);


    /** computes cell-centered level-operator divergence of
        cell-centered vector field u; if uCrse != NULL, does coarse-fine
        boundary conditions for u -- if  quadInterp == true, uses quadratic
        coarse-fine boundary conditions, otherwise, use extrap BC's */
    static void levelDivergenceCC (LevelData<FArrayBox>&       a_div,
                                   LevelData<FArrayBox>&       a_u,
                                   const LevelData<FArrayBox>* a_uCrsePtr,
                                   const bool                  a_quadInterp,
                                   const LevelGeometry&        a_levGeo,
                                   const Real                  a_time = BOGUS_TIME,
                                   const BC_type*              a_fluxBC = NULL);

    /** computes cell-centered level-operator divergence of
        cell-centered vector field u; if uCrse != NULL, does coarse-fine
        boundary conditions for u -- if  quadInterp == true, uses quadratic
        coarse-fine boundary conditions, otherwise, use extrap BC's.  This
        one takes a pre-allocated QuadCFInterp object as an argument */
    static void levelDivergenceCC (LevelData<FArrayBox>&       a_div,
                                   LevelData<FArrayBox>&       a_u,
                                   const LevelData<FArrayBox>* a_uCrsePtr,
                                   const bool                  a_quadInterp,
                                   MappedQuadCFInterp&         a_cfInterp,
                                   const LevelGeometry&        a_levGeo,
                                   const Real                  a_time = BOGUS_TIME,
                                   const BC_type*              a_fluxBC = NULL);


    // Composite divergence functions ------------------------------------------

    /** computes composite cell-centered divergence of edge-centered
        vector field uEdge; if finer level data uEdgeFine exists, use
        flux-matching condition to compute divergence along coarse-fine
        interface. */
    static void compDivergenceMAC (LevelData<FArrayBox>&          a_div,
                                   LevelData<FluxBox>&            a_uEdge,
                                   const LevelData<FluxBox>*      a_uEdgeFinePtr,
                                   const LevelGeometry&           a_levGeo,
                                   const Real                     a_time = BOGUS_TIME,
                                   const BC_type*                 a_fluxBC = NULL);

    /** just like normal compDivergenceMAC, but pass in a predefined
        flux register (more efficient) */
    static void compDivergenceMAC (LevelData<FArrayBox>&     a_div,
                                   LevelData<FluxBox>&       a_uEdge,
                                   const LevelData<FluxBox>* a_uEdgeFinePtr,
                                   MappedLevelFluxRegister*  a_fluxRegPtr,
                                   const LevelGeometry&      a_levGeo,
                                   const Real                a_time = BOGUS_TIME,
                                   const BC_type*            a_fluxBC = NULL);

    /** computes cell-centered composite-operator divergence
        of cell-centered vector field u; uses same coarse-fine BC's as
        levelDivergenceCC -- if uFine != NULL, also does flux-matching
        BC with finer level uFine, using same coarse-fine BC's as used
        for coarse level BC. */
    static void compDivergenceCC (LevelData<FArrayBox>& a_div,
                                  LevelData<FArrayBox>& a_u,
                                  const LevelData<FArrayBox>* a_uCrsePtr,
                                  LevelData<FArrayBox>* a_uFinePtr,
                                  const bool            a_quadInterp,
                                  const LevelGeometry&  a_levGeo,
                                  const Real            a_time = BOGUS_TIME,
                                  const BC_type*        a_fluxBC = NULL);

    /** computes cell-centered composite-operator divergence
        of cell-centered vector field u; uses same coarse-fine BC's as
        levelDivergenceCC -- if uFine != NULL, also does flux-matching
        BC with finer level uFine, using same coarse-fine BC's as used
        for coarse level BC. In this one, pre-allocated LevelFluxRegister
        and QuadCFInterp objects are passed in as an argument (saves
        allocation inside fn.)*/
    static void compDivergenceCC (LevelData<FArrayBox>&       a_div,
                                  LevelData<FArrayBox>&       a_u,
                                  const LevelData<FArrayBox>* a_uCrsePtr,
                                  LevelData<FArrayBox>*       a_uFinePtr,
                                  const bool                  a_quadInterp,
                                  MappedLevelFluxRegister*    a_fluxRegFinePtr,
                                  MappedQuadCFInterp&         a_cfInterpCrse,
                                  MappedQuadCFInterp&         a_cfInterpFine,
                                  const LevelGeometry&        a_levGeo,
                                  const Real                  a_time = BOGUS_TIME,
                                  const BC_type*              a_fluxBC = NULL);
};


#endif //!__DIVERGENCE_H__INCLUDED__
