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
#ifndef __PHYSBCUTIL_H__INCLUDED__
#define __PHYSBCUTIL_H__INCLUDED__

#include "EllipticBCInterface.H"
#include "LevelGeometry.H"
#include "LevelData.H"
#include "Tuple.H"


// -----------------------------------------------------------------------------
// Encapsulates all physical BC's for mapped AMRINS solves.
//
// This class is designed to contain all the information needed to specify all
// physical BCs needed by the AMR code. This includes BC's for velocity,
// pressures, and scalars. While basic solid-wall BCs are coded in to this
// class, the idea is that for a more complicated set of physical boundary
// conditions (inflow, etc), a problem-dependent class will be derived from this
// one which will override the generic BC's in this class.
// -----------------------------------------------------------------------------
class PhysBCUtil
{
public:
    // Constructor
    PhysBCUtil ();

    // Destructor
    virtual ~PhysBCUtil ();

    // Derived classes should be thier own factories.
    virtual PhysBCUtil* newPhysBCUtil () const = 0;

    // Reads BC data from file
    virtual void define ();


    // ICs / background fields -------------------------------------------------
    // These functions have defaults that just set the FABs to zero and should
    // be overriden as needed.

    // Fills a FAB with the initial velocity.
    // Locations are in mapped space, but components are Cartesian.
    virtual void setVelIC (FArrayBox&           a_velFAB,
                           const int            a_velComp,
                           const LevelGeometry& a_levGeo,
                           const DataIndex&     a_di) const;

    // Fills a FAB with the initial scalars
    virtual void setScalarIC (FArrayBox&           a_scalarFAB,
                              const int            a_scalarComp,
                              const LevelGeometry& a_levGeo,
                              const DataIndex&     a_di) const;

    // Fills a FAB with the background scalar
    virtual void setBackgroundScalar (FArrayBox&           a_scalarFAB,
                                      const int            a_scalarComp,
                                      const LevelGeometry& a_levGeo,
                                      const DataIndex&     a_di,
                                      const Real           a_time) const;


    // Utility functions -------------------------------------------------------
    // Most of these functions just call thier FAB versions and should not
    // need to be overridden.

    // Sets the velocity ICs over an entire level
    virtual void setVelIC (LevelData<FArrayBox>& a_vel,
                           const int             a_velComp,
                           const LevelGeometry&  a_levGeo) const;

    // Sets the scalar ICs over an entire level
    virtual void setScalarIC (LevelData<FArrayBox>& a_scalar,
                              const int             a_scalarComp,
                              const LevelGeometry&  a_levGeo) const;

    // Sets the background scalars over an entire level
    // If you plan to use a background stratification, overriding this function
    // is all you need to do.
    virtual void setBackgroundScalar (LevelData<FArrayBox>& a_scalar,
                                      const int             a_scalarComp,
                                      const Real            a_time,
                                      const LevelGeometry&  a_levGeo) const;

    // Adds the background scalar to an existing LevelData
    virtual void addBackgroundScalar (LevelData<FArrayBox>& a_scalar,
                                      const int             a_scalarComp,
                                      const Real            a_time,
                                      const LevelGeometry&  a_levGeo) const;

    // Subtracts the background scalar from an existing LevelData
    virtual void subtractBackgroundScalar (LevelData<FArrayBox>& a_scalar,
                                           const int             a_scalarComp,
                                           const Real            a_time,
                                           const LevelGeometry&  a_levGeo) const;

    // Computes the Brunt–Väisälä frequency on a single grid.
    // All ins and outs must be CC.
    // a_bFAB is the background buoyancy field.
    // a_dXidz needs SpaceDim comps...(dXi/dz, dNu/dz, dZeta/dz).
    // This function is a bit raw, but dXidz is expensive to compute and should
    // be cached by the user.
    virtual void computeNSq (FArrayBox&       a_NsqFAB,
                             const FArrayBox& a_bFAB,
                             const FArrayBox& a_dXidzFAB,
                             const Box&       a_destBox,
                             const RealVect&  a_dx) const;

    // Computes the Brunt–Väisälä frequency over an entire level.
    // a_Nsq must be CC.
    virtual void computeNSq (LevelData<FArrayBox>& a_Nsq,
                             const LevelGeometry&  a_levGeo,
                             const Real            a_time) const;

    // Simple accessor
    inline static bool useBackgroundScalar ();

    // This is in case the BC's have an effect on the timestep.
    // Pass in currently computed dt, along with the cfl and dx. If the effect
    // of the BCs requires a decreased timestep, then the newly reduced timestep
    // is returned. In the default case, this just returns a_dt back; however,
    // derived classes may actually have an effect.
    virtual void computeBoundaryDt (Real&                a_dt,
                                    const Real           a_cfl,
                                    const LevelGeometry& a_levGeo) const;


    // Sponge layer functions --------------------------------------------------

    // Sets the sponge ramp profile. This is a mapping [0,1] -> [0,1].
    inline static bool useSpongeLayer ();

    // Default is a linear increase as a_ratio approaches unity.
    inline virtual Real spongeLayerRamp (Real a_ratio) const;

    // Splits the domain into its sponge layers and the interior. The locations of
    // the splitting are given as face indices. If the domain is not split, these
    // indices will lie outside of the domain.
    static void computeSpongeRegions (Tuple<Box, 2>&       a_spongeBox,
                                      Tuple<int, 2>&       a_splitFaceIndex,
                                      Box&                 a_interior,
                                      const LevelGeometry& a_levGeo,
                                      const int            a_dir);

    // Sets a_srcTerm = (target - state) / (layer time scale).
    // To fill velocity sponge layer sources, just use a_comp's default value.
    virtual void fillSpongeLayerSrcTerm (LevelData<FArrayBox>&       a_srcTerm,
                                         const LevelData<FArrayBox>& a_state,
                                         const Real                  a_time,
                                         const Real                  a_dt,
                                         const LevelGeometry&        a_levGeo,
                                         const int                   a_comp = -1);

    // Sets the Cartesian-based target velocity for the sponge layer.
    // By default, this function throws an error.
    virtual void fillVelSpongeLayerTarget (FArrayBox&           a_target,
                                           const int            a_velComp,
                                           const int            a_spongeDir,
                                           const Side::LoHiSide a_spongeSide,
                                           const LevelGeometry& a_levGeo,
                                           const DataIndex&     a_di,
                                           const Real           a_time);

    // Sets the target values for the scalar sponge layer. If we are using a
    // background scalar, then this function set the perturbation to zero.
    // Otherwise, an error is thrown and this function will need to be overridden.
    virtual void fillScalarSpongeLayerTarget (FArrayBox&           a_target,
                                              const int            a_scalarComp,
                                              const int            a_spongeDir,
                                              const Side::LoHiSide a_spongeSide,
                                              const LevelGeometry& a_levGeo,
                                              const DataIndex&     a_di,
                                              const Real           a_time);

    // Velocity BCs ------------------------------------------------------------

    // Pre-projection velocity BC.
    virtual Tuple<BCMethodHolder, SpaceDim> uStarFuncBC (bool a_isViscous) const;

    // Sets ghosts needed to calculate the viscous source term nu.L[vel]
    virtual Tuple<BCMethodHolder, SpaceDim> viscousSourceFuncBC () const;

    // Used in single-compoent velocity TGA solves
    virtual BCMethodHolder viscousSolveFuncBC (int a_dir) const;

    // returns single-component BC for viscous refluxing solves
    virtual BCMethodHolder viscousRefluxBC (int a_dir) const;

    // Sets BCs on a generic viscous velocity field (Used by fill funcs)
    virtual Tuple<BCMethodHolder, SpaceDim> viscousVelFuncBC () const;

    // Sets BCs for tracing (Always inviscid. Used by fill funcs)
    virtual Tuple<BCMethodHolder, SpaceDim> tracingVelFuncBC () const;

    // Sets BCs on CC (u.Grad)u term
    virtual Tuple<BCMethodHolder, SpaceDim> uDelUFuncBC (bool a_isViscous) const;

    // Sets BC for FC advecting velocity
    virtual Tuple<BCMethodHolder, SpaceDim> advectingVelFuncBC (bool a_isViscous) const;

    // sets BCs on velocity when computing vorticity
    virtual Tuple<BCMethodHolder, SpaceDim> vortFuncBC (bool a_isViscous) const;

    // Sets BCs on FC velocity in the Riemann solver.
    // NOTE: This needs to set BCs on the velocity in the Cartesian basis.
    virtual Tuple<BCMethodHolder,SpaceDim> velRiemannBC (int a_velComp, bool a_isViscous) const;

    // Sets CC boundary slopes (undivided differences) on velocity in the tracing scheme.
    // NOTE: This needs to set BCs on the velocity in the Cartesian basis.
    virtual Tuple<BCMethodHolder,SpaceDim> velSlopeBC (int a_velComp, bool a_isViscous) const;


    // Scalar BCs --------------------------------------------------------------

    // Sets ghosts needed to calculate the diffusive term nu.L[scalar]
    virtual BCMethodHolder diffusiveSourceFuncBC () const;

    /// used in scalar TGA solves
    virtual BCMethodHolder diffusiveSolveFuncBC () const;

    /// used in implicit refluxing solves
    virtual BCMethodHolder scalarRefluxSolveBC (int a_scalarType) const;

    /// added by petermc, 5 nov 2007
    virtual BCMethodHolder scalarTraceFuncBC (int a_scalarType) const;

    // Sets BCs on FC scalars in the Riemann solver.
    virtual Tuple<BCMethodHolder,SpaceDim> scalarRiemannBC (int a_scalarType) const;

    // Sets CC boundary slopes (undivided differences) on scalars in the tracing scheme.
    virtual Tuple<BCMethodHolder,SpaceDim> scalarSlopeBC (int scalarType) const;


    // Freestream preservation BCs ---------------------------------------------

    /// volume discrepancy BC's
    virtual BCMethodHolder lambdaFuncBC () const;

    // added by petermc, 30 oct 2007
    virtual BCMethodHolder FreestreamCorrFuncBC () const;

    // added by petermc, 30 oct 2007
    virtual BCMethodHolder gradELambdaFuncBC () const;

    // Sets BCs on FC scalars in the Riemann solver.
    virtual Tuple<BCMethodHolder,SpaceDim> lambdaRiemannBC () const;

    // Sets CC boundary slopes (undivided differences) on lambda in the tracing scheme.
    virtual Tuple<BCMethodHolder,SpaceDim> lambdaSlopeBC () const;


    // Pressure BCs ------------------------------------------------------------

    // Used in levelMacProject solver
    virtual BCMethodHolder MacPressureFuncBC () const;

    // Used to calculate Grad[phi] in MAC projection
    virtual BCMethodHolder gradMacPressureFuncBC () const;

    // Used in LevelProject solver
    virtual BCMethodHolder LevelPressureFuncBC () const;

    // Used to calculate CCGrad[Pi] in LevelProject
    virtual BCMethodHolder gradPiFuncBC () const;

    // Used in sync projection solver
    virtual BCMethodHolder SyncProjFuncBC () const;

    // Used to calculate CCGrad[eSync] in sync projection
    virtual BCMethodHolder gradESyncFuncBC () const;


    // Miscellaneous BCs -------------------------------------------------------

    // Sets post-regrid smoothing BCs.
    virtual BCMethodHolder smoothingSolverBC () const;

    // BCs for the streamfunction solver.
    virtual BCMethodHolder streamSolverBC (int comp) const;


protected:
    // The basic BC functions --------------------------------------------------

    /// Sets physical BCs on velocities
    virtual BCMethodHolder basicVelFuncBC (int a_veldir, bool a_isViscous) const;

    /// Sets physical BCs on a generic passive scalar
    virtual BCMethodHolder basicScalarFuncBC () const;

    /// Sets physical BCs on pressures (used by the Poisson solvers)
    virtual BCMethodHolder basicPressureFuncBC (bool a_isHomogeneous) const;

    /// Sets physical BCs on pressures before taking gradients
    virtual BCMethodHolder basicGradPressureFuncBC () const;


    // Member variables and utility functions ----------------------------------

    // Have we read params from the input file yet?
    static bool s_staticMembersSet;

    static RealVect s_domLength;
    static bool s_useBackgroundScalar;

    static RealVect s_tidalU0;
    static Real s_tidalOmega;
    static bool s_doTidalFlow;

    static bool s_useSpongeLayer;
    static Real s_spongeWidth[CH_SPACEDIM][2];
    static Real s_spongeDtMult[CH_SPACEDIM][2];
};



// -----------------------------------------------------------------------------
// Simple accessor
// -----------------------------------------------------------------------------
inline bool PhysBCUtil::useBackgroundScalar ()
{
    return s_useBackgroundScalar;
}


// -----------------------------------------------------------------------------
// Simple accessor
// -----------------------------------------------------------------------------
inline bool PhysBCUtil::useSpongeLayer ()
{
    return s_useSpongeLayer;
}


// -----------------------------------------------------------------------------
// Default is a linear increase as a_ratio approaches unity.
// -----------------------------------------------------------------------------
inline Real PhysBCUtil::spongeLayerRamp (Real a_ratio) const
{
    return a_ratio;
}


#endif //!__PHYSBCUTIL_H__INCLUDED__

