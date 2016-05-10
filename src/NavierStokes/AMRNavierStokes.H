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
#ifndef _AMRNAVIERSTOKES_H_
#define _AMRNAVIERSTOKES_H_

#include "MappedAMRLevel.H"
#include "MappedCoarseAverage.H"
#include "MappedQuadCFInterp.H"
#include "MappedAMRPoissonOp.H"
#include "VelocityAMRPoissonOp.H"
#include "PhysBCUtil.H"
#include "BiCGStabSolver.H"
#include "MappedBaseLevelHeatSolver.H"

#include "LevelMACProjector.H"
#include "LevelCCProjector.H"

#include "CornerCopier.H"
#include "AdvectUtil.H"
#include "AlteredMetric.H"  // Perturbed metric for the BV solver
#include "StressMetric.H"   // Viscous stress tensor
#include "Debug.H"


#ifdef CH_USE_DOUBLE
#   define TIME_EPS 1.0e-10
#else
#   define TIME_EPS 1.0e-5
#endif

// Small parameter used in defining TGA constants for viscous solves
#define TGA_EPS 1.0e-12


class AMRNavierStokes: public MappedAMRLevel
{
    friend class AMRNavierStokesFactory;

public:
    // Default constructor
    AMRNavierStokes ();

    // unused...
    AMRNavierStokes (MappedAMRLevel* a_coarser_level_ptr,
                     const Box&      a_prob_domain,
                     int             a_level,
                     const IntVect&  a_ref_ratio);

    // unused...
    AMRNavierStokes (MappedAMRLevel*      a_coarser_level_ptr,
                     const ProblemDomain& a_prob_domain,
                     int                  a_level,
                     const IntVect&       a_ref_ratio);

    // Destructor
    virtual ~AMRNavierStokes ();

    // unused...
    virtual MappedAMRLevel* makeAMRLevel (MappedAMRLevel* a_coarser_level_ptr,
                                          const Box&      a_problem_domain,
                                          int             a_level,
                                          const IntVect&  a_ref_ratio) const;

    // unused...
    virtual MappedAMRLevel* makeAMRLevel (MappedAMRLevel*      a_coarser_level_ptr,
                                          const ProblemDomain& a_problem_domain,
                                          int                  a_level,
                                          const IntVect&       a_ref_ratio) const;

    // unused...
    virtual void define (MappedAMRLevel* a_coarse_level_ptr,
                         const Box&      a_problem_domain,
                         int             a_level,
                         const IntVect&  a_ref_ratio);

    // Full virtual constructor
    virtual void define (MappedAMRLevel*      a_coarse_level_ptr,
                         const ProblemDomain& a_problem_domain,
                         int                  a_level,
                         const IntVect&       a_ref_ratio);

    // advance by one timestep (returns maximum allowable timestep)
    virtual Real advance ();

    // things to do after a timestep
    virtual void postTimeStep ();

    // create tags
    virtual void tagCells (IntVectSet& a_tags);

    // create tags at initialization
    virtual void tagCellsInit (IntVectSet& a_tags);

    // perform any pre-regridding ops -- lBase is the finest unchanged level
    virtual void preRegrid (int a_lBase, const Vector<Vector<Box> >& a_newGrids);

    // regrid
    virtual void regrid (const Vector<Box>& a_new_grids);

    // perform any post-regridding ops -- lBase is the finest unchanged level
    virtual void postRegrid (int a_lBase);

    // initialize grids
    virtual void initialGrid (const Vector<Box>& a_new_grids);

    // Performs operations required after the grid has been defined but before
    // data initialization. This will also be called after readCheckpointLevel
    // during a restart procedure with argument a_restart set to true.
    virtual void postInitialGrid (const bool a_restart);

    // initialize data
    virtual void initialData ();

    // things do do after initialization
    virtual void postInitialize ();

    // compute dt
    virtual Real computeDt ();

    // compute dt with initial data
    virtual Real computeInitialDt ();


#ifdef CH_USE_HDF5
    // AMRNavierStokesIO.cpp ---------------------------------------------------

    // write checkpoint header
    virtual void writeCheckpointHeader (HDF5Handle& a_handle) const;

    // write checkpoint data for this level
    virtual void writeCheckpointLevel (HDF5Handle& a_handle) const;

    // read checkpoint header
    virtual void readCheckpointHeader (HDF5Handle& a_handle);

    // read checkpoint data for this level
    virtual void readCheckpointLevel (HDF5Handle& a_handle);

    // write plotfile header
    virtual void writePlotHeader (HDF5Handle& a_handle) const;

    // write plotfile data for this level
    virtual void writePlotLevel (HDF5Handle& a_handle) const;

    // Calculate number of components in plotfiles
    virtual int numPlotComps () const;

    // Stuffs the plotfile data
    virtual void getPlotData (LevelData<FArrayBox>& a_plot_data) const;
#endif


protected:
    // AMRNavierStokesAdvance.cpp ----------------------------------------------

    // Increments the coarse and fine registers
    virtual void updateVelFluxRegister (LevelData<FluxBox>& a_flux,
                                        const Real          a_scale);

    // Increments the coarse and fine registers
    virtual void updateLambdaFluxRegister (LevelData<FluxBox>& a_flux,
                                           const Real          a_scale);

    // Increments the coarse and fine registers
    virtual void updateScalarFluxRegister (LevelData<FluxBox>& a_flux,
                                           const Real          a_scale,
                                           const int           a_scalComp);


    // AMRNavierStokesAdvancePPM.cpp -------------------------------------------
    virtual void PPMTimeStep (const Real a_old_time,
                              const Real a_dt,
                              const bool a_updatePassiveScalars,
                              const bool a_doLevelProj);

    // Compute the half-time, VD corrected, projected, FC advecting velocity.
    virtual void computeAdvectingVelocities (LevelData<FluxBox>&   a_adv_vel,
                                             LevelData<FArrayBox>& a_old_vel,
                                             const Real            a_old_time,
                                             const Real            a_dt);

    // Calculate d(lambda)/dt = RHS centered at the time of a_old_vel + a_dt/2
    virtual void getNewLambda (LevelData<FArrayBox>&       a_rhs,
                               LevelData<FArrayBox>&       a_new_lambda,
                               LevelData<FArrayBox>&       a_old_lambda,
                               const LevelData<FArrayBox>& a_old_vel,
                               const LevelData<FluxBox>&   a_adv_vel,
                               const Real                  a_old_time,
                               const Real                  a_dt,
                               const Real                  a_FRscale);

    // Calculate d(scalar)/dt = RHS centered at the time of a_old_vel + a_dt/2
    // Since diffusive scalars are updated implicitly, this function actually
    // calculates a_new_scalar, then sets a_rhs = (a_new_scalar - a_old_scalar) / dt
    //
    // WARNING: At some point, I just assumed that a_comp = 0. If you plan to use
    // a scalar other than buoyancy, plan to debug!
    virtual void getNewScalar (LevelData<FArrayBox>&       a_rhs,
                               LevelData<FArrayBox>&       a_new_scalar,
                               LevelData<FArrayBox>&       a_old_scalar,
                               const LevelData<FArrayBox>& a_old_vel,
                               const LevelData<FluxBox>&   a_adv_vel,
                               const Real                  a_old_time,
                               const Real                  a_dt,
                               const Real                  a_FRscale,
                               const int                   a_comp);

    // Calculate d(velocity)/dt = RHS centered at the time of a_old_vel + a_dt/2
    virtual void getNewVelocity (LevelData<FArrayBox>&       a_rhs,
                                 LevelData<FArrayBox>&       a_new_vel,
                                 LevelData<FArrayBox>&       a_old_vel,
                                 const LevelData<FluxBox>&   a_adv_vel,
                                 const Real                  a_old_time,
                                 const Real                  a_dt,
                                 const Real                  a_FRscale);

    // Performs a CC (approximate) level projection of a_vel
    virtual void doCCProjection (LevelData<FArrayBox>& a_vel,
                                 const Real            a_new_time,
                                 const Real            a_dt,
                                 const bool            a_doProj = true);

    // Uses the CC velocity at a_oldTime to predict each component at all FCs at
    // time a_oldTime + 0.5 * a_dt. The exchange on the result is up to the caller.
    // This function also does not perform a projection or VD correction.
    // the input and output velocities are all in the mapped basis.
    virtual void predictVelocities (LevelData<FluxBox>&       a_predVel,
                                    LevelData<FArrayBox>&     a_oldVel,
                                    const LevelData<FluxBox>& a_advVel,
                                    const Real                a_oldTime,
                                    const Real                a_dt);


    // AMRNavierStokesAdvanceIGPPM.cpp -----------------------------------------

    virtual void PPMIGTimeStep (const Real a_old_time,
                                const Real a_dt,
                                const bool a_updatePassiveScalars,
                                const bool a_doLevelProj);

    // Implicitly handles the gravity forcing and projection.
    // The old* inputs should be the values at t^n.
    // The new* inputs should be the updated values from the TGA solver.
    virtual void doCCIGProjection (LevelData<FArrayBox>&       a_newVel,
                                   LevelData<FArrayBox>&       a_newB,
                                   const LevelData<FArrayBox>& a_oldVel,
                                   const LevelData<FArrayBox>& a_oldB,
                                   const LevelData<FluxBox>&   a_advVel,
                                   const Real                  a_newTime,
                                   const Real                  a_dt,
                                   const bool                  a_doProj = true);

    // AMRNavierStokesAdvanceRK3.cpp -------------------------------------------

    // Computes the source terms for a Runge-Kutta update.
    // This assumes only one scalar (not counting lambda).
    virtual void RK3TimeStep (const Real a_oldTime,
                              const Real a_dt,
                              const bool a_updatePassiveScalars,
                              const bool a_doLevelProj);

    // Computes:
    // Su = advective source + external force. Returned in the Cartesian basis.
    // Sb = advective source + external force.
    virtual void computeMOLSources (LevelData<FArrayBox>& a_Su,
                                    LevelData<FArrayBox>& a_Sb,
                                    LevelData<FArrayBox>& a_u,
                                    LevelData<FArrayBox>& a_b,
                                    const Real            a_stateTime,
                                    const int             a_stage);

    // Performs the TGA solve.
    // [I - h*L/2]q = [I + h*L/2]q + h*(beta*newS + zeta*oldS)
    // a_time is the time a_newSu is at.
    // a_stage can be 1, 2, or 3.
    virtual void updateState (LevelData<FArrayBox>&       a_u,
                              LevelData<FArrayBox>&       a_b,
                              Real&                       a_stateTime,
                              const LevelData<FArrayBox>& a_newSu,
                              const LevelData<FArrayBox>& a_newSb,
                              const LevelData<FArrayBox>& a_oldSu,
                              const LevelData<FArrayBox>& a_oldSb,
                              const int                   a_stage);

    // Erases the data in oldPtr and swaps the pointers.
    virtual void swapPointers (LevelData<FArrayBox>*& a_newPtr,
                               LevelData<FArrayBox>*& a_oldPtr);


    // AMRNavierStokesInit.cpp -------------------------------------------------

    // Set up data structures associated with this level
    void levelSetup (const DisjointBoxLayout& grids);

    // Solve for the internal wave speed due to stratification.
    void initializeInternalWaveSpeed ();

    // Define the viscous and diffusive solvers
    void defineViscousMGSolver (const DisjointBoxLayout& a_grids,
                                const DisjointBoxLayout* a_crseGridsPtr,
                                const IntVect&           a_refCrse);

    // Manages pressure initialization after init or regrid
    void initializeGlobalPressure ();

public:
    // AMRNavierStokesUtil.cpp -------------------------------------------------

    // Encapsulates the finer-level-pointer casting
    AMRNavierStokes* fineNSPtr () const;

    // Encapsulates the coarser-level pointer casting
    AMRNavierStokes* crseNSPtr () const;

    // Returns a pointer to the finest level
    AMRNavierStokes* finestNSPtr ();

    // Returns a pointer to the finest level (const version)
    const AMRNavierStokes* finestNSPtr () const;

    // Returns a pointer to the coarsest level (const version)
    const AMRNavierStokes* coarsestNSPtr () const;

    // Is this this finest level?
    bool finestLevel () const;

    // Is this an empty level?
    bool isEmpty () const;

protected:
    // Set whether this is the finest level or not
    void finestLevel (bool a_finest_level);

    // If s_set_bogus_values, assigns s_bogus_value to  vel, lambda, and scalars.
    void setAllBogus ();

public:
    //
    LevelData<FArrayBox>& newVel ();

    //
    const LevelData<FArrayBox>& newVel () const;

    //
    LevelData<FArrayBox>& oldVel ();

    //
    const LevelData<FArrayBox>& oldVel () const;

    //
    LevelData<FArrayBox>& newLambda ();

    //
    const LevelData<FArrayBox>& newLambda () const;

    //
    LevelData<FArrayBox>& oldLambda ();

    //
    const LevelData<FArrayBox>& oldLambda () const;

    //
    LevelData<FArrayBox>& newScal (const int a_comp);

    //
    const LevelData<FArrayBox>& newScal (const int a_comp) const;

    //
    LevelData<FArrayBox>& oldScal (const int a_comp);

    //
    const LevelData<FArrayBox>& oldScal (const int a_comp) const;

protected:
    // Collects the new velocities from each level at and above a_startNSPtr.
    // Lower levels are set to NULL. If a_startPtr is NULL, all data is retrieved.
    virtual Vector<LevelData<FArrayBox>*> gatherNewVelStartingWith (AMRNavierStokes* a_startNSPtr);

    // Collects the sync pressures from each level at and above this level-1.
    // Lower levels are set to NULL.
    virtual Vector<LevelData<FArrayBox>*> gatherSyncPressure ();

    // Sets m_syncPressureState on this and all higher levels.
    virtual void setSyncPressureStates (int a_state);

    // Collects eLambda from each level at and above this level-1.
    // Lower levels are set to NULL.
    virtual Vector<LevelData<FArrayBox>*> gatherELambda ();

    // Sets m_eLambdaState on this and all higher levels.
    virtual void setELambdaStates (int a_state);

    // Sets m_gradELambdaState on this and all higher levels.
    virtual void setGradELambdaStates (int a_state);

    // Move new-time state into old-time state -- also advances time by dt
    void swapOldAndNewStates ();

    // Exchange old and new states, resets time to a_time
    void resetStates (const Real a_time);


    // AMRNavierStokesFill.cpp -------------------------------------------------

    // Sets all ghosts on velocity. (Does not exchange.)
    void setGhostsVelocity (LevelData<FArrayBox>& a_vel,
                            VelBCHolder&          a_velBC,
                            const Real            a_time) const;

    // Sets all ghosts on lambda. (Does not exchange.)
    void setGhostsLambda (LevelData<FArrayBox>& a_lambda,
                          BCMethodHolder&       a_lambdaBC,
                          const Real            a_time) const;

    // Sets all ghosts on a scalar. (Does not exchange.)
    void setGhostsScalar (LevelData<FArrayBox>& a_scal,
                          BCMethodHolder&       a_scalBC,
                          const Real            a_time,
                          const int             a_comp,
                          const bool            a_onlyAtValid = false) const;

    // Returns velocity at time a_time.
    // This does a linear interpolation in time between old and new time velocities,
    // but does not fill ghosts or interpolate in space from a coarser level.
    void velocity (LevelData<FArrayBox>& a_vel,
                   const Real            a_time) const;

    // Returns lambda at time a_time.
    // This does a linear interpolation in time between old and new time lambdas,
    // but does not fill ghosts or interpolate in space from a coarser level.
    void lambda (LevelData<FArrayBox>& a_lambda,
                 const Real            a_time) const;

    // Returns scalar at time a_time.
    // This does a linear interpolation in time between old and new time scalars,
    // but does not fill ghosts or interpolate in space from a coarser level.
    void scalar (LevelData<FArrayBox>& a_scal,
                 const Real            a_time,
                 const int             a_comp) const;

    // Fill grown velocity holder using piecewise-constant interp
    void fillVelocity (LevelData<FArrayBox>& a_vel,
                       const Real            a_time) const;

    // Resize and fill temporary lambda for advection tracing computation.
    // This fills both interior and boundary cells.
    void fillLambda (LevelData<FArrayBox>& a_lambda,
                     const Real            a_time) const;

    // Resize and fill temporary scalar field for advection tracing computation.
    void fillScalars (LevelData<FArrayBox>& a_scal,
                      const Real            a_time,
                      const int             a_comp,
                      const bool            a_addBackground = false) const;

    // Computes the viscous source term, nu.L[u].
    // All inputs and outputs are in the Cartesian basis.
    // This assumes the ghosts of a_cartVel have already been set.
    void fillViscousSource (LevelData<FArrayBox>&       a_viscSource,
                            const LevelData<FArrayBox>& a_cartVel,
                            const Real                  a_time);

    // Returns the gravitational source term in mapped coordinates.
    // This is -g*zhat in Cartesian. a_gravSource needs SpaceDim comps.
    void fillGravSource (LevelData<FArrayBox>& a_gravSource,
                         const Real            a_time,
                         const bool            a_addBackground = false) const;

    // Returns the tidal forcing source term.
    // This is U0*omega*cos(omega*t) in the Cartesian basis.
    // a_tidalSource needs SpaceDim comps.
    void fillTidalSource (LevelData<FArrayBox>& a_tidalSource,
                          const Real            a_time,
                          const Real            a_dt) const;

    // Returns the internal wave speed
    void fillInternalWaveSpeed (LevelData<FArrayBox>& a_c0i) const;

    // Returns the vertical structure function.
    void fillVerticalStructure (LevelData<FArrayBox>& a_phi0) const;


    // AMRNavierStokesRegrid.cpp -----------------------------------------------

    // Assigns procs to the new set of boxes.
    DisjointBoxLayout loadBalance (const Vector<Box>& a_grids);

    // Constructs the RHS of the smooth operations to be performed after regridding.
    void setupPostRegridSmoothing (int a_lBase);

    // Smooths composite fields after regridding
    void doPostRegridSmoothing (int lbase);

    // Defines operator factory used in smoothing during regridding
    void defineRegridAMROp (MappedAMRPoissonOpFactory&       a_factory,
                            const int&                       a_lBase,
                            const Real                       a_viscCoeff);


    // AMRNavierStokesDerive.cpp -----------------------------------------------

    // Compute Grad[macPressure]. This is really just a wrapper around
    // m_macProjector.computeGrad to help keep the code neat. This also avoids
    // using a possibly undefined projector in the compressible case.
    void gradMACPressure (LevelData<FluxBox>& a_gradPressure,
                          const Real          a_scale);

    // Compute Grad[ccPressure]. This is really just a wrapper around
    // m_ccProjector.computeGrad to help keep the code neat. This also avoids
    // using a possibly undefined projector in the compressible case.
    void gradCCPressure (LevelData<FArrayBox>& a_gradPressure,
                         const Real            a_scale);

    // Compute Laplacian(velocity), set physical BCs, and exchange
    void computeLapVel (LevelData<FArrayBox>&       a_lapVel,
                        const LevelData<FArrayBox>& a_vel,
                        const LevelData<FArrayBox>* a_crseVelPtr,
                        const CFRegion*             a_CFRegionPtr = NULL,
                        const Real                  a_time = BOGUS_TIME);

    // Compute Laplacian(scalars)
    void computeLapScal (LevelData<FArrayBox>&       a_lapScal,
                         LevelData<FArrayBox>&       a_scal,
                         const LevelData<FArrayBox>* a_crseScalPtr,
                         const CFRegion*             a_CFRegionPtr = NULL);

    // compute vorticity
    // NOTE: This returns two-form components, not a vector!
    void computeVorticity (LevelData<FArrayBox>& a_vorticity) const;

    // Compute the stream function on this level.
    void computeStreamFunction (LevelData<FArrayBox>& a_stream) const;

    // Computes the energy integral
    Real totalEnergy () const;


    // AMRNavierStokesSync.cpp -------------------------------------------------

    // Do implicit scalar refluxing
    void doImplicitScalarReflux ();

    // Do implicit momentum refluxing (used in sync ops)
    void doImplicitMomentumReflux (const Vector<LevelData<FArrayBox>*>& a_compVel);

    // Compute the volume discrepancy correction.
    // Set a_init to true during initialization VD solves and false for sync solves.
    void computeVDCorrection (const bool a_init);

    // Diagnostic functions
    void syncBottomDiagnostics (const Vector<LevelData<FArrayBox>*>& a_compVel);
    void syncSingleGridDiagnostics ();
    void syncTermDiagnostics ();


    // AMRNavierStokesBasics.cpp -----------------------------------------------

    // Read parameters from parmParse database
    void readParameters ();

protected:

    // Sets physical IBCs
    PhysBCUtil* m_physBCPtr;

    // Contains this level's mapping info.
    LevelGeometry* m_levGeoPtr;

    // velocity vector at old time
    LevelData<FArrayBox>* m_vel_old_ptr;

    // velocity vector at new time
    LevelData<FArrayBox>* m_vel_new_ptr;

    // freestream preservation qty at old time
    LevelData<FArrayBox>* m_lambda_old_ptr;

    // freestream preservation qty at new time
    LevelData<FArrayBox>* m_lambda_new_ptr;

    // scalars at old time
    Vector<LevelData<FArrayBox>*> m_scal_old;

    // scalars at new time
    Vector<LevelData<FArrayBox>*> m_scal_new;

    // number of scalars (not including lambda)
    static int s_num_scal_comps;

    // names of velocity components
    static const char* s_vel_names[CH_SPACEDIM];

    // names of scalars
    static Vector<std::string> s_scal_names;

    // The level MAC projector's pressure.
    // NOTE: This is not the level's dynamical pressure -- the CC projector
    // produces that. We only keep this around so we don't have to project
    // the advecting velocity every time it is needed.
    LevelData<FArrayBox> m_macPressure;

    // The level CC projector's pressure.
    // The complete dynamical pressure is this plus the sync pressure.
    LevelData<FArrayBox> m_ccPressure;

    // This is just for debugging purposes.
    struct CCPressureState {
        enum {
            UNDEFINED = -2, // Not allocated/defined.
            BOGUS     = -1, // Set to s_bogus_value for debugging purposes.
            ZERO      =  0, // Set to zero prior to pressure initialization.
            VALID     =  1  // Contains valid pressure data.
        };
    };
    int m_ccPressureState;

    // The AMR sync pressure correction.
    // The complete dynamical pressure is this plus the CC pressure.
    LevelData<FArrayBox> m_syncPressure;

    // This is just for debugging purposes.
    struct SyncPressureState {
        enum {
            UNDEFINED = -2, // Not allocated/defined.
            BOGUS     = -1, // Set to s_bogus_value for debugging purposes.
            ZERO      =  0, // Set to zero prior to projecting.
            INIT      =  1, // Contains remnants of the postInitialize() initialization projection.
            REGRID    =  2, // Contains remnants of the postRegrid() intitialization projections.
            SYNC      =  3  // Contains valid pressure data.
        };
    };
    int m_syncPressureState;

    // The volume discrepancy correction
    LevelData<FArrayBox> m_eLambda;
    struct ELambdaState {
        enum {
            UNDEFINED = -2,
            BOGUS     = -1,
            ZERO      =  0,
            VALID     =  1
        };
    };
    int m_eLambdaState;

    LevelData<FluxBox> m_gradELambda;
    struct GradELambdaState {
        enum {
            UNDEFINED = -2,
            BOGUS     = -1,
            ZERO      =  0,
            VALID     =  1
        };
    };
    int m_gradELambdaState;

    // Have we initialized static variables from parmParse?
    static bool s_ppInit;

    // Size of physical domain
    static RealVect s_domLength;

    // debugging option -- initially set arrays to s_bogus_value
    static const bool s_set_bogus_values = DEBUG_BOOL;

    // debugging option -- value to which arrays are initially set
    static Real s_bogus_value;

    // do we specify initial grids in a gridFile?
    static bool s_specifyInitialGrids;

    // gridfile where initial grids are spec'd, if relevant
    static string s_initialGridFile;


    // amount by which to grow tagged regions before calling meshrefine
    static int s_tags_grow;

    // Tag on vorticity.
    // Criterion: |vort| > s_magvort_tag_quota * max(vort)
    // Set to zero to turn off.
    static Real s_magvort_tag_quota;

    // New vorticity tagging strategy.
    // Criterion: |vort[dir]*dA[dir]| > s_vort_tag_tol[dir]
    //            where dir={2} in 2D and dir={0,1,2} in 3D.
    // Set all comps to zero vector to turn off.
    static Tuple<Real,3> s_vort_tag_tol;

    // Tag on velocity differences.
    // Criterion: |UD velocity[dir]| > s_vel_tag_tol
    //            where UD is the max undivided difference operator.
    // Set to zero to turn off.
    static Real s_vel_tag_tol;

    // Tag on buoyancy differences.
    // Criterion: |UD buoyancy| > s_buoyancy_tag_tol
    //            where UD is the max undivided difference operator.
    // Set to zero to turn off.
    static Real s_buoyancy_tag_tol;

    // Tag on pressure differences.
    // Criterion: |UD pressure| > s_pressure_tag_tol
    //            where UD is the max undivided difference operator.
    // Set to zero to turn off.
    static Real s_pressure_tag_tol;

    // Vertically extrudes all tags when set to true.
    // This should guarantee that each *level* will vertically span the domain.
    // This DOES NOT guarantee the same for each *patch* itself.
    static bool s_vert_extrude_tags;

    // Is this fluid incompressible? (false turns off the projector)
    static bool s_isIncompressible;

    // Number of initial (approximate) velocity projections
    static int s_initial_projection_iters;

    // Number of sync (approximate) velocity projections
    static int s_sync_projection_iters;

    // Number of passes to initialize pressures after regridding
    static int s_initial_pressure_iters;

    // Number of level (approximate) velocity projections
    static int s_level_projection_iters;

    // Do antidiffusive smoothing after regridding?
    static bool s_smooth_after_regrid;

    // Antidiffusive smoothing coefficient (multiplies nu*dtCrse)
    static Real s_regrid_smoothing_coeff;

    // Internal indicator whether post-regrid smoothing setup has been done
    bool m_regrid_smoothing_done;

    /// How should we perform refluxing?
    static bool s_advective_momentum_reflux;
    static bool s_diffusive_momentum_reflux;
    static bool s_implicit_momentum_reflux;
    static bool s_advective_scalar_reflux;
    static bool s_diffusive_scalar_reflux;
    static bool s_implicit_scalar_reflux;
    static bool s_advective_lambda_reflux;

    /// apply freestream preservation correction?
    static bool s_applyFreestreamCorrection;
    static Real s_etaLambda;

    // Sets the form of the nonlinear velocity forcing term.
    // Div[u u], u.Grad[u], etc...
    static int s_nonlinearDifferencingForm;

    // With gravity off, a zero initial velocity better stay zero!
    static int s_gravityMethod;
    static Real s_gravityTheta;

    // The tidal parameters (set to zero for no tide)
    static Real s_tidalOmega;
    static RealVect s_tidalU0;

    // The Coriolis parameter
    static Real s_coriolisF;

    // -- I/O options --

    // Should we write to standard out?
    static bool s_write_stdout;

    // include divergence in plotfile?
    static bool s_write_divergence;

    // include lambda in plotfile?
    static bool s_write_lambda;

    // include freestream preservation correcion in plotfile?
    static bool s_write_grad_eLambda;

    // include pressure in plotfile?
    static bool s_write_pressure;

    // include vorticity in plotfile?
    static bool s_write_vorticity;

    // include streamfunction in plotfile?
    static bool s_write_streamfunction;

    // include scalars in plotfile?
    static bool s_write_scalars;

    // include scalars-background in plotfile?
    static bool s_write_scalarsMinusBackground;

    // include processor distribution in plotfile?
    static bool s_write_proc_ids;

    // include level assignments in plotfile?
    static bool s_write_level_ids;

    // include geometric displacements from Cartesian in plotfile?
    static bool s_write_displacement;

    // write out grids and processor assignments
    static bool s_write_grids;

    // Should we write the Jacobian, metric, etc?
    static bool s_write_geometry;

    // The total (composite) energy
    static Real s_totalEnergy;

    // factor to shrink initial timestep
    static Real s_init_shrink;

    // maximum allowable timestep
    static Real s_max_dt;
    static Real s_max_dt_grow;

    // Fixed timestep (if desired)
    static Real s_prescribedDt;

    static int s_updateScheme;              // Finite volume method or RK3?

    // Velocity advection scheme parameters
    static int  s_normalPredOrderVel;       // CTU, PLM, or PPM
    static bool s_useFourthOrderSlopesVel;  // Helps remove the effects of diffusive errors
    static bool s_useLimitingVel;           // Use a van-Leer limiter?
    static bool s_useHighOrderLimiterVel;   // Only used in PPM method
    static bool s_useUpwindingVel;

    // Scalar advection scheme parameters
    static int  s_normalPredOrderScal;      // CTU, PLM, or PPM
    static bool s_useFourthOrderSlopesScal; // Helps remove the effects of diffusive errors
    static bool s_useLimitingScal;          // Use a van-Leer limiter?
    static bool s_useHighOrderLimiterScal;  // Only used in PPM method
    static bool s_useUpwindingScal;

    // How should we limit dt?
    static bool s_limitDtViaViscosity;
    static bool s_limitDtViaDiffusion;
    static bool s_limitDtViaPressureGradient;
    static bool s_limitDtViaInternalWaveSpeed;

    // Should MappedCoarseAverage consider varying cell sizes?
    static const bool s_considerCellSizes;

    // true if this is the finest extant level
    bool m_finest_level;

    // true if this level is an empty level (no grids)
    bool m_is_empty;

    // cfl number
    Real m_cfl;

    // Exchange copiers for hyperbolic tracing scheme
    IntVect m_tracingGhosts;
    Copier m_tracingExCopier;
    CornerCopier m_tracingExCornerCopier;

    // Exchange copier for any field on this level with one ghost.
    Copier m_oneGhostExCopier;

    // averaging from fine to coarse level for velocity
    MappedCoarseAverage m_coarse_average;

    // averaging from fine to coarse level for scalars
    MappedCoarseAverage m_coarse_average_scal;

    // The MAC level projector
    LevelMACProjector m_macProjector;

    // The CC level projector
    LevelCCProjector m_ccProjector;

    // The metric used by the implicit BV projector
    AlteredMetric m_alteredMetric;

#   ifdef USE_STRESSMETRIC
        // Stress tensor
        StressMetric m_stressMetric;
#   endif

    // Advection utilities
    MappedAdvectionUtil m_advectUtilVel;
    MappedAdvectionUtil m_advectUtilScal;
    MappedAdvectionUtil m_advectUtilLambda;

    // momentum flux register
    MappedLevelFluxRegister m_vel_flux_reg;

    // freestream preservation qty flux register
    MappedLevelFluxRegister m_lambda_flux_reg;

    // conserved scalar flux register
    Vector<MappedLevelFluxRegister*> m_scal_fluxreg_ptrs;

    // does quadratic coarse-fine interpolation for velocities
    MappedQuadCFInterp m_velCFInterp;

    // MappedAMRPoissonOp for velocities
    VelocityAMRPoissonOp m_velocityAMRPoissonOp;

    // MappedAMRPoissonOp for advected scalars (for diffusion)
    MappedAMRPoissonOp m_scalarsAMRPoissonOp;

    // Bottom solver for viscous/diffusive solves.
    BiCGStabSolver<LevelData<FArrayBox> > m_viscousBottomSolver;

    // levelsolver for viscous solves for velocity
    Tuple< RefCountedPtr< MappedAMRMultiGrid<LevelData<FArrayBox> > >, SpaceDim> m_viscAMRMGPtrs;

    // levelsolver for advectDiffuseScalar solves
    Vector< RefCountedPtr< MappedAMRMultiGrid<LevelData<FArrayBox> > > > m_diffAMRMGPtrs;

    // TGA solver for computeUstar():  LevelTGA has no default constructor
    Tuple< RefCountedPtr<MappedBaseLevelHeatSolver>, SpaceDim> m_viscSolverPtrs;

    // TGA solver for advectDiffuseScalar():  LevelTGA has no default constructor
    Vector< RefCountedPtr<MappedBaseLevelHeatSolver> > m_diffSolverPtrs;

    // Selects HeatSolverSchemes for the viscous and diffusive solves.
    static int s_viscSolverScheme;
    static int s_diffSolverScheme;

    // Viscosity coefficient
    static Real s_nu;

    // Diffusion coefficients for scalar advection-diffusion
    static Vector<Real> s_scal_coeffs;

    // Solver parameters
    static Real s_AMRMG_eps;                        // Solver tolerance
    static int s_AMRMG_num_smooth_down;             // MG pre smoothing iters
    static int s_AMRMG_num_smooth_up;               // MG post smoothing iters
    static int s_AMRMG_num_smooth_bottom;           // MG bottom smoothing iters
    static int s_AMRMG_num_smooth_precond;          // MG preconditioner smoothing iters
    static int s_AMRMG_numMG;                       // 1=V-cycle, 2=W-cycle, etc...
    static int s_AMRMG_imin;                        // Min number of V-cycles
    static int s_AMRMG_imax;                        // Max number of V-cycles
    static Real s_AMRMG_hang;                       // Defaults to 1e-15
    static Real s_AMRMG_normThresh;                 // Defaults to 1e-30
    static int s_AMRMG_maxDepth;                    // Max MG depth (-1 for as deep as possible)
    static int s_AMRMG_verbosity;                   //
    static int s_AMRMG_relaxMode;
    static int s_AMRMG_precondMode;

    static Real s_bottom_eps;                       // Solver tolerance
    static Real s_bottom_reps;                      // Solver relative tolerance
    static int s_bottom_imax;                       // Max BiCGStab iterations
    static int s_bottom_numRestarts;                // Max BiCGStab restarts
    static Real s_bottom_hang;                      //
    static Real s_bottom_small;                     //
    static int s_bottom_normType;                   //
    static int s_bottom_verbosity;                  //

    // Viscous solver overrides
    static Real s_viscous_AMRMG_eps;                        // Solver tolerance
    static int s_viscous_AMRMG_num_smooth_down;             // MG pre smoothing iters
    static int s_viscous_AMRMG_num_smooth_up;               // MG post smoothing iters
    static int s_viscous_AMRMG_num_smooth_bottom;           // MG bottom smoothing iters
    static int s_viscous_AMRMG_num_smooth_precond;          // MG preconditioner smoothing iters
    static int s_viscous_AMRMG_numMG;                       // 1=V-cycle, 2=W-cycle, etc...
    static int s_viscous_AMRMG_imin;                        // Min number of V-cycles
    static int s_viscous_AMRMG_imax;                        // Max number of V-cycles
    static Real s_viscous_AMRMG_hang;                       // Defaults to 1e-15
    static Real s_viscous_AMRMG_normThresh;                 // Defaults to 1e-30
    static int s_viscous_AMRMG_maxDepth;                    // Max MG depth (-1 for as deep as possible)
    static int s_viscous_AMRMG_relaxMode;
    static int s_viscous_AMRMG_precondMode;
    static int s_viscous_AMRMG_verbosity;                   //

    static Real s_viscous_bottom_eps;                       // Solver tolerance
    static Real s_viscous_bottom_reps;                      // Solver relative tolerance
    static int s_viscous_bottom_imax;                       // Max BiCGStab iterations
    static int s_viscous_bottom_numRestarts;                // Max BiCGStab restarts
    static Real s_viscous_bottom_hang;                      //
    static Real s_viscous_bottom_small;                     //
    static int s_viscous_bottom_normType;                   //
    static int s_viscous_bottom_verbosity;                  //

    // The current step number (of the coarsest level)
    static int s_step_number;

    // Cache info for updateDt.
    LevelData<FArrayBox>* m_c0iPtr;

    // This is used to project fields.
    LevelData<FArrayBox>* m_phi0Ptr;
};


#endif
