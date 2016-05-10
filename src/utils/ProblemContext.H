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
#ifndef __Parameters_H__INCLUDED__
#define __Parameters_H__INCLUDED__

#include "IntVect.H"
#include "ProblemDomain.H"
#include "RealVect.H"
#include "Tuple.H"
#include <string>
class PhysBCUtil;
class GeoSourceInterface;


// -----------------------------------------------------------------------------
// Reads the problem parameters from input.* and creates a packages the data
// into a context class.
// -----------------------------------------------------------------------------
class ProblemContext
{
private:
    // This is a singleton class -- private constructor, static pointer!
    ProblemContext ();
    static ProblemContext* s_singletonPtr;

public:
    // Destructor
    ~ProblemContext();

    // This returns the single ProblemContext object.
    static const ProblemContext* getInstance ();

    // This deletes the static pointer in hopes valgrind will take notice.
    static void freeMemory ();

private:
    // The amr.* parameters
    void readAMR ();
    static bool s_AMRParamsRead;

public:
    IntVect nx;
    IntVect nx_offset;
    bool isPeriodic[CH_SPACEDIM];
    ProblemDomain domain;
    RealVect domainLength;
    Real max_dt_grow;
    Real fixed_dt;
    Real stopTime;
    Real cfl;
    bool useSubcycling;
    int maxsteps;
    int verbosity;
    int max_level;
    int numlevels;
    Vector<int> regrid_intervals;
    Vector<IntVect> refRatios;
    int block_factor;
    int bufferSize;
    Real fill_ratio;
    Vector<int> splitDirs;
    IntVect maxGridSize;
    IntVect maxBaseGridSize;
    bool isRestart;
    std::string restart_file;
    bool hasPredefinedGrids;
    Vector<Vector<Box> > predefinedGrids;

    int tags_grow;
    Real magvort_tag_quota;
    Tuple<Real,3> vort_tag_tol;
    Real vel_tag_tol;
    Real buoyancy_tag_tol;
    Real pressure_tag_tol;
    bool vert_extrude_tags;

    bool write_stdout;
    Real init_dt_multiplier;
    Real max_dt;
    bool limitDtViaViscosity;
    bool limitDtViaDiffusion;
    bool limitDtViaPressureGradient;
    bool limitDtViaInternalWaveSpeed;
    Real bogus_value;
    bool smooth_after_regrid;
    Real regrid_smoothing_coeff;
    bool advective_momentum_reflux;
    bool diffusive_momentum_reflux;
    bool implicit_momentum_reflux;
    bool advective_scalar_reflux;
    bool diffusive_scalar_reflux;
    bool implicit_scalar_reflux;
    bool advective_lambda_reflux;

    struct NonlinearDifferencingForm {
        enum {
            NONE         = -1,
            CONSERVATIVE =  0,
            ADVECTIVE    =  1,
            _NUM_NONLINEAR_DIFFERENCING_FORMS
        };
    };
    int nonlinearDifferencingForm;

    struct GravityMethod {
        enum {
            NONE     = 0,
            EXPLICIT = 1,
            IMPLICIT = 2,
            _NUM_GRAVITYMETHODS
        };
    };
    int gravityMethod;
    Real gravityTheta;
    Real coriolisF;

    struct HeatSolverScheme {
        enum {
            BACKWARD_EULER = 0,
            CRANK_NICOLSON = 1,
            TGA            = 2,
            EXPLICIT       = 3,
            _NUM_HEATSOLVERS
        };
    };
    int viscSolverScheme;
    int diffSolverScheme;

    int num_scal_comps;
    Vector<Real> scal_coeffs;

    Real nu;

private:
    // The geometry.* parameters
    void readGeometry ();
    static bool s_geometryParamsRead;

public:
    // Use this to request a new GeoSourceInterface cast from the
    // appropriate child.
    GeoSourceInterface* newGeoSourceInterface () const;

    struct CoordMap {
        enum {
            UNDEFINED        =-1,
            CARTESIAN        = 0,
            TWISTED          = 1,
            BEAMGENERATOR    = 2,
            //VERTBDRYSTRETCH  = 3, // Removed - no longer available.
            CYLINDRICAL      = 4,
            LEDGE            = 5,
            // MASSBAY          = 6, // Removed - no longer available.
            // NEWBEAMGENERATOR = 7,
            DEMMAP           = 8,
            _NUM_COORD_MAPS
        };
    };
    int coordMap;

    RealVect pert;

    // Specific to LedgeMap
    int ledgeMapTransitionOrder;
    Real ledgeMapHl;
    Real ledgeMapHr;
    Real ledgeMapXl;
    Real ledgeMapXr;

    // Specific to BeamGeneratorMap
    Real beamGenMapAlpha;

    // Specific to DEMMap
    std::string demFile;
    int interpOrder;

private:
    // The plot.* parameters
    void readPlot ();
    static bool s_plotParamsRead;

public:
    int plot_interval;
    Real plot_period;
    std::string plot_prefix;
    int checkpoint_interval;
    std::string check_prefix;

    bool write_divergence;
    bool write_lambda;
    bool write_grad_eLambda;
    bool write_pressure;
    bool write_vorticity;
    bool write_streamfunction;
    bool write_scalars;
    bool write_scalarsMinusBackground;
    bool write_proc_ids;
    bool write_level_ids;
    bool write_grids;
    bool write_displacement;
    bool write_geometry;

private:
    // The ibc.* parameters.
    void readIBC ();
    static bool s_IBCParamsRead;

public:
    // Use this to request a new (undefined) PhysBCUtil cast from the
    // appropriate child.
    PhysBCUtil* newPhysBCUtil () const;

    struct ProblemType {
        enum {
            ADVECTION_TEST  = 0,
            LOCK_EXCHANGE   = 1,
            BEAM_GENERATION = 2,
            INTERNAL_WAVE   = 3,
            TAYLOR_GREEN    = 4,
            VORTEX_STREET   = 5,
            HORIZ_CONV      = 6,
            SOLITARYWAVE    = 7,
            DJL             = 8,
            _NUM_PROBLEM_TYPES
        };
    };
    int problem;

    bool useBackgroundScalar;
    bool useSpongeLayer;
    Real spongeWidth[CH_SPACEDIM][2];
    Real spongeDtMult[CH_SPACEDIM][2];
    Real tidalOmega;
    RealVect tidalU0;
    RealVect inflowVel;
    bool specifyInitialGrids;
    std::string initialGridFile;

    Real solitaryWave_H2;
    Real solitaryWave_xcenter;
    Real solitaryWave_amp;
    Real solitaryWave_rho0;
    Real solitaryWave_drho;
    Real solitaryWave_delta0;

private:
    // The advection.* parameters
    void readAdvection ();
    static bool s_advectionParamsRead;

public:
    struct UpdateScheme {
        enum {
            FiniteVolume = 0,
            RK3 = 1,
            _NUM_UPDATE_SCHEMES
        };
    };
    int updateScheme;

    int normalPredOrderVel;
    bool useFourthOrderSlopesVel;
    bool useLimitingVel;
    bool useHighOrderLimiterVel;
    bool useUpwindingVel;

    int normalPredOrderScal;
    bool useFourthOrderSlopesScal;
    bool useLimitingScal;
    bool useHighOrderLimiterScal;
    bool useUpwindingScal;

private:
    // The solver parameters. This includes the projector, viscous solver, etc.
    void readSolver ();
    static bool s_solverParamsRead;

public:
    struct PrecondMode {
        enum {
            None = -1,
            DiagRelax = 0,
            DiagLineRelax = 1,
            _NUM_PRECOND_MODES
        };
    };

    struct RelaxMode {
        enum {
            NORELAX     = -1,
            JACOBI      =  0,
            LEVEL_GSRB  =  1,
            LOOSE_GSRB  =  2,
            LINE_GSRB   =  3,
            NUM_RELAX_MODES
        };
    };

    Real AMRMG_eps;                        // Solver tolerance
    int AMRMG_num_smooth_down;             // MG pre smoothing iters
    int AMRMG_num_smooth_up;               // MG post smoothing iters
    int AMRMG_num_smooth_bottom;           // MG bottom smoothing iters
    int AMRMG_num_smooth_precond;          // MG preconditioner smoothing iters
    int AMRMG_numMG;                       // 1=V-cycle, 2=W-cycle, etc...
    int AMRMG_imin;                        // Min number of V-cycles
    int AMRMG_imax;                        // Max number of V-cycles
    Real AMRMG_hang;                       // Defaults to 1e-15
    Real AMRMG_normThresh;                 // Defaults to 1e-30
    int AMRMG_maxDepth;                    // Max MG depth (-1 for as deep as possible)
    int AMRMG_verbosity;                   //
    int AMRMG_relaxMode;
    int AMRMG_precondMode;

    Real bottom_eps;                       // Solver tolerance
    Real bottom_reps;                      // Solver relative tolerance
    int bottom_imax;                       // Max BiCGStab iterations
    int bottom_numRestarts;                // Max BiCGStab restarts
    Real bottom_hang;                      //
    Real bottom_small;                     //
    int bottom_normType;                   //
    int bottom_verbosity;                  //

    Real viscous_AMRMG_eps;                        // Solver tolerance
    int viscous_AMRMG_num_smooth_down;             // MG pre smoothing iters
    int viscous_AMRMG_num_smooth_up;               // MG post smoothing iters
    int viscous_AMRMG_num_smooth_bottom;           // MG bottom smoothing iters
    int viscous_AMRMG_num_smooth_precond;          // MG preconditioner smoothing iters
    int viscous_AMRMG_numMG;                       // 1=V-cycle, 2=W-cycle, etc...
    int viscous_AMRMG_imin;                        // Min number of V-cycles
    int viscous_AMRMG_imax;                        // Max number of V-cycles
    Real viscous_AMRMG_hang;                       // Defaults to 1e-15
    Real viscous_AMRMG_normThresh;                 // Defaults to 1e-30
    int viscous_AMRMG_maxDepth;                    // Max MG depth (-1 for as deep as possible)
    int viscous_AMRMG_verbosity;                   //
    int viscous_AMRMG_relaxMode;
    int viscous_AMRMG_precondMode;

    Real viscous_bottom_eps;                       // Solver tolerance
    Real viscous_bottom_reps;                      // Solver relative tolerance
    int viscous_bottom_imax;                       // Max BiCGStab iterations
    int viscous_bottom_numRestarts;                // Max BiCGStab restarts
    Real viscous_bottom_hang;                      //
    Real viscous_bottom_small;                     //
    int viscous_bottom_normType;                   //
    int viscous_bottom_verbosity;                  //

private:
    // The projection.* parameters
    void readProjection ();
    static bool s_projectionParamsRead;

public:
    bool isIncompressible;
    int initial_projection_iters;
    int initial_pressure_iters;
    int level_projection_iters;
    bool doSyncProjection;
    int sync_projection_iters;
    bool applyVDCorrection;
    Real etaLambda;
    bool applySyncCorrection;

    Real MACprojection_AMRMG_eps;                        // Solver tolerance
    int MACprojection_AMRMG_num_smooth_down;             // MG pre smoothing iters
    int MACprojection_AMRMG_num_smooth_up;               // MG post smoothing iters
    int MACprojection_AMRMG_num_smooth_bottom;           // MG bottom smoothing iters
    int MACprojection_AMRMG_num_smooth_precond;          // MG preconditioner smoothing iters
    int MACprojection_AMRMG_numMG;                       // 1=V-cycle, 2=W-cycle, etc...
    int MACprojection_AMRMG_imin;                        // Min number of V-cycles
    int MACprojection_AMRMG_imax;                        // Max number of V-cycles
    Real MACprojection_AMRMG_hang;                       // Defaults to 1e-15
    Real MACprojection_AMRMG_normThresh;                 // Defaults to 1e-30
    int MACprojection_AMRMG_maxDepth;                    // Max MG depth (-1 for as deep as possible)
    int MACprojection_AMRMG_relaxMode;                   // 0 = custom, 1 = Jacobi, 2 = looseGSRB, 3 = levelGSRB, 4 = overlapGSRB
    int MACprojection_AMRMG_precondMode;
    int MACprojection_AMRMG_verbosity;                   //

    Real MACprojection_bottom_eps;                       // Solver tolerance
    Real MACprojection_bottom_reps;                      // Solver relative tolerance
    int MACprojection_bottom_imax;                       // Max BiCGStab iterations
    int MACprojection_bottom_numRestarts;                // Max BiCGStab restarts
    Real MACprojection_bottom_hang;                      //
    Real MACprojection_bottom_small;                     //
    int MACprojection_bottom_normType;                   //
    int MACprojection_bottom_verbosity;                  //

    // Approximate projection solver overrides (For CC velocities)
    Real CCprojection_AMRMG_eps;                        // Solver tolerance
    int CCprojection_AMRMG_num_smooth_down;             // MG pre smoothing iters
    int CCprojection_AMRMG_num_smooth_up;               // MG post smoothing iters
    int CCprojection_AMRMG_num_smooth_bottom;           // MG bottom smoothing iters
    int CCprojection_AMRMG_num_smooth_precond;          // MG preconditioner smoothing iters
    int CCprojection_AMRMG_numMG;                       // 1=V-cycle, 2=W-cycle, etc...
    int CCprojection_AMRMG_imin;                        // Min number of V-cycles
    int CCprojection_AMRMG_imax;                        // Max number of V-cycles
    Real CCprojection_AMRMG_hang;                       // Defaults to 1e-15
    Real CCprojection_AMRMG_normThresh;                 // Defaults to 1e-30
    int CCprojection_AMRMG_maxDepth;                    // Max MG depth (-1 for as deep as possible)
    int CCprojection_AMRMG_relaxMode;                   // 0 = custom, 1 = Jacobi, 2 = looseGSRB, 3 = levelGSRB, 4 = overlapGSRB
    int CCprojection_AMRMG_precondMode;
    int CCprojection_AMRMG_verbosity;                   //

    Real CCprojection_bottom_eps;                       // Solver tolerance
    Real CCprojection_bottom_reps;                      // Solver relative tolerance
    int CCprojection_bottom_imax;                       // Max BiCGStab iterations
    int CCprojection_bottom_numRestarts;                // Max BiCGStab restarts
    Real CCprojection_bottom_hang;                      //
    Real CCprojection_bottom_small;                     //
    int CCprojection_bottom_normType;                   //
    int CCprojection_bottom_verbosity;                  //

    // Sync (approximate) projection solver overrides (For composite CC velocities)
    Real syncProjection_AMRMG_eps;                        // Solver tolerance
    int syncProjection_AMRMG_num_smooth_down;             // MG pre smoothing iters
    int syncProjection_AMRMG_num_smooth_up;               // MG post smoothing iters
    int syncProjection_AMRMG_num_smooth_bottom;           // MG bottom smoothing iters
    int syncProjection_AMRMG_num_smooth_precond;          // MG preconditioner smoothing iters
    int syncProjection_AMRMG_numMG;                       // 1=V-cycle, 2=W-cycle, etc...
    int syncProjection_AMRMG_imin;                        // Min number of V-cycles
    int syncProjection_AMRMG_imax;                        // Max number of V-cycles
    Real syncProjection_AMRMG_hang;                       // Defaults to 1e-15
    Real syncProjection_AMRMG_normThresh;                 // Defaults to 1e-30
    int syncProjection_AMRMG_maxDepth;                    // Max MG depth (-1 for as deep as possible)
    int syncProjection_AMRMG_relaxMode;                   // 0 = custom, 1 = Jacobi, 2 = looseGSRB, 3 = levelGSRB, 4 = overlapGSRB
    int syncProjection_AMRMG_precondMode;
    int syncProjection_AMRMG_verbosity;                   //

    Real syncProjection_bottom_eps;                       // Solver tolerance
    Real syncProjection_bottom_reps;                      // Solver relative tolerance
    int syncProjection_bottom_imax;                       // Max BiCGStab iterations
    int syncProjection_bottom_numRestarts;                // Max BiCGStab restarts
    Real syncProjection_bottom_hang;                      //
    Real syncProjection_bottom_small;                     //
    int syncProjection_bottom_normType;                   //
    int syncProjection_bottom_verbosity;                  //

    // Volume discrepancy (freestream preservation) solver overrides
    Real VD_AMRMG_eps;                        // Solver tolerance
    int VD_AMRMG_num_smooth_down;             // MG pre smoothing iters
    int VD_AMRMG_num_smooth_up;               // MG post smoothing iters
    int VD_AMRMG_num_smooth_bottom;           // MG bottom smoothing iters
    int VD_AMRMG_num_smooth_precond;          // MG preconditioner smoothing iters
    int VD_AMRMG_numMG;                       // 1=V-cycle, 2=W-cycle, etc...
    int VD_AMRMG_imin;                        // Min number of V-cycles
    int VD_AMRMG_imax;                        // Max number of V-cycles
    Real VD_AMRMG_hang;                       // Defaults to 1e-15
    Real VD_AMRMG_normThresh;                 // Defaults to 1e-30
    int VD_AMRMG_maxDepth;                    // Max MG depth (-1 for as deep as possible)
    int VD_AMRMG_relaxMode;                   // 0 = custom, 1 = Jacobi, 2 = looseGSRB, 3 = levelGSRB, 4 = overlapGSRB
    int VD_AMRMG_precondMode;
    int VD_AMRMG_verbosity;                   //

    Real VD_bottom_eps;                       // Solver tolerance
    Real VD_bottom_reps;                      // Solver relative tolerance
    int VD_bottom_imax;                       // Max BiCGStab iterations
    int VD_bottom_numRestarts;                // Max BiCGStab restarts
    Real VD_bottom_hang;                      //
    Real VD_bottom_small;                     //
    int VD_bottom_normType;                   //
    int VD_bottom_verbosity;                  //
};


#endif //!__Parameters_H__INCLUDED__
