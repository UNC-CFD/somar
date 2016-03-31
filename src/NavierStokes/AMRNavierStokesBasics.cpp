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
#include "AMRNavierStokes.H"
#include "Constants.H"
#include "ProblemContext.H"


// Initialize static variables here

// These are constants
int AMRNavierStokes::s_num_scal_comps = 1;            // number of scalars (not including lambda)
const bool AMRNavierStokes::s_considerCellSizes = false;

// These will be set by input file
bool AMRNavierStokes::s_ppInit = false;
RealVect AMRNavierStokes::s_domLength = RealVect::Zero;
Real AMRNavierStokes::s_bogus_value = quietNAN;

bool AMRNavierStokes::s_specifyInitialGrids = false;
string AMRNavierStokes::s_initialGridFile;

int  AMRNavierStokes::s_tags_grow = 0;
Real AMRNavierStokes::s_magvort_tag_quota = 0.0;
Tuple<Real,3> AMRNavierStokes::s_vort_tag_tol;
Real AMRNavierStokes::s_vel_tag_tol = 0.0;
Real AMRNavierStokes::s_buoyancy_tag_tol = 0.0;
Real AMRNavierStokes::s_pressure_tag_tol = 0.0;
bool AMRNavierStokes::s_vert_extrude_tags = false;

Real AMRNavierStokes::s_init_shrink = 1.0;
Real AMRNavierStokes::s_max_dt = 1.0e8;
Real AMRNavierStokes::s_prescribedDt = -1.0;
Real AMRNavierStokes::s_max_dt_grow = 1.0e8;

bool AMRNavierStokes::s_limitDtViaViscosity = true;
bool AMRNavierStokes::s_limitDtViaDiffusion = true;
bool AMRNavierStokes::s_limitDtViaPressureGradient = false;
bool AMRNavierStokes::s_limitDtViaInternalWaveSpeed = false;

Real AMRNavierStokes::s_tidalOmega = 0.0;
RealVect AMRNavierStokes::s_tidalU0 = RealVect::Zero;

bool AMRNavierStokes::s_advective_momentum_reflux = false;
bool AMRNavierStokes::s_diffusive_momentum_reflux = false;
bool AMRNavierStokes::s_implicit_momentum_reflux = false;
bool AMRNavierStokes::s_advective_scalar_reflux = false;
bool AMRNavierStokes::s_diffusive_scalar_reflux = false;
bool AMRNavierStokes::s_implicit_scalar_reflux = false;
bool AMRNavierStokes::s_advective_lambda_reflux = false;

bool AMRNavierStokes::s_isIncompressible = true;
int  AMRNavierStokes::s_initial_projection_iters = 1;
int  AMRNavierStokes::s_initial_pressure_iters = 1;
int  AMRNavierStokes::s_level_projection_iters = 1;
int  AMRNavierStokes::s_sync_projection_iters = 1;
bool AMRNavierStokes::s_smooth_after_regrid = false;
Real AMRNavierStokes::s_regrid_smoothing_coeff = 4.0;
int  AMRNavierStokes::s_step_number = 0;

bool AMRNavierStokes::s_applyFreestreamCorrection = false;
Real AMRNavierStokes::s_etaLambda = 0.0;

int AMRNavierStokes::s_updateScheme = ProblemContext::UpdateScheme::FiniteVolume;  // Finite volume method or RK3?

int  AMRNavierStokes::s_normalPredOrderVel = MappedAdvectionUtil::PPM_NORMAL_PRED;
bool AMRNavierStokes::s_useFourthOrderSlopesVel = true;    // This helps quite a bit!
bool AMRNavierStokes::s_useLimitingVel = false;            // With 4th order slopes, this only turns off some of the limiters.
bool AMRNavierStokes::s_useHighOrderLimiterVel = false;    // Unnecessary for low mach flows.
bool AMRNavierStokes::s_useUpwindingVel = true;

int  AMRNavierStokes::s_normalPredOrderScal = MappedAdvectionUtil::PPM_NORMAL_PRED;
bool AMRNavierStokes::s_useFourthOrderSlopesScal = true;    // This helps quite a bit!
bool AMRNavierStokes::s_useLimitingScal = true;             // With 4th order slopes, this only turns off some of the limiters.
bool AMRNavierStokes::s_useHighOrderLimiterScal = false;    // Unnecessary for low mach flows.
bool AMRNavierStokes::s_useUpwindingScal = true;

int AMRNavierStokes::s_viscSolverScheme = ProblemContext::HeatSolverScheme::TGA;
int AMRNavierStokes::s_diffSolverScheme = ProblemContext::HeatSolverScheme::TGA;
Real AMRNavierStokes::s_nu = 0.0;
Vector<Real> AMRNavierStokes::s_scal_coeffs;

// Solver parameters
Real AMRNavierStokes::s_AMRMG_eps = 1e-6;                   // Solver tolerance
int AMRNavierStokes::s_AMRMG_num_smooth_down = 2;           // MG pre smoothing iters
int AMRNavierStokes::s_AMRMG_num_smooth_up = 2;             // MG post smoothing iters
int AMRNavierStokes::s_AMRMG_num_smooth_bottom = 2;         // MG bottom smoothing iters
int AMRNavierStokes::s_AMRMG_num_smooth_precond = 2;        // MG preconditioner smoothing iters
int AMRNavierStokes::s_AMRMG_numMG = 1;                     // 1=V-cycle, 2=W-cycle, etc...
int AMRNavierStokes::s_AMRMG_imin = 5;                      // Min number of V-cycles
int AMRNavierStokes::s_AMRMG_imax = 20;                     // Max number of V-cycles
Real AMRNavierStokes::s_AMRMG_hang = 1e-15;                 // Defaults to 1e-15
Real AMRNavierStokes::s_AMRMG_normThresh = 1e-30;           // Defaults to 1e-30
int AMRNavierStokes::s_AMRMG_maxDepth = -1;                 // Max MG depth (-1 for as deep as possible)
int AMRNavierStokes::s_AMRMG_verbosity = AMRNavierStokes::s_verbosity;
int AMRNavierStokes::s_AMRMG_relaxMode = ProblemContext::RelaxMode::LEVEL_GSRB;
int AMRNavierStokes::s_AMRMG_precondMode = ProblemContext::RelaxMode::LEVEL_GSRB;

Real AMRNavierStokes::s_bottom_eps = 1e-6;                  // Bottom solver tolerance
Real AMRNavierStokes::s_bottom_reps = 1e-12;                // Bottom solver relative tolerance
int AMRNavierStokes::s_bottom_imax = 80;                    // Max BiCGStab iterations
int AMRNavierStokes::s_bottom_numRestarts = 5;              // Max BiCGStab restarts
Real AMRNavierStokes::s_bottom_hang = 1e-15;                //
Real AMRNavierStokes::s_bottom_small = 1e-30;               //
int AMRNavierStokes::s_bottom_normType = 2;
int AMRNavierStokes::s_bottom_verbosity = AMRNavierStokes::s_verbosity;

// Viscous solver overrides
Real AMRNavierStokes::s_viscous_AMRMG_eps = 1e-6;                   // Solver tolerance
int AMRNavierStokes::s_viscous_AMRMG_num_smooth_down = 2;           // MG pre smoothing iters
int AMRNavierStokes::s_viscous_AMRMG_num_smooth_up = 2;             // MG post smoothing iters
int AMRNavierStokes::s_viscous_AMRMG_num_smooth_bottom = 2;         // MG bottom smoothing iters
int AMRNavierStokes::s_viscous_AMRMG_num_smooth_precond = 2;        // MG preconditioner smoothing iters
int AMRNavierStokes::s_viscous_AMRMG_numMG = 1;                     // 1=V-cycle, 2=W-cycle, etc...
int AMRNavierStokes::s_viscous_AMRMG_imin = 5;                      // Min number of V-cycles
int AMRNavierStokes::s_viscous_AMRMG_imax = 20;                     // Max number of V-cycles
Real AMRNavierStokes::s_viscous_AMRMG_hang = 1e-15;                 // Defaults to 1e-15
Real AMRNavierStokes::s_viscous_AMRMG_normThresh = 1e-30;           // Defaults to 1e-30
int AMRNavierStokes::s_viscous_AMRMG_maxDepth = -1;                 // Max MG depth (-1 for as deep as possible)
int AMRNavierStokes::s_viscous_AMRMG_relaxMode = ProblemContext::RelaxMode::LEVEL_GSRB;
int AMRNavierStokes::s_viscous_AMRMG_precondMode = ProblemContext::RelaxMode::LEVEL_GSRB;
int AMRNavierStokes::s_viscous_AMRMG_verbosity = AMRNavierStokes::s_verbosity;

Real AMRNavierStokes::s_viscous_bottom_eps = 1e-6;                  // Bottom solver tolerance
Real AMRNavierStokes::s_viscous_bottom_reps = 1e-12;                // Bottom solver relative tolerance
int AMRNavierStokes::s_viscous_bottom_imax = 80;                    // Max BiCGStab iterations
int AMRNavierStokes::s_viscous_bottom_numRestarts = 5;              // Max BiCGStab restarts
Real AMRNavierStokes::s_viscous_bottom_hang = 1e-15;                //
Real AMRNavierStokes::s_viscous_bottom_small = 1e-30;               //
int AMRNavierStokes::s_viscous_bottom_normType = 2;                 //
int AMRNavierStokes::s_viscous_bottom_verbosity = AMRNavierStokes::s_verbosity;

int AMRNavierStokes::s_gravityMethod = ProblemContext::GravityMethod::EXPLICIT;
Real AMRNavierStokes::s_gravityTheta = 0.6;
int AMRNavierStokes::s_nonlinearDifferencingForm = ProblemContext::NonlinearDifferencingForm::CONSERVATIVE;

bool AMRNavierStokes::s_write_stdout = true;

// initialize plotfile options
bool AMRNavierStokes::s_write_scalars = true;
bool AMRNavierStokes::s_write_scalarsMinusBackground = true;
bool AMRNavierStokes::s_write_pressure = true;
bool AMRNavierStokes::s_write_divergence = false;
bool AMRNavierStokes::s_write_lambda = false;
bool AMRNavierStokes::s_write_grad_eLambda = false;
bool AMRNavierStokes::s_write_vorticity = false;
bool AMRNavierStokes::s_write_streamfunction = false;
bool AMRNavierStokes::s_write_proc_ids = false;
bool AMRNavierStokes::s_write_level_ids = false;
bool AMRNavierStokes::s_write_grids = false;
bool AMRNavierStokes::s_write_displacement = true;
bool AMRNavierStokes::s_write_geometry = false;

// Total (composite) energy
Real AMRNavierStokes::s_totalEnergy = 1e8;

// names of state variables
const char* AMRNavierStokes::s_vel_names[CH_SPACEDIM] =
{
    "x-velocity",
#if CH_SPACEDIM >= 2
    "y-velocity",
#endif
#if CH_SPACEDIM >= 3
    "z-velocity"
#endif
};
Vector<std::string> AMRNavierStokes::s_scal_names(0);


// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
AMRNavierStokes::AMRNavierStokes ()
: m_finest_level(false),
  m_is_empty(true),
  m_cfl(0.8)
{
    if (s_verbosity >= 5) {
        pout() << "AMRNavierStokes default constructor" << endl;
    }

    m_levGeoPtr = NULL;
    m_physBCPtr = NULL;

    m_vel_new_ptr = NULL;
    m_vel_old_ptr = NULL;
    m_lambda_new_ptr = NULL;
    m_lambda_old_ptr = NULL;
    m_scal_new.resize(0);
    m_scal_old.resize(0);

    m_scal_fluxreg_ptrs.resize(0);

    m_c0iPtr = NULL;
    m_phi0Ptr = NULL;

    m_regrid_smoothing_done = false;
}


// -----------------------------------------------------------------------------
// Destructor
// NOTE: Do not delete m_physBCPtr. The factory owns it.
// -----------------------------------------------------------------------------
AMRNavierStokes::~AMRNavierStokes ()
{
    if (s_verbosity >= 5) {
        pout () << "AMRNavierStokes destructor" << endl;
    }

    if (m_vel_new_ptr != NULL) {
        delete m_vel_new_ptr;
        m_vel_new_ptr = NULL;
    }

    if (m_vel_old_ptr != NULL) {
        delete m_vel_old_ptr;
        m_vel_old_ptr = NULL;
    }

    if (m_lambda_new_ptr != NULL) {
        delete m_lambda_new_ptr;
        m_lambda_new_ptr = NULL;
    }

    if (m_lambda_old_ptr != NULL) {
        delete m_lambda_old_ptr;
        m_lambda_old_ptr = NULL;
    }

    // loop over scalars and delete
    int nScalComp = m_scal_new.size();
    for (int comp = 0; comp < nScalComp; ++comp) {
        if (m_scal_new[comp] != NULL) {
            delete m_scal_new[comp];
            m_scal_new[comp] = NULL;
        }
    }

    nScalComp = m_scal_old.size();
    for (int comp = 0; comp < nScalComp; ++comp) {
        if (m_scal_old[comp] != NULL) {
            delete m_scal_old[comp];
            m_scal_old[comp] = NULL;
        }
    }

    if (m_levGeoPtr != NULL) {
        delete m_levGeoPtr;
        m_levGeoPtr = NULL;
    }

    nScalComp = m_scal_fluxreg_ptrs.size();
    for (int comp = 0; comp < nScalComp; ++comp) {
        if (m_scal_fluxreg_ptrs[comp] != NULL) {
            delete m_scal_fluxreg_ptrs[comp];
            m_scal_fluxreg_ptrs[comp] = NULL;
        }
    }

    if (m_c0iPtr != NULL) {
        delete m_c0iPtr;
        m_c0iPtr = NULL;
    }

    if (m_phi0Ptr != NULL) {
        delete m_phi0Ptr;
        m_phi0Ptr = NULL;
    }
}


// -----------------------------------------------------------------------------
// Define for starting from scratch
// -----------------------------------------------------------------------------
void AMRNavierStokes::define (MappedAMRLevel*       a_coarser_level_ptr,
                              const ProblemDomain&  a_problem_domain,
                              int                   a_level,
                              const IntVect&        a_ref_ratio)
{
    if (s_verbosity >= 5) {
        pout () << "AMRNavierStokes::define " << a_level << endl;
    }

    // Call the generic AMRLevel define
    MappedAMRLevel::define(a_coarser_level_ptr,
                           a_problem_domain,
                           a_level,
                           a_ref_ratio);

    // Read stuff from parmParse database
    if (!s_ppInit) {
        this->readParameters();
    }

    // Create this level's levGeo object
    m_levGeoPtr = new LevelGeometry;

    // To define the new levGeo, we will need these.
    // (The coarser levGeo ptr will be reset in initialGrid and regrid though.)
    const RealVect m_dx = s_domLength / RealVect(m_problem_domain.size());
    LevelGeometry* crseLevGeoPtr = NULL;

    if (a_coarser_level_ptr != NULL) {
        // Set coarser level's finer_level_ptr to point to this
        a_coarser_level_ptr->finerLevelPtr(static_cast<MappedAMRLevel*>(this));

        // Cast the coarser level ptr into a NS ptr.
        AMRNavierStokes* crseNSPtr = dynamic_cast<AMRNavierStokes*>(a_coarser_level_ptr);
        if (crseNSPtr == NULL) {
            MayDay::Error ("AMRNavierStokes::define: a_coarser_level_ptr is not castable to AMRNavierStokes*");
        }

        // Grab data from the coarser level
        crseLevGeoPtr = crseNSPtr->m_levGeoPtr;
        m_cfl = crseNSPtr->m_cfl;
    }

    // Now finish defining this level.
    m_levGeoPtr->define(m_dx, crseLevGeoPtr);
}


// -----------------------------------------------------------------------------
// readParameters
// -----------------------------------------------------------------------------
void AMRNavierStokes::readParameters ()
{
    if (s_verbosity >= 5) {
        pout () << "AMRNavierStokes::readParameters " << endl;
    }

    // Name the scalars
    s_scal_names.resize(s_num_scal_comps);
    for (int comp = 0; comp < s_num_scal_comps; ++comp) {
        ostringstream name;
        name << "scalar_" << comp;
        s_scal_names[comp] = name.str();
    }

    const ProblemContext* ctx = ProblemContext::getInstance();

    D_TERM(s_domLength[0] = ctx->domainLength[0];,
           s_domLength[1] = ctx->domainLength[1];,
           s_domLength[2] = ctx->domainLength[2];)
    s_specifyInitialGrids = ctx->specifyInitialGrids;
    if (s_specifyInitialGrids) {
        s_initialGridFile = ctx->initialGridFile;
    }
    s_tags_grow = ctx->tags_grow;
    CH_assert(s_tags_grow >= 0);

    s_magvort_tag_quota = ctx->magvort_tag_quota;
    s_vort_tag_tol = ctx->vort_tag_tol;
    s_vel_tag_tol = ctx->vel_tag_tol;
    s_buoyancy_tag_tol = ctx->buoyancy_tag_tol;
    s_pressure_tag_tol = ctx->pressure_tag_tol;
    s_vert_extrude_tags = ctx->vert_extrude_tags;

    s_write_stdout = ctx->write_stdout;
    s_init_shrink = ctx->init_dt_multiplier;
    s_max_dt = ctx->max_dt;
    s_prescribedDt = ctx->fixed_dt;
    s_max_dt_grow = ctx->max_dt_grow;

    s_limitDtViaViscosity = ctx->limitDtViaViscosity;
    s_limitDtViaDiffusion = ctx->limitDtViaDiffusion;
    s_limitDtViaPressureGradient = ctx->limitDtViaPressureGradient;
    s_limitDtViaInternalWaveSpeed = ctx->limitDtViaInternalWaveSpeed;
    s_tidalOmega = ctx->tidalOmega;
    s_tidalU0 = ctx->tidalU0;

    s_bogus_value = ctx->bogus_value;
    s_smooth_after_regrid = ctx->smooth_after_regrid;
    s_regrid_smoothing_coeff = ctx->regrid_smoothing_coeff;
    s_advective_momentum_reflux = ctx->advective_momentum_reflux;
    s_diffusive_momentum_reflux = ctx->diffusive_momentum_reflux;
    s_implicit_momentum_reflux = ctx->implicit_momentum_reflux;
    s_advective_scalar_reflux = ctx->advective_scalar_reflux;
    s_diffusive_scalar_reflux = ctx->diffusive_scalar_reflux;
    s_implicit_scalar_reflux = ctx->implicit_scalar_reflux;
    s_advective_lambda_reflux = ctx->advective_lambda_reflux;

    s_gravityMethod = int(ctx->gravityMethod);
    if (s_gravityMethod == ProblemContext::GravityMethod::IMPLICIT) {
        s_gravityTheta = ctx->gravityTheta;
    }

    s_nonlinearDifferencingForm = ctx->nonlinearDifferencingForm;

    s_viscSolverScheme = ctx->viscSolverScheme;
    s_nu = ctx->nu;

    s_diffSolverScheme = ctx->diffSolverScheme;
    s_num_scal_comps = ctx->num_scal_comps;
    s_scal_coeffs = ctx->scal_coeffs;

    // Advection settings
    s_updateScheme = ctx->updateScheme;

    s_normalPredOrderVel = ctx->normalPredOrderVel;
    s_useFourthOrderSlopesVel = ctx->useFourthOrderSlopesVel;
    s_useLimitingVel = ctx->useLimitingVel;
    s_useHighOrderLimiterVel = ctx->useHighOrderLimiterVel;
    s_useUpwindingVel = ctx->useUpwindingVel;

    s_normalPredOrderScal = ctx->normalPredOrderScal;
    s_useFourthOrderSlopesScal = ctx->useFourthOrderSlopesScal;
    s_useLimitingScal = ctx->useLimitingScal;
    s_useHighOrderLimiterScal = ctx->useHighOrderLimiterScal;
    s_useUpwindingScal = ctx->useUpwindingScal;

    // Solver parameters
    s_AMRMG_eps = ctx->AMRMG_eps;
    s_AMRMG_num_smooth_down = ctx->AMRMG_num_smooth_down;
    s_AMRMG_num_smooth_up = ctx->AMRMG_num_smooth_up;
    s_AMRMG_num_smooth_bottom = ctx->AMRMG_num_smooth_bottom;
    s_AMRMG_num_smooth_precond = ctx->AMRMG_num_smooth_precond;
    s_AMRMG_numMG = ctx->AMRMG_numMG;
    s_AMRMG_imin = ctx->AMRMG_imin;
    s_AMRMG_imax = ctx->AMRMG_imax;
    s_AMRMG_hang = ctx->AMRMG_hang;
    s_AMRMG_normThresh = ctx->AMRMG_normThresh;
    s_AMRMG_maxDepth = ctx->AMRMG_maxDepth;
    s_AMRMG_verbosity = ctx->AMRMG_verbosity;
    s_AMRMG_relaxMode = ctx->AMRMG_relaxMode;
    s_AMRMG_precondMode = ctx->AMRMG_precondMode;

    s_bottom_eps = ctx->bottom_eps;
    s_bottom_reps = ctx->bottom_reps;
    s_bottom_imax = ctx->bottom_imax;
    s_bottom_numRestarts = ctx->bottom_numRestarts;
    s_bottom_hang = ctx->bottom_hang;
    s_bottom_small = ctx->bottom_small;
    s_bottom_normType = ctx->bottom_normType;
    s_bottom_verbosity = ctx->bottom_verbosity;

    // Viscous solver overrides
    s_viscous_AMRMG_eps = ctx->viscous_AMRMG_eps;
    s_viscous_AMRMG_num_smooth_down = ctx->viscous_AMRMG_num_smooth_down;
    s_viscous_AMRMG_num_smooth_up = ctx->viscous_AMRMG_num_smooth_up;
    s_viscous_AMRMG_num_smooth_bottom = ctx->viscous_AMRMG_num_smooth_bottom;
    s_viscous_AMRMG_num_smooth_precond = ctx->viscous_AMRMG_num_smooth_precond;
    s_viscous_AMRMG_numMG = ctx->viscous_AMRMG_numMG;
    s_viscous_AMRMG_imin = ctx->viscous_AMRMG_imin;
    s_viscous_AMRMG_imax = ctx->viscous_AMRMG_imax;
    s_viscous_AMRMG_hang = ctx->viscous_AMRMG_hang;
    s_viscous_AMRMG_normThresh = ctx->viscous_AMRMG_normThresh;
    s_viscous_AMRMG_maxDepth = ctx->viscous_AMRMG_maxDepth;
    s_viscous_AMRMG_verbosity = ctx->viscous_AMRMG_verbosity;
    s_viscous_AMRMG_relaxMode = ctx->viscous_AMRMG_relaxMode;
    s_viscous_AMRMG_precondMode = ctx->viscous_AMRMG_precondMode;

    s_viscous_bottom_eps = ctx->viscous_bottom_eps;
    s_viscous_bottom_reps = ctx->viscous_bottom_reps;
    s_viscous_bottom_imax = ctx->viscous_bottom_imax;
    s_viscous_bottom_numRestarts = ctx->viscous_bottom_numRestarts;
    s_viscous_bottom_hang = ctx->viscous_bottom_hang;
    s_viscous_bottom_small = ctx->viscous_bottom_small;
    s_viscous_bottom_normType = ctx->viscous_bottom_normType;
    s_viscous_bottom_verbosity = ctx->viscous_bottom_verbosity;

    // Projection stuff
    s_isIncompressible = ctx->isIncompressible;
    s_initial_projection_iters = ctx->initial_projection_iters;
    s_initial_pressure_iters = ctx->initial_pressure_iters;
    s_level_projection_iters = ctx->level_projection_iters;
    s_sync_projection_iters = ctx->sync_projection_iters;
    s_applyFreestreamCorrection = ctx->applyVDCorrection;
    s_etaLambda = ctx->etaLambda;

    // Plot stuff
    s_write_divergence = ctx->write_divergence;
    s_write_lambda = ctx->write_lambda;
    s_write_grad_eLambda = ctx->write_grad_eLambda;
    s_write_pressure = ctx->write_pressure;
    s_write_vorticity = ctx->write_vorticity;
    s_write_streamfunction = ctx->write_streamfunction;
    s_write_scalars = ctx->write_scalars;
    s_write_scalarsMinusBackground = ctx->write_scalarsMinusBackground;
    s_write_proc_ids = ctx->write_proc_ids;
    s_write_level_ids = ctx->write_level_ids;
    s_write_grids = ctx->write_grids;
    s_write_displacement = ctx->write_displacement;
    s_write_geometry = ctx->write_geometry;

    // set flag to indicate that we've done this
    s_ppInit = true;
}
