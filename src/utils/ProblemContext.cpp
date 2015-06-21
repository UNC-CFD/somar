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
#include "ProblemContext.H"
#include "ParmParse.H"
#include "Constants.H"
#include <cstdio>
#include <sstream>

#include "CartesianMap.H"
#include "TwistedMap.H"
#include "BeamGeneratorMap.H"
#include "NewBeamGeneratorMap.H"
#include "CylindricalMap.H"
#include "LedgeMap.H"

#include "DEMMap.H"
#include "AdvectionTestBCUtil.H"
#include "LockExchangeBCUtil.H"
#include "BeamGenerationBCUtil.H"
#include "InternalWaveBCUtil.H"
#include "TaylorGreenBCUtil.H"
#include "VortexStreetBCUtil.H"
#include "HorizConvBCUtil.H"
#include "SolitaryWaveBCUtil.H"
#include "DJLBCUtil.H"


ProblemContext* ProblemContext::s_singletonPtr = NULL;
bool ProblemContext::s_AMRParamsRead = false;
bool ProblemContext::s_geometryParamsRead = false;
bool ProblemContext::s_plotParamsRead = false;
bool ProblemContext::s_IBCParamsRead = false;
bool ProblemContext::s_advectionParamsRead = false;
bool ProblemContext::s_solverParamsRead = false;
bool ProblemContext::s_projectionParamsRead = false;


// -----------------------------------------------------------------------------
// Constructor
// This will only be called once. It reads the input file.
// -----------------------------------------------------------------------------
ProblemContext::ProblemContext ()
{
    this->readAMR();
    this->readGeometry();
    this->readPlot();
    this->readIBC();
    this->readAdvection();
    this->readSolver();
    this->readProjection();
    this->readPlot();
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
ProblemContext::~ProblemContext ()
{;}


// -----------------------------------------------------------------------------
// This returns the single ProblemContext object.
// -----------------------------------------------------------------------------
const ProblemContext* ProblemContext::getInstance ()
{
    if (s_singletonPtr == NULL) {
        s_singletonPtr = new ProblemContext();
    }
    return s_singletonPtr;
}


// -----------------------------------------------------------------------------
// This deletes the static pointer in hopes valgrind will take notice.
// -----------------------------------------------------------------------------
void ProblemContext::freeMemory ()
{
    delete s_singletonPtr;
    s_singletonPtr = NULL;
}


// -----------------------------------------------------------------------------
// Read the amr.* inputs.
// -----------------------------------------------------------------------------
void ProblemContext::readAMR ()
{
    if (s_AMRParamsRead) return;
    s_AMRParamsRead = true;

    ParmParse ppAMR("amr");
    pout() << "ProblemContext::readAMR:" << endl;

    // Temp vars
    Vector<int> vint(SpaceDim);
    Vector<Real> vreal(SpaceDim);

    ppAMR.getarr("nx", vint, 0, SpaceDim);
    nx = IntVect(vint);
    pout() << "\tnx = " << nx << endl;

    vint = Vector<int>(SpaceDim, 0);
    ppAMR.queryarr("nx_offset", vint, 0, SpaceDim);
    nx_offset = IntVect(vint);
    pout() << "\tnx_offset = " << nx_offset << endl;

    ppAMR.getarr("length", vreal, 0, SpaceDim);
    domainLength = RealVect(vreal);
    pout() << "\tdomainLength = " << domainLength << endl;

    RealVect dx = domainLength / RealVect(nx);
    pout() << "\tcoarse level dx = " << dx << endl;

    Real lepticity = Min(dx[0], dx[SpaceDim-2]) / domainLength[SpaceDim-1];
    pout() << "leptic ratio = " << lepticity << endl;
    pout() << "inverse square leptic ratio = " << pow(lepticity, -2.0) << endl;


    vreal = Vector<Real>(SpaceDim, 0.0);
    if (ppAMR.queryarr("offset", vreal, 0, SpaceDim)) {
        D_TERM(nx_offset[0] = int(vreal[0] / dx[0]);,
               nx_offset[1] = int(vreal[1] / dx[1]);,
               nx_offset[2] = int(vreal[2] / dx[2]);)
        pout() << "\toffset (x0) = " << RealVect(vreal) << endl;
        pout() << "\tnx_offset = " << nx_offset << endl;
    }

    vint = Vector<int>(SpaceDim, 0);
    ppAMR.queryarr("isPeriodic", vint, 0, SpaceDim);
    isPeriodic[0] = (vint[0] != 0);
    pout() << "\tisPeriodic = (" << vint[0];
    for(int dir = 1; dir < SpaceDim; ++dir) {
        isPeriodic[dir] = (vint[dir] != 0);
        pout() << ", " << vint[dir];
    }
    pout() << ")" << endl;

    domain.define(nx_offset, nx_offset + nx - 1, isPeriodic);
    pout() << "\tdomain = " << domain << endl;

    max_dt_grow = 1.5;
    ppAMR.query("max_dt_grow", max_dt_grow);
    if (max_dt_grow <= 0.0) max_dt_grow = 1.0e8;
    pout() << "\tmax_dt_grow = " << max_dt_grow << endl;

    fixed_dt = -1.0;
    if (ppAMR.query("fixed_dt", fixed_dt)) {
        pout() << "\tfixed_dt = " << fixed_dt << endl;
    } else {
        pout() << "\tfixed_dt = ---" << endl;
    }

    ppAMR.get("final", stopTime);
    pout() << "\tstopTime = " << stopTime << endl;

    ppAMR.get("cfl", cfl);
    pout() << "\tcfl = " << cfl << endl;

    useSubcycling = true;
    ppAMR.query("useSubcycling", useSubcycling);
    pout() << "\tuseSubcycling = " << (useSubcycling? "true": "false") << endl;

    ppAMR.get("maxsteps", maxsteps);
    pout() << "\tnstop = " << maxsteps << endl;

    verbosity = 0;
    ppAMR.query("verbosity", verbosity);
    pout() << "\tverbosity = " << verbosity << std::endl;
    CH_assert (verbosity >= 0);

    ppAMR.get("maxlevel", max_level);
    pout() << "\tmax_level = " << max_level << endl;
    numlevels = max_level + 1;
    CH_assert(numlevels > 0);

    const int num_read_levels = Max(max_level,1);
    Vector<int> vintLevels(num_read_levels, 1);

    regrid_intervals = vintLevels;
    if (max_level > 0) {
        ppAMR.getarr("regrid_intervals", regrid_intervals, 0, num_read_levels);
    }
    pout() << "\tregrid_intervals = " << regrid_intervals << endl;

    Vector<int> levRefRatio(SpaceDim, 1);
    bool defaultSet = ppAMR.queryarr("refratio", levRefRatio, 0, SpaceDim);

    refRatios = Vector<IntVect>(numlevels, IntVect(levRefRatio));
    for (int lev = 0; lev < numlevels-1; ++lev) {
        std::ostringstream str;
        str << "refratio_lev" << lev;
        bool levSet = ppAMR.queryarr(str.str().c_str(), levRefRatio, 0, SpaceDim);

        if (!defaultSet && !levSet) {
            MayDay::Error("Must set a default ref ratio or set each level's ref ratio individually.");
        }

        if (levSet) {
            refRatios[lev] = IntVect(D_DECL(levRefRatio[0], levRefRatio[1], levRefRatio[2]));
        }

        pout() << "\tlevel " << lev << " ref ratio = " << refRatios[lev] << endl;
    }
    refRatios[max_level] = IntVect::Unit; // This is needed by Chombo.

    block_factor = 8;
    ppAMR.query("block_factor", block_factor);
    pout() << "\tblock_factor = " << block_factor << endl;

    bufferSize = 1;
    ppAMR.query("grid_buffer_size", bufferSize);
    pout() << "\tbufferSize = " << bufferSize << std::endl;

    fill_ratio = 0.80;
    ppAMR.query("fill_ratio", fill_ratio);
    pout() << "\tfill_ratio = " << fill_ratio << std::endl;

    splitDirs = Vector<int>(SpaceDim, 1);
    ppAMR.queryarr("splitDirs", splitDirs, 0, SpaceDim);
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (splitDirs[dir] == 0) {
            pout() << "\tGrids will not be split in direction " << dir << "." << endl;
        } else {
            splitDirs[dir] = 1;
        }
    }

    const int nproc = numProc();
    maxGridSize = IntVect::Zero;
    IntVect userDefined = IntVect::Zero;
    {
        Vector<int> _maxGridSize(SpaceDim,0);
        ppAMR.queryarr("max_grid_size", _maxGridSize, 0, SpaceDim);

        for (int dir = 0; dir < SpaceDim; ++dir) {
            // Begin with the blocking factor. Never let the grid sizes
            // become smaller than this value.
            maxGridSize[dir] = block_factor;

            if (_maxGridSize[dir] != 0) {
                // Always use the user's values.
                userDefined[dir] = 1;
                maxGridSize[dir] = Max(_maxGridSize[dir], block_factor);
                if (maxGridSize[dir] != _maxGridSize[dir]) {
                    ostringstream msg;
                    msg << "maxGridSize[" << dir << "] defaulted to the blocking factor "
                        << "(" << block_factor << ") because input value "
                        << "(" << _maxGridSize[dir] << ") is too small";
                    MayDay::Warning(msg.str().c_str());
                }

            } else {
                // No default value -- compute an appropriate maxGridSize.
                if (splitDirs[dir] == 0) {
                    // We do not want to split in this direction.
                    maxGridSize[dir] = nx[dir];

                } else if (nx[dir] % nproc == 0) {
                    // Just divide side length over procs.
                    maxGridSize[dir] = Max(nx[dir] / nproc, block_factor);

                } else {
                    // We're out of options. Maybe the next direction will
                    // be able to be split.
                    maxGridSize[dir] = nx[dir];
                }
            }

            CH_assert(maxGridSize[dir] >= 4); // For Poisson stencil stuff.
            CH_assert(maxGridSize[dir] >= block_factor);
        }
    }
    pout() << "\tmax_grid_size = " << maxGridSize << endl;

    maxBaseGridSize = IntVect::Zero;
    {
        Vector<int> _maxBaseGridSize(SpaceDim,0);
        ppAMR.queryarr("max_base_grid_size", _maxBaseGridSize, 0, SpaceDim);
        for (int dir = 0; dir < SpaceDim; ++dir) {
            // Begin with the blocking factor. Never let the grid sizes
            // become smaller than this value.
            maxBaseGridSize[dir] = block_factor;

            if (_maxBaseGridSize[dir] != 0) {
                // Always use the user's values.
                maxBaseGridSize[dir] = Max(_maxBaseGridSize[dir], block_factor);
                if (maxBaseGridSize[dir] != _maxBaseGridSize[dir]) {
                    ostringstream msg;
                    msg << "maxBaseGridSize[" << dir << "] defaulted to the blocking factor "
                        << "(" << block_factor << ") because input value "
                        << "(" << _maxBaseGridSize[dir] << ") is too small";
                    MayDay::Warning(msg.str().c_str());
                }

            } else if (userDefined[dir] != 0) {
                // Secondary default -- the user's maxGridSize
                maxBaseGridSize[dir] = Max(maxGridSize[dir], block_factor);
                if (maxBaseGridSize[dir] != maxGridSize[dir]) {
                    ostringstream msg;
                    msg << "maxBaseGridSize[" << dir << "] tried to default to  "
                        << "maxGridSize[" << dir << "] = " << maxGridSize[dir]
                        << " which is less than the blocking factor = " << block_factor;
                    MayDay::Error(msg.str().c_str());
                }

            } else {
                // No default value -- compute an appropriate maxBaseGridSize.
                if (splitDirs[dir] == 0) {
                    // We do not want to split in this direction.
                    maxBaseGridSize[dir] = nx[dir];

                } else if (nx[dir] % nproc == 0) {
                    // Just divide side length over procs.
                    maxBaseGridSize[dir] = Max(nx[dir] / nproc, block_factor);

                } else {
                    // We're out of options. Maybe the next direction will
                    // be able to be split.
                    maxBaseGridSize[dir] = nx[dir];
                }
            }

            CH_assert(maxBaseGridSize[dir] >= 4); // For Poisson stencil stuff.
            CH_assert(maxBaseGridSize[dir] >= block_factor);
        }
    }
    pout() << "\tmax_base_grid_size = " << maxBaseGridSize << endl;

    // TODO: This should be in readPlot!
    isRestart = ppAMR.contains("restart_file");
    if (isRestart) {
        ppAMR.get("restart_file", restart_file);
    }
    pout() << "\trestart_file = " << restart_file << endl;

    hasPredefinedGrids = ppAMR.contains("gridfile");
    if (hasPredefinedGrids) {
        std::string gridfile;
        ppAMR.get("gridfile", gridfile);
        pout() << "\tgridFile = " << gridfile << endl;

        MayDay::Error("Predefined grid code not complete in amrins.cpp");
        // readAMRGrids(amrGrids,
        //              gridfile.c_str(),
        //              domain,
        //              maxBaseGridSize,
        //              maxGridSize);
    }

    specifyInitialGrids = false;
    ppAMR.query("specifyInitialGrids", specifyInitialGrids);
    pout() << "\tspecifyInitialGrids = " << specifyInitialGrids << endl;

    if (specifyInitialGrids) {
        ppAMR.get("initialGridFile", initialGridFile);
        pout() << "\tinitialGridFile = " << initialGridFile << endl;
    }

    magvort_tag_quota = 0.0;
    ppAMR.query("magvort_tag_quota", magvort_tag_quota);
    pout() << "\tmagvort_tag_quota = " << magvort_tag_quota << endl;

    // Deprecation warning
    if (ppAMR.contains("vorticity_tag_factor")) {
        if (ppAMR.contains("magvort_tag_quota")) {
            MayDay::Warning("vorticity_tag_factor is deprecated. Found magvort_tag_quota. Using that instead.");
        } else {
            MayDay::Warning("vorticity_tag_factor is deprecated. Re-interpreting as magvort_tag_quota.");
            magvort_tag_quota = 0.0;
            ppAMR.query("vorticity_tag_factor", magvort_tag_quota);
            pout() << "\tmagvort_tag_quota = " << magvort_tag_quota << endl;
        }
    }

    Vector<Real> _threeRV(3, 0.0);
    ppAMR.queryarr("vort_tag_tol", _threeRV, 0, 3);
    for (int dir = 0; dir < 3; ++dir) {
        vort_tag_tol[dir] = _threeRV[dir];
        CH_assert(vort_tag_tol[dir] >= 0.0);
    }
    pout() << "\tvort_tag_tol = "
           << vort_tag_tol[0] << ", "
           << vort_tag_tol[1] << ", "
           << vort_tag_tol[2] << endl;

    vel_tag_tol = 0.0;
    ppAMR.query("vel_tag_tol", vel_tag_tol);
    pout() << "\tvel_tag_tol = " << vel_tag_tol << endl;

    buoyancy_tag_tol = 0.0;
    ppAMR.query("buoyancy_tag_tol", buoyancy_tag_tol);
    pout() << "\tbuoyancy_tag_tol = " << buoyancy_tag_tol << endl;

    // Deprecation warning
    if (ppAMR.contains("density_tag_tol")) {
        if (ppAMR.contains("buoyancy_tag_tol")) {
            MayDay::Warning("density_tag_tol is deprecated. Found buoyancy_tag_tol. Using that instead.");
        } else {
            MayDay::Warning("density_tag_tol is deprecated. Re-interpreting as buoyancy_tag_tol.");
            buoyancy_tag_tol = 0.0;
            ppAMR.query("density_tag_tol", buoyancy_tag_tol);
            pout() << "\tbuoyancy_tag_tol = " << buoyancy_tag_tol << endl;
        }
    }

    pressure_tag_tol = 0.0;
    ppAMR.query("pressure_tag_tol", pressure_tag_tol);
    pout() << "\tpressure_tag_tol = " << pressure_tag_tol << endl;

    vert_extrude_tags = false;
    ppAMR.query("vert_extrude_tags", vert_extrude_tags);
    pout() << "\tvert_extrude_tags = " << vert_extrude_tags << endl;

    tags_grow = 0;
    ppAMR.query("tags_grow", tags_grow);
    pout() << "\ttags_grow = " << tags_grow << endl;

    write_stdout = true;
    ppAMR.query("write_stdout", write_stdout);
    pout() << "\twrite_stdout = " << write_stdout << endl;

    init_dt_multiplier = 1.0;
    ppAMR.query("init_dt_multiplier", init_dt_multiplier);
    pout() << "\tinit_dt_multiplier = " << init_dt_multiplier << endl;

    max_dt = -1.0;
    ppAMR.query("max_dt", max_dt);
    if (max_dt <= 0.0) max_dt = 1.0e8;
    pout() << "\tmax_dt = " << max_dt << endl;

    limitDtViaViscosity = false;
    ppAMR.query("limitDtViaViscosity", limitDtViaViscosity);
    pout() << "\tlimitDtViaViscosity = " << limitDtViaViscosity << endl;

    limitDtViaDiffusion = false;
    ppAMR.query("limitDtViaDiffusion", limitDtViaDiffusion);
    pout() << "\tlimitDtViaDiffusion = " << limitDtViaDiffusion << endl;

    limitDtViaPressureGradient = true;
    ppAMR.query("limitDtViaPressureGradient", limitDtViaPressureGradient);
    pout() << "\tlimitDtViaPressureGradient = " << limitDtViaPressureGradient << endl;

    limitDtViaInternalWaveSpeed = false;
    {
        int useBackgroundScalar = 0;
        ParmParse("ibc").query("useBackgroundScalar", useBackgroundScalar);
        if (useBackgroundScalar) {
            limitDtViaInternalWaveSpeed = true;
        }
    }
    ppAMR.query("limitDtViaInternalWaveSpeed", limitDtViaInternalWaveSpeed);
    pout() << "\tlimitDtViaInternalWaveSpeed = " << limitDtViaInternalWaveSpeed << endl;

    bogus_value = quietNAN;
    ppAMR.query("bogus_val", bogus_value);
    pout() << "\tbogus_val = " << bogus_value << endl;

    smooth_after_regrid = false;
    ppAMR.query("smooth_after_regrid", smooth_after_regrid);
    pout() << "\tsmooth_after_regrid = " << smooth_after_regrid << endl;

    regrid_smoothing_coeff = 0.0;
    if (smooth_after_regrid) {
        ppAMR.get("postRegrid_smoothing_coeff", regrid_smoothing_coeff);
        pout() << "\tpostRegrid_smoothing_coeff = " << regrid_smoothing_coeff << endl;
    }

    advective_momentum_reflux = false;
    ppAMR.query("advective_momentum_reflux", advective_momentum_reflux);
    pout() << "\tadvective_momentum_reflux = " << advective_momentum_reflux << endl;

    diffusive_momentum_reflux = false;
    ppAMR.query("diffusive_momentum_reflux", diffusive_momentum_reflux);
    pout() << "\tdiffusive_momentum_reflux = " << diffusive_momentum_reflux << endl;

    implicit_momentum_reflux = false;
    ppAMR.query("implicit_momentum_reflux", implicit_momentum_reflux);
    pout() << "\timplicit_momentum_reflux = " << implicit_momentum_reflux << endl;

    advective_scalar_reflux = false;
    ppAMR.query("advective_scalar_reflux", advective_scalar_reflux);
    pout() << "\tadvective_scalar_reflux = " << advective_scalar_reflux << endl;

    diffusive_scalar_reflux = false;
    ppAMR.query("diffusive_scalar_reflux", diffusive_scalar_reflux);
    pout() << "\tdiffusive_scalar_reflux = " << diffusive_scalar_reflux << endl;

    implicit_scalar_reflux = false;
    ppAMR.query("implicit_scalar_reflux", implicit_scalar_reflux);
    pout() << "\timplicit_scalar_reflux = " << implicit_scalar_reflux << endl;

    advective_lambda_reflux = false;
    ppAMR.query("advective_lambda_reflux", advective_lambda_reflux);
    pout() << "\tadvective_lambda_reflux = " << advective_lambda_reflux << endl;

    gravityMethod = ProblemContext::GravityMethod::EXPLICIT;
    ppAMR.query("gravityMethod", gravityMethod);
    CH_assert(0 <= gravityMethod);
    CH_assert(gravityMethod < ProblemContext::GravityMethod::_NUM_GRAVITYMETHODS);
    pout() << "\tgravityMethod = " << gravityMethod << endl;

    if (gravityMethod == ProblemContext::GravityMethod::IMPLICIT) {
        ppAMR.query("gravityTheta", gravityTheta);
        pout() << "\tgravityTheta = " << gravityTheta << endl;
    }

    viscSolverScheme = HeatSolverScheme::CRANK_NICOLSON;
    ppAMR.query("viscous_solver_type", viscSolverScheme);
    pout() << "\tviscSolverScheme = " << viscSolverScheme << endl;

    diffSolverScheme = HeatSolverScheme::CRANK_NICOLSON;
    ppAMR.query("diffusive_solver_type", diffSolverScheme);
    pout() << "\tdiffSolverScheme = " << diffSolverScheme << endl;

    num_scal_comps = 1;
    ppAMR.query("num_scal_comps", num_scal_comps);
    pout() << "\tnum_scal_comps = " << num_scal_comps << endl;

    if (num_scal_comps > 0) {
        scal_coeffs.resize(num_scal_comps);
        ppAMR.getarr("scal_diffusion_coeffs", scal_coeffs, 0, num_scal_comps);
        for (int idx = 0; idx < num_scal_comps; ++idx) {
            pout() << "\tscal_coeffs[" << idx << "] = " << scal_coeffs[idx] << endl;
        }
    }

    ppAMR.get("viscosity", nu);
    pout() << "\tnu = " << nu << endl;

    pout() << endl;
}


// -----------------------------------------------------------------------------
// The geometry.* parameters
// -----------------------------------------------------------------------------
void ProblemContext::readGeometry ()
{
    if (s_geometryParamsRead) return;
    s_geometryParamsRead = true;

    ParmParse ppGeo("geometry");
    pout() << "ProblemContext::readGeometry:" << endl;

    ppGeo.get("coordMap", coordMap);
    pout() << "\tcoordMap = " << coordMap << endl;
    CH_assert(0 <= coordMap && coordMap < CoordMap::_NUM_COORD_MAPS);

    switch (coordMap) {
    case CoordMap::TWISTED:
        {
            Vector<Real> _pert(CH_SPACEDIM, 0.0);
            ppGeo.getarr("perturbations", _pert, 0, SpaceDim);
            pert = RealVect(_pert);
            pout() << "\tperturbations = " << pert;
        }
        break;
    case CoordMap::LEDGE:
        {
            ppGeo.get("ledgeMapTransitionOrder", ledgeMapTransitionOrder);
            pout() << "\tledgeMapTransitionOrder = " << ledgeMapTransitionOrder << endl;

            ppGeo.get("ledgeMapHl", ledgeMapHl);
            pout() << "\tledgeMapHl = " << ledgeMapHl << endl;

            ppGeo.get("ledgeMapHr", ledgeMapHr);
            pout() << "\tledgeMapHr = " << ledgeMapHr << endl;

            ppGeo.get("ledgeMapXl", ledgeMapXl);
            pout() << "\tledgeMapXl = " << ledgeMapXl << endl;

            ppGeo.get("ledgeMapXr", ledgeMapXr);
            pout() << "\tledgeMapXr = " << ledgeMapXr << endl;
        }
        break;
    case CoordMap::BEAMGENERATOR:
        {
            ppGeo.get("alpha", beamGenMapAlpha);
            pout() << "\talpha = " << beamGenMapAlpha << endl;
            beamGenMapAlpha *= PI / 180.0;
        }
        break;
    case CoordMap::NEWBEAMGENERATOR:
        {
            ppGeo.get("alpha", beamGenMapAlpha);
            pout() << "\talpha = " << beamGenMapAlpha << endl;
            beamGenMapAlpha *= PI / 180.0;
        }
        break;
    case CoordMap::DEMMAP:
        {
            ppGeo.get("DemFile", demFile);
            pout() << "\tDEM file = " << demFile << endl;

            ppGeo.get("Interpolation_Order", interpOrder);
            pout() << "\tInterpolation_Order = " << interpOrder << endl;
        }
        break;
    }

    pout() << endl;
}


// -----------------------------------------------------------------------------
// Use this to request a new GeoSourceInterface cast from the
// appropriate child.
// -----------------------------------------------------------------------------
GeoSourceInterface* ProblemContext::newGeoSourceInterface () const
{
    CH_assert(s_geometryParamsRead);

    GeoSourceInterface* geoSourcePtr = NULL;

    switch (coordMap) {
    case ProblemContext::CoordMap::UNDEFINED:
        MayDay::Error("ProblemContext::newGeoSourceInterface received "
                      "coordMap = UNDEFINED");
        break;
    case ProblemContext::CoordMap::CARTESIAN:
        geoSourcePtr = new CartesianMap;
        break;
    case ProblemContext::CoordMap::TWISTED:
        geoSourcePtr = new TwistedMap;
        break;
    case ProblemContext::CoordMap::BEAMGENERATOR:
        geoSourcePtr = new BeamGeneratorMap;
        break;
    case ProblemContext::CoordMap::CYLINDRICAL:
        geoSourcePtr = new CylindricalMap;
        break;
    case ProblemContext::CoordMap::LEDGE:
        geoSourcePtr = new LedgeMap;
        break;
    case ProblemContext::CoordMap::NEWBEAMGENERATOR:
        geoSourcePtr = new NewBeamGeneratorMap;
        break;
    case ProblemContext::CoordMap::DEMMAP:
        geoSourcePtr = new DEMMap;
        break;
    default:
        MayDay::Error("ProblemContext::newGeoSourceInterface received "
                      "an invalid s_coordMap");
    }

    return geoSourcePtr;
}


// -----------------------------------------------------------------------------
// Read the plot.* inputs.
// -----------------------------------------------------------------------------
void ProblemContext::readPlot ()
{
    if (s_plotParamsRead) return;
    s_plotParamsRead = true;

    // The levgeo is needed by write_displacement.
    this->readGeometry();

    ParmParse ppPlot("plot");
    pout() << "ProblemContext::readPlot:" << endl;

    bool plotScheduled = false;
    plot_interval = -1;
    plot_period = -1.0;
    plot_prefix = std::string("plot_");

    if (ppPlot.query("plot_interval", plot_interval)) {
        plotScheduled = true;
        pout() << "\tplot_interval = " << plot_interval << endl;
    }
    if (ppPlot.query("plot_period", plot_period)) {
        plotScheduled = true;
        pout() << "\tplot_period = " << plot_period << endl;
    }
    if (plotScheduled) {
        ppPlot.query("plot_prefix", plot_prefix);
        pout() << "\tplot_prefix = " << plot_prefix << std::endl;
    } else {
        MayDay::Warning("No plots scheduled");
    }

    bool checkpointScheduled = false;
    checkpoint_interval = -1;
    check_prefix = std::string("chkpt_");

    if (ppPlot.query("checkpoint_interval", checkpoint_interval)) {
        checkpointScheduled = true;
        pout() << "\tcheckpoint_interval = " << checkpoint_interval << std::endl;
    }
    if (checkpointScheduled) {
        ppPlot.query("checkpoint_prefix", check_prefix);
        pout() << "\tcheckpoint_prefix = " << check_prefix << std::endl;
    } else {
        MayDay::Warning("No checkpoints scheduled");
    }

    if (!plotScheduled && !checkpointScheduled) {
        MayDay::Error("You must schedule a plot or a checkpoint");
    }

    write_divergence = false;
    ppPlot.query("writeDivergence", write_divergence);
    pout() << "\twrite_divergence = " << write_divergence << endl;

    write_lambda = false;
    ppPlot.query("writeLambda", write_lambda);
    pout() << "\twrite_lambda = " << write_lambda << endl;

    write_grad_eLambda = false;
    ppPlot.query("writeGradELambda", write_grad_eLambda);
    pout() << "\twrite_grad_eLambda = " << write_grad_eLambda << endl;

    write_pressure = true;
    ppPlot.query("writePressure", write_pressure);
    pout() << "\twrite_pressure = " << write_pressure << endl;

    write_vorticity = false;
    ppPlot.query("writeVorticity", write_vorticity);
    pout() << "\twrite_vorticity = " << write_vorticity << endl;

    write_streamfunction = false;
    ppPlot.query("writeStreamfunction", write_streamfunction);
    pout() << "\twrite_streamfunction = " << write_streamfunction << endl;

    write_scalars = true;
    ppPlot.query("writeScalars", write_scalars);
    pout() << "\twrite_scalars = " << write_scalars << endl;

    write_scalarsMinusBackground = false;
    {
        int useBackgroundScalar = 0;
        ParmParse("ibc").query("useBackgroundScalar", useBackgroundScalar);
        if (useBackgroundScalar) {
            write_scalarsMinusBackground = true;
        }
    }
    ppPlot.query("writeScalarsMinusBackground", write_scalarsMinusBackground);
    pout() << "\twrite_scalarsMinusBackground = " << write_scalarsMinusBackground << endl;

    write_proc_ids = false;
    ppPlot.query("writeProcIDs", write_proc_ids);
    pout() << "\twrite_proc_ids = " << write_proc_ids << endl;

    write_level_ids = false;
    ppPlot.query("writeLevelIDs", write_level_ids);
    pout() << "\twrite_level_ids = " << write_level_ids << endl;

    write_grids = false;
    ppPlot.query("writeGrids", write_grids);
    pout() << "\twrite_grids = " << write_grids << endl;

    write_displacement = true;
    ppPlot.query("writeDisplacement", write_displacement);
    pout() << "\twrite_displacement = " << write_displacement << endl;

    write_geometry = false;
    ppPlot.query("writeGeometry", write_geometry);
    pout() << "\twrite_geometry = " << write_geometry << endl;

    pout() << endl;
}


// -----------------------------------------------------------------------------
// The ibc.* parameters.
// -----------------------------------------------------------------------------
void ProblemContext::readIBC ()
{
    if (s_IBCParamsRead) return;
    s_IBCParamsRead = true;

    // We need to know about the domain length.
    this->readAMR();

    ParmParse ppIBC("ibc");
    pout() << "ProblemContext::readIBC:" << endl;

    ppIBC.get("problem", problem);
    pout() << "\tproblem = " << problem << endl;
    CH_assert(0 <= problem);
    CH_assert(problem < ProblemType::_NUM_PROBLEM_TYPES);

    ppIBC.get("useBackgroundScalar", useBackgroundScalar);
    pout() << "\tuseBackgroundScalar = " << useBackgroundScalar << endl;

    // Should we use a sponge layer?
    useSpongeLayer = false;
    ppIBC.query("useSpongeLayer", useSpongeLayer);
    pout() << "\tuseSpongeLayer = " << useSpongeLayer << endl;

    if (useSpongeLayer) {
        Vector<Real> vreal;

        // Lo side widths
        vreal = Vector<Real>(SpaceDim, 0.0);
        if (ppIBC.queryarr("spongeWidthLo", vreal, 0, SpaceDim)) {
            D_TERM(spongeWidth[0][0] = vreal[0];,
                   spongeWidth[1][0] = vreal[1];,
                   spongeWidth[2][0] = vreal[2];)
        } else if (ppIBC.queryarr("spongeWidthFracLo", vreal, 0, SpaceDim)) {
            D_TERM(spongeWidth[0][0] = vreal[0] * domainLength[0];,
                   spongeWidth[1][0] = vreal[1] * domainLength[1];,
                   spongeWidth[2][0] = vreal[2] * domainLength[2];)
        } else {
            MayDay::Error("You must specify ibc.spongeWidthLo or ibc.spongeWidthFracLo");
        }
        pout() << "\tspongeWidth (Lo side) = "
        D_TERM(<< "("  << spongeWidth[0][0],
               << ", " << spongeWidth[1][0],
               << ", " << spongeWidth[2][0])
               << ")" << endl;

        // Hi side widths
        vreal = Vector<Real>(SpaceDim, 0.0);
        if (ppIBC.queryarr("spongeWidthHi", vreal, 0, SpaceDim)) {
            D_TERM(spongeWidth[0][1] = vreal[0];,
                   spongeWidth[1][1] = vreal[1];,
                   spongeWidth[2][1] = vreal[2];)
        } else if (ppIBC.queryarr("spongeWidthFracHi", vreal, 0, SpaceDim)) {
            D_TERM(spongeWidth[0][1] = vreal[0] * domainLength[0];,
                   spongeWidth[1][1] = vreal[1] * domainLength[1];,
                   spongeWidth[2][1] = vreal[2] * domainLength[2];)
        } else {
            MayDay::Error("You must specify ibc.spongeWidthHi or ibc.spongeWidthFracHi");
        }
        pout() << "\tspongeWidth (Hi side) = "
        D_TERM(<< "("  << spongeWidth[0][1],
               << ", " << spongeWidth[1][1],
               << ", " << spongeWidth[2][1])
               << ")" << endl;

        // Lo side multipliers
        vreal = Vector<Real>(SpaceDim, 0.0);
        ppIBC.getarr("spongeDtMultLo", vreal, 0, SpaceDim);
        D_TERM(spongeDtMult[0][0] = vreal[0];,
               spongeDtMult[1][0] = vreal[1];,
               spongeDtMult[2][0] = vreal[2];)
        pout() << "\tspongeDtMult (Lo side) = "
        D_TERM(<< "("  << spongeDtMult[0][0],
               << ", " << spongeDtMult[1][0],
               << ", " << spongeDtMult[2][0])
               << ")" << endl;

        // Hi side multipliers
        vreal = Vector<Real>(SpaceDim, 0.0);
        ppIBC.getarr("spongeDtMultHi", vreal, 0, SpaceDim);
        D_TERM(spongeDtMult[0][1] = vreal[0];,
               spongeDtMult[1][1] = vreal[1];,
               spongeDtMult[2][1] = vreal[2];)
        pout() << "\tspongeDtMult (Hi side) = "
        D_TERM(<< "("  << spongeDtMult[0][1],
               << ", " << spongeDtMult[1][1],
               << ", " << spongeDtMult[2][1])
               << ")" << endl;
    }

    tidalOmega = 0.0;
    ppIBC.query("tidalOmega", tidalOmega);
    pout() << "\ttidalOmega = " << tidalOmega << endl;

    tidalU0 = 0.0;
    ppIBC.query("tidalU0", tidalU0);
    pout() << "\ttidalU0 = " << tidalU0 << endl;


    // Add problem-specific stuff here ----------------------
    switch (problem) {
    case ProblemType::VORTEX_STREET:
        {
            Vector<Real> vreal(SpaceDim, 0.0);
            ppIBC.getarr("inflowVel", vreal, 0, SpaceDim);
            inflowVel = RealVect(vreal);
            pout() << "\tinflowVel = " << inflowVel << endl;
        }
        break;
    case ProblemType::HORIZ_CONV:
        // TODO
        break;
    case ProblemType::SOLITARYWAVE:
        {
            ppIBC.get("solitaryWave_xcenter", solitaryWave_xcenter);
            pout() << "\tsolitaryWave_xcenter = " << solitaryWave_xcenter << endl;

            ppIBC.get("solitaryWave_H2", solitaryWave_H2);
            pout() << "\tsolitaryWave_H2 = " << solitaryWave_H2 << endl;

            ppIBC.get("solitaryWave_amp", solitaryWave_amp);
            pout() << "\tsolitaryWave_amp = " << solitaryWave_amp << endl;

            solitaryWave_rho0 = 0.0;
            ppIBC.query("solitaryWave_rho0", solitaryWave_rho0);
            pout() << "\tsolitaryWave_rho0 = " << solitaryWave_rho0 << endl;

            ppIBC.get("solitaryWave_drho", solitaryWave_drho);
            pout() << "\tsolitaryWave_drho = " << solitaryWave_drho << endl;

            ppIBC.get("solitaryWave_delta0", solitaryWave_delta0);
            pout() << "\tsolitaryWave_delta0 = " << solitaryWave_delta0 << endl;
        }
        break;
    }
    // End problem-specific stuff ---------------------------

    pout() << endl;
}


// -----------------------------------------------------------------------------
// Use this to request a new (undefined) PhysBCUtil cast from the
// appropriate child.
// -----------------------------------------------------------------------------
PhysBCUtil* ProblemContext::newPhysBCUtil () const
{
    CH_assert(s_IBCParamsRead);

    PhysBCUtil* physBCPtr = NULL;

    switch (problem) {
    case ProblemType::ADVECTION_TEST:
        physBCPtr = new AdvectionTestBCUtil;
        break;
    case ProblemType::LOCK_EXCHANGE:
        physBCPtr = new LockExchangeBCUtil;
        break;
    case ProblemType::BEAM_GENERATION:
        physBCPtr = new BeamGenerationBCUtil;
        break;
    case ProblemType::INTERNAL_WAVE:
        physBCPtr = new InternalWaveBCUtil;
        break;
    case ProblemType::TAYLOR_GREEN:
        physBCPtr = new TaylorGreenBCUtil;
        break;
    case ProblemType::VORTEX_STREET:
        physBCPtr = new VortexStreetBCUtil;
        break;
    case ProblemType::HORIZ_CONV:
        physBCPtr = new HorizConvBCUtil;
        break;
    case ProblemType::SOLITARYWAVE:
        physBCPtr = new SolitaryWaveBCUtil;
        break;
    case ProblemType::DJL:
        physBCPtr = new DJLBCUtil;
        break;
    default:
        // Undefined problem
        MayDay::Error("Bad problem type");
    }

    return physBCPtr;
}


// -----------------------------------------------------------------------------
// The advection.* parameters
// -----------------------------------------------------------------------------
void ProblemContext::readAdvection ()
{
    if (s_advectionParamsRead) return;
    s_advectionParamsRead = true;

    ParmParse ppAdvection("advection");
    pout() << "ProblemContext::readAdvection:" << endl;

    updateScheme = UpdateScheme::FiniteVolume;
    ppAdvection.query("updateScheme", updateScheme);
    pout() << "\tupdateScheme = " << updateScheme << endl;


    // Velocity...
    normalPredOrderVel = 2;
    ppAdvection.query("normalPredOrderVel", normalPredOrderVel);
    pout() << "\tnormalPredOrderVel = " << normalPredOrderVel << endl;

    useFourthOrderSlopesVel = true;
    ppAdvection.query("useFourthOrderSlopesVel", useFourthOrderSlopesVel);
    pout() << "\tuseFourthOrderSlopesVel = " << useFourthOrderSlopesVel << endl;

    useLimitingVel = false;
    ppAdvection.query("useLimitingVel", useLimitingVel);
    pout() << "\tuseLimitingVel = " << useLimitingVel << endl;

    useHighOrderLimiterVel = false;
    ppAdvection.query("useHighOrderLimiterVel", useHighOrderLimiterVel);
    pout() << "\tuseHighOrderLimiterVel = " << useHighOrderLimiterVel << endl;

    useUpwindingVel = true;
    ppAdvection.query("useUpwindingVel", useUpwindingVel);
    pout() << "\tuseUpwindingVel = " << useUpwindingVel << endl;

    nonlinearDifferencingForm = NonlinearDifferencingForm::CONSERVATIVE;
    ppAdvection.query("nonlinearDifferencingForm", nonlinearDifferencingForm);
    pout() << "\tnonlinearDifferencingForm = " << nonlinearDifferencingForm << endl;


    // Scalars...
    normalPredOrderScal = 2;
    ppAdvection.query("normalPredOrderScal", normalPredOrderScal);
    pout() << "\tnormalPredOrderScal = " << normalPredOrderScal << endl;

    useFourthOrderSlopesScal = true;
    ppAdvection.query("useFourthOrderSlopesScal", useFourthOrderSlopesScal);
    pout() << "\tuseFourthOrderSlopesScal = " << useFourthOrderSlopesScal << endl;

    useLimitingScal = true;
    ppAdvection.query("useLimitingScal", useLimitingScal);
    pout() << "\tuseLimitingScal = " << useLimitingScal << endl;

    useHighOrderLimiterScal = true;
    ppAdvection.query("useHighOrderLimiterScal", useHighOrderLimiterScal);
    pout() << "\tuseHighOrderLimiterScal = " << useHighOrderLimiterScal << endl;

    useUpwindingScal = true;
    ppAdvection.query("useUpwindingScal", useUpwindingScal);
    pout() << "\tuseUpwindingScal = " << useUpwindingScal << endl;

    pout() << endl;
}


// -----------------------------------------------------------------------------
// The solver parameters. This includes the projector, viscous solver, etc.
// -----------------------------------------------------------------------------
void ProblemContext::readSolver ()
{
    if (s_solverParamsRead) return;
    s_solverParamsRead = true;

    pout() << "ProblemContext::readSolver:" << endl;

    {   // Solver settings
        ParmParse ppAMRMG("AMRMG");

        AMRMG_eps = 1e-6;
        ppAMRMG.query("eps", AMRMG_eps);
        pout() << "\tAMRMG.eps = " << AMRMG_eps << endl;

        AMRMG_num_smooth_down = 2;
        ppAMRMG.query("num_smooth_down", AMRMG_num_smooth_down);
        pout() << "\tAMRMG.num_smooth_down = " << AMRMG_num_smooth_down << endl;

        AMRMG_num_smooth_up = 2;
        ppAMRMG.query("num_smooth_up", AMRMG_num_smooth_up);
        pout() << "\tAMRMG.num_smooth_up = " << AMRMG_num_smooth_up << endl;

        AMRMG_num_smooth_bottom = 2;
        ppAMRMG.query("num_smooth_bottom", AMRMG_num_smooth_bottom);
        pout() << "\tAMRMG.num_smooth_bottom = " << AMRMG_num_smooth_bottom << endl;

        AMRMG_num_smooth_precond = 2;
        ppAMRMG.query("num_smooth_precond", AMRMG_num_smooth_precond);
        pout() << "\tAMRMG.num_smooth_precond = " << AMRMG_num_smooth_precond << endl;

        AMRMG_numMG = 1;
        ppAMRMG.query("numMG", AMRMG_numMG);
        pout() << "\tAMRMG.numMG = " << AMRMG_numMG << endl;

        AMRMG_imin = 5;
        ppAMRMG.query("imin", AMRMG_imin);
        pout() << "\tAMRMG.imin = " << AMRMG_imin << endl;

        AMRMG_imax = 20;
        ppAMRMG.query("imax", AMRMG_imax);
        pout() << "\tAMRMG.imax = " << AMRMG_imax << endl;

        AMRMG_hang = 1e-15;
        ppAMRMG.query("hang", AMRMG_hang);
        pout() << "\tAMRMG.hang = " << AMRMG_hang << endl;

        AMRMG_normThresh = 1e-30;
        ppAMRMG.query("normThresh", AMRMG_normThresh);
        pout() << "\tAMRMG.normThresh = " << AMRMG_normThresh << endl;

        AMRMG_maxDepth = -1;
        ppAMRMG.query("maxDepth", AMRMG_maxDepth);
        pout() << "\tAMRMG.maxDepth = " << AMRMG_maxDepth << endl;

        AMRMG_verbosity = verbosity;
        ppAMRMG.query("verbosity", AMRMG_verbosity);
        pout() << "\tAMRMG.verbosity = " << AMRMG_verbosity << endl;

        AMRMG_relaxMode = RelaxMode::LEVEL_GSRB;
        ppAMRMG.query("relax_mode", AMRMG_relaxMode);
        pout() << "\tAMRMG.relaxMode = " << AMRMG_relaxMode << endl;

        AMRMG_precondMode = PrecondMode::DiagRelax;
        ppAMRMG.query("precond_mode", AMRMG_precondMode);
        pout() << "\tAMRMG.precondMode = " << AMRMG_precondMode << endl;

        ParmParse ppBottom("bottom");

        bottom_eps = 1e-6;
        ppBottom.query("eps", bottom_eps);
        pout() << "\tbottom.eps = " << bottom_eps << endl;

        bottom_reps = 1e-12;
        ppBottom.query("reps", bottom_reps);
        pout() << "\tbottom.reps = " << bottom_reps << endl;

        bottom_imax = 80;
        ppBottom.query("imax", bottom_imax);
        pout() << "\tbottom.imax = " << bottom_imax << endl;

        bottom_numRestarts = 5;
        ppBottom.query("numRestarts", bottom_numRestarts);
        pout() << "\tbottom.numRestarts = " << bottom_numRestarts << endl;

        bottom_hang = 1e-15;
        ppBottom.query("hang", bottom_hang);
        pout() << "\tbottom.hang = " << bottom_hang << endl;

        bottom_small = 1e-30;
        ppBottom.query("small", bottom_small);
        pout() << "\tbottom.small = " << bottom_small << endl;

        bottom_normType = 2;
        ppBottom.query("normType", bottom_normType);
        pout() << "\tbottom.normType = " << bottom_normType << endl;

        bottom_verbosity = verbosity;
        ppBottom.query("verbosity", bottom_verbosity);
        pout() << "\tbottom.verbosity = " << bottom_verbosity << endl;

        pout() << endl;
    }

    {   // Viscous solver overrides
        ParmParse ppAMRMG("viscous_AMRMG");

        viscous_AMRMG_eps = AMRMG_eps;
        ppAMRMG.query("eps", viscous_AMRMG_eps);
        pout() << "\tviscous_AMRMG.eps = " << viscous_AMRMG_eps << endl;

        viscous_AMRMG_num_smooth_down = AMRMG_num_smooth_down;
        ppAMRMG.query("num_smooth_down", viscous_AMRMG_num_smooth_down);
        pout() << "\tviscous_AMRMG.num_smooth_down = " << viscous_AMRMG_num_smooth_down << endl;

        viscous_AMRMG_num_smooth_up = AMRMG_num_smooth_up;
        ppAMRMG.query("num_smooth_up", viscous_AMRMG_num_smooth_up);
        pout() << "\tviscous_AMRMG.num_smooth_up = " << viscous_AMRMG_num_smooth_up << endl;

        viscous_AMRMG_num_smooth_bottom = AMRMG_num_smooth_bottom;
        ppAMRMG.query("num_smooth_bottom", viscous_AMRMG_num_smooth_bottom);
        pout() << "\tviscous_AMRMG.num_smooth_bottom = " << viscous_AMRMG_num_smooth_bottom << endl;

        viscous_AMRMG_num_smooth_precond = AMRMG_num_smooth_precond;
        ppAMRMG.query("num_smooth_precond", viscous_AMRMG_num_smooth_precond);
        pout() << "\tviscous_AMRMG.num_smooth_precond = " << viscous_AMRMG_num_smooth_precond << endl;

        viscous_AMRMG_numMG = AMRMG_numMG;
        ppAMRMG.query("numMG", viscous_AMRMG_numMG);
        pout() << "\tviscous_AMRMG.numMG = " << viscous_AMRMG_numMG << endl;

        viscous_AMRMG_imin = AMRMG_imin;
        ppAMRMG.query("imin", viscous_AMRMG_imin);
        pout() << "\tviscous_AMRMG.imin = " << viscous_AMRMG_imin << endl;

        viscous_AMRMG_imax = AMRMG_imax;
        ppAMRMG.query("imax", viscous_AMRMG_imax);
        pout() << "\tviscous_AMRMG.imax = " << viscous_AMRMG_imax << endl;

        viscous_AMRMG_hang = AMRMG_hang;
        ppAMRMG.query("hang", viscous_AMRMG_hang);
        pout() << "\tviscous_AMRMG.hang = " << viscous_AMRMG_hang << endl;

        viscous_AMRMG_normThresh = AMRMG_normThresh;
        ppAMRMG.query("normThresh", viscous_AMRMG_normThresh);
        pout() << "\tviscous_AMRMG.normThresh = " << viscous_AMRMG_normThresh << endl;

        viscous_AMRMG_maxDepth = AMRMG_maxDepth;
        ppAMRMG.query("maxDepth", viscous_AMRMG_maxDepth);
        pout() << "\tviscous_AMRMG.maxDepth = " << viscous_AMRMG_maxDepth << endl;

        viscous_AMRMG_verbosity = AMRMG_verbosity;
        ppAMRMG.query("verbosity", viscous_AMRMG_verbosity);
        pout() << "\tviscous_AMRMG.verbosity = " << viscous_AMRMG_verbosity << endl;

        viscous_AMRMG_relaxMode = AMRMG_relaxMode;
        ppAMRMG.query("relax_mode", viscous_AMRMG_relaxMode);
        pout() << "\tviscous_AMRMG.relax_mode = " << viscous_AMRMG_relaxMode << endl;

        viscous_AMRMG_precondMode = AMRMG_precondMode;
        ppAMRMG.query("precond_mode", viscous_AMRMG_precondMode);
        pout() << "\tviscous_AMRMG.precondMode = " << viscous_AMRMG_precondMode << endl;


        ParmParse ppBottom("viscous_bottom");

        viscous_bottom_eps = bottom_eps;
        ppBottom.query("eps", viscous_bottom_eps);
        pout() << "\tviscous_bottom.eps = " << viscous_bottom_eps << endl;

        viscous_bottom_reps = bottom_reps;
        ppBottom.query("reps", viscous_bottom_reps);
        pout() << "\tviscous_bottom.reps = " << viscous_bottom_reps << endl;

        viscous_bottom_imax = bottom_imax;
        ppBottom.query("imax", viscous_bottom_imax);
        pout() << "\tviscous_bottom.imax = " << viscous_bottom_imax << endl;

        viscous_bottom_numRestarts = bottom_numRestarts;
        ppBottom.query("numRestarts", viscous_bottom_numRestarts);
        pout() << "\tviscous_bottom.numRestarts = " << viscous_bottom_numRestarts << endl;

        viscous_bottom_hang = bottom_hang;
        ppBottom.query("hang", viscous_bottom_hang);
        pout() << "\tviscous_bottom.hang = " << viscous_bottom_hang << endl;

        viscous_bottom_small = bottom_small;
        ppBottom.query("small", viscous_bottom_small);
        pout() << "\tviscous_bottom.small = " << viscous_bottom_small << endl;

        viscous_bottom_normType = bottom_normType;
        ppBottom.query("normType", viscous_bottom_normType);
        pout() << "\tviscous_bottom.normType = " << viscous_bottom_normType << endl;

        viscous_bottom_verbosity = bottom_verbosity;
        ppBottom.query("verbosity", viscous_bottom_verbosity);
        pout() << "\tviscous_bottom.verbosity = " << viscous_bottom_verbosity << endl;

        pout() << endl;
    }
}


// -----------------------------------------------------------------------------
// The projection.* parameters
// -----------------------------------------------------------------------------
void ProblemContext::readProjection ()
{
    if (s_projectionParamsRead) return;
    s_projectionParamsRead = true;

    // We use the base solver settings as defaults.
    this->readSolver();

    pout() << "ProblemContext::readProjection:" << endl;
    ParmParse ppProj("projection");

    isIncompressible = true;
    ppProj.query("isIncompressible", isIncompressible);
    pout() << "\tisIncompressible = " << isIncompressible << endl;

    if (isIncompressible) {
        initial_projection_iters = 1;
        ppProj.query("initial_projection_iters", initial_projection_iters);
        pout() << "\tinitial_projection_iters = " << initial_projection_iters << endl;

        initial_pressure_iters = 1;
        ppProj.query("initial_pressure_iters", initial_pressure_iters);
        pout() << "\tinitial_pressure_iters = " << initial_pressure_iters << endl;

        level_projection_iters = 1;
        ppProj.query("level_projection_iters", level_projection_iters);
        pout() << "\tlevel_projection_iters = " << level_projection_iters << endl;

        doSyncProjection = true;
        ppProj.query("doSyncProjection", doSyncProjection);
        pout() << "\tdoSyncProjection = " << doSyncProjection << endl;

        sync_projection_iters = 1;
        if (doSyncProjection) {
            ppProj.query("sync_projection_iters", sync_projection_iters);
            pout() << "\tsync_projection_iters = " << sync_projection_iters << endl;
        } else {
            sync_projection_iters = 0;
        }

        applyVDCorrection = false;
        ppProj.query("applyVDCorrection", applyVDCorrection);
        pout() << "\tapplyVDCorrection = " << applyVDCorrection << endl;

        etaLambda = 0.0;
        if (applyVDCorrection) {
            if (ppProj.query("eta", etaLambda)) {
                MayDay::Warning("projection.eta is deprecated. Re-interpreting as etaLambda");
            } else {
                ppProj.get("etaLambda", etaLambda);
            }
            pout() << "\tetaLambda = " << etaLambda << endl;
        } else {
            if (ppProj.query("eta", etaLambda)) {
                MayDay::Warning("projection.eta is deprecated. Re-interpreting as etaLambda");
            } else {
                ppProj.query("etaLambda", etaLambda);
            }
            pout() << "\tetaLambda = " << etaLambda << endl;
        }

        applySyncCorrection = true;
        ppProj.query("applySyncCorrection", applySyncCorrection);
        pout() << "\tapplySyncCorrection = " << applySyncCorrection << endl;

    } else {
        initial_projection_iters = 0;
        pout() << "\tinitial_projection_iters = " << initial_projection_iters << endl;

        initial_pressure_iters = 0;
        pout() << "\tinitial_pressure_iters = " << initial_pressure_iters << endl;

        level_projection_iters = 0;
        pout() << "\tlevel_projection_iters = " << level_projection_iters << endl;

        doSyncProjection = false;
        pout() << "\tdoSyncProjection = " << doSyncProjection << endl;

        sync_projection_iters = 0;
        pout() << "\tsync_projection_iters = " << sync_projection_iters << endl;

        applyVDCorrection = false;
        etaLambda = 0.0;
        pout() << "\tapplyVDCorrection = " << applyVDCorrection << endl;

        applySyncCorrection = false;
        pout() << "\tapplySyncCorrection = " << applySyncCorrection << endl;
    }

    // Debug. It looks like some input files use applyFreestreamCorrection. Oops.
    if (ppProj.contains("applyFreestreamCorrection")) {
        MayDay::Error("applyFreestreamCorrection is outdated. Use applyVDCorrection");
    }
    pout() << endl;


    // Exact projection solver overrides (For FC advecting velocities)
    if (isIncompressible) {
        ParmParse ppAMRMG("MACprojection_AMRMG");

        MACprojection_AMRMG_eps = AMRMG_eps;
        ppAMRMG.query("eps", MACprojection_AMRMG_eps);
        pout() << "\tMACprojection_AMRMG.eps = " << MACprojection_AMRMG_eps << endl;

        MACprojection_AMRMG_num_smooth_down = AMRMG_num_smooth_down;
        ppAMRMG.query("num_smooth_down", MACprojection_AMRMG_num_smooth_down);
        pout() << "\tMACprojection_AMRMG.num_smooth_down = " << MACprojection_AMRMG_num_smooth_down << endl;

        MACprojection_AMRMG_num_smooth_up = AMRMG_num_smooth_up;
        ppAMRMG.query("num_smooth_up", MACprojection_AMRMG_num_smooth_up);
        pout() << "\tMACprojection_AMRMG.num_smooth_up = " << MACprojection_AMRMG_num_smooth_up << endl;

        MACprojection_AMRMG_num_smooth_bottom = AMRMG_num_smooth_bottom;
        ppAMRMG.query("num_smooth_bottom", MACprojection_AMRMG_num_smooth_bottom);
        pout() << "\tMACprojection_AMRMG.num_smooth_bottom = " << MACprojection_AMRMG_num_smooth_bottom << endl;

        MACprojection_AMRMG_num_smooth_precond = AMRMG_num_smooth_precond;
        ppAMRMG.query("num_smooth_precond", MACprojection_AMRMG_num_smooth_precond);
        pout() << "\tMACprojection_AMRMG.num_smooth_precond = " << MACprojection_AMRMG_num_smooth_precond << endl;

        MACprojection_AMRMG_numMG = AMRMG_numMG;
        ppAMRMG.query("numMG", MACprojection_AMRMG_numMG);
        pout() << "\tMACprojection_AMRMG.numMG = " << MACprojection_AMRMG_numMG << endl;

        MACprojection_AMRMG_imin = AMRMG_imin;
        ppAMRMG.query("imin", MACprojection_AMRMG_imin);
        pout() << "\tMACprojection_AMRMG.imin = " << MACprojection_AMRMG_imin << endl;

        MACprojection_AMRMG_imax = AMRMG_imax;
        ppAMRMG.query("imax", MACprojection_AMRMG_imax);
        pout() << "\tMACprojection_AMRMG.imax = " << MACprojection_AMRMG_imax << endl;

        MACprojection_AMRMG_hang = AMRMG_hang;
        ppAMRMG.query("hang", MACprojection_AMRMG_hang);
        pout() << "\tMACprojection_AMRMG.hang = " << MACprojection_AMRMG_hang << endl;

        MACprojection_AMRMG_normThresh = AMRMG_normThresh;
        ppAMRMG.query("normThresh", MACprojection_AMRMG_normThresh);
        pout() << "\tMACprojection_AMRMG.normThresh = " << MACprojection_AMRMG_normThresh << endl;

        MACprojection_AMRMG_maxDepth = AMRMG_maxDepth;
        ppAMRMG.query("maxDepth", MACprojection_AMRMG_maxDepth);
        pout() << "\tMACprojection_AMRMG.maxDepth = " << MACprojection_AMRMG_maxDepth << endl;

        MACprojection_AMRMG_relaxMode = AMRMG_relaxMode;
        ppAMRMG.query("relax_mode", MACprojection_AMRMG_relaxMode);
        pout() << "\tMACprojection_AMRMG.relax_mode = " << MACprojection_AMRMG_relaxMode << endl;

        MACprojection_AMRMG_precondMode = AMRMG_precondMode;
        ppAMRMG.query("precond_mode", MACprojection_AMRMG_precondMode);
        pout() << "\tMACprojection_AMRMG.precondMode = " << MACprojection_AMRMG_precondMode << endl;

        MACprojection_AMRMG_verbosity = AMRMG_verbosity;
        ppAMRMG.query("verbosity", MACprojection_AMRMG_verbosity);
        pout() << "\tMACprojection_AMRMG.verbosity = " << MACprojection_AMRMG_verbosity << endl;


        ParmParse ppBottom("MACprojection_bottom");

        MACprojection_bottom_eps = bottom_eps;
        ppBottom.query("eps", MACprojection_bottom_eps);
        pout() << "\tMACprojection_bottom.eps = " << MACprojection_bottom_eps << endl;

        MACprojection_bottom_reps = bottom_reps;
        ppBottom.query("reps", MACprojection_bottom_reps);
        pout() << "\tMACprojection_bottom.reps = " << MACprojection_bottom_reps << endl;

        MACprojection_bottom_imax = bottom_imax;
        ppBottom.query("imax", MACprojection_bottom_imax);
        pout() << "\tMACprojection_bottom.imax = " << MACprojection_bottom_imax << endl;

        MACprojection_bottom_numRestarts = bottom_numRestarts;
        ppBottom.query("numRestarts", MACprojection_bottom_numRestarts);
        pout() << "\tMACprojection_bottom.numRestarts = " << MACprojection_bottom_numRestarts << endl;

        MACprojection_bottom_hang = bottom_hang;
        ppBottom.query("hang", MACprojection_bottom_hang);
        pout() << "\tMACprojection_bottom.hang = " << MACprojection_bottom_hang << endl;

        MACprojection_bottom_small = bottom_small;
        ppBottom.query("small", MACprojection_bottom_small);
        pout() << "\tMACprojection_bottom.small = " << MACprojection_bottom_small << endl;

        MACprojection_bottom_normType = bottom_normType;
        ppBottom.query("normType", MACprojection_bottom_normType);
        pout() << "\tMACprojection_bottom.normType = " << MACprojection_bottom_normType << endl;

        MACprojection_bottom_verbosity = bottom_verbosity;
        ppBottom.query("verbosity", MACprojection_bottom_verbosity);
        pout() << "\tMACprojection_bottom.verbosity = " << MACprojection_bottom_verbosity << endl;

        pout() << endl;
    }

    // Approximate projection overrides (For CC velocities)
    if (isIncompressible) {
        ParmParse ppAMRMG("CCprojection_AMRMG");

        CCprojection_AMRMG_eps = AMRMG_eps;
        ppAMRMG.query("eps", CCprojection_AMRMG_eps);
        pout() << "\tCCprojection_AMRMG.eps = " << CCprojection_AMRMG_eps << endl;

        CCprojection_AMRMG_num_smooth_down = AMRMG_num_smooth_down;
        ppAMRMG.query("num_smooth_down", CCprojection_AMRMG_num_smooth_down);
        pout() << "\tCCprojection_AMRMG.num_smooth_down = " << CCprojection_AMRMG_num_smooth_down << endl;

        CCprojection_AMRMG_num_smooth_up = AMRMG_num_smooth_up;
        ppAMRMG.query("num_smooth_up", CCprojection_AMRMG_num_smooth_up);
        pout() << "\tCCprojection_AMRMG.num_smooth_up = " << CCprojection_AMRMG_num_smooth_up << endl;

        CCprojection_AMRMG_num_smooth_bottom = AMRMG_num_smooth_bottom;
        ppAMRMG.query("num_smooth_bottom", CCprojection_AMRMG_num_smooth_bottom);
        pout() << "\tCCprojection_AMRMG.num_smooth_bottom = " << CCprojection_AMRMG_num_smooth_bottom << endl;

        CCprojection_AMRMG_num_smooth_precond = AMRMG_num_smooth_precond;
        ppAMRMG.query("num_smooth_precond", CCprojection_AMRMG_num_smooth_precond);
        pout() << "\tCCprojection_AMRMG.num_smooth_precond = " << CCprojection_AMRMG_num_smooth_precond << endl;

        CCprojection_AMRMG_numMG = AMRMG_numMG;
        ppAMRMG.query("numMG", CCprojection_AMRMG_numMG);
        pout() << "\tCCprojection_AMRMG.numMG = " << CCprojection_AMRMG_numMG << endl;

        CCprojection_AMRMG_imin = AMRMG_imin;
        ppAMRMG.query("imin", CCprojection_AMRMG_imin);
        pout() << "\tCCprojection_AMRMG.imin = " << CCprojection_AMRMG_imin << endl;

        CCprojection_AMRMG_imax = AMRMG_imax;
        ppAMRMG.query("imax", CCprojection_AMRMG_imax);
        pout() << "\tCCprojection_AMRMG.imax = " << CCprojection_AMRMG_imax << endl;

        CCprojection_AMRMG_hang = AMRMG_hang;
        ppAMRMG.query("hang", CCprojection_AMRMG_hang);
        pout() << "\tCCprojection_AMRMG.hang = " << CCprojection_AMRMG_hang << endl;

        CCprojection_AMRMG_normThresh = AMRMG_normThresh;
        ppAMRMG.query("normThresh", CCprojection_AMRMG_normThresh);
        pout() << "\tCCprojection_AMRMG.normThresh = " << CCprojection_AMRMG_normThresh << endl;

        CCprojection_AMRMG_maxDepth = AMRMG_maxDepth;
        ppAMRMG.query("maxDepth", CCprojection_AMRMG_maxDepth);
        pout() << "\tCCprojection_AMRMG.maxDepth = " << CCprojection_AMRMG_maxDepth << endl;

        CCprojection_AMRMG_relaxMode = AMRMG_relaxMode;
        ppAMRMG.query("relax_mode", CCprojection_AMRMG_relaxMode);
        pout() << "\tCCprojection_AMRMG.relax_mode = " << CCprojection_AMRMG_relaxMode << endl;

        CCprojection_AMRMG_precondMode = AMRMG_precondMode;
        ppAMRMG.query("precond_mode", CCprojection_AMRMG_precondMode);
        pout() << "\tCCprojection_AMRMG.precondMode = " << CCprojection_AMRMG_precondMode << endl;

        CCprojection_AMRMG_verbosity = AMRMG_verbosity;
        ppAMRMG.query("verbosity", CCprojection_AMRMG_verbosity);
        pout() << "\tCCprojection_AMRMG.verbosity = " << CCprojection_AMRMG_verbosity << endl;


        ParmParse ppBottom("CCprojection_bottom");

        CCprojection_bottom_eps = bottom_eps;
        ppBottom.query("eps", CCprojection_bottom_eps);
        pout() << "\tCCprojection_bottom.eps = " << CCprojection_bottom_eps << endl;

        CCprojection_bottom_reps = bottom_reps;
        ppBottom.query("reps", CCprojection_bottom_reps);
        pout() << "\tCCprojection_bottom.reps = " << CCprojection_bottom_reps << endl;

        CCprojection_bottom_imax = bottom_imax;
        ppBottom.query("imax", CCprojection_bottom_imax);
        pout() << "\tCCprojection_bottom.imax = " << CCprojection_bottom_imax << endl;

        CCprojection_bottom_numRestarts = bottom_numRestarts;
        ppBottom.query("numRestarts", CCprojection_bottom_numRestarts);
        pout() << "\tCCprojection_bottom.numRestarts = " << CCprojection_bottom_numRestarts << endl;

        CCprojection_bottom_hang = bottom_hang;
        ppBottom.query("hang", CCprojection_bottom_hang);
        pout() << "\tCCprojection_bottom.hang = " << CCprojection_bottom_hang << endl;

        CCprojection_bottom_small = bottom_small;
        ppBottom.query("small", CCprojection_bottom_small);
        pout() << "\tCCprojection_bottom.small = " << CCprojection_bottom_small << endl;

        CCprojection_bottom_normType = bottom_normType;
        ppBottom.query("normType", CCprojection_bottom_normType);
        pout() << "\tCCprojection_bottom.normType = " << CCprojection_bottom_normType << endl;

        CCprojection_bottom_verbosity = bottom_verbosity;
        ppBottom.query("verbosity", CCprojection_bottom_verbosity);
        pout() << "\tCCprojection_bottom.verbosity = " << CCprojection_bottom_verbosity << endl;

        pout() << endl;
    }

    // Sync (approximate) projection solver overrides (For composite CC velocities)
    if (isIncompressible) {
        ParmParse ppAMRMG("syncProjection_AMRMG");

        syncProjection_AMRMG_eps = AMRMG_eps;
        ppAMRMG.query("eps", syncProjection_AMRMG_eps);
        pout() << "\tsyncProjection_AMRMG.eps = " << syncProjection_AMRMG_eps << endl;

        syncProjection_AMRMG_num_smooth_down = AMRMG_num_smooth_down;
        ppAMRMG.query("num_smooth_down", syncProjection_AMRMG_num_smooth_down);
        pout() << "\tsyncProjection_AMRMG.num_smooth_down = " << syncProjection_AMRMG_num_smooth_down << endl;

        syncProjection_AMRMG_num_smooth_up = AMRMG_num_smooth_up;
        ppAMRMG.query("num_smooth_up", syncProjection_AMRMG_num_smooth_up);
        pout() << "\tsyncProjection_AMRMG.num_smooth_up = " << syncProjection_AMRMG_num_smooth_up << endl;

        syncProjection_AMRMG_num_smooth_bottom = AMRMG_num_smooth_bottom;
        ppAMRMG.query("num_smooth_bottom", syncProjection_AMRMG_num_smooth_bottom);
        pout() << "\tsyncProjection_AMRMG.num_smooth_bottom = " << syncProjection_AMRMG_num_smooth_bottom << endl;

        syncProjection_AMRMG_num_smooth_precond = AMRMG_num_smooth_precond;
        ppAMRMG.query("num_smooth_precond", syncProjection_AMRMG_num_smooth_precond);
        pout() << "\tsyncProjection_AMRMG.num_smooth_precond = " << syncProjection_AMRMG_num_smooth_precond << endl;

        syncProjection_AMRMG_numMG = AMRMG_numMG;
        ppAMRMG.query("numMG", syncProjection_AMRMG_numMG);
        pout() << "\tsyncProjection_AMRMG.numMG = " << syncProjection_AMRMG_numMG << endl;

        syncProjection_AMRMG_imin = AMRMG_imin;
        ppAMRMG.query("imin", syncProjection_AMRMG_imin);
        pout() << "\tsyncProjection_AMRMG.imin = " << syncProjection_AMRMG_imin << endl;

        syncProjection_AMRMG_imax = AMRMG_imax;
        ppAMRMG.query("imax", syncProjection_AMRMG_imax);
        pout() << "\tsyncProjection_AMRMG.imax = " << syncProjection_AMRMG_imax << endl;

        syncProjection_AMRMG_hang = AMRMG_hang;
        ppAMRMG.query("hang", syncProjection_AMRMG_hang);
        pout() << "\tsyncProjection_AMRMG.hang = " << syncProjection_AMRMG_hang << endl;

        syncProjection_AMRMG_normThresh = AMRMG_normThresh;
        ppAMRMG.query("normThresh", syncProjection_AMRMG_normThresh);
        pout() << "\tsyncProjection_AMRMG.normThresh = " << syncProjection_AMRMG_normThresh << endl;

        syncProjection_AMRMG_maxDepth = AMRMG_maxDepth;
        ppAMRMG.query("maxDepth", syncProjection_AMRMG_maxDepth);
        pout() << "\tsyncProjection_AMRMG.maxDepth = " << syncProjection_AMRMG_maxDepth << endl;

        syncProjection_AMRMG_relaxMode = AMRMG_relaxMode;
        ppAMRMG.query("relax_mode", syncProjection_AMRMG_relaxMode);
        pout() << "\tsyncProjection_AMRMG.relax_mode = " << syncProjection_AMRMG_relaxMode << endl;

        syncProjection_AMRMG_precondMode = AMRMG_precondMode;
        ppAMRMG.query("precond_mode", syncProjection_AMRMG_precondMode);
        pout() << "\tsyncProjection_AMRMG.precondMode = " << syncProjection_AMRMG_precondMode << endl;

        syncProjection_AMRMG_verbosity = AMRMG_verbosity;
        ppAMRMG.query("verbosity", syncProjection_AMRMG_verbosity);
        pout() << "\tsyncProjection_AMRMG.verbosity = " << syncProjection_AMRMG_verbosity << endl;


        ParmParse ppBottom("syncProjection_bottom");

        syncProjection_bottom_eps = bottom_eps;
        ppBottom.query("eps", syncProjection_bottom_eps);
        pout() << "\tsyncProjection_bottom.eps = " << syncProjection_bottom_eps << endl;

        syncProjection_bottom_reps = bottom_reps;
        ppBottom.query("reps", syncProjection_bottom_reps);
        pout() << "\tsyncProjection_bottom.reps = " << syncProjection_bottom_reps << endl;

        syncProjection_bottom_imax = bottom_imax;
        ppBottom.query("imax", syncProjection_bottom_imax);
        pout() << "\tsyncProjection_bottom.imax = " << syncProjection_bottom_imax << endl;

        syncProjection_bottom_numRestarts = bottom_numRestarts;
        ppBottom.query("numRestarts", syncProjection_bottom_numRestarts);
        pout() << "\tsyncProjection_bottom.numRestarts = " << syncProjection_bottom_numRestarts << endl;

        syncProjection_bottom_hang = bottom_hang;
        ppBottom.query("hang", syncProjection_bottom_hang);
        pout() << "\tsyncProjection_bottom.hang = " << syncProjection_bottom_hang << endl;

        syncProjection_bottom_small = bottom_small;
        ppBottom.query("small", syncProjection_bottom_small);
        pout() << "\tsyncProjection_bottom.small = " << syncProjection_bottom_small << endl;

        syncProjection_bottom_normType = bottom_normType;
        ppBottom.query("normType", syncProjection_bottom_normType);
        pout() << "\tsyncProjection_bottom.normType = " << syncProjection_bottom_normType << endl;

        syncProjection_bottom_verbosity = bottom_verbosity;
        ppBottom.query("verbosity", syncProjection_bottom_verbosity);
        pout() << "\tsyncProjection_bottom.verbosity = " << syncProjection_bottom_verbosity << endl;

        pout() << endl;
    }

    // Volume discrepancy (freestream preservation) solver overrides
    if (etaLambda > 0.0) {
        ParmParse ppAMRMG("VD_AMRMG");

        VD_AMRMG_eps = AMRMG_eps;
        ppAMRMG.query("eps", VD_AMRMG_eps);
        pout() << "\tVD_AMRMG.eps = " << VD_AMRMG_eps << endl;

        VD_AMRMG_num_smooth_down = AMRMG_num_smooth_down;
        ppAMRMG.query("num_smooth_down", VD_AMRMG_num_smooth_down);
        pout() << "\tVD_AMRMG.num_smooth_down = " << VD_AMRMG_num_smooth_down << endl;

        VD_AMRMG_num_smooth_up = AMRMG_num_smooth_up;
        ppAMRMG.query("num_smooth_up", VD_AMRMG_num_smooth_up);
        pout() << "\tVD_AMRMG.num_smooth_up = " << VD_AMRMG_num_smooth_up << endl;

        VD_AMRMG_num_smooth_bottom = AMRMG_num_smooth_bottom;
        ppAMRMG.query("num_smooth_bottom", VD_AMRMG_num_smooth_bottom);
        pout() << "\tVD_AMRMG.num_smooth_bottom = " << VD_AMRMG_num_smooth_bottom << endl;

        VD_AMRMG_num_smooth_precond = AMRMG_num_smooth_precond;
        ppAMRMG.query("num_smooth_precond", VD_AMRMG_num_smooth_precond);
        pout() << "\tVD_AMRMG.num_smooth_precond = " << VD_AMRMG_num_smooth_precond << endl;

        VD_AMRMG_numMG = AMRMG_numMG;
        ppAMRMG.query("numMG", VD_AMRMG_numMG);
        pout() << "\tVD_AMRMG.numMG = " << VD_AMRMG_numMG << endl;

        VD_AMRMG_imin = AMRMG_imin;
        ppAMRMG.query("imin", VD_AMRMG_imin);
        pout() << "\tVD_AMRMG.imin = " << VD_AMRMG_imin << endl;

        VD_AMRMG_imax = AMRMG_imax;
        ppAMRMG.query("imax", VD_AMRMG_imax);
        pout() << "\tVD_AMRMG.imax = " << VD_AMRMG_imax << endl;

        VD_AMRMG_hang = AMRMG_hang;
        ppAMRMG.query("hang", VD_AMRMG_hang);
        pout() << "\tVD_AMRMG.hang = " << VD_AMRMG_hang << endl;

        VD_AMRMG_normThresh = AMRMG_normThresh;
        ppAMRMG.query("normThresh", VD_AMRMG_normThresh);
        pout() << "\tVD_AMRMG.normThresh = " << VD_AMRMG_normThresh << endl;

        VD_AMRMG_maxDepth = AMRMG_maxDepth;
        ppAMRMG.query("maxDepth", VD_AMRMG_maxDepth);
        pout() << "\tVD_AMRMG.maxDepth = " << VD_AMRMG_maxDepth << endl;

        VD_AMRMG_relaxMode = AMRMG_relaxMode;
        ppAMRMG.query("relax_mode", VD_AMRMG_relaxMode);
        pout() << "\tVD_AMRMG.relax_mode = " << VD_AMRMG_relaxMode << endl;

        VD_AMRMG_precondMode = AMRMG_precondMode;
        ppAMRMG.query("precond_mode", VD_AMRMG_precondMode);
        pout() << "\tVD_AMRMG.precondMode = " << VD_AMRMG_precondMode << endl;

        VD_AMRMG_verbosity = AMRMG_verbosity;
        ppAMRMG.query("verbosity", VD_AMRMG_verbosity);
        pout() << "\tVD_AMRMG.verbosity = " << VD_AMRMG_verbosity << endl;


        ParmParse ppBottom("VD_bottom");

        VD_bottom_eps = bottom_eps;
        ppBottom.query("eps", VD_bottom_eps);
        pout() << "\tVD_bottom.eps = " << VD_bottom_eps << endl;

        VD_bottom_reps = bottom_reps;
        ppBottom.query("reps", VD_bottom_reps);
        pout() << "\tVD_bottom.reps = " << VD_bottom_reps << endl;

        VD_bottom_imax = bottom_imax;
        ppBottom.query("imax", VD_bottom_imax);
        pout() << "\tVD_bottom.imax = " << VD_bottom_imax << endl;

        VD_bottom_numRestarts = bottom_numRestarts;
        ppBottom.query("numRestarts", VD_bottom_numRestarts);
        pout() << "\tVD_bottom.numRestarts = " << VD_bottom_numRestarts << endl;

        VD_bottom_hang = bottom_hang;
        ppBottom.query("hang", VD_bottom_hang);
        pout() << "\tVD_bottom.hang = " << VD_bottom_hang << endl;

        VD_bottom_small = bottom_small;
        ppBottom.query("small", VD_bottom_small);
        pout() << "\tVD_bottom.small = " << VD_bottom_small << endl;

        VD_bottom_normType = bottom_normType;
        ppBottom.query("normType", VD_bottom_normType);
        pout() << "\tVD_bottom.normType = " << VD_bottom_normType << endl;

        VD_bottom_verbosity = bottom_verbosity;
        ppBottom.query("verbosity", VD_bottom_verbosity);
        pout() << "\tVD_bottom.verbosity = " << VD_bottom_verbosity << endl;

        pout() << endl;
    }
}
