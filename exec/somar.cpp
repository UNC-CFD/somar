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


// This is the driver of the Navier-Stokes solver. It contains the main() C++
// function. This file was not written to be user-friendly and should not need
// to be altered. -ES


// Standard headers
#include <iostream>
#include <sys/utsname.h>
using std::cout;
using std::endl;

#ifdef CH_MPI
#   include "mpi.h"
#endif

// Chombo headers
#include "SPMD.H"
#include "SPACE.H"
#include "parstream.H"
#include "ParmParse.H"
#include "memusage.H"

#ifndef NDEBUG
#   include "CH_Attach.H"
#endif

#ifdef CH_USE_MEMORY_TRACKING
#   include "Pool.H"
#endif

#ifndef CH_NTIMER
#   include "OldTimer.H"
#endif

// Project headers
#include "Printing.H"
#include "LevelGeometry.H"
#include "AMRNavierStokesFactory.H"
#include "LepticAMR.H"
#include "FORT_PROTO.H"
#include "ProblemContext.H"


// Function prototypes
void nsrun ();

// #define TRAP_FPE  //(should be off by default)
#ifdef TRAP_FPE
    // Previous versions of glibc require the following code:
#   include "parstream.H"
    extern "C"
    {
#       include <fpu_control.h>
    }

    /* IM: Invalid operation mask
     * DM: Denormalized operand mask
     * ZM: Zero-divide mask
     * OM: Overflow mask
     * UM: Underflow mask
     * PM: Precision (inexact result) mask */
    static void __attribute__ ((constructor)) trapfpe(void)
    {
        pout() << " Turning on floating-point traps! " << std::endl;
        //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_UM);
        //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_ZM);
        //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_OM | _FPU_MASK_UM);
        fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_UM);
        //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM);
        //fpu_control_t cw = _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_DM | _FPU_MASK_UM);
        //fpu_control_t cw = _FPU_DEFAULT;
        _FPU_SETCW(cw);
        /* On x86, this expands to: */
        /* unsigned int cw = 0x037f & ~(0x01 | 0x04 | 0x08); */
        /* __asm__ ("fldcw %0" : : "m" (*&cw));              */
    }
#endif


// -----------------------------------------------------------------------------
// Setup and shutdown routines. This is the main()
// function, but the real work is done in mainLoop().
// -----------------------------------------------------------------------------
int main(int argc, char* argv[])
{
#   ifdef CH_MPI
        MPI_Init(&argc, &argv);
#   endif

    // Reset the stdout color
    cout << color::none << std::flush;

    // Open the input file
    char* in_file = argv[1];
    ParmParse pp(argc-2, argv+2, NULL, in_file);
    pout() << "\nReading from file: " << in_file << std::endl;

#   ifdef CH_MPI
#       ifdef CH_AIX
            H5dont_atexit();
#       endif
#   endif

#   ifdef TRAP_FPE
        trapfpe();
#   endif

#   ifndef CH_NTIMER
        // Start yer timers.
        OldTimer Everything;
        Everything.start();
#   endif

#   ifdef CH_MPI
        pout() << "Using MPI." << endl;

        int chomboRank = procID();
        int chomboNumRanks = numProc();
        pout() << "\nchomboRank = " << chomboRank
               << "\nchomboNumRanks = " << chomboNumRanks
               << endl;
#   else
        pout() << "Not using MPI." << endl;
#   endif

    // Make sure the user knows if we are in debug mode.
#   ifndef NDEBUG
        pout() << "Proc #" << procID()
               << "/" << numProc()-1
               << ": *** DEBUG MODE ***"
               << endl;
#   else
        pout() << "Proc #" << procID()
               << "/" << numProc()-1
               << ": *** RELEASE MODE ***"
               << endl;
#   endif

    pout() << "SpaceDim = " << SpaceDim << endl;

#   ifdef  CH_USE_MEMORY_TRACKING
        pout() << "Using memory tracking." << endl;
#   endif

    // Get the host name
    utsname hostInfo;
    int errCode = uname(&hostInfo);
    if (errCode != 0) {
        pout() << "\nerrCode = " << errCode << " recieved from uname function." << endl;
    } else {
        pout() << "\nHost info:"
               << "\n  sysname    = " << hostInfo.sysname
               << "\n  nodename   = " << hostInfo.nodename
               << "\n  release    = " << hostInfo.release
               << "\n  version    = " << hostInfo.version
               << "\n  machine    = " << hostInfo.machine
               << endl;
#       ifdef _GNU_SOURCE
            pout() << "  domainname = " << hostInfo.domainname << endl;
#       endif
    }

#   ifndef NDEBUG
        // Only register the debugger on certain systems.
        bool useDebugger = false;
        useDebugger |= std::string(hostInfo.nodename) == std::string("iceman");
        useDebugger |= (std::string(hostInfo.nodename) == std::string("scorpion") && numProc() == 1);
        if (useDebugger) {
            pout() << "Registering debugger..." << endl;
            registerDebugger();
        }
#   endif
    pout() << endl;

    // Setup AMR and run the simulation
    nsrun();
    cout << flush;

#   ifndef CH_NTIMER
        // Stop timing.
        Everything.stop();

        Real end_memory = get_memory_usage_from_OS();
        std::string end_time = formatTime(Everything.wc_time());

        pout()  << "\nEverything done.\n"
                << "mem usage: " << end_memory << "MB\n"
                << "elapsed time: " << end_time << endl;
        tout(0) << "Elapsed time = " << end_time << "." << endl;
#   endif

#   ifdef CH_USE_MEMORY_TRACKING
        dumpmemoryatexit();
#   endif

#   ifndef CH_NTIMER
        CH_TIMER_REPORT();
#   endif

#   ifdef CH_MPI
        MPI_Finalize();
#   endif

    return 0;
}


// -----------------------------------------------------------------------------
void nsrun ()
{
#ifndef CH_NTIMER
    // Time the setup routines
    OldTimer setupTmr;
    setupTmr.start();
#endif

    // This is probably the first request for a ProblemContext object.
    // After this, the input file will be read and finished with.
    const ProblemContext* ctx = ProblemContext::getInstance();

    // Create AMRFactory object
    AMRNavierStokesFactory amrns_fact;

    // Create AMR object
    LepticAMR thisAMR;
    thisAMR.maxGridSize(ctx->maxGridSize);
    thisAMR.maxBaseGridSize(ctx->maxBaseGridSize);
    thisAMR.define(ctx->max_level, ctx->refRatios, ctx->domain, &amrns_fact);
    thisAMR.verbosity(ctx->verbosity);
    MappedAMRLevel::verbosity(ctx->verbosity);

    thisAMR.useSubcyclingInTime(ctx->useSubcycling);

    thisAMR.plotInterval(ctx->plot_interval);
    thisAMR.plotPeriod(ctx->plot_period);
    thisAMR.plotPrefix(ctx->plot_prefix);
    thisAMR.checkpointInterval(ctx->checkpoint_interval);
    thisAMR.checkpointPrefix(ctx->check_prefix);
    thisAMR.gridBufferSize(ctx->bufferSize);

    thisAMR.maxGridSize(ctx->maxGridSize);
    thisAMR.maxBaseGridSize(ctx->maxBaseGridSize);
    thisAMR.splitDirs(ctx->splitDirs);
    thisAMR.fillRatio(ctx->fill_ratio);
    thisAMR.blockFactor(ctx->block_factor);
    thisAMR.regridIntervals(ctx->regrid_intervals);

    if (ctx->fixed_dt > 0) {
        thisAMR.fixedDt(ctx->fixed_dt);
    }
    thisAMR.maxDtGrow(ctx->max_dt_grow);

    if (ctx->isRestart) {
        // Initialize from restart file
#ifdef CH_USE_HDF5
        HDF5Handle handle(ctx->restart_file, HDF5Handle::OPEN_RDONLY);
        thisAMR.setupForRestart(handle);
        handle.close();
#else
        MayDay::Error("AMRNavierStokes restart only defined with HDF5");
#endif
    } else {
        // New run
        if (ctx->hasPredefinedGrids) {
            thisAMR.setupForFixedHierarchyRun(ctx->predefinedGrids, 1);
        } else {
            thisAMR.setupForNewAMRRun();
        }
    }

#ifndef CH_NTIMER
    setupTmr.stop();
    pout() << "Total setup time = " << formatTime(setupTmr.wc_time()) << endl;
#endif

    // Run
    thisAMR.run(ctx->stopTime, ctx->maxsteps);

    // Conclude
    thisAMR.conclude();

    // Free statically allocated memory.
    LevelGeometry::staticUndefine();
    LepticMeshRefine::deleteBuffer();
    ProblemContext::freeMemory();
}

