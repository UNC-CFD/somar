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
#include "AMRLESMeta.H"
#include "Printing.H"
#include "LevelGeometry.H"
#include "AMRNavierStokesFactory.H"
#include "LepticAMR.H"
#include "FORT_PROTO.H"
#include "ProblemContext.H"


// Function prototypes
void setupMPIComms (const Real a_lesProcFrac);
void testMPIComms ();
void nsrun ();

// The LES code is not for public use.
// #define USE_LES
#ifdef USE_LES
extern "C" {
    void FORTRAN_NAME(DIABLO,diablo) (
        int*, int*, int*, int*
    );
}
#endif //USE_LES

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

    ParmParse ppLES("les");

    Real lesProcFrac = 0.0;
#   ifdef USE_LES
        ppLES.query("proc_frac", lesProcFrac);
        pout() << "\tles.proc_frac = " << lesProcFrac << endl;
#   endif

#   ifndef CH_MPI
        if (lesProcFrac != 0.0) {
            MayDay::Warning("Cannot run LES without MPI. Continuing les.proc_frac = 0.0");
            lesProcFrac = 0.0;
        }
#   endif

    // This creates groups and comms for communication among amr and les.
    // Call this even if we are not using MPI.
    setupMPIComms(lesProcFrac);
    testMPIComms();

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

    // BEGIN: Chombo-only code.
    if (AMRLESMeta::amrSize > 0) {
        if (AMRLESMeta::isGroupMember(AMRLESMeta::amrGroup)) {
#           ifdef CH_MPI
                pout() << "Using MPI." << endl;

                int chomboRank = procID();
                int chomboNumRanks = numProc();
                pout() << "\nchomboRank = " << chomboRank
                       << "\nchomboNumRanks = " << chomboNumRanks
                       << endl;
#           else
                pout() << "Not using MPI." << endl;
#           endif

            // Make sure the user knows if we are in debug mode.
#           ifndef NDEBUG
                pout() << "Proc #" << AMRLESMeta::worldRank
                       << "/" << AMRLESMeta::worldSize-1
                       << ": *** DEBUG MODE ***"
                       << endl;
#           else
                pout() << "Proc #" << AMRLESMeta::worldRank
                       << "/" << AMRLESMeta::worldSize-1
                       << ": *** RELEASE MODE ***"
                       << endl;
#           endif

            pout() << "SpaceDim = " << SpaceDim << endl;

#           ifdef  CH_USE_MEMORY_TRACKING
                pout() << "Using memory tracking." << endl;
#           endif

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
#               ifdef _GNU_SOURCE
                    pout() << "  domainname = " << hostInfo.domainname << endl;
#               endif
            }

#           ifndef NDEBUG
                // Only register the debugger on certain systems.
                bool useDebugger = false;
                useDebugger |= std::string(hostInfo.nodename) == std::string("iceman");
                useDebugger |= (std::string(hostInfo.nodename) == std::string("scorpion") && numProc() == 1);
                if (useDebugger) {
                    pout() << "Registering debugger..." << endl;
                    registerDebugger();
                }
#           endif
            pout() << endl;

            // Setup AMR and run the simulation
            nsrun();

            // MPI_Barrier(AMRLESMeta::amrComm);
            printf("worldRank %d: AMR code finished.\n", AMRLESMeta::worldRank);
        }
    }
    // END: Chombo-only code.

#ifdef USE_LES
    // BEGIN: LES-only code.
    if (AMRLESMeta::lesSize > 0) {
        if (AMRLESMeta::isGroupMember(AMRLESMeta::lesGroup)) {

            // Register the debugger?
#           ifndef NDEBUG
                if (false) {
                    registerDebugger();
                }
#           endif

             FORTRAN_NAME(DIABLO,diablo)(
                 &AMRLESMeta::lesComm,
                 &AMRLESMeta::interComm,
                 &AMRLESMeta::amr2lesLeader,
                 &AMRLESMeta::les2amrLeader
             );

            // MPI_Barrier(AMRLESMeta::lesComm);
            printf("worldRank %d: LES code finished.\n", AMRLESMeta::worldRank);
        }
    }
    // End LES-only code
#endif //USE_LES

    cout << flush;
    MPI_Barrier(MPI_COMM_WORLD);

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
// This creates groups and comms for communication among amr and les.
// -----------------------------------------------------------------------------
void setupMPIComms (const Real a_lesProcFrac)
{
    using namespace AMRLESMeta;

    // Gather world info
    worldComm = MPI_COMM_WORLD;
    worldGroup = getCommGroup(MPI_COMM_WORLD);
    worldSize = getGroupSize(worldGroup);
    worldRank = getGroupRank(worldGroup);


    // Compute the sizes of the groups.
    amrSize = int(round( Real(worldSize)*(1.0-a_lesProcFrac) ));
    if (amrSize < 0 || worldSize < amrSize) {
        MayDay::Error("setupMPIComms produced an amrSize that is out of range");
    }

    lesSize = worldSize - amrSize;
    if (lesSize < 0 || worldSize < lesSize) {
        MayDay::Error("setupMPIComms produced an lesSize that is out of range");
    }

    // Now, we have to split worldComm.
    int membershipKey = ((worldRank < amrSize)? 0: 1);
    MPI_Comm splitComm;
    if (MPI_Comm_split(worldComm, membershipKey, worldRank, &splitComm) != MPI_SUCCESS) {
        MayDay::Error("setupMPIComms could not split MPI_COMM_WORLD");
    }
    if (worldRank < amrSize) {
        amrComm = splitComm;
        amrGroup = getCommGroup(amrComm);
        amrRank = getGroupRank(amrGroup);

        lesComm = MPI_COMM_NULL;
        lesGroup = MPI_GROUP_NULL;
        lesRank = MPI_UNDEFINED_RANK;

    } else {
        lesComm = splitComm;
        lesGroup = getCommGroup(lesComm);
        lesRank = getGroupRank(lesGroup);

        amrComm = MPI_COMM_NULL;
        amrGroup = MPI_GROUP_NULL;
        amrRank = MPI_UNDEFINED_RANK;
    }

    // Set Chombo's communicator
    if (isGroupMember(amrGroup)) {
        Chombo_MPI::comm = amrComm;
    }

    // Create the intercommunicator for the AMR and LES groups.
    if (amrSize > 0 && lesSize > 0) {
        amr2lesLeader = 0; // This is the rank in amrComm.
        les2amrLeader = 0; // This is the rank in lesComm.
        int tag = 435;

        // Use a dedicated peer communicator
        if (MPI_Comm_dup(MPI_COMM_WORLD, &amrlesPeerComm) != MPI_SUCCESS) {
            MayDay::Error("setupMPIComms had an error creating amrlesPeerComm");
        }

        if (isGroupMember(amrGroup)) {
            int remoteLeader = amrSize; // The local group specifies the worldRank.
            if (MPI_Intercomm_create(splitComm, amr2lesLeader, amrlesPeerComm, remoteLeader, tag, &interComm) != MPI_SUCCESS) {
                MayDay::Error("setupMPIComms had an error creating the AMR to LES intercommunicator");
            }

            // Collective communications never want a local rank, they want a special macro.
            amr2lesLeader = ((amrRank == amr2lesLeader)? MPI_ROOT: MPI_PROC_NULL);

        } else {
            int remoteLeader = 0; // The remote group specifies the remoteRank  (remote = les).
            if (MPI_Intercomm_create(splitComm, amr2lesLeader, amrlesPeerComm, remoteLeader, tag, &interComm) != MPI_SUCCESS) {
                MayDay::Error("setupMPIComms had an error creating the AMR to LES intercommunicator");
            }

            // Collective communications never want a local rank, they want a special macro.
            les2amrLeader = ((lesRank == les2amrLeader)? MPI_ROOT: MPI_PROC_NULL);
        }
    } else {
        amr2lesLeader = MPI_UNDEFINED_RANK;
        amrlesPeerComm = MPI_COMM_NULL;
        interComm = MPI_COMM_NULL;
    }

    // Name the communicators to help debugging along.
    if (isGroupMember(amrGroup)) {
        MPI_Comm_set_name(amrComm, "amrComm");
    }
    if (isGroupMember(lesGroup)) {
        MPI_Comm_set_name(lesComm, "lesComm");
    }
    if (amrSize > 0 && lesSize > 0) {
        MPI_Comm_set_name(interComm, "interComm");
    }
}


// -----------------------------------------------------------------------------
// Passes magic values among the AMR and LES groups to test the
// intercommunicator.
// -----------------------------------------------------------------------------
void testMPIComms () {
    using namespace AMRLESMeta;

    // Do we have an intercommunicator?
    if (amrSize <= 0 || lesSize <= 0) return;

    // Make sure Chombo is using amrComm.
    if (isGroupMember(amrGroup)) {
        if (Chombo_MPI::comm != amrComm) {
            MayDay::Error("testMPIComms found that Chombo is not using amrComm");
        }
        if (procID() != amrRank) {
            MayDay::Error("testMPIComms found that Chombo is not using amrComm to compute ranks");
        }
        if (numProc() != amrSize) {
            MayDay::Error("testMPIComms found that Chombo is not using amrComm to compute comm size");
        }
    }

    // AMR to LES intercommunication
    {
        const double targetMagic = 12345.09876;

        double localMagic = -1.0;
        if (amr2lesLeader == MPI_ROOT) {
            localMagic = targetMagic;
        }

        if (MPI_Bcast(&localMagic, 1, MPI_DOUBLE, amr2lesLeader, interComm) != MPI_SUCCESS) {
            MayDay::Error("testMPIComms could not broadcast a magic number from AMR to LES");
        }

        if (isGroupMember(lesGroup)) {
            if (localMagic != targetMagic) {
                MayDay::Error("testMPIComms received an invalid magic number while testing AMR to LES intercommunication");
            }
        } else {
            if (amr2lesLeader == MPI_ROOT) {
                if (localMagic != targetMagic) {
                    MayDay::Error("testMPIComms clobbered the magic number while testing AMR to LES intercommunication");
                }
            } else {
                if (localMagic != -1.0) {
                    MayDay::Error("testMPIComms broadcasted the magic number to too many ranks during AMR to LES intercommunication");
                }
            }
        }
    }

    // AMR to LES intercommunication
    {
        const double targetMagic = 92753.41608;

        double localMagic = -1.0;
        if (les2amrLeader == MPI_ROOT) {
            localMagic = targetMagic;
        }

        if (MPI_Bcast(&localMagic, 1, MPI_DOUBLE, les2amrLeader, interComm) != MPI_SUCCESS) {
            MayDay::Error("testMPIComms could not broadcast a magic number from LES to AMR");
        }

        if (isGroupMember(amrGroup)) {
            if (localMagic != targetMagic) {
                MayDay::Error("testMPIComms received an invalid magic number while testing LES to AMR intercommunication");
            }
        } else {
            if (les2amrLeader == MPI_ROOT) {
                if (localMagic != targetMagic) {
                    MayDay::Error("testMPIComms clobbered the magic number while testing LES to AMR intercommunication");
                }
            } else {
                if (localMagic != -1.0) {
                    MayDay::Error("testMPIComms broadcasted the magic number to too many ranks during LES to AMR intercommunication");
                }
            }
        }
    }
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

