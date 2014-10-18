/*******************************************************************************
 *    FILE: AMRLESMeta.H
 *    DESC: Maintains static information about the different codes.
 *    DATE: Fri 24 Jan 2014 04:34:52 PM EST
 *    MAIL: EdSantilli@physics.unc.edu
 ******************************************************************************/

#ifndef __AMRLESMeta_H__INCLUDED__
#define __AMRLESMeta_H__INCLUDED__

#include "mpi.h"
#include <vector>
#include <utility>
#include "REAL.H"


// This is not defined on all mpi versions
#ifndef MPI_UNDEFINED_RANK
#   define MPI_UNDEFINED_RANK MPI_UNDEFINED
#endif


namespace AMRLESMeta
{
    // Tags..
    static const int TAG_AMR2LES_MAGIC = 93001;


    // Communicators and such...

    // World intracommunication
    extern int worldRank;
    extern int worldSize;
    extern MPI_Comm worldComm;
    extern MPI_Group worldGroup;

    // AMR intracommunication
    extern int amrRank;
    extern int amrSize;
    extern MPI_Comm amrComm;
    extern MPI_Group amrGroup;

    // LES intracommunication
    extern int lesRank;
    extern int lesSize;
    extern MPI_Comm lesComm;
    extern MPI_Group lesGroup;

    // AMR/LES intercommunication
    extern int amr2lesLeader;
    extern int les2amrLeader;
    extern MPI_Comm interComm;
    extern MPI_Comm amrlesPeerComm;


    // Helper functions...

    // Gets the group from a communicator.
    MPI_Group getCommGroup (const MPI_Comm& a_comm);

    // Just like numProc, but acts on any group instead of the amrComm.
    int getGroupSize (const MPI_Group& a_group);

    // Just like procID, but returns MPI_UNDEFINED_RANK if the caller
    // is not in the group.
    int getGroupRank (const MPI_Group& a_group);

    // Is the caller a group member?
    inline bool isGroupMember (const MPI_Group& a_group);

    // Finds the procs that are in the union of group1 and group2.
    // Suppose we find that 3 procs are in the union. The return values will
    // look like this:
    //   vector[0] = pair<proc 0's group1 rank, proc 0's group2 rank>
    //   vector[1] = pair<proc 1's group1 rank, proc 1's group2 rank>
    //   vector[2] = pair<proc 2's group1 rank, proc 2's group2 rank>
    // The vector elements 0, 1, and 2 are not chosen in any particular order.
    std::vector<std::pair<int,int> >
    getRanksInGroupIntersection (const MPI_Group& a_group1, const MPI_Group& a_group2);
};


// -----------------------------------------------------------------------------
// Is the caller a group member?
// -----------------------------------------------------------------------------
inline bool AMRLESMeta::isGroupMember (const MPI_Group& a_group)
{
    return (getGroupRank(a_group) != MPI_UNDEFINED);
}


#endif //!__AMRLESMeta_H__INCLUDED__
