#include "AMRLESMeta.H"
#include "MayDay.H"
#include <sstream>
#include "SPMD.H"


// World intracommunication
int       AMRLESMeta::worldRank  = MPI_UNDEFINED_RANK;
int       AMRLESMeta::worldSize  = MPI_UNDEFINED;
MPI_Comm  AMRLESMeta::worldComm  = MPI_COMM_NULL;
MPI_Group AMRLESMeta::worldGroup = MPI_GROUP_NULL;

// AMR intracommunication
int       AMRLESMeta::amrRank  = MPI_UNDEFINED_RANK;
int       AMRLESMeta::amrSize  = MPI_UNDEFINED;
MPI_Comm  AMRLESMeta::amrComm  = MPI_COMM_NULL;
MPI_Group AMRLESMeta::amrGroup = MPI_GROUP_NULL;

// LES intracommunication
int       AMRLESMeta::lesRank  = MPI_UNDEFINED_RANK;
int       AMRLESMeta::lesSize  = MPI_UNDEFINED;
MPI_Comm  AMRLESMeta::lesComm  = MPI_COMM_NULL;
MPI_Group AMRLESMeta::lesGroup = MPI_GROUP_NULL;

// AMR/LES intercommunication
int AMRLESMeta::amr2lesLeader = MPI_UNDEFINED_RANK;
int AMRLESMeta::les2amrLeader = MPI_UNDEFINED_RANK;
MPI_Comm AMRLESMeta::interComm = MPI_COMM_NULL;
MPI_Comm AMRLESMeta::amrlesPeerComm = MPI_COMM_NULL;



// -----------------------------------------------------------------------------
// Gets the group from a communicator.
// -----------------------------------------------------------------------------
MPI_Group AMRLESMeta::getCommGroup (const MPI_Comm& a_comm)
{
    if (a_comm == MPI_COMM_NULL) {
        return MPI_GROUP_NULL;
    }

    MPI_Group group;
    int ierr = MPI_Comm_group(a_comm, &group);
    if (ierr != MPI_SUCCESS) {
        std::ostringstream errmsg;
        errmsg << "MPI_Comm_group failed. Error " << ierr << std::endl;
        MayDay::Error(errmsg.str().c_str());
    }

    return group;
}


// -----------------------------------------------------------------------------
// Just like numProc, but acts on any group instead of the amrComm.
// -----------------------------------------------------------------------------
int AMRLESMeta::getGroupSize (const MPI_Group& a_group)
{
    if (a_group == MPI_GROUP_NULL) {
        return MPI_UNDEFINED;
    }

    int np;
    int ierr = MPI_Group_size(a_group, &np);
    if (ierr != MPI_SUCCESS) {
        std::ostringstream errmsg;
        errmsg << "MPI_Group_size failed. Error " << ierr << std::endl;
        MayDay::Error(errmsg.str().c_str());
    }

    return np;
}


// -----------------------------------------------------------------------------
// Just like procID, but returns MPI_UNDEFINED_RANK if the caller
// is not in the group.
// -----------------------------------------------------------------------------
int AMRLESMeta::getGroupRank (const MPI_Group& a_group)
{
    if (a_group == MPI_GROUP_NULL) {
        return MPI_UNDEFINED_RANK;
    }

    int rank;
    int ierr = MPI_Group_rank(a_group, &rank);
    if (ierr != MPI_SUCCESS) {
        std::ostringstream errmsg;
        errmsg << "MPI_Group_rank failed. Error " << ierr << std::endl;
        MayDay::Error(errmsg.str().c_str());
    }

    return rank;
}


// -----------------------------------------------------------------------------
// Finds the procs that are in the union of group1 and group2.
// Suppose we find that 3 procs are in the union. The return values will
// look like this:
//   vector[0] = pair<proc 0's group1 rank, proc 0's group2 rank>
//   vector[1] = pair<proc 1's group1 rank, proc 1's group2 rank>
//   vector[2] = pair<proc 2's group1 rank, proc 2's group2 rank>
// The vector elements 0, 1, and 2 are not chosen in any particular order.
// -----------------------------------------------------------------------------
std::vector<std::pair<int,int> >
AMRLESMeta::getRanksInGroupIntersection (const MPI_Group& a_group1, const MPI_Group& a_group2)
{
    if (a_group1 == MPI_GROUP_NULL) {
        MayDay::Error("AMRLESMeta::getRanksInGroupIntersection was passed a NULL group1");
    }
    if (a_group2 == MPI_GROUP_NULL) {
        MayDay::Error("AMRLESMeta::getRanksInGroupIntersection was passed a NULL group2");
    }

    using std::vector;
    using std::pair;
    using std::make_pair;

    vector<pair<int,int> > retvec(0);

    // Is this proc a member of both groups?
    int groupRank[2];
    groupRank[0] = getGroupRank(a_group1);
    groupRank[1] = getGroupRank(a_group2);
    if (groupRank[0] == MPI_UNDEFINED || groupRank[1] == MPI_UNDEFINED) {
        groupRank[0] = MPI_UNDEFINED;
        groupRank[1] = MPI_UNDEFINED;
    }

    // All-to-all communication of results.
    int rbufCount = getGroupSize(worldGroup);
    int* rbuf = new int[2*rbufCount];
    MPI_Allgather(&groupRank, 2, MPI_INT, rbuf, 2, MPI_INT, worldComm);

    // Iterate through the results and collect results.
    for (int i = 0; i < 2*rbufCount; i += 2) {
        if (rbuf[i] == MPI_UNDEFINED) {
            if (rbuf[i+1] != MPI_UNDEFINED) {
                MayDay::Error("getRanksInGroupIntersection received a non-matching pair");
            }
        } else {
            if (rbuf[i+1] == MPI_UNDEFINED) {
                MayDay::Error("getRanksInGroupIntersection received a non-matching pair");
            } else {
                retvec.push_back(make_pair(rbuf[i],rbuf[i+1]));
            }
        }
    }

    // Free memory and exit.
    delete[] rbuf;
    rbuf = NULL;

    return retvec;
}
