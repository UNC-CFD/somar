/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Copyright (C) 2014 Edward Santilli & Alberto Scotti
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
#include "EllipticBCInterface.H"
#include "FluxBox.H"

static const Real BOGUS_TIME = -1.0e300;


// -----------------------------------------------------------------------------
// Default constructor.
// This object will not be useable until BCMethods are added.
// -----------------------------------------------------------------------------
BCMethodHolder::BCMethodHolder ()
{

    m_BCGhostFuncPtrs.resize(0);
    m_BCFluxFuncPtrs.resize(0);
    m_BCGhostObjPtrs.resize(0);
    m_BCFluxObjPtrs.resize(0);
}


// -----------------------------------------------------------------------------
// Add a BCGhostFunc method.
// -----------------------------------------------------------------------------
void BCMethodHolder::addBCMethod (BCGhostFunc a_methodPtr)
{

    for (int idx = 0; idx < m_BCGhostFuncPtrs.size(); ++idx) {
        if (m_BCGhostFuncPtrs[idx] == a_methodPtr)
            MayDay::Error("BCMethodHolder received a redundant BCGhostFunc");
    }
    m_BCGhostFuncPtrs.push_back(a_methodPtr);
}


// -----------------------------------------------------------------------------
// Add a BCFluxFunc method.
// -----------------------------------------------------------------------------
void BCMethodHolder::addBCMethod (BCFluxFunc a_methodPtr)
{

    for (int idx = 0; idx < m_BCFluxFuncPtrs.size(); ++idx) {
        if (m_BCFluxFuncPtrs[idx] == a_methodPtr)
            MayDay::Error("BCMethodHolder received a redundant BCFluxFunc");
    }
    m_BCFluxFuncPtrs.push_back(a_methodPtr);
}


// -----------------------------------------------------------------------------
// Add a BCGhostClass method.
// -----------------------------------------------------------------------------
void BCMethodHolder::addBCMethod (RefCountedPtr<BCGhostClass>& a_methodPtr)
{

    for (int idx = 0; idx < m_BCGhostObjPtrs.size(); ++idx) {
        if (m_BCGhostObjPtrs[idx] == a_methodPtr)
            MayDay::Error("BCMethodHolder received a redundant BCGhostClass");
    }
    m_BCGhostObjPtrs.push_back(a_methodPtr);

    m_ghostDescriptor |= a_methodPtr->getBCDescriptor();
}


// -----------------------------------------------------------------------------
// Add a BCFluxClass method.
// -----------------------------------------------------------------------------
void BCMethodHolder::addBCMethod (RefCountedPtr<BCFluxClass>& a_methodPtr)
{

    for (int idx = 0; idx < m_BCFluxObjPtrs.size(); ++idx) {
        if (m_BCFluxObjPtrs[idx] == a_methodPtr)
            MayDay::Error("BCMethodHolder received a redundant BCFluxClass");
    }
    m_BCFluxObjPtrs.push_back(a_methodPtr);

    m_fluxDescriptor |= a_methodPtr->getBCDescriptor();
}


// -----------------------------------------------------------------------------
// Sets ghost cells
// -----------------------------------------------------------------------------
const BCDescriptor& BCMethodHolder::setGhosts (FArrayBox&           a_state,
                                               const FArrayBox*     a_extrapPtr,
                                               const Box&           a_valid,
                                               const ProblemDomain& a_domain,
                                               const RealVect&      a_dx,
                                               const DataIndex&     a_index,
                                               const FluxBox*       a_JgupPtr,
                                               bool                 a_homogeneous,
                                               Real                 a_time,
                                               const Interval&      a_interval) const
{
    CH_assert(this->isDefined());

    int numFuncs = m_BCGhostFuncPtrs.size();
    int numObjs = m_BCGhostObjPtrs.size();

    for (int idx = 0; idx < numFuncs; ++idx)
        m_BCGhostFuncPtrs[idx] (a_state, a_extrapPtr, a_valid, a_domain, a_dx,
                                a_index, a_JgupPtr, a_homogeneous, a_time, a_interval);

    for (int idx = 0; idx < numObjs; ++idx)
        m_BCGhostObjPtrs[idx]->operator() (a_state, a_extrapPtr, a_valid, a_domain, a_dx,
                                           a_index, a_JgupPtr, a_homogeneous, a_time, a_interval);

    return m_ghostDescriptor;
}


// -----------------------------------------------------------------------------
// Sets fluxes in one dir
// -----------------------------------------------------------------------------
const BCDirDescriptor& BCMethodHolder::setFluxes (FArrayBox&           a_state,
                                                  const FArrayBox*     a_extrapPtr,
                                                  const Box&           a_valid,
                                                  const ProblemDomain& a_domain,
                                                  const RealVect&      a_dx,
                                                  const DataIndex&     a_index,
                                                  const FluxBox*       a_JgupPtr,
                                                  int                  a_dir,
                                                  bool                 a_homogeneous,
                                                  Real                 a_time,
                                                  const Interval&      a_interval) const
{
    CH_assert(this->isDefined());

    int numFuncs = m_BCFluxFuncPtrs.size();
    int numObjs = m_BCFluxObjPtrs.size();

    for (int idx = 0; idx < numFuncs; ++idx)
        m_BCFluxFuncPtrs[idx] (a_state, a_extrapPtr, a_valid, a_domain, a_dx,
                               a_index, a_JgupPtr, a_dir,
                               a_homogeneous, a_time, a_interval);

    for (int idx = 0; idx < numObjs; ++idx)
        m_BCFluxObjPtrs[idx]->operator() (a_state, a_extrapPtr, a_valid, a_domain, a_dx,
                                          a_index, a_JgupPtr, a_dir,
                                          a_homogeneous, a_time, a_interval);

    return m_fluxDescriptor[a_dir];
}


// -----------------------------------------------------------------------------
// Sets fluxes in all dirs
// -----------------------------------------------------------------------------
const BCDescriptor& BCMethodHolder::setFluxes (FluxBox&             a_state,
                                               const FArrayBox*     a_extrapPtr,
                                               const Box&           a_valid,
                                               const ProblemDomain& a_domain,
                                               const RealVect&      a_dx,
                                               const DataIndex&     a_index,
                                               const FluxBox*       a_JgupPtr,
                                               bool                 a_homogeneous,
                                               Real                 a_time,
                                               const Interval&      a_interval) const
{
    CH_assert(this->isDefined());

    int numFuncs = m_BCFluxFuncPtrs.size();
    int numObjs = m_BCFluxObjPtrs.size();

    for (int idx = 0; idx < numFuncs; ++idx) {
        D_TERM(
            m_BCFluxFuncPtrs[idx] (a_state[0], a_extrapPtr, a_valid, a_domain,
                                   a_dx, a_index, a_JgupPtr, 0,
                                   a_homogeneous, a_time, a_interval);,
            m_BCFluxFuncPtrs[idx] (a_state[1], a_extrapPtr, a_valid, a_domain,
                                   a_dx, a_index, a_JgupPtr, 1,
                                   a_homogeneous, a_time, a_interval);,
            m_BCFluxFuncPtrs[idx] (a_state[2], a_extrapPtr, a_valid, a_domain,
                                   a_dx, a_index, a_JgupPtr, 2,
                                   a_homogeneous, a_time, a_interval);
        )
    }

    for (int idx = 0; idx < numObjs; ++idx) {
        D_TERM(
            m_BCFluxObjPtrs[idx]->operator() (a_state[0], a_extrapPtr, a_valid, a_domain,
                                              a_dx, a_index, a_JgupPtr, 0,
                                              a_homogeneous, a_time, a_interval);,
            m_BCFluxObjPtrs[idx]->operator() (a_state[1], a_extrapPtr, a_valid, a_domain,
                                              a_dx, a_index, a_JgupPtr, 1,
                                              a_homogeneous, a_time, a_interval);,
            m_BCFluxObjPtrs[idx]->operator() (a_state[2], a_extrapPtr, a_valid, a_domain,
                                              a_dx, a_index, a_JgupPtr, 2,
                                              a_homogeneous, a_time, a_interval);
        )
    }

    return m_fluxDescriptor;
}


// -----------------------------------------------------------------------------
// Search boundaries for non-periodic, non-Neum BCs.
// -----------------------------------------------------------------------------
bool BCMethodHolder::hasNullSpace (const DisjointBoxLayout& a_grids,
                                   const CFRegion&          a_cfregion,
                                   const IntVect&           a_activeDirs) const
{
    CH_assert(this->isDefined());

    // Let each processor search its grids.
    bool localVal = this->hasNullSpaceNoComm(a_grids, a_cfregion, a_activeDirs);

    // Gather/broadcast the overall results. ('and'ing all local results...)
    // MPI_C_BOOL seems to be problematic, so I'm just using ints.
#ifdef CH_MPI
    int localInt = (localVal? 1: 0);
    int globalInt;

    int result = MPI_Allreduce(&localInt, &globalInt, 1, MPI_INT, MPI_BAND, Chombo_MPI::comm);
    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in hasNullSpace");
    }

    bool globalVal = (globalInt != 0);
#else
    bool globalVal = localVal;
#endif

    // We are done.
    return globalVal;
}


// -----------------------------------------------------------------------------
// Returns true if all BCs on a level are periodic or Neumann.
// This does not perform communication.
// -----------------------------------------------------------------------------
bool BCMethodHolder::hasNullSpaceNoComm (const DisjointBoxLayout& a_grids,
                                         const CFRegion&          a_cfregion,
                                         const IntVect&           a_activeDirs) const
{
    CH_assert(this->isDefined());

    // This function hereby swears not to alter this structure.
    CFRegion& cfregionRef = (CFRegion&)a_cfregion;

    const BCDescriptor& ghostDesc = this->getGhostDescriptor();
    const BCDescriptor& fluxDesc = this->getFluxDescriptor();

    const ProblemDomain& domain = a_grids.physDomain();
    const Box& domBox = domain.domainBox();
    DataIterator dit = a_grids.dataIterator();

    for (int dir = 0; dir < SpaceDim; ++dir) {
        // If this is not an active dir, then ignore it altogether.
        if (a_activeDirs[dir] == 0) continue;

        // Periodic dirs need avgs removed.
        if (domain.isPeriodic(dir)) continue;

        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide& iside = sit();

            // If we get here, we know we don't have periodic BCs in this dir.
            // If we also don't have Neum BCs, then we can stop searching.
            // Warning: This check is a bit tricky. You had better be very
            // careful if you decide to change things around with DeMorgan's Law
            // and whatnot. Remember, this is just a time saver. You can always
            // skip this and perform the full checks below.
            if (!(ghostDesc[dir][iside] == BCType::Neum || fluxDesc[dir][iside] == BCType::Neum)) {
                return false;
            }

            for (dit.reset(); dit.ok(); ++dit) {
                const Box& valid = a_grids[dit];

                // We are at a dir and side that is critical. If the CFRegion is not
                // empty, then we are dealing with Diri BCs and we do not want to
                // remove the avg.
                if (   (iside == Side::Lo && !cfregionRef.loCFIVS(dit(), dir).isEmpty()) ||
                       (iside == Side::Hi && !cfregionRef.hiCFIVS(dit(), dir).isEmpty())   ) {
                    return false;
                }

                // If this grid abuts the physical boundary and the BC is not
                // Neum or periodic, then we do not want to remove the avg.
                int thisDesc = ghostDesc.stencil(valid, domain, dir, iside);
                if (thisDesc != BCType::None && thisDesc != BCType::Neum) {
                    return false;
                }

                thisDesc = fluxDesc.stencil(valid, domain, dir, iside);
                if (thisDesc != BCType::None && thisDesc != BCType::Neum) {
                    return false;
                }

            } // end loop over grids (dit)
        } // end loop over sides (sit)
    } // end loop over dirs (dir)

    // We did not find a non-periodic, non-Neum BC.
    // Remove the average.
    return true;
}

