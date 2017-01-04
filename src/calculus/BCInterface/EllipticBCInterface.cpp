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
#include "EllipticBCInterface.H"
#include "FluxBox.H"

const Real BOGUS_TIME = -1.0e300;


// -----------------------------------------------------------------------------
// Default constructor.
// This object will not be usable until BCMethods are added.
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

