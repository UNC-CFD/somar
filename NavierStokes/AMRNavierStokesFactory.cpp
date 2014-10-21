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
#include "AMRNavierStokesFactory.H"
#include "AMRNavierStokes.H"
#include "ProblemContext.H"

#include "AdvectionTestBCUtil.H"
#include "LockExchangeBCUtil.H"
#include "BeamGenerationBCUtil.H"
#include "InternalWaveBCUtil.H"
#include "TaylorGreenBCUtil.H"
#include "VortexStreetBCUtil.H"
#include "HorizConvBCUtil.H"
#include "SolitaryWaveBCUtil.H"
#include "DJLBCUtil.H"


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
AMRNavierStokesFactory::AMRNavierStokesFactory ()
{
    const ProblemContext* ctx = ProblemContext::getInstance();

    m_cfl = ctx->cfl;

    // Determine the correct PhysBCUtil to pass to the AMRLevels.
    switch (ctx->problem) {
    case ProblemContext::ProblemType::ADVECTION_TEST:
        m_physBCPtr = new AdvectionTestBCUtil;
        break;
    case ProblemContext::ProblemType::LOCK_EXCHANGE:
        m_physBCPtr = new LockExchangeBCUtil;
        break;
    case ProblemContext::ProblemType::BEAM_GENERATION:
        m_physBCPtr = new BeamGenerationBCUtil;
        break;
    case ProblemContext::ProblemType::INTERNAL_WAVE:
        m_physBCPtr = new InternalWaveBCUtil;
        break;
    case ProblemContext::ProblemType::TAYLOR_GREEN:
        m_physBCPtr = new TaylorGreenBCUtil;
        break;
    case ProblemContext::ProblemType::VORTEX_STREET:
        m_physBCPtr = new VortexStreetBCUtil;
        break;
    case ProblemContext::ProblemType::HORIZ_CONV:
        m_physBCPtr = new HorizConvBCUtil;
        break;
    case ProblemContext::ProblemType::SOLITARYWAVE:
        m_physBCPtr = new SolitaryWaveBCUtil;
        break;
    case ProblemContext::ProblemType::DJL:
        m_physBCPtr = new DJLBCUtil;
        break;
    default:
        // Undefined problem
        MayDay::Error("Bad problem type");
    }

    m_physBCPtr->define();
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
AMRNavierStokesFactory::~AMRNavierStokesFactory ()
{
    if (m_physBCPtr != NULL) {
        delete m_physBCPtr;
        m_physBCPtr = NULL;
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void AMRNavierStokesFactory::setPhysBC (const PhysBCUtil& a_physBC)
{
    MayDay::Error("AMRNavierStokesFactory::setPhysBC: Do we really need this?");

    // Remember to delete the old pointer...
    if (m_physBCPtr != NULL) {
        delete m_physBCPtr;
        m_physBCPtr = NULL;
    }

    // ...before replacing it.
    m_physBCPtr = a_physBC.newPhysBCUtil();
}


// -----------------------------------------------------------------------------
// Virtual override of level factory
// -----------------------------------------------------------------------------
MappedAMRLevel* AMRNavierStokesFactory::new_amrlevel () const
{
    CH_assert (m_physBCPtr != NULL);

    // Create the new AMRNavierStokes level
    AMRNavierStokes* amrns_ptr = new AMRNavierStokes();

    // Set its parameters
    amrns_ptr->m_physBCPtr = m_physBCPtr;
    amrns_ptr->m_cfl = m_cfl;

    // Return an AMRLevel slice of the new object.
    return (static_cast<MappedAMRLevel*>(amrns_ptr));
}
