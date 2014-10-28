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


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
AMRNavierStokesFactory::AMRNavierStokesFactory ()
{
    const ProblemContext* ctx = ProblemContext::getInstance();

    m_cfl = ctx->cfl;

    // Create a new PhysBCUtil for this level.
    m_physBCPtr = ctx->newPhysBCUtil();
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
