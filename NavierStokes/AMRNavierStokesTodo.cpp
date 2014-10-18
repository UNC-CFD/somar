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
#include "AMRNavierStokes.H"
#include "Debug.H"

// This file contains unused/stubbed-out function definitions.


// -----------------------------------------------------------------------------
// Unused function
// -----------------------------------------------------------------------------
AMRNavierStokes::AMRNavierStokes (MappedAMRLevel* a_coarser_level_ptr,
                                  const Box&      a_prob_domain,
                                  int             a_level,
                                  const IntVect&  a_ref_ratio)
{
    TODO();
}


// -----------------------------------------------------------------------------
// Unused function
// -----------------------------------------------------------------------------
AMRNavierStokes::AMRNavierStokes (MappedAMRLevel*      a_coarser_level_ptr,
                                  const ProblemDomain& a_prob_domain,
                                  int                  a_level,
                                  const IntVect&       a_ref_ratio)
{
    TODO();
}


// -----------------------------------------------------------------------------
// Unused function
// -----------------------------------------------------------------------------
MappedAMRLevel* AMRNavierStokes::makeAMRLevel (MappedAMRLevel* a_coarser_level_ptr,
                                               const Box&      a_problem_domain,
                                               int             a_level,
                                               const IntVect&  a_ref_ratio) const
{
    TODO();
    return NULL;
}


// -----------------------------------------------------------------------------
// Unused function
// -----------------------------------------------------------------------------
MappedAMRLevel* AMRNavierStokes::makeAMRLevel (MappedAMRLevel*      a_coarser_level_ptr,
                                               const ProblemDomain& a_problem_domain,
                                               int                  a_level,
                                               const IntVect&       a_ref_ratio) const
{
    TODO();
    return NULL;
}


// -----------------------------------------------------------------------------
// Unused function
// -----------------------------------------------------------------------------
void AMRNavierStokes::define (MappedAMRLevel* a_coarse_level_ptr,
                              const Box&      a_problem_domain,
                              int             a_level,
                              const IntVect&  a_ref_ratio)
{
    TODO();
}
