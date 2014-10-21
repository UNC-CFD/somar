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
#include "BeamGeneratorMap.H"
#include "BeamGeneratorMapF_F.H"
#include "Subspace.H"
#include "ProblemContext.H"

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
BeamGeneratorMap::BeamGeneratorMap ()
: BathymetricBaseMap()
{
    const ProblemContext* ctx = ProblemContext::getInstance();
    m_alpha = ctx->beamGenMapAlpha;
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
BeamGeneratorMap::~BeamGeneratorMap ()
{;}


// -----------------------------------------------------------------------------
// Must return the name of the coordinate mapping
// -----------------------------------------------------------------------------
const char* BeamGeneratorMap::getCoorMapName () const
{
    return "BeamGeneratorMap";
}


// -----------------------------------------------------------------------------
// Must return whether or not this metric is diagonal
// -----------------------------------------------------------------------------
bool BeamGeneratorMap::isDiagonal () const
{
    return false;
}


// -----------------------------------------------------------------------------
// Fills a NodeFAB with the bathymetric data. a_dest must be flat in the
// vertical. Upon return, each point in the horizontal (Xi,Eta) of a_dest
// will contain the (positive) local elevation from the sea floor.
// NOTE: This vertical distance is measured in a straight line perpendicular
// the the surface. We are measuring this distance along the Cartesian
// vertical coordinate line, not the mapped vertical coordinate line.
// -----------------------------------------------------------------------------
void BeamGeneratorMap::fill_bathymetry (FArrayBox&       a_dest,
                                        const int        a_destComp,
                                        const FArrayBox& a_cartPos,
                                        const RealVect&  a_dXi) const
{
    const Box& destBox = a_dest.box();
    const IntVect destBoxType = destBox.type();

    // The holder needs to be flat and nodal in the vertical.
    CH_assert(destBox == horizontalDataBox(destBox));
    CH_assert(destBoxType[SpaceDim-1] == 1);

    FORT_FILL_BEAMGENERATORMAPBATHYMETRY (
        CHF_FRA1(a_dest, a_destComp),
        CHF_BOX(destBox),
        CHF_CONST_INTVECT(destBoxType),
        CHF_CONST_FRA(a_cartPos),
        CHF_CONST_REALVECT(a_dXi),
        CHF_CONST_REALVECT(m_L),
        CHF_CONST_REAL(m_alpha));
}
