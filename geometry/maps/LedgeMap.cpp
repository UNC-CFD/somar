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
#include "LedgeMap.H"
#include "ProblemContext.H"
#include "Subspace.H"
#include "BoxIterator.H"
#include "Debug.H"


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
LedgeMap::LedgeMap ()
: BathymetricBaseMap()
{
    const ProblemContext* ctx = ProblemContext::getInstance();
    m_transitionOrder = ctx->ledgeMapTransitionOrder;
    m_hl = ctx->ledgeMapHl;
    m_hr = ctx->ledgeMapHr;
    m_xl = ctx->ledgeMapXl;
    m_xr = ctx->ledgeMapXr;

    Real dh = m_hr - m_hl;
    Real dx = m_xr - m_xl;
    Real invdx3 = pow(dx, -3.0);

    if (m_transitionOrder == 1) {
        // Compute the linear coefficients.
        m_coeff0 = m_hr - m_xr*dh/dx;
        m_coeff1 = dh/dx;
        m_coeff2 = 0.0;
        m_coeff3 = 0.0;
    } else if (m_transitionOrder == 3) {
        // Compute the cubic coefficients.
        m_coeff0 = m_hr + dh*(3.0*m_xl-m_xr)*m_xr*m_xr*invdx3;
        m_coeff1 = -6.0*dh*m_xl*m_xr*invdx3;
        m_coeff2 = 3.0*dh*(m_xl+m_xr)*invdx3;
        m_coeff3 = -2.0*dh*invdx3;
    } else {
        MayDay::Error("LedgeMap::m_transitionOrder must be 1 or 3");
    }
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
LedgeMap::~LedgeMap ()
{;}


// -----------------------------------------------------------------------------
// Must return the name of the coordinate mapping
// -----------------------------------------------------------------------------
const char* LedgeMap::getCoorMapName () const
{
    return "LedgeMap";
}


// -----------------------------------------------------------------------------
// Must return whether or not this metric is diagonal
// -----------------------------------------------------------------------------
bool LedgeMap::isDiagonal () const
{
    return false;
}


// -----------------------------------------------------------------------------
// Fills a NodeFAB with the bathymetric data. a_dest must be flat in the
// vertical. Upon return, each point in the horizontal (Xi,Eta) of a_dest
// will contain the (positive) local depth.
// NOTE: This vertical distance is measured in a straight line perpendicular
// to the surface. We are measuring this distance along the Cartesian
// vertical coordinate line, not the mapped vertical coordinate line.
// -----------------------------------------------------------------------------
void LedgeMap::fill_bathymetry (FArrayBox&       a_dest,
                                const int        a_destComp,
                                const FArrayBox& a_cartPos,
                                const RealVect&  a_dXi) const
{
    const Box& destBox = a_dest.box();
    const IntVect destBoxType = destBox.type();

    // The holder needs to be flat and nodal in the vertical.
    CH_assert(destBox == horizontalDataBox(destBox));
    CH_assert(destBoxType[SpaceDim-1] == 1);

#if CH_SPACEDIM == 3
    // Alberto's Gaussian bump...
    // TODO: Make this it's own class.
    Real x, y, R2;
    BoxIterator bit(destBox);
    for (bit.reset(); bit.ok(); ++bit){
        const IntVect& iv = bit();
        x = a_cartPos(iv,0) - m_xl;
        y = a_cartPos(iv,1) - m_xl;
        R2 = x*x + y*y;
        a_dest(iv,a_destComp) = (m_hl-m_hr)*exp(-R2/(m_xr*m_xr)) + m_hr;
    }

#else
    Real x;
    BoxIterator bit(destBox);

    switch (m_transitionOrder) {
    case 1: // Linear transition region
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& iv = bit();
            x = a_cartPos(iv, 0);

            if (x < m_xl) {
                a_dest(iv,a_destComp) = m_hl;
            } else if (x > m_xr) {
                a_dest(iv,a_destComp) = m_hr;
            } else {
                a_dest(iv,a_destComp) = m_coeff0 + x*m_coeff1;
            }
        }
        break;

    case 3: // Cubic transition region
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& iv = bit();
            x = a_cartPos(iv, 0);

            if (x < m_xl) {
                a_dest(iv,a_destComp) = m_hl;
            } else if (x > m_xr) {
                a_dest(iv,a_destComp) = m_hr;
            } else {
                a_dest(iv,a_destComp) = m_coeff0 + x*(m_coeff1 + x*(m_coeff2 + x*m_coeff3));
            }
        }
        break;

    default:
        MayDay::Error("LedgeMap::fill_bathymetry received an invalid m_transitionOrder");
    }
#endif
}
