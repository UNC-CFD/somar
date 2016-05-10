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
#include "AlteredMetric.H"
#include "LevelGeometry.H"
#include "PhysBCUtil.H"
#include "BoxIterator.H"
#include "CellToEdge.H"
#include "StratUtils.H"


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
AlteredMetric::AlteredMetric ()
: m_geoSourcePtr(NULL),
  m_physBCPtr(NULL)
{;}


// -----------------------------------------------------------------------------
// Full constructor
// -----------------------------------------------------------------------------
AlteredMetric::AlteredMetric (const GeoSourceInterface* a_geoSourcePtr,
                              const PhysBCUtil*         a_physBCPtr,
                              const Real                a_dtTheta,
                              const Real                a_coriolisF)
{
    this->define(a_geoSourcePtr, a_physBCPtr, a_dtTheta, a_coriolisF);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
AlteredMetric::~AlteredMetric ()
{
    m_geoSourcePtr = NULL;
    m_physBCPtr = NULL;
}


// -----------------------------------------------------------------------------
// Full constructor
// -----------------------------------------------------------------------------
void AlteredMetric::define (const GeoSourceInterface* a_geoSourcePtr,
                            const PhysBCUtil*         a_physBCPtr,
                            const Real                a_dtTheta,
                            const Real                a_coriolisF)
{
    m_geoSourcePtr = a_geoSourcePtr;
    m_physBCPtr = a_physBCPtr;
    m_dtTheta = a_dtTheta;
    m_coriolisF = a_coriolisF;
}


// -----------------------------------------------------------------------------
// Set dest = J*(gup^{ij} - (omega^2 / (1 + omega^2)) * dXi^i/dz * dXi^j/dz)
// where omega = (dt*N*theta)^2
// -----------------------------------------------------------------------------
void AlteredMetric::fill_Jgup (FArrayBox&       a_dest,
                               const int        a_destComp,
                               const int        a_mu,
                               const int        a_nu,
                               const RealVect&  a_dXi,
                               const Real       a_scale) const
{
    // Create work space
    FArrayBox destAlias(Interval(a_destComp, a_destComp), a_dest);
    const Box& destBox = destAlias.box();
    const IntVect& destBoxType = destBox.type();
    FArrayBox tmpFAB(destBox, 1);

    CH_assert(destBoxType == BASISV(a_mu)); // For now.

    // Fill destAlias with Nsq.
    {
        // To fill bbar, we need a simple LevGeo to pass to physBCUtil.
        // Only getDx and getGeoSource should be called from in there,
        // so this define doesn't need to do much and should be fast.
        LevelGeometry simplelevGeo(a_dXi);

        // We will need some CC space large enough to enclose destBox.
        Box CCBox = destBox;
        CCBox.grow(1);
        CCBox.enclosedCells();

        // Fill an FC holder with bbar
        FluxBox bbarFB(CCBox, 1);
        const Real dummyTime = -1.0e300;
        const DataIndex dummyDi;
        D_TERM(m_physBCPtr->setBackgroundScalar(bbarFB[0], 0, simplelevGeo, dummyDi, dummyTime);,
               m_physBCPtr->setBackgroundScalar(bbarFB[1], 0, simplelevGeo, dummyDi, dummyTime);,
               m_physBCPtr->setBackgroundScalar(bbarFB[2], 0, simplelevGeo, dummyDi, dummyTime);)

        // Compute Nsq
        FArrayBox NsqFAB(CCBox, 1);
        computeBVFreq(NsqFAB, bbarFB, CCBox, simplelevGeo, 2.0);
        CellToEdge(NsqFAB, destAlias, a_mu);
    }

    // * (dt*theta)^2
    destAlias *= m_dtTheta*m_dtTheta;

    // = -omega^2/(1+omega^2)
    tmpFAB.copy(destAlias);
    tmpFAB += 1.0;
    destAlias /= tmpFAB;
    destAlias *= -1.0;

    // = +ftilde^2/(1+ftilde^2)
    const Real ftilde = m_coriolisF * m_dtTheta;
    const Real ftildesq = ftilde*ftilde;
    const Real invfCoeff = 1.0 / (1.0 + ftildesq);
    destAlias += ftildesq * invfCoeff;

    // * dXi^mu/dz
    m_geoSourcePtr->fill_dXidx(tmpFAB, 0, a_mu, SpaceDim-1, a_dXi);
    destAlias *= tmpFAB;

    // * dXi^nu/dz
    if (a_mu != a_nu) {
        m_geoSourcePtr->fill_dXidx(tmpFAB, 0, a_nu, SpaceDim-1, a_dXi);
    }
    destAlias *= tmpFAB;

    // + ftilde*invfCoeff*(Horizontal Jacobian)
    if (a_mu != a_nu) {
        FArrayBox ix(destBox, 1);
        FArrayBox jy(destBox, 1);
        FArrayBox iy(destBox, 1);
        FArrayBox jx(destBox, 1);

        m_geoSourcePtr->fill_dXidx(ix, 0, a_mu, 0, a_dXi);
        m_geoSourcePtr->fill_dXidx(jy, 0, a_nu, 1, a_dXi);
        m_geoSourcePtr->fill_dXidx(iy, 0, a_mu, 1, a_dXi);
        m_geoSourcePtr->fill_dXidx(jx, 0, a_nu, 0, a_dXi);

        BoxIterator bit(destBox);
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& cc = bit();

            destAlias(cc) += ftilde * invfCoeff
                           * (ix(cc)*jy(cc) - iy(cc)*jx(cc));
        }
    }

    // + g^{mu,nu} * invfCoeff
    m_geoSourcePtr->fill_gup(tmpFAB, 0, a_mu, a_nu, a_dXi);
    tmpFAB *= invfCoeff;
    destAlias += tmpFAB;

    // * J*a_scale
    m_geoSourcePtr->fill_J(tmpFAB, 0, a_dXi, a_scale);
    destAlias *= tmpFAB;
}
