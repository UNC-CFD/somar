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
                              const Real                a_dtTheta)
{
    this->define(a_geoSourcePtr, a_physBCPtr, a_dtTheta);
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
                            const Real                a_dtTheta)
{
    m_geoSourcePtr = a_geoSourcePtr;
    m_physBCPtr = a_physBCPtr;
    m_theta = a_dtTheta;
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

        // Fill dXidz
        FArrayBox dXidzFAB(CCBox, SpaceDim);
        D_TERM(m_geoSourcePtr->fill_dXidx(dXidzFAB, 0, 0, SpaceDim-1, a_dXi);,
               m_geoSourcePtr->fill_dXidx(dXidzFAB, 1, 1, SpaceDim-1, a_dXi);,
               m_geoSourcePtr->fill_dXidx(dXidzFAB, 2, 2, SpaceDim-1, a_dXi);)

        // Compute Nsq
        FArrayBox NsqFAB(CCBox, 1);
        BoxIterator bit(CCBox);
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& cc = bit();
            NsqFAB(cc) = D_TERM(-dXidzFAB(cc,0) * (bbarFB[0](cc+BASISV(0)) - bbarFB[0](cc)) / a_dXi[0],
                                -dXidzFAB(cc,1) * (bbarFB[1](cc+BASISV(1)) - bbarFB[1](cc)) / a_dXi[1],
                                -dXidzFAB(cc,2) * (bbarFB[2](cc+BASISV(2)) - bbarFB[2](cc)) / a_dXi[2]);
        }

        // Send to destAlias
        if (destBoxType == IntVect::Zero) {
            destAlias.copy(NsqFAB);
        } else {
            D_TERM(int FCdir = 0;,
                   if (destBoxType[1] == 1) FCdir = 1;,
                   if (destBoxType[2] == 1) FCdir = 2;)

            CellToEdge(NsqFAB, destAlias, FCdir);
        }
    }

    // * (dt*theta)^2
    destAlias *= m_theta*m_theta;

    // = -omega^2/(1+omega^2)
    tmpFAB.copy(destAlias);
    tmpFAB += 1.0;
    destAlias /= tmpFAB;
    destAlias *= -1.0;

    // * dXi^mu/dz
    m_geoSourcePtr->fill_dXidx(tmpFAB, 0, a_mu, SpaceDim-1, a_dXi);
    destAlias *= tmpFAB;

    // * dXi^nu/dz
    if (a_mu != a_nu) {
        m_geoSourcePtr->fill_dXidx(tmpFAB, 0, a_nu, SpaceDim-1, a_dXi);
    }
    destAlias *= tmpFAB;

    // + g^{mu,nu}
    m_geoSourcePtr->fill_gup(tmpFAB, 0, a_mu, a_nu, a_dXi);
    destAlias += tmpFAB;

    // * J*a_scale
    m_geoSourcePtr->fill_J(tmpFAB, 0, a_dXi, a_scale);
    destAlias *= tmpFAB;
}
