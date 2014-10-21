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
#include "TwistedMap.H"
#include "TwistedMapF_F.H"
#include "ProblemContext.H"



// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
TwistedMap::TwistedMap ()
{
    const ProblemContext* ctx = ProblemContext::getInstance();
    m_L = ctx->domainLength;
    m_pert = ctx->pert;

    // 0 = Original, 1+eps*sin version (JCOMP 2952-2976)
    // 1 = New, 1+cos*cos*sin version.
    m_twistType = 1;
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
TwistedMap::~TwistedMap ()
{;}


// -----------------------------------------------------------------------------
// 1. Must return the name of the coordinate mapping
// -----------------------------------------------------------------------------
const char* TwistedMap::getCoorMapName () const
{
    static char name[20];
    sprintf(name, "Twisted (type %d)", m_twistType);
    return name;
}


// -----------------------------------------------------------------------------
// 2. Must return whether or not this metric is diagonal
// -----------------------------------------------------------------------------
bool TwistedMap::isDiagonal () const
{
    return false;
}


// -----------------------------------------------------------------------------
// 3. Must return whether or not this metric is uniform
// -----------------------------------------------------------------------------
bool TwistedMap::isUniform () const
{
    return false;
}


// -----------------------------------------------------------------------------
// 4. Must fill a mapped box with Cartesian locations.
// Typically, only a_dXi[a_mu] will be used, I've included all dXi comps
// just in case we ever need them for some unforseen reason.
// -----------------------------------------------------------------------------
void TwistedMap::fill_physCoor (FArrayBox&      a_dest,
                                const int       a_destComp,
                                const int       a_mu,
                                const RealVect& a_dXi) const
{
    CH_TIME("TwistedMap::fill_physCoor");

    // Sanity checks
    CH_assert(0 <= a_destComp);
    CH_assert(a_destComp < a_dest.nComp());
    CH_assert(0 <= a_mu);
    CH_assert(a_mu < SpaceDim);

    const Box& destBox = a_dest.box();
    const IntVect& destBoxType = destBox.type();
    const Real& pert = m_pert[a_mu];

    if (m_twistType == 0) {
        FORT_TWISTED0_FILL_PHYSCOOR(
            CHF_FRA1(a_dest,a_destComp),
            CHF_CONST_INT(a_mu),
            CHF_CONST_REALVECT(a_dXi),
            CHF_CONST_REAL(pert),
            CHF_BOX(destBox),
            CHF_CONST_INTVECT(destBoxType));
    } else {
        FORT_TWISTED1_FILL_PHYSCOOR(
            CHF_FRA1(a_dest,a_destComp),
            CHF_CONST_INT(a_mu),
            CHF_CONST_REALVECT(a_dXi),
            CHF_CONST_REALVECT(m_L),
            CHF_CONST_REAL(pert),
            CHF_BOX(destBox),
            CHF_CONST_INTVECT(destBoxType));
    }
}


// -----------------------------------------------------------------------------
// Fills a mapped box with Cartesian locations (a_dest must have SpaceDim comps)
// Mild time saver.
// -----------------------------------------------------------------------------
void TwistedMap::fill_physCoor (FArrayBox&      a_dest,
                                const RealVect& a_dXi,
                                const RealVect  a_scale) const
{
    CH_TIME("TwistedMap::fill_physCoor (all comps)");

    CH_assert(a_dest.nComp() == SpaceDim);

    const Box& destBox = a_dest.box();
    const IntVect& destBoxType = destBox.type();

    if (m_twistType == 0) {
        FORT_TWISTED0_FILL_PHYSCOOR_ALL_COMPS(
            CHF_FRA(a_dest),
            CHF_CONST_REALVECT(a_dXi),
            CHF_CONST_REALVECT(m_pert),
            CHF_CONST_REALVECT(a_scale),
            CHF_BOX(destBox),
            CHF_CONST_INTVECT(destBoxType));
    } else {
        GeoSourceInterface::fill_physCoor(a_dest, a_dXi, a_scale);
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the Jacobian matrix elements d[x^mu] / d[xi^nu].
// Major speedup!!!
// -----------------------------------------------------------------------------
void TwistedMap::fill_dxdXi (FArrayBox&      a_dest,
                             const int       a_destComp,
                             const int       a_mu,
                             const int       a_nu,
                             const RealVect& a_dXi,
                             const Real      a_scale) const
{
    CH_TIME("TwistedMap::fill_dxdXi");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));
    CH_assert((0 <= a_nu) && (a_nu < SpaceDim));

    const Box& destBox = a_dest.box();
    const IntVect& destBoxType = destBox.type();
    const Real& pert = m_pert[a_mu];

    CH_assert(destBoxType.product() == 0 || destBoxType.sum() == 1);

    if (m_twistType == 0) {
        if (a_mu == a_nu) {
            a_dest.setVal(a_scale, a_destComp);
        } else {
            FORT_TWISTED0_FILL_DXDXI(
                CHF_FRA1(a_dest,a_destComp),
                CHF_CONST_INT(a_mu),
                CHF_CONST_INT(a_nu),
                CHF_CONST_REALVECT(a_dXi),
                CHF_CONST_REAL(pert),
                CHF_CONST_REAL(a_scale),
                CHF_BOX(destBox),
                CHF_CONST_INTVECT(destBoxType));
        }
    } else {
        GeoSourceInterface::fill_dxdXi(a_dest,
                                       a_destComp,
                                       a_mu,
                                       a_nu,
                                       a_dXi,
                                       a_scale);
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with J = det[Jacobian]
// Major speedup!
// Without the analytic version of this function, machine-epsilon sized
// noise is introduced.
// -----------------------------------------------------------------------------
void TwistedMap::fill_J (FArrayBox&      a_dest,
                         const int       a_destComp,
                         const RealVect& a_dXi,
                         const Real      a_scale) const
{
    CH_TIME("TwistedMap::fill_J");

    // Sanity check
    CH_assert(a_dest.interval().contains(a_destComp));

    const Box& destBox = a_dest.box();
    const IntVect& destBoxType = destBox.type();

    if (m_twistType == 0) {
        FORT_TWISTED0_FILL_J(
            CHF_FRA1(a_dest,a_destComp),
            CHF_CONST_REALVECT(a_dXi),
            CHF_CONST_REALVECT(m_pert),
            CHF_CONST_REAL(a_scale),
            CHF_BOX(destBox),
            CHF_CONST_INTVECT(destBoxType));
    } else {
        GeoSourceInterface::fill_J(a_dest,
                                   a_destComp,
                                   a_dXi,
                                   a_scale);
    }
}

