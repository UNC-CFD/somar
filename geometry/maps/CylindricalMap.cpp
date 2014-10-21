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
#include "CylindricalMap.H"
#include "CylindricalMapF_F.H"


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
CylindricalMap::CylindricalMap ()
{;}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
CylindricalMap::~CylindricalMap ()
{;}


// -----------------------------------------------------------------------------
// 1. Must return the name of the coordinate mapping
// -----------------------------------------------------------------------------
const char* CylindricalMap::getCoorMapName () const
{
    return "Cylindrical";
}


// -----------------------------------------------------------------------------
// 2. Must return whether or not this metric is diagonal
// -----------------------------------------------------------------------------
bool CylindricalMap::isDiagonal () const
{
    return true;
}


// -----------------------------------------------------------------------------
// 3. Must return whether or not this metric is uniform
// -----------------------------------------------------------------------------
bool CylindricalMap::isUniform () const
{
    return false;
}


// -----------------------------------------------------------------------------
// 4. Must fill a mapped box with Cartesian locations.
// -----------------------------------------------------------------------------
void CylindricalMap::fill_physCoor (FArrayBox&      a_dest,
                                    const int       a_destComp,
                                    const int       a_mu,
                                    const RealVect& a_dXi) const
{
    CH_TIME("CylindricalMap::fill_physCoor");

    // Sanity checks
    CH_assert(0 <= a_destComp);
    CH_assert(a_destComp < a_dest.nComp());
    CH_assert(0 <= a_mu);
    CH_assert(a_mu < SpaceDim);

    const Box& destBox = a_dest.box();
    const IntVect& destBoxType = destBox.type();

    FORT_CYLINDRICAL_FILL_PHYSCOOR (
        CHF_FRA1(a_dest,a_destComp),
        CHF_CONST_INT(a_mu),
        CHF_CONST_REALVECT(a_dXi),
        CHF_BOX(destBox),
        CHF_CONST_INTVECT(destBoxType));
}


// -----------------------------------------------------------------------------
// Fills a mapped box with Cartesian locations (a_dest must have SpaceDim comps)
// Mild time saver.
// -----------------------------------------------------------------------------
void CylindricalMap::fill_physCoor (FArrayBox&      a_dest,
                                    const RealVect& a_dXi,
                                    const RealVect  a_scale) const
{
    CH_TIME("CylindricalMap::fill_physCoor (all comps)");

    const Box& destBox = a_dest.box();
    const IntVect& destBoxType = destBox.type();

    CH_assert(a_dest.nComp() == SpaceDim);

    FORT_CYLINDRICAL_FILL_PHYSCOOR_ALL_COMPS (
        CHF_FRA(a_dest),
        CHF_CONST_REALVECT(a_dXi),
        CHF_CONST_REALVECT(a_scale),
        CHF_BOX(destBox),
        CHF_CONST_INTVECT(destBoxType));
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the Jacobian matrix elements d[x^mu] / d[xi^nu].
// Major speedup!!!
// -----------------------------------------------------------------------------
void CylindricalMap::fill_dxdXi (FArrayBox&      a_dest,
                                 const int       a_destComp,
                                 const int       a_mu,
                                 const int       a_nu,
                                 const RealVect& a_dXi,
                                 const Real      a_scale) const
{
    CH_TIME("CylindricalMap::fill_dxdXi");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));
    CH_assert((0 <= a_nu) && (a_nu < SpaceDim));

    const Box& destBox = a_dest.box();
    const IntVect& destBoxType = destBox.type();

    CH_assert(destBoxType.product() == 0 || destBoxType.sum() == 1);

    if (a_mu == 2 || a_nu == 2) {
        if (a_mu != a_nu) {
            a_dest.setVal(0.0, a_destComp);
        } else {
            a_dest.setVal(a_scale, a_destComp);
        }
    } else {
        FORT_CYLINDRICAL_FILL_DXDXI (
            CHF_FRA1(a_dest,a_destComp),
            CHF_CONST_INT(a_mu),
            CHF_CONST_INT(a_nu),
            CHF_CONST_REALVECT(a_dXi),
            CHF_CONST_REAL(a_scale),
            CHF_BOX(destBox),
            CHF_CONST_INTVECT(destBoxType));
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with J = det[Jacobian]
// Major speedup!
// Without the analytic version of this function, machine-epsilon sized
// noise is introduced.
// -----------------------------------------------------------------------------
void CylindricalMap::fill_J (FArrayBox&      a_dest,
                             const int       a_destComp,
                             const RealVect& a_dXi,
                             const Real      a_scale) const
{
    CH_TIME("CylindricalMap::fill_J");

    // Sanity check
    CH_assert(a_dest.interval().contains(a_destComp));

    const Box& destBox = a_dest.box();
    const int destBoxType0 = destBox.type()[0];
    const Real dXi0 = a_dXi[0];

    FORT_CYLINDRICAL_FILL_J (
        CHF_FRA1(a_dest,a_destComp),
        CHF_CONST_REAL(dXi0),
        CHF_CONST_REAL(a_scale),
        CHF_BOX(destBox),
        CHF_CONST_INT(destBoxType0));
}

