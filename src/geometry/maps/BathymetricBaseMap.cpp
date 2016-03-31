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
#include "BathymetricBaseMap.H"
#include "BathymetricBaseMapF_F.H"
#include "ProblemContext.H"
#include "NodeInterpF_F.H"
#include "ConvertFABF_F.H"
#include "Subspace.H"
#include "Debug.H"


// -----------------------------------------------------------------------------
// Construtor
// -----------------------------------------------------------------------------
BathymetricBaseMap::BathymetricBaseMap ()
{
    const ProblemContext* ctx = ProblemContext::getInstance();
    m_L = ctx->domainLength;
    m_lev0DXi = m_L / RealVect(ctx->nx);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
BathymetricBaseMap::~BathymetricBaseMap ()
{;}


// -----------------------------------------------------------------------------
// Fills a mapped box with Cartesian locations.
// This new code tries to generate a grid on AMR level 0, then perform
// a piecewise linear interpolation of the nodal points. This is an attempt to
// reduce boundary face mismatches among different AMR levels.
// -----------------------------------------------------------------------------
void BathymetricBaseMap::fill_physCoor (FArrayBox&      a_dest,
                                        const int       a_destComp,
                                        const int       a_mu,
                                        const RealVect& a_dXi) const
{
    CH_TIME("BathymetricBaseMap::fill_physCoor");

    // Sanity checks
    CH_assert(0 <= a_destComp);
    CH_assert(a_destComp < a_dest.nComp());
    CH_assert(0 <= a_mu);
    CH_assert(a_mu < SpaceDim);

    // Define the computation regions.
    const Box& destBox = a_dest.box();
    const IntVect& destBoxType = destBox.type();

    // Calculate the nodal positions, x(Xi).
    if (!isDiagonal()) {
        if (a_mu < SpaceDim-1) {
            // Horizontal components are identity mapped or simply stretched.
            const int destBoxDirType = destBoxType[a_mu];
            const Real dXi = a_dXi[a_mu];
            const Real domLength = m_L[a_mu];
            FORT_FILL_BATHYHORIZMAP(
                CHF_FRA1(a_dest, a_destComp),
                CHF_BOX(destBox),
                CHF_CONST_INT(destBoxDirType),
                CHF_CONST_INT(a_mu),
                CHF_CONST_REAL(domLength),
                CHF_CONST_REAL(dXi));

        } else {
            // Vertical components are streched by the local elevation,
            // which is a function of the horizontal position.

            // Collect data and create computation regions.
            const IntVect& nodeBoxType = IntVect::Unit;
            const Box nodeBox(surroundingNodes(destBox));
            const Box bottomNodeBox(horizontalDataBox(nodeBox));

            // Get node-centered horizontal positions.
            FArrayBox seaFloorFAB(bottomNodeBox, SpaceDim);
            for (int dir = 0; dir < SpaceDim-1; ++dir) {
                this->fill_physCoor(seaFloorFAB, dir, dir, a_dXi);
            }

            // Get the node-centered local depths.
            this->fill_bathymetry(seaFloorFAB, SpaceDim-1, seaFloorFAB, a_dXi);

            // Calculate the z-coordinate at each node.
            FArrayBox nodeZFAB(nodeBox, 1);
            FORT_FILL_BATHYVERTMAP(
                CHF_FRA1(nodeZFAB, 0),
                CHF_BOX(nodeBox),
                CHF_CONST_INTVECT(nodeBoxType),
                CHF_CONST_FRA1(seaFloorFAB, CH_SPACEDIM-1),
                CHF_CONST_REALVECT(m_L),
                CHF_CONST_REALVECT(a_dXi));

            // Convert data to destBoxType.
            FORT_CONVERTFAB(
                CHF_FRA1(a_dest, a_destComp),
                CHF_BOX(destBox),
                CHF_CONST_INTVECT(destBoxType),
                CHF_CONST_FRA1(nodeZFAB, 0),
                CHF_CONST_INTVECT(nodeBoxType));

            // // If this trips, we should have just set a_dest directly.
            // CH_assert(destBoxType != nodeBoxType);
        }
    } else {
        // TODO: Need to solve a least-squares problem.
        MayDay::Error("BathymetricBaseMap::fill_physCoor: Cannot solve least-squares problem for orthogonal coordinates yet");
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the Jacobian matrix elements d[x^mu] / d[xi^nu].
// -----------------------------------------------------------------------------
void BathymetricBaseMap::fill_dxdXi (FArrayBox&      a_dest,
                                     const int       a_destComp,
                                     const int       a_mu,
                                     const int       a_nu,
                                     const RealVect& a_dXi,
                                     const Real      a_scale) const
{
    // If this coordinate map is the result of a least-squares problem, Then
    // the Jacobian matrix is complicated and needs to be computed numerically.
    if (isDiagonal()) {
        GeoSourceInterface::fill_dxdXi(a_dest, a_destComp, a_mu, a_nu, a_dXi, a_scale);
        return;
    }

    if (a_mu != SpaceDim-1) {
        if (a_nu != SpaceDim-1) {
            // (H,H) block
            if (a_mu == a_nu) {
                const Box& destBox = a_dest.box();
                const IntVect& destBoxType = destBox.type();
                FORT_FILL_BATHYDXDXI(
                    CHF_FRA1(a_dest,a_destComp),
                    CHF_BOX(destBox),
                    CHF_CONST_INTVECT(destBoxType),
                    CHF_CONST_INT(a_nu),
                    CHF_CONST_REALVECT(m_L),
                    CHF_CONST_REALVECT(a_dXi),
                    CHF_CONST_REAL(a_scale));

            } else {
                a_dest.setVal(0.0, a_destComp);
            }
        } else {
            // (H,V) block
            // x is not a function of Zeta. This element is zero.
            a_dest.setVal(0.0, a_destComp);
        }
    } else {
        if (a_nu != SpaceDim-1) {
            // (V,H) block

            // Collect data and create computation regions.
            const Box& destBox = a_dest.box();
            const IntVect& destBoxType = destBox.type();
            const IntVect& edgeBoxType = IntVect::Unit - BASISV(a_nu);

            const IntVect growIV = (destBoxType[a_nu] == 1)? BASISV(a_nu): IntVect::Zero;
            const Box nodeBox(surroundingNodes(destBox).grow(growIV));
            const Box bottomNodeBox(horizontalDataBox(nodeBox));

            // Get the horizontal Cartesian positions
            FArrayBox cartPosFAB(bottomNodeBox, SpaceDim-1);
            for (int dir = 0; dir < SpaceDim-1; ++dir) {
                this->fill_physCoor(cartPosFAB, dir, dir, a_dXi);
            }

            // Get the node-centered bathymetric data.
            FArrayBox seaFloorFAB(bottomNodeBox, 1);
            this->fill_bathymetry(seaFloorFAB, 0, cartPosFAB, a_dXi);

            // Compute the Jacobian element at edges in the a_nu-direction.
            if (destBoxType == edgeBoxType) {
                FORT_FILL_BATHYDZDXI(
                    CHF_FRA1(a_dest,a_destComp),
                    CHF_BOX(destBox),
                    CHF_CONST_INT(a_nu),
                    CHF_CONST_FRA1(seaFloorFAB,0),
                    CHF_CONST_REALVECT(m_L),
                    CHF_CONST_REALVECT(a_dXi),
                    CHF_CONST_REAL(a_scale));
            } else {
                const Box edgeBox(enclosedCells(nodeBox, a_nu));
                FArrayBox edgeElemFAB(edgeBox, 1);
                FORT_FILL_BATHYDZDXI(
                    CHF_FRA1(edgeElemFAB,0),
                    CHF_BOX(edgeBox),
                    CHF_CONST_INT(a_nu),
                    CHF_CONST_FRA1(seaFloorFAB,0),
                    CHF_CONST_REALVECT(m_L),
                    CHF_CONST_REALVECT(a_dXi),
                    CHF_CONST_REAL(a_scale));

                // Convert to the desired centering
                FORT_CONVERTFAB(
                    CHF_FRA1(a_dest, a_destComp),
                    CHF_BOX(destBox),
                    CHF_CONST_INTVECT(destBoxType),
                    CHF_CONST_FRA1(edgeElemFAB, 0),
                    CHF_CONST_INTVECT(edgeBoxType));
            }

        } else {
            // (V,V) block

            // Collect data and create computation regions.
            const Box& destBox = a_dest.box();
            const IntVect destBoxType = destBox.type();
            const IntVect nodeBoxType = IntVect::Unit;

            const Box nodeBox(surroundingNodes(destBox));
            const Box bottomNodeBox(horizontalDataBox(nodeBox));

            // Get the horizontal Cartesian positions
            FArrayBox cartPosFAB(bottomNodeBox, SpaceDim-1);
            for (int dir = 0; dir < SpaceDim-1; ++dir) {
                this->fill_physCoor(cartPosFAB, dir, dir, a_dXi);
            }

            // Get the node-centered bathymetric data.
            FArrayBox seaFloorFAB(bottomNodeBox, 1);
            this->fill_bathymetry(seaFloorFAB, 0, cartPosFAB, a_dXi);

            // Calculate the Jacobian element at the nodes
            if (destBoxType == nodeBoxType) {
                FORT_FILL_BATHYDZDZETA(
                    CHF_FRA1(a_dest, a_destComp),
                    CHF_BOX(destBox),
                    CHF_CONST_FRA1(seaFloorFAB, 0),
                    CHF_CONST_REALVECT(m_L),
                    CHF_CONST_REALVECT(a_dXi),
                    CHF_CONST_REAL(a_scale));
            } else {
                FArrayBox nodeVals(nodeBox, 1);
                FORT_FILL_BATHYDZDZETA(
                    CHF_FRA1(nodeVals, 0),
                    CHF_BOX(nodeBox),
                    CHF_CONST_FRA1(seaFloorFAB, 0),
                    CHF_CONST_REALVECT(m_L),
                    CHF_CONST_REALVECT(a_dXi),
                    CHF_CONST_REAL(a_scale));

                // Convert to the desired centering
                FORT_CONVERTFAB(
                    CHF_FRA1(a_dest, a_destComp),
                    CHF_BOX(destBox),
                    CHF_CONST_INTVECT(destBoxType),
                    CHF_CONST_FRA1(nodeVals, 0),
                    CHF_CONST_INTVECT(nodeBoxType));
            }
        } // end if a_nu is horizontal or vertical
    } // end if a_mu is horizontal or vertical
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with J = det[Jacobian] = (dx/dXi) * (dy/dNu) * (dz/dZeta)
// -----------------------------------------------------------------------------
void BathymetricBaseMap::fill_J (FArrayBox&      a_dest,
                                 const int       a_destComp,
                                 const RealVect& a_dXi,
                                 const Real      a_scale) const
{
    if (isDiagonal()) {
        // This coordinate map is the result of a least-squares problem. The
        // Jacobian is complicated and needs to be computed numerically.
        GeoSourceInterface::fill_J(a_dest, a_destComp, a_dXi, a_scale);
    } else {
        // GeoSourceInterface::fill_J(a_dest, a_destComp, a_dXi, a_scale);

        // Initialize with dz/dZeta * scale
        this->fill_dxdXi(a_dest, a_destComp, SpaceDim-1, SpaceDim-1, a_dXi, a_scale);

        // Compute dx/dXi and dy/dNu, and accumulate product
        const Box destBox = a_dest.box();
        FArrayBox dxdXiFAB(destBox, 1);

        this->fill_dxdXi(dxdXiFAB, 0, 0, 0, a_dXi, 1.0);
        a_dest.mult(dxdXiFAB, destBox, destBox, 0, a_destComp, 1);

#if CH_SPACEDIM > 2
        this->fill_dxdXi(dxdXiFAB, 0, 1, 1, a_dXi, 1.0);
        a_dest.mult(dxdXiFAB, destBox, destBox, 0, a_destComp, 1);
#endif
    }
}

