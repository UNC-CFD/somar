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
#include "BoxIterator.H"
#include "MayDay.H"
#include "CellToEdge.H"
#include "GeoSourceInterface.H"
#include "GeoSourceInterfaceF_F.H"
#include "LevelGeometry.H"


// -----------------------------------------------------------------------------
// Default destructor
// -----------------------------------------------------------------------------
GeoSourceInterface::~GeoSourceInterface ()
{;}


// -----------------------------------------------------------------------------
// Fills a mapped box with Cartesian locations (a_dest must have SpaceDim comps)
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_physCoor (FArrayBox&      a_dest,
                                        const RealVect& a_dXi,
                                        const RealVect  a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_physCoor (all comps)");

    for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
        this->fill_physCoor(a_dest,
                            dir,
                            dir,
                            a_dXi);

        if (a_scale[dir] != 1.0) {
            a_dest.mult(a_scale[dir], dir, 1);
        }
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the Jacobian matrix elements d[x^mu] / d[Xi^nu].
// This is a speed bottleneck!!!
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_dxdXi (FArrayBox&      a_dest,
                                     const int       a_destComp,
                                     const int       a_mu,
                                     const int       a_nu,
                                     const RealVect& a_dXi,
                                     const Real      a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_dxdXi");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));
    CH_assert((0 <= a_nu) && (a_nu < SpaceDim));

    // Stagger the xmu FAB wrt the dest FAB.
    const Box& destBox = a_dest.box();
    const int destBoxType = destBox.type()[a_nu];

    Box xmuBox = destBox;
    if (destBoxType == 0) {
        xmuBox.surroundingNodes(a_nu);
    } else {
        xmuBox.grow(a_nu, 1);
        xmuBox.enclosedCells(a_nu);
    }

    // Fill a FAB with the mapping x^{mu}(Xi)
    FArrayBox xmu(xmuBox, 1);
    this->fill_physCoor(xmu, 0, a_mu, a_dXi);

    // Differentiate the mapping function and apply the scaling
    Real scaleOnDXi = a_scale / a_dXi[a_nu];

    if (destBoxType == 0) {
        FORT_SIMPLECCDERIV(
            CHF_FRA1(a_dest,a_destComp),
            CHF_CONST_FRA1(xmu,0),
            CHF_CONST_REAL(scaleOnDXi),
            CHF_CONST_INT(a_nu),
            CHF_BOX(destBox));
    } else {
        FORT_SIMPLEFCDERIV(
            CHF_FRA1(a_dest,a_destComp),
            CHF_CONST_FRA1(xmu,0),
            CHF_CONST_REAL(scaleOnDXi),
            CHF_CONST_INT(a_nu),
            CHF_BOX(destBox));
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with J = det[Jacobian]
// This is a speed bottleneck!!!
// TODO: try calling dxdXi and using its results. This way, the speed
// bottleneck may be circumvented by rewriting one function instead of two.
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_J (FArrayBox&      a_dest,
                                 const int       a_destComp,
                                 const RealVect& a_dXi,
                                 const Real      a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_J");

    // Sanity check
    CH_assert(a_dest.interval().contains(a_destComp));

    Box destBox = a_dest.box();
    if (destBox.type() == IntVect::Zero) {
        // Calculate the pre-mapped coordinates
        FluxBox x(destBox, SpaceDim);
        for (int nu = 0; nu < CH_SPACEDIM; ++nu) {
            this->fill_physCoor(x[nu], a_dXi);
        }

        // Calculate J (unscaled)
#       if CH_SPACEDIM == 2
            FORT_DEFAULT_FILL_J_2D(
                CHF_FRA1(a_dest,a_destComp),
                CHF_CONST_FRA(x[0]),
                CHF_CONST_FRA(x[1]),
                CHF_BOX(destBox));
#       elif CH_SPACEDIM == 3
            FORT_DEFAULT_FILL_J_3D(
                CHF_FRA1(a_dest,a_destComp),
                CHF_CONST_FRA(x[0]),
                CHF_CONST_FRA(x[1]),
                CHF_CONST_FRA(x[2]),
                CHF_BOX(destBox));
#       else
#           error Bad SPACEDIM
#       endif
    } else {
        // This can only handle destBoxes that are FC in at most 1 dir.
        CH_assert(destBox.type().sum() == 1);

        // Find the FC direction
        D_TERM(
        int FCdir = 0;,
        if (destBox.type()[1] != 0) FCdir = 1;,
        else if (destBox.type()[2] != 0) FCdir = 2;)

        // We will need a CC version of J that we can average to FC.
        destBox.grow(FCdir, 1);
        destBox.enclosedCells(FCdir);

        // Calculate the pre-mapped coordinates
        FluxBox x(destBox, SpaceDim);
        for (int nu = 0; nu < CH_SPACEDIM; ++nu) {
            this->fill_physCoor(x[nu], a_dXi);
        }

        // Calculate the CC version of J (unscaled)
        FArrayBox CCJ(destBox, 1);
#       if CH_SPACEDIM == 2
            FORT_DEFAULT_FILL_J_2D(
                CHF_FRA1(CCJ,0),
                CHF_CONST_FRA(x[0]),
                CHF_CONST_FRA(x[1]),
                CHF_BOX(destBox));
#       elif CH_SPACEDIM == 3
            FORT_DEFAULT_FILL_J_3D(
                CHF_FRA1(CCJ,0),
                CHF_CONST_FRA(x[0]),
                CHF_CONST_FRA(x[1]),
                CHF_CONST_FRA(x[2]),
                CHF_BOX(destBox));
#       else
#           error Bad SPACEDIM
#       endif

        // Average to an FC J
        CellToEdge(CCJ, 0, a_dest, a_destComp, FCdir);
    }

    // Scale the result
    a_dest.mult(a_scale / a_dXi.product(), a_destComp, 1);
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the inverse Jacobian matrix elements d[xi^mu] / d[x^nu].
// TODO: Make this faster by requesting detJ.
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_dXidx (FArrayBox&       a_dest,
                                     const int        a_destComp,
                                     const int        a_mu,
                                     const int        a_nu,
                                     const RealVect&  a_dXi,
                                     const Real       a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_dXidx");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));
    CH_assert((0 <= a_nu) && (a_nu < SpaceDim));

    // CH_assert(SpaceDim == 3);   // TODO: Create a 2D version of this function
    const Box& destBox = a_dest.box();

#if CH_SPACEDIM == 2
    // Get the corresponding Jacobian matrix element
    int mu1 = (a_nu + 1) % SpaceDim;
    int nu1 = (a_mu + 1) % SpaceDim;
    this->fill_dxdXi(a_dest, a_destComp, mu1, nu1, a_dXi);

    // Divide by Det[J]
    FArrayBox detJ(destBox, 1);
    this->fill_J(detJ, 0, a_dXi);
    a_dest.divide(detJ, 0, a_destComp, 1);    // Changed to divide from mult on May 7, 2013.

    // Flip sign of anti-diagonal elements and apply the scaling.
    if ((a_mu + a_nu) % 2 == 1) {
        a_dest.mult(-a_scale, a_destComp, 1);
    } else if (a_scale != 1.0) {
        a_dest.mult(a_scale, a_destComp, 1);
    }

#elif CH_SPACEDIM == 3
    // Start with zero
    a_dest.setVal(0.0, a_destComp);

    // Fill a FAB with the Jacobian matrix elements for the first term
    FArrayBox dxdXi(destBox, 2);

    int mu1 = (a_nu + 1) % SpaceDim;
    int mu2 = (a_nu + 2) % SpaceDim;

    int nu1 = (a_mu + 1) % SpaceDim;
    int nu2 = (a_mu + 2) % SpaceDim;

    this->fill_dxdXi(dxdXi, 0, mu1, nu1, a_dXi);
    this->fill_dxdXi(dxdXi, 1, mu2, nu2, a_dXi);

    // Add on first term
    FORT_ADDPROD2(
        CHF_FRA1(a_dest,a_destComp),
        CHF_CONST_FRA1(dxdXi,0),
        CHF_CONST_FRA1(dxdXi,1),
        CHF_BOX(destBox));

    // Fill a FAB with the Jacobian matrix elements for the second term
    this->fill_dxdXi(dxdXi, 0, mu1, nu2, a_dXi);
    this->fill_dxdXi(dxdXi, 1, mu2, nu1, a_dXi);

    // Subtract second term
    FORT_SUBPROD2(
        CHF_FRA1(a_dest,a_destComp),
        CHF_CONST_FRA1(dxdXi,0),
        CHF_CONST_FRA1(dxdXi,1),
        CHF_BOX(destBox));

    // Divide by Det[J]
    FArrayBox detJ(Interval(0,0), dxdXi);
    this->fill_J(detJ, 0, a_dXi);
    a_dest.divide(detJ, 0, a_destComp, 1);  // Changed to divide from mult on May 7, 2013.

    // Apply the scaling
    if (a_scale != 1.0) {
        a_dest.mult(a_scale, a_destComp, 1);
    }

#else
#   error Bad Spacedim
#endif
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with 1/J
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_Jinv (FArrayBox&       a_dest,
                                    const int        a_destComp,
                                    const RealVect&  a_dXi,
                                    const Real       a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_Jinv");

    // Sanity check
    CH_assert(a_dest.interval().contains(a_destComp));

    // Calculate or copy J
    this->fill_J(a_dest, a_destComp, a_dXi);

    // Invert and apply scaling
    a_dest.invert(a_scale);
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the covariant metric elements
// gdn_{mu,nu} = Sum over rho [ dx^{rho}/dXi^{mu} * dx^{rho}/dXi^{nu} ]
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_gdn (FArrayBox&      a_dest,
                                   const int       a_destComp,
                                   const int       a_mu,
                                   const int       a_nu,
                                   const RealVect& a_dXi,
                                   const Real      a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_gdn");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));
    CH_assert((0 <= a_nu) && (a_nu < SpaceDim));

    // Start with zero
    a_dest.setVal(0.0, a_destComp);

    // If the metric is diagonal (orthogonal) and we are evaluating
    // an off-diagonal element, then we are done.
    if (this->isDiagonal() && (a_mu != a_nu)) return;

    // Create space for Jacobian matrix elements
    const Box& destBox = a_dest.box();
    FArrayBox dxdXi(destBox, 2);

    for (int rho = 0; rho < CH_SPACEDIM; ++rho) {
        // Calculate dx^{rho}/dXi^{mu} and dx^{rho}/dXi^{nu}
        this->fill_dxdXi(dxdXi, 0, rho, a_mu, a_dXi);
        this->fill_dxdXi(dxdXi, 1, rho, a_nu, a_dXi);

        // Add product to dest
        FORT_ADDPROD2(
            CHF_FRA1(a_dest,a_destComp),
            CHF_CONST_FRA1(dxdXi,0),
            CHF_CONST_FRA1(dxdXi,1),
            CHF_BOX(destBox));
    }

    // Apply the scaling
    if (a_scale != 1.0) {
        a_dest.mult(a_scale, a_destComp, 1);
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the contravariant metric elements
// gup^{mu,nu} = Sum over rho [ dXi^{mu}/dx^{rho} * dXi^{nu}/dx^{rho} ]
// TODO: Make this faster by requesting detJ.
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_gup (FArrayBox&       a_dest,
                                   const int        a_destComp,
                                   const int        a_mu,
                                   const int        a_nu,
                                   const RealVect&  a_dXi,
                                   const Real       a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_gup");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));
    CH_assert((0 <= a_nu) && (a_nu < SpaceDim));

    // Start with zero
    a_dest.setVal(0.0, a_destComp);

    // If the metric is diagonal (orthogonal) and we are evaluating
    // an off-diagonal element, then we are done.
    if (this->isDiagonal() && (a_mu != a_nu)) return;

    // Create space for Jacobian matrix elements
    const Box& destBox = a_dest.box();
    FArrayBox dXidx(destBox, 2);

    for (int rho = 0; rho < CH_SPACEDIM; ++rho) {
        // Calculate dXi^{mu}/dx^{rho} and dXi^{nu}/dx^{rho}.
        this->fill_dXidx(dXidx, 0, a_mu, rho, a_dXi);
        this->fill_dXidx(dXidx, 1, a_nu, rho, a_dXi);

        // Add product to dest
        FORT_ADDPROD2(
            CHF_FRA1(a_dest,a_destComp),
            CHF_CONST_FRA1(dXidx,0),
            CHF_CONST_FRA1(dXidx,1),
            CHF_BOX(destBox));
    }

    // Apply the scaling
    if (a_scale != 1.0) {
        a_dest.mult(a_scale, a_destComp, 1);
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with detJ * gup
// TODO: Make this faster by requesting detJ.
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_Jgup (FArrayBox&       a_dest,
                                    const int        a_destComp,
                                    const int        a_mu,
                                    const int        a_nu,
                                    const RealVect&  a_dXi,
                                    const Real       a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_Jgup");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));
    CH_assert((0 <= a_nu) && (a_nu < SpaceDim));

    // If the metric is diagonal (orthogonal) and we are evaluating
    // an off-diagonal element, then fill a_dest with zeros and return.
    if (this->isDiagonal() && (a_mu != a_nu)) {
        a_dest.setVal(0.0, a_destComp);
        return;
    }

    // Calculate gup
    this->fill_gup(a_dest, a_destComp, a_mu, a_nu, a_dXi);

    // Multiply by detJ
    FArrayBox detJ(a_dest.box(), 1);
    this->fill_J(detJ, 0, a_dXi);
    a_dest.mult(detJ, 0, a_destComp, 1);

    // Apply the scaling
    if (a_scale != 1.0) {
        a_dest.mult(a_scale, a_destComp, 1);
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the connection elements.
// Note that although this function is provided and works, it is never used
// by our NS algorithm.
// -----------------------------------------------------------------------------
void GeoSourceInterface::fill_Gamma (FArrayBox&      a_dest,
                                     const int       a_destComp,
                                     const int       a_up,
                                     const int       a_dn1,
                                     const int       a_dn2,
                                     const RealVect& a_dXi,
                                     const Real      a_scale) const
{
    CH_TIME("GeoSourceInterface::fill_Gamma");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_up) && (a_up < SpaceDim));
    CH_assert((0 <= a_dn1) && (a_dn1 < SpaceDim));
    CH_assert((0 <= a_dn2) && (a_dn2 < SpaceDim));

    const int symComps = SpaceDim * (SpaceDim+1) / 2;
    const Box destBox = a_dest.box();

    // Calculate gup
    FArrayBox gup(destBox, SpaceDim);
    for (int mu = 0; mu < CH_SPACEDIM; ++mu) {
        this->fill_gup(gup, mu, a_up, mu, a_dXi);
    }

    // Calculate gdn
    Box gdnBox = destBox;
    gdnBox.grow(IntVect::Unit);
    FArrayBox gdn(gdnBox, symComps);
    for (int mu = 0; mu < CH_SPACEDIM; ++mu) {
        for (int nu = mu; nu < CH_SPACEDIM; ++nu) {
            int comp = LevelGeometry::symTensorCompCC(mu, nu);
            this->fill_gdn(gdn, comp, mu, nu, a_dXi);
        }
    }

    // Calculate the gamma entry
    FORT_DEFAULT_FILL_GAMMA(
        CHF_FRA1(a_dest,a_destComp),
        CHF_CONST_FRA(gup),
        CHF_CONST_FRA(gdn),
        CHF_CONST_INT(a_up),
        CHF_CONST_INT(a_dn1),
        CHF_CONST_INT(a_dn2),
        CHF_CONST_REALVECT(a_dXi),
        CHF_CONST_REAL(a_scale),
        CHF_BOX(destBox));
}
