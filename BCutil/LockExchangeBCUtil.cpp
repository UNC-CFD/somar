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
#include "LockExchangeBCUtil.H"
#include "EllipticBCUtils.H"


// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
LockExchangeBCUtil::LockExchangeBCUtil ()
{;}


// -----------------------------------------------------------------------------
// Default destructor
// -----------------------------------------------------------------------------
LockExchangeBCUtil::~LockExchangeBCUtil ()
{;}


// -----------------------------------------------------------------------------
// This object is its own factory
// -----------------------------------------------------------------------------
PhysBCUtil* LockExchangeBCUtil::newPhysBCUtil () const
{
    PhysBCUtil* newBCPtr = new LockExchangeBCUtil();
    return newBCPtr;
}


// ************************ ICs / background fields ****************************



#include "CONSTANTS.H"
#include "BoxIterator.H"
// -----------------------------------------------------------------------------
// Fills a FAB with the initial scalars
// -----------------------------------------------------------------------------
void LockExchangeBCUtil::setScalarIC (FArrayBox&           a_scalarFAB,
                                      const int            a_scalarComp,
                                      const LevelGeometry& a_levGeo,
                                      const DataIndex&     a_di) const
{
    CH_assert(a_scalarFAB.nComp() == 1);

    if (a_scalarComp == 0) {
        // Gather geometric info
        const Box domBox = a_levGeo.getDomain().domainBox();
        const Box dataBox = a_scalarFAB.box();
        const RealVect& dx = a_levGeo.getDx();
        const IntVect ex = BASISV(0);

        // Compute Cartesian cell coordinates
        FArrayBox xposFAB(surroundingNodes(dataBox,0), 1);
        LevelGeometry::getGeoSourcePtr()->fill_physCoor(xposFAB, 0, 0, dx);

        FArrayBox yposFAB(dataBox, 1);
        LevelGeometry::getGeoSourcePtr()->fill_physCoor(yposFAB, 0, 1, dx);

        // Define IC parameters

        // Option #1: Interface at center of domain.
        // const Real x0 = Real(domBox.smallEnd(0)) * dx[0];
        // const Real xhalf = x0 + 0.5 * a_levGeo.getDomainLength(0);

        // Option #2: Interface at x = 0.
        const Real xhalf = 0.0;

        const Real pertA = 0.025;
        const Real pertK = 2.0 * Pi / a_levGeo.getDomainLength(1);
        const Real smoothingFactor = 2.0;

        const Real bmin = 0.0;
        const Real bmax = 1.0;

        // Loop over databox and set buoyancy IC.
        BoxIterator bit(dataBox);
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& cc = bit();
            Real xl = xposFAB(cc   ,0);
            Real xr = xposFAB(cc+ex,0);
            Real y  = yposFAB(cc   ,0);
            Real ifx = xhalf + pertA * sin(pertK * y);

            if (xr < ifx) {
                // Left values
                a_scalarFAB(cc,a_scalarComp) = bmin;
            } else if (ifx < xl) {
                // Right values
                a_scalarFAB(cc,a_scalarComp) = bmax;
            } else {
                // Interface is inside this cell. Use smoothing.
                Real frac = (ifx - xr) / (xl - xr);
                frac = 2.0 * frac - 1.0;
                frac = tanh(smoothingFactor * frac);
                a_scalarFAB(cc,a_scalarComp) = bmin + bmax * 0.5 * (frac + 1.0);
            }
        }

    } else {
        MayDay::Error("LockExchangeBCUtil::setScalarIC received a_scalarComp > 0");
    }
}


// -----------------------------------------------------------------------------
// basicScalarFuncBC   (Extrapolate BCs)
// Sets physical BCs on a generic passive scalar.
// Chombo uses 1st order extrap
// -----------------------------------------------------------------------------
BCMethodHolder LockExchangeBCUtil::basicScalarFuncBC () const
{
    BCMethodHolder holder;

    // Solid wall and outflow
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticExtrapBCGhostClass(1, // extrap order
                                       IntVect::Unit,
                                       IntVect::Unit)
    );
    holder.addBCMethod(BCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// basicVelFuncBC
// Sets physical BCs on velocities.
// -----------------------------------------------------------------------------
BCMethodHolder LockExchangeBCUtil::basicVelFuncBC (int a_veldir, bool a_isViscous) const
{
    bool isViscous = a_isViscous; // Set to false for free-slip BCs.

    BCMethodHolder holder;

    RefCountedPtr<BCGhostClass> velBCPtr = RefCountedPtr<BCGhostClass>(
        new BasicVelocityBCGhostClass(0.0,      //s_inflowVel,
                                      -1,       //s_inflowDir,
                                      Side::Lo, //s_inflowSide,
                                      -1,       //s_outflowDir,
                                      Side::Lo, //s_outflowSide,
                                      a_veldir,
                                      isViscous)
    );
    holder.addBCMethod(velBCPtr);

    return holder;
}
