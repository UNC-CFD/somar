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
#include "VortexStreetBCUtil.H"
#include "ProblemContext.H"
#include "EllipticBCUtils.H"
#include "BoxIterator.H"


bool VortexStreetBCUtil::s_paramsRead = false;
RealVect VortexStreetBCUtil::s_inflowVel = RealVect::Zero;


// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
VortexStreetBCUtil::VortexStreetBCUtil ()
{
    if (s_paramsRead) return;

    const ProblemContext* ctx = ProblemContext::getInstance();
    s_inflowVel = ctx->inflowVel;

    s_paramsRead = true;
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
VortexStreetBCUtil::~VortexStreetBCUtil ()
{;}


// -----------------------------------------------------------------------------
// Factory
// -----------------------------------------------------------------------------
PhysBCUtil* VortexStreetBCUtil::newPhysBCUtil () const
{
    PhysBCUtil* newBCPtr = new VortexStreetBCUtil();
    return newBCPtr;
}


// ************************ ICs / background fields ****************************

// -----------------------------------------------------------------------------
// Fills a FAB with the initial velocity
// Locations are in mapped space, but components are Cartesian.
// -----------------------------------------------------------------------------
void VortexStreetBCUtil::setVelIC (FArrayBox&           a_velFAB,
                                   const int            a_velComp,
                                   const LevelGeometry& a_levGeo,
                                   const DataIndex&     a_di) const
{
    CH_assert(a_velFAB.nComp() == SpaceDim);
    CH_assert(0 <= a_velComp);
    CH_assert(a_velComp < SpaceDim);

    a_velFAB.setVal(s_inflowVel[a_velComp], a_velComp);

    // Add perturbations do initiate the vortices
    const Real pertMag = 1e-2;
    BoxIterator bit(a_velFAB.box());
    for (bit.reset(); bit.ok(); ++bit) {
        const IntVect& cc = bit();
        a_velFAB(cc, a_velComp) += pertMag * Real(rand()-RAND_MAX/2) / Real(RAND_MAX);
    }
}


// -----------------------------------------------------------------------------
// Sets the Cartesian-based target velocity for the sponge layer.
// -----------------------------------------------------------------------------
void VortexStreetBCUtil::fillVelSpongeLayerTarget (FArrayBox&           a_target,
                                                   const int            a_velComp,
                                                   const int            a_spongeDir,
                                                   const Side::LoHiSide a_spongeSide,
                                                   const LevelGeometry& a_levGeo,
                                                   const DataIndex&     a_di,
                                                   const Real           a_time)
{
    // Sanity checks
    CH_assert(useSpongeLayer());
    CH_assert(a_target.nComp() == 1);
    CH_assert(0 <= a_spongeDir);
    CH_assert(a_spongeDir < SpaceDim);

    a_target.setVal(s_inflowVel[a_velComp]);
}


// -----------------------------------------------------------------------------
// Sets the target values for the scalar sponge layer. If we are using a
// background scalar, then this function set the perturbation to zero.
// Otherwise, an error is thrown and this function will need to be overridden.
// -----------------------------------------------------------------------------
void VortexStreetBCUtil::fillScalarSpongeLayerTarget (FArrayBox&           a_target,
                                                      const int            a_scalarComp,
                                                      const int            a_spongeDir,
                                                      const Side::LoHiSide a_spongeSide,
                                                      const LevelGeometry& a_levGeo,
                                                      const DataIndex&     a_di,
                                                      const Real           a_time)
{
    // Sanity checks
    CH_assert(useSpongeLayer());
    CH_assert(a_target.nComp() == 1);
    CH_assert(0 <= a_spongeDir);
    CH_assert(a_spongeDir < SpaceDim);

    const char* msg = "VortexStreetBCUtil::fillScalarSpongeLayerTarget called";
    MayDay::Error(msg);
}


// -----------------------------------------------------------------------------
// basicVelFuncBC
// Sets physical BCs on velocities.
// -----------------------------------------------------------------------------
BCMethodHolder VortexStreetBCUtil::basicVelFuncBC (int a_veldir, bool a_isViscous) const
{
    // return PhysBCUtil::basicVelFuncBC(a_veldir, a_isViscous);
    BCMethodHolder holder;

    //            Diri 0
    //  Diri 0 |==========| Extrap
    //            Diri 0

    // Low order extrap in inflow/outflow (sponged) directions
    const IntVect inflowUnit = BASISV(0);
    const int extrapOrder = 0;

    // Extrap
    RefCountedPtr<BCGhostClass> horizBCPtr(
        new EllipticExtrapBCGhostClass(extrapOrder,
                                       IntVect::Zero,
                                       inflowUnit)
    );
    holder.addBCMethod(horizBCPtr);

    // No-slip
    RefCountedPtr<BCGhostClass> hiVertBCPtr = RefCountedPtr<BCGhostClass>(
        new BasicVelocityBCGhostClass(0.0,             // inflowVel
                                      -1,              // inflowDir
                                      Side::Lo,        // inflowSide
                                      -1,              // outflowDir
                                      Side::Hi,        // outflowSide
                                      a_veldir,
                                      a_isViscous,     // isViscous
                                      IntVect::Unit,
                                      IntVect::Unit - inflowUnit)
    );
    holder.addBCMethod(hiVertBCPtr);

    return holder;
}
