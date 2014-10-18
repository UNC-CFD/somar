/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Copyright (C) 2014 Edward Santilli & Alberto Scotti
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
#include "InternalWaveBCUtil.H"
#include "BoxIterator.H"
#include "Debug.H"


// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
InternalWaveBCUtil::InternalWaveBCUtil ()
{;}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
InternalWaveBCUtil::~InternalWaveBCUtil ()
{;}


// -----------------------------------------------------------------------------
// Factory
// -----------------------------------------------------------------------------
PhysBCUtil* InternalWaveBCUtil::newPhysBCUtil () const
{
    PhysBCUtil* newBCPtr = new InternalWaveBCUtil();
    return newBCPtr;
}


// -----------------------------------------------------------------------------
// This is in case the BC's have an effect on the timestep.
// Pass in currently computed dt, along with the cfl and dx. If the effect
// of the BCs requires a decreased timestep, then the newly reduced timestep
// is returned. In the default case, this just returns a_dt back; however,
// derived classes may actually have an effect.
// -----------------------------------------------------------------------------
void InternalWaveBCUtil::computeBoundaryDt (Real&                a_dt,
                                            const Real           a_cfl,
                                            const LevelGeometry& a_levGeo) const
{
    // Do nothing.
}


// -----------------------------------------------------------------------------
// Fills a FAB with the initial scalars
// -----------------------------------------------------------------------------
void InternalWaveBCUtil::setScalarIC (FArrayBox&           a_scalarFAB,
                                      const int            a_scalarComp,
                                      const LevelGeometry& a_levGeo,
                                      const DataIndex&     a_di) const
{
    // Sanity checks
    CH_assert(a_scalarFAB.nComp() == 1);
    CH_assert(a_scalarComp == 0);

    // Set the background density. This uses the setBackgroundScalar
    // function. This way, if we decide to change the background buoyancy,
    // we only need to rewrite one piece of code!
    {
        const bool useBkgdSave = s_useBackgroundScalar;
        s_useBackgroundScalar = true;
        this->setBackgroundScalar(a_scalarFAB, a_scalarComp, a_levGeo, a_di, 0.0);
        s_useBackgroundScalar = useBkgdSave;
    }

    // Set the mixed region.
#   define USE_AS_VERSION
#   ifdef USE_AS_VERSION
    {
        const Box& domBox = a_levGeo.getDomain().domainBox();
        int chopPoint = domBox.bigEnd(0) - domBox.size(0) / 8;

        Box mixedBox_x = domBox;
        mixedBox_x = mixedBox_x.chop(0, chopPoint);
        chopPoint = domBox.bigEnd(CH_SPACEDIM-1) - domBox.size(CH_SPACEDIM-1) / 2;

        Box mixedBox_z = domBox;
        mixedBox_z = mixedBox_z.chop(CH_SPACEDIM-1, chopPoint);
        Box mixedBox = (mixedBox_z & mixedBox_x) & a_scalarFAB.box();
        if (!mixedBox.isEmpty()) {
            a_scalarFAB.setVal(-0.5, mixedBox, a_scalarComp);
        }
    }
#   else // Use ES version
    {
        const Box& domBox = a_levGeo.getDomain().domainBox();
        int chopPoint = domBox.bigEnd(0) - domBox.size(0) / 64;
        Box mixedBox = domBox;
        mixedBox = mixedBox.chop(0, chopPoint);
        mixedBox &= a_scalarFAB.box();

        if (!mixedBox.isEmpty()) {
            a_scalarFAB.setVal(0.0, mixedBox, a_scalarComp);
        }
    }
#   endif
}


// -----------------------------------------------------------------------------
// Fills a FAB with the background scalar
// -----------------------------------------------------------------------------
void InternalWaveBCUtil::setBackgroundScalar (FArrayBox&           a_scalarFAB,
                                              const int            a_scalarComp,
                                              const LevelGeometry& a_levGeo,
                                              const DataIndex&     a_di,
                                              const Real           a_time) const
{
    CH_assert(a_scalarFAB.nComp() == 1);

    if (s_useBackgroundScalar && a_scalarComp == 0) {
        // Get Cartesian coordinates at each point of a_scalarFAB.
        FArrayBox posFAB(a_scalarFAB.box(), 1);
        const RealVect& dx = a_levGeo.getDx();
        a_levGeo.getGeoSourcePtr()->fill_physCoor(posFAB, 0, SpaceDim-1, dx);

        // Loop over each point of a_scalarFAB and set to background value.
        BoxIterator bit(a_scalarFAB.box());
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& iv = bit();

            Real pos[CH_SPACEDIM];
            pos[CH_SPACEDIM-1] = posFAB(iv);

            int dirDummy;
            Side::LoHiSide sideDummy;
            Real derivScaleDummy = 1.0e300;
            Real value[1];
            this->bscalBCValues(pos, &dirDummy, &sideDummy, value, derivScaleDummy, a_time);
            a_scalarFAB(iv,a_scalarComp) = value[0];
        }
    } else {
        // We do not have a background scalar.
        a_scalarFAB.setVal(0.0);
    }
}


// -----------------------------------------------------------------------------
// Sets a_value to the boundary/background state value at a_pos.
// -----------------------------------------------------------------------------
void InternalWaveBCUtil::bscalBCValues (Real*           a_pos,
                                        int*            a_dir,
                                        Side::LoHiSide* a_side,
                                        Real*           a_value,
                                        Real            a_derivScale,
                                        Real            a_time)
{
     const Real& H = s_domLength[CH_SPACEDIM-1];
     const Real drho = 0.5;
     const Real z_pycnocline = .75;
     const Real Thickness = .1;
     a_value[0] = -drho * tanh((a_pos[CH_SPACEDIM-1]/H - z_pycnocline) / Thickness);
 }

