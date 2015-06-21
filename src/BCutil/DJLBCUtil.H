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
#ifndef __DJLBCUTIL_H__INCLUDED__
#define __DJLBCUTIL_H__INCLUDED__

#include "PhysBCUtil.H"
#include "Constants.H"


// -----------------------------------------------------------------------------
// IC/BC info to read the DJL soln from file.
// -----------------------------------------------------------------------------
class DJLBCUtil: public PhysBCUtil
{
public:

    // Default constructor
    DJLBCUtil ();

    // Default destructor
    virtual ~DJLBCUtil ();

    // This object is its own factory
    virtual PhysBCUtil* newPhysBCUtil () const;

    // Fills a FAB with the initial velocity.
    // Locations are in mapped space, but components are Cartesian.
    virtual void setVelIC (FArrayBox&           a_velFAB,
                           const int            a_velComp,
                           const LevelGeometry& a_levGeo,
                           const DataIndex&     a_di) const;

    // Fills a FAB with the initial scalars
    virtual void setScalarIC (FArrayBox&           a_scalarFAB,
                              const int            a_scalarComp,
                              const LevelGeometry& a_levGeo,
                              const DataIndex&     a_di) const;

    // Fills a FAB with the background scalar
    virtual void setBackgroundScalar (FArrayBox&           a_scalarFAB,
                                      const int            a_scalarComp,
                                      const LevelGeometry& a_levGeo,
                                      const DataIndex&     a_di,
                                      const Real           a_time) const;

    // Sets the sponge ramp profile. This is a mapping [0,1] -> [0,1].
    inline virtual Real spongeLayerRamp (Real a_ratio) const;

    // Sets the Cartesian-based target velocity for the sponge layer.
    // By default, this function throws an error.
    virtual void fillVelSpongeLayerTarget (FArrayBox&           a_target,
                                           const int            a_velComp,
                                           const int            a_spongeDir,
                                           const Side::LoHiSide a_spongeSide,
                                           const LevelGeometry& a_levGeo,
                                           const DataIndex&     a_di,
                                           const Real           a_time);

protected:
    // Reads a_eta from the DJLIC_[a_nx]x[a_nz].bin file.
    // a_eta[i] is a Vector<Real> containing eta(x) at z[i].
    // (a_nx, a_nz) are the number of cell centers in the domain, not nodes!
    // Returns c.
    static Real readDJLICFile (Vector<Vector<Real> >& a_eta,
                               const int              a_nx,
                               const int              a_nz);

    // Computes the u DJL solution over a horizontal slice.
    // u = c*eta_z (but this assumes eta is already scaled by c)
    static void fill_uDJL (Vector<Real>&       a_uDJL,    // CC
                           const Vector<Real>& a_etaTop,  // NC
                           const Vector<Real>& a_etaBot,  // NC
                           const Real          a_dz);

    // Computes the w DJL solution over a horizontal slice.
    // w = -c*eta_x (but this assumes eta is already scaled by c)
    static void fill_wDJL (Vector<Real>&       a_wDJL,    // CC
                           const Vector<Real>& a_etaTop,  // NC
                           const Vector<Real>& a_etaBot,  // NC
                           const Real          a_dx);

    // Computes the b DJL solution over a horizontal slice.
    // b(x,z) = eta(z-eta(x)) (but this assumes eta is scaled by c)
    // Returns the CC background buoyancy at this slice.
    static Real fill_bDJL (Vector<Real>&       a_bDJL, // CC
                           const Vector<Real>& a_eta,  // CC
                           const Real          a_c,    // long wave speed
                           const Real          a_z);   // slice location

    // Converts a NC horizontal slice of data in a vector to CC.
    static void convertSliceNC2CC (Vector<Real>&       a_cc,
                                   const Vector<Real>& a_ncTop,
                                   const Vector<Real>& a_ncBot);

    // Envelope function for extrusion
    static Real extrusionEnvelope (const Real a_yprime);

    // // Sets physical BCs on velocities
    // virtual BCMethodHolder basicVelFuncBC (int a_veldir, bool a_isViscous) const;

    // // Sets physical BCs on a generic passive scalar
    // virtual BCMethodHolder basicScalarFuncBC () const;

    // /// Sets physical BCs on pressures (used by the Poisson solvers)
    // virtual BCMethodHolder basicPressureFuncBC (bool a_isHomogeneous) const;

    // Simply sets a_value to the background density at any given point, which
    // is all we need for the boundary conditions.
    // This function conforms to the EllipticBCValueFunc typedef.
    static void bscalBCValues (Real*           a_pos,
                               int*            a_dir,
                               Side::LoHiSide* a_side,
                               Real*           a_value,
                               Real            a_derivScale,
                               Real            a_time);

    // Member variables
    static RealVect s_L;

    static const Real s_d;        // Thickness of pycnocline
    static const Real s_z0;       // Location of pycnocline

    static const Real s_envSlope; // Envelope slope
    static const Real s_envWidth; // Envelope width
    static const Real s_offsetx;  // Pushes the IC into the domain
    static const Real s_offsety;  // Pushes the IC into the domain
    static const Real s_rotAngle;

    static Real s_rho_bottom;
    static Real s_rho_top;
    static Real s_rho_scale;
};



// -----------------------------------------------------------------------------
// Sets the sponge ramp profile. This is a mapping [0,1] -> [0,1].
// -----------------------------------------------------------------------------
inline Real DJLBCUtil::spongeLayerRamp (Real a_ratio) const
{
    return a_ratio*a_ratio;
}



#endif //!__DJLBCUTIL_H__INCLUDED__
