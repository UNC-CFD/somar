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
#include "DJLBCUtil.H"
#include "ProblemContext.H"
#include "BoxIterator.H"
#include <iostream>
#include <fstream>
#include <sstream>
#include "ConvertFABF_F.H"
#include "EllipticBCUtils.H"

#include "AMRIO.H"
#include "CubicSpline.H"
#include "Debug.H"


// Static variable definitions...
RealVect DJLBCUtil::s_L = RealVect::Zero;

Real DJLBCUtil::s_rho_bottom = quietNAN;
Real DJLBCUtil::s_rho_top = quietNAN;
Real DJLBCUtil::s_rho_scale = quietNAN;

const Real DJLBCUtil::s_d = 0.1;        // Thickness of pycnocline
const Real DJLBCUtil::s_z0 = 0.8;       // Location of pycnocline

// // For non-oblique 3D problem
// const Real DJLBCUtil::s_envSlope = 1.0 / 6.4;   // Envelope slope
// const Real DJLBCUtil::s_envWidth = 192.0;       // Envelope width
// const Real DJLBCUtil::s_offsetx = 0.0;          // Pushes the IC into the domain
// const Real DJLBCUtil::s_offsety = 128.0;        // Pushes the IC into the domain
// const Real DJLBCUtil::s_rotAngle = 0.0 * Pi/180.0;

// For oblique 3D problem
const Real DJLBCUtil::s_envSlope = 1.0 / 6.4;   // Envelope slope
const Real DJLBCUtil::s_envWidth = 192.0;       // Envelope width
const Real DJLBCUtil::s_offsetx = 128.0;        // Pushes the IC into the domain
const Real DJLBCUtil::s_offsety = 128.0;        // Pushes the IC into the domain
const Real DJLBCUtil::s_rotAngle = 45.0 * Pi/180.0;



// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
DJLBCUtil::DJLBCUtil ()
{
    static bool paramsRead = false;
    if (!paramsRead) {
        const ProblemContext* ctx = ProblemContext::getInstance();

        s_L = ctx->domainLength;

        // TODO: Read all static parameters from input file.

        s_rho_top    = 0.5 * (1.0 - tanh((1.0 - s_z0) / s_d));
        s_rho_bottom = 0.5 * (1.0 - tanh((0.0 - s_z0) / s_d));
        s_rho_scale = 1.0 / (s_rho_bottom - s_rho_top);
        paramsRead = true;
    }
}


// -----------------------------------------------------------------------------
// Default destructor
// -----------------------------------------------------------------------------
DJLBCUtil::~DJLBCUtil ()
{;}


// -----------------------------------------------------------------------------
// This object is its own factory
// -----------------------------------------------------------------------------
PhysBCUtil* DJLBCUtil::newPhysBCUtil () const
{
    PhysBCUtil* newBCPtr = new DJLBCUtil();
    return newBCPtr;
}


// -----------------------------------------------------------------------------
// Fills a FAB with the initial velocity
// Locations are in mapped space, but components are Cartesian.
// -----------------------------------------------------------------------------
void DJLBCUtil::setVelIC (FArrayBox&           a_velFAB,
                          const int            a_velComp,
                          const LevelGeometry& a_levGeo,
                          const DataIndex&     a_di) const
{
    // Sanity checks
    CH_assert(a_velFAB.nComp() == SpaceDim);
    CH_assert(0 <= a_velComp);
    CH_assert(a_velComp < SpaceDim);
    CH_assert(a_velFAB.box().type() == IntVect::Zero);

    // Gather domain data
    const Box domBox = a_levGeo.getDomain().domainBox();
    const Box valid = a_levGeo.getBoxes()[a_di] & a_velFAB.box();
    const RealVect physDx = a_levGeo.getDx();
    const IntVect& Nx = domBox.size();
    const RealVect L = a_levGeo.getDomainLength();

    const Real sinA = sin(s_rotAngle);
    const Real cosA = cos(s_rotAngle);

    // Read eta from file.
    Vector<Vector<Real> > eta(Nx[SpaceDim-1]+1, Vector<Real>(Nx[0]+1, 0.0));
    readDJLICFile(eta, Nx[0], Nx[SpaceDim-1]);

    // Compute locations of cell centers.
    Vector<Real> x(Nx[0]);
    for (int i = 0; i < Nx[0]; ++i) {
        x[i] = (Real(domBox.smallEnd(0) + i) + 0.5) * physDx[0];
    }
    const Real minX = x[0];
    const Real maxX = x[Nx[0]-1];


    // Loop over horizontal slices
    IntVect cc = valid.smallEnd();
    int k = cc[SpaceDim-1] - domBox.smallEnd(SpaceDim-1);
    for (; cc[SpaceDim-1] <= valid.bigEnd(SpaceDim-1); ++cc[SpaceDim-1], ++k) {

        // Construct CC DJL velocity
        const Real thisZ = (Real(cc[SpaceDim-1]) + 0.5) * physDx[SpaceDim-1];
        Vector<Real> velDJL;
        if (a_velComp == SpaceDim-1) {
            fill_wDJL(velDJL, eta[k+1], eta[k], physDx[0]);
        } else {
            fill_uDJL(velDJL, eta[k+1], eta[k], physDx[SpaceDim-1]);
        }

#if CH_SPACEDIM == 2
        // Loop over this horizontal slice and fill FAB.
        // No rotations needed.
        cc[0] = valid.smallEnd(0);
        int i = cc[0] - domBox.smallEnd(0);
        for (; cc[0] <= valid.bigEnd(0); ++cc[0], ++i) {
            a_velFAB(cc, a_velComp) = velDJL[i];
        }

#else // CH_SPACEDIM == 3
        // Construct splines of DJL velocity
        CubicSpline velDJLSpline;
        velDJLSpline.solve(velDJL, x);

        // Loop over this horizontal slice
        cc[0] = valid.smallEnd(0);
        int i = cc[0] - domBox.smallEnd(0);
        for (; cc[0] <= valid.bigEnd(0); ++cc[0], ++i) {

            cc[1] = valid.smallEnd(1);
            int j = cc[1] - domBox.smallEnd(1);
            for (; cc[1] <= valid.bigEnd(1); ++cc[1], ++j) {

                // Compute this cell's location.
                Real thisX = (Real(cc[0]) + 0.5) * physDx[0] - s_offsetx;
                Real thisY = (Real(cc[1]) + 0.5) * physDx[1] - s_offsety;

                // Compute the location of the "source" DJL solution
                // and our distance away from the center of extrusion.
                Real xprime =  thisX*cosA + thisY*sinA;
                Real yprime = -thisX*sinA + thisY*cosA;

                // Interpolate the DJL solution at the source location.
                Real srcVel = ((minX <= xprime && xprime <= maxX)? velDJLSpline.interp(xprime): 0.0);

                // Now, rotate the source solution into position
                if (a_velComp == 0) {
                    srcVel *= cosA;
                } else if (a_velComp == 1) {
                    srcVel *= sinA;
                }

                // Compute envelope at this distance from the center of extrusion.
                Real env = extrusionEnvelope(yprime);

                // Set 3D field.
                a_velFAB(cc, a_velComp) = env * srcVel;

            } // end loop over y (cc[1] and j)
        } // end loop over x (cc[0] and i)
#endif // CH_SPACEDIM == 2 or 3
    } // end loop over z (cc[2] and k)
}


// -----------------------------------------------------------------------------
// Fills a FAB with the initial scalars
// -----------------------------------------------------------------------------
void DJLBCUtil::setScalarIC (FArrayBox&           a_scalarFAB,
                             const int            a_scalarComp,
                             const LevelGeometry& a_levGeo,
                             const DataIndex&     a_di) const
{
    // Sanity checks
    CH_assert(a_scalarFAB.nComp() == 1);
    CH_assert(a_scalarComp == 0);
    CH_assert(a_scalarFAB.box().type() == IntVect::Zero);

    // Gather domain data
    const Box domBox = a_levGeo.getDomain().domainBox();
    const Box valid = a_levGeo.getBoxes()[a_di] & a_scalarFAB.box();
    const RealVect physDx = a_levGeo.getDx();
    const IntVect& Nx = domBox.size();
    const RealVect L = a_levGeo.getDomainLength();

    const Real sinA = sin(s_rotAngle);
    const Real cosA = cos(s_rotAngle);

    // Read eta from file.
    Vector<Vector<Real> > eta(Nx[SpaceDim-1]+1, Vector<Real>(Nx[0]+1, 0.0));
    const Real c = readDJLICFile(eta, Nx[0], Nx[SpaceDim-1]);

    // Compute locations of cell centers.
    Vector<Real> x(Nx[0]);
    for (int i = 0; i < Nx[0]; ++i) {
        x[i] = (Real(domBox.smallEnd(0) + i) + 0.5) * physDx[0];
    }
    const Real minX = x[0];
    const Real maxX = x[Nx[0]-1];


    // Loop over horizontal slices
    IntVect cc = valid.smallEnd();
    int k = cc[SpaceDim-1] - domBox.smallEnd(SpaceDim-1);
    for (; cc[SpaceDim-1] <= valid.bigEnd(SpaceDim-1); ++cc[SpaceDim-1], ++k) {

        // Construct CC DJL buoyancy
        const Real thisZ = (Real(cc[SpaceDim-1]) + 0.5) * physDx[SpaceDim-1];
        Vector<Real> ccEta, bDJL;
        convertSliceNC2CC(ccEta, eta[k+1], eta[k]);
        const Real bgScalar = fill_bDJL(bDJL, ccEta, c, thisZ);


#if CH_SPACEDIM == 2
        // Loop over this horizontal slice and fill FAB.
        // No rotations needed.
        cc[0] = valid.smallEnd(0);
        int i = cc[0] - domBox.smallEnd(0);
        for (; cc[0] <= valid.bigEnd(0); ++cc[0], ++i) {
            a_scalarFAB(cc,0) = bDJL[i];
        }

#else // CH_SPACEDIM == 3
        // Construct splines of DJL buoyancy
        CubicSpline bDJLSpline;
        bDJLSpline.solve(bDJL, x);

        // Loop over this horizontal slice
        cc[0] = valid.smallEnd(0);
        for (; cc[0] <= valid.bigEnd(0); ++cc[0]) {

            cc[1] = valid.smallEnd(1);
            for (; cc[1] <= valid.bigEnd(1); ++cc[1]) {

                // Compute this cell's location.
                Real thisX = (Real(cc[0]) + 0.5) * physDx[0] - s_offsetx;
                Real thisY = (Real(cc[1]) + 0.5) * physDx[1] - s_offsety;

                // Compute the location of the "source" DJL solution
                // and our distance away from the center of extrusion.
                Real xprime =  thisX*cosA + thisY*sinA;
                Real yprime = -thisX*sinA + thisY*cosA;

                // Interpolate the DJL solution at the source location.
                Real srcB = ((minX <= xprime && xprime <= maxX)? bDJLSpline.interp(xprime): bgScalar);

                // Compute envelope at this distance from the center of extrusion.
                Real env = extrusionEnvelope(yprime);

                // Set 3D field.
                a_scalarFAB(cc,0) = env*srcB + (1.0-env)*bgScalar;

            } // end loop over y (cc[1])
        } // end loop over x (cc[0])
#endif // CH_SPACEDIM == 2 or 3
    } // end loop over z (cc[2] and k)
}


// -----------------------------------------------------------------------------
// Fills a FAB with the background scalar
// -----------------------------------------------------------------------------
void DJLBCUtil::setBackgroundScalar (FArrayBox&           a_scalarFAB,
                                     const int            a_scalarComp,
                                     const LevelGeometry& a_levGeo,
                                     const DataIndex&     a_di,
                                     const Real           a_time) const
{
    CH_assert(a_scalarFAB.nComp() == 1);

    if (s_useBackgroundScalar && a_scalarComp == 0) {
        // Gather Cartesian coordinates.
        FArrayBox posFAB(a_scalarFAB.box(), 1);
        const RealVect& dx = a_levGeo.getDx();
        a_levGeo.getGeoSourcePtr()->fill_physCoor(posFAB, 0, SpaceDim-1, dx);

        // Loop over a_scalarFAB and set background scalar values.
        // We get the values from the bscalBCValues function.
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
        // This scalar does not have a background state.
        a_scalarFAB.setVal(0.0, a_scalarComp);
    }
}


// -----------------------------------------------------------------------------
// Simply sets a_value to the background density at any given point, which
// is all we need for the boundary conditions.
// This function conforms to the EllipticBCValueFunc typedef.
// -----------------------------------------------------------------------------
void DJLBCUtil::bscalBCValues (Real*           a_pos,
                               int*            a_dir,
                               Side::LoHiSide* a_side,
                               Real*           a_value,
                               Real            a_derivScale,
                               Real            a_time)
{
    Real z = a_pos[CH_SPACEDIM-1];
    Real rho = 0.5 * (1.0 - tanh((z - s_z0) / s_d));
    a_value[0] = (rho - s_rho_top) * s_rho_scale;
}


// -----------------------------------------------------------------------------
// Sets the (vector scaled) target velocity for the sponge layer. By default,
// this function persuades the velocity field to approach its inflow value.
// -----------------------------------------------------------------------------
void DJLBCUtil::fillVelSpongeLayerTarget (FArrayBox&           a_target,
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

    // Assume the boundaries are in the far-field where not much is happening.
    if (a_spongeDir < SpaceDim-1) {
        a_target.setVal(0.0);
    } else {
        MayDay::Error("DJLBCUtil::fillVelSpongeLayerTarget "
                      "can only set a sponge target when a_spongeDir < SpaceDim-1");
    }
}


// -----------------------------------------------------------------------------
// Static utility
// Reads a_eta from the DJLIC_[a_nx]x[a_nz].bin file.
// a_eta[i] is a Vector<Real> containing eta(x) at z[i].
// (a_nx, a_nz) are the number of cell centers in the domain, not nodes!
// Returns c.
// -----------------------------------------------------------------------------
Real DJLBCUtil::readDJLICFile (Vector<Vector<Real> >& a_eta,
                               const int              a_nx,
                               const int              a_nz)
{
    const ProblemContext* ctx = ProblemContext::getInstance();
    const int L = int(ctx->domainLength[0]);
    const int H = int(ctx->domainLength[SpaceDim-1]);

    char infileName[100];
    sprintf(infileName, "DJLIC_%don%dx%don%d.bin", L, a_nx, H, a_nz);
    pout() << "infileName = " << infileName << endl;
    std::ifstream infile;
    infile.open(infileName, ios::in | ios::binary);

    if (!infile.is_open()) {
        std::ostringstream errmsg;
        errmsg << "Could not open " << infileName;
        MayDay::Error(errmsg.str().c_str());
    }

    infile.seekg(0, ios::beg);

    // nmax
    double nmax = 0.0;
    infile.seekg(4, ios::cur);
    infile.read((char*)&nmax, sizeof(double));
    infile.seekg(4, ios::cur);

    // c
    double c = 0.0;
    infile.seekg(4, ios::cur);
    infile.read((char*)&c, sizeof(double));
    pout() << "c = " << c << endl;
    infile.seekg(4, ios::cur);

    // x
    Vector<double> x(a_nx+1, 0.0);
    infile.seekg(4, ios::cur);
    infile.read((char*)&x[0], sizeof(double)*x.size());
    infile.seekg(4, ios::cur);

    Real fileDx = x[1] - x[0];
    Real thisDx = s_L[0] / a_nx;
    if (abs(fileDx-thisDx) > 1.0e-9) {
        pout() << "fileDx = " << fileDx << "\tthisDx = " << thisDx << endl;
        MayDay::Error("dx is not properly set");
    }

    // z
    Vector<double> z(a_nz+1, 0.0);
    infile.seekg(4, ios::cur);
    infile.read((char*)&z[0], sizeof(double)*z.size());
    infile.seekg(4, ios::cur);

    Real fileDz = z[1] - z[0];
    Real thisDz = 1.0 / a_nz;
    if (abs(fileDz-thisDz) > 1.0e-9) {
        pout() << "fileDz = " << fileDz << "\tthisDz = " << thisDz << endl;
        MayDay::Error("dz is not properly set");
    }

    // eta
    CH_assert(a_eta.size() >= a_nz+1);
    for (int k = 0; k <= a_nz; ++k) {
        Vector<double> dblVec(a_nx+1, 0.0);
        infile.seekg(4, ios::cur);
        infile.read((char*)&dblVec[0], sizeof(double)*(a_nx+1));
        infile.seekg(4, ios::cur);

        CH_assert(a_eta[k].size() >= a_nx+1);
        for (int i = 0; i <= a_nx; ++i) {
            a_eta[k][i] = dblVec[i];
        }
    }

    infile.close();

    return ((Real)c);
}


// -----------------------------------------------------------------------------
// Static utility
// Computes the u DJL solution over a horizontal slice.
// u = c*eta_z (but this assumes eta is already scaled by c)
// -----------------------------------------------------------------------------
void DJLBCUtil::fill_uDJL (Vector<Real>&       a_uDJL,    // CC
                           const Vector<Real>& a_etaTop,  // NC
                           const Vector<Real>& a_etaBot,  // NC
                           const Real          a_dz)
{
    // Allocate solution
    const int Nx = a_etaTop.size() - 1;
    CH_assert(a_etaBot.size() == Nx+1);
    a_uDJL.resize(Nx);

    // Compute solution
    const Real scale = 0.5 / a_dz;
    for (int i = 0; i < Nx; ++i) {
        Real detar = a_etaTop[i+1] - a_etaBot[i+1];
        Real detal = a_etaTop[i  ] - a_etaBot[i  ];
        a_uDJL[i] = (detar + detal) * scale;
    }
}


// -----------------------------------------------------------------------------
// Static utility
// Computes the w DJL solution over a horizontal slice.
// w = -c*eta_x (but this assumes eta is already scaled by c)
// -----------------------------------------------------------------------------
void DJLBCUtil::fill_wDJL (Vector<Real>&       a_wDJL,    // CC
                           const Vector<Real>& a_etaTop,  // NC
                           const Vector<Real>& a_etaBot,  // NC
                           const Real          a_dx)
{
    // Allocate solution
    const int Nx = a_etaTop.size() - 1;
    CH_assert(a_etaBot.size() == Nx+1);
    a_wDJL.resize(Nx);

    // Compute solution
    const Real scale = -0.5 / a_dx;
    for (int i = 0; i < Nx; ++i) {
        Real detat = a_etaTop[i+1] - a_etaTop[i  ];
        Real detab = a_etaBot[i+1] - a_etaBot[i  ];
        a_wDJL[i] = (detat + detab) * scale;
    }
}


// -----------------------------------------------------------------------------
// Static utility
// Computes the b DJL solution over a horizontal slice.
// b(x,z) = eta(z-eta(x)) (but this assumes eta is scaled by c)
// Returns the CC background buoyancy at this slice.
// -----------------------------------------------------------------------------
Real DJLBCUtil::fill_bDJL (Vector<Real>&       a_bDJL, // CC
                           const Vector<Real>& a_eta,  // CC
                           const Real          a_c,    // long wave speed
                           const Real          a_z)    // slice location
{
    static const Real s_z0 = 0.8; // TEMP!!!
    static const Real s_d = 0.1;  // TEMP!!!

    // Allocate solution
    const int Nx = a_eta.size();
    a_bDJL.resize(Nx);

    // Compute background solution
    Real rho = 0.5 * (1.0 - tanh((a_z - s_z0) / s_d));
    Real bgScalar = (rho - s_rho_top) * s_rho_scale;

    // Compute the DJL solution
    for (int i = 0; i < Nx; ++i) {
        Real ztilde = a_z - (a_eta[i] / a_c);
        rho = 0.5 * (1.0 - tanh((ztilde - s_z0) / s_d));
        a_bDJL[i] = (rho - s_rho_top) * s_rho_scale;
    }

    return bgScalar;
}


// -----------------------------------------------------------------------------
// Static utility
// Converts a NC horizontal slice of data in a vector to CC.
// -----------------------------------------------------------------------------
void DJLBCUtil::convertSliceNC2CC (Vector<Real>&       a_cc,
                                   const Vector<Real>& a_ncTop,
                                   const Vector<Real>& a_ncBot)
{
    // Allocate
    const int Nx = a_ncTop.size() - 1;
    CH_assert(a_ncBot.size() == Nx+1);
    a_cc.resize(Nx);

    // Convert
    const Real scale = 0.25;
    for (int i = 0; i < Nx; ++i) {
        Real topSum = a_ncTop[i  ] + a_ncTop[i+1];
        Real botSum = a_ncBot[i  ] + a_ncBot[i+1];
        a_cc[i] = scale * (topSum + botSum);
    }
}


// -----------------------------------------------------------------------------
// Static utility
// Envelope function for extrusion
// -----------------------------------------------------------------------------
Real DJLBCUtil::extrusionEnvelope (const Real a_yprime)
{
    return 0.5*(  tanh(s_envSlope*(a_yprime+0.5*s_envWidth))
                - tanh(s_envSlope*(a_yprime-0.5*s_envWidth))  );
}
