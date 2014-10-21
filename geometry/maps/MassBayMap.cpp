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

// NOTE: This is a work in progress.
// MassBayMap is not for public use yet.

#include "MassBayMap.H"
#include "LevelData.H"
#include "AMRIO.H"
#include "CH_HDF5.H"
#include "Debug.H"
#include "Subspace.H"
#include "Constants.H"
// #include "BathymetricBaseMap.H"
// #include "BathymetricBaseMapF_F.H"
#include "ProblemContext.H"
// #include "NodeInterpF_F.H"
// #include "ConvertFABF_F.H"
#include "CubicSpline.H"
#include "HermiteInterp.H"
#include "EllipticBCUtils.H"
#include "BiCGStabSolver.H"
#include "MappedAMRPoissonOp.H"
#include "ConvertFABF_F.H"
#include "ExtrapolationUtils.H"
#include "LepticMeshRefine.H"
#include "NodeAMRIO.H"
#include "FCDataFactory.H"


bool MassBayMap::s_fileIsRead = false;
RefCountedPtr<LevelData<FArrayBox> > MassBayMap::s_depthPtr;

// -----------------------------------------------------------------------------
// Construtor
// -----------------------------------------------------------------------------
MassBayMap::MassBayMap ()
: BathymetricBaseMap()
{
    if (!s_fileIsRead) {
        this->createBaseMap();

        // CH_assert(s_depthPtr.isNull());
        // s_depthPtr = RefCountedPtr<LevelData<NodeFArrayBox> >(new LevelData<NodeFArrayBox>);

        // std::string filename("../geometry/maps/MassBayDepths.hdf5");
        // Vector<DisjointBoxLayout> vectGrids;
        // Vector<LevelData<NodeFArrayBox>*> vectData(1, &*s_depthPtr);
        // Vector<std::string> vectNames;
        // Box domBox;
        // Real dx, dt, time;
        // Vector<int> refRatio;
        // int numLevels;
        // bool setGhost = true;

        // int err = ReadAMRHierarchyHDF5(filename,
        //                                vectGrids,
        //                                vectData,
        //                                vectNames,
        //                                domBox,
        //                                dx,
        //                                dt,
        //                                time,
        //                                refRatio,
        //                                numLevels,
        //                                setGhost);
        // if (err) {
        //     ostringstream msg;
        //     msg << "Could not read " << filename << ". Error code = " << err;
        //     MayDay::Error(msg.str().c_str());
        // }

        s_fileIsRead = true;
    }
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
MassBayMap::~MassBayMap ()
{;}


// -----------------------------------------------------------------------------
// Must return the name of the coordinate mapping
// -----------------------------------------------------------------------------
const char* MassBayMap::getCoorMapName () const
{
    return "MassBayMap";
}


// -----------------------------------------------------------------------------
// Must return whether or not this metric is diagonal
// -----------------------------------------------------------------------------
bool MassBayMap::isDiagonal () const
{
    return false;
}


// -----------------------------------------------------------------------------
// Fills a NodeFAB with the bathymetric data. a_dest must be flat in the
// vertical. Upon return, each point in the horizontal (Xi,Eta) of a_dest
// will contain the (positive) local depth.
// NOTE: This vertical distance is measured in a straight line perpendicular
// the the surface. We are measuring this distance along the Cartesian
// vertical coordinate line, not the mapped vertical coordinate line.
// -----------------------------------------------------------------------------
void MassBayMap::fill_bathymetry (FArrayBox&       a_dest,
                                  const int        a_destComp,
                                  const FArrayBox& a_cartPos,
                                  const RealVect&  a_dXi) const
{
    const Box& destBox = a_dest.box();
    const IntVect destBoxType = destBox.type();

    // The holder needs to be flat and nodal in the vertical.
    CH_assert(destBox == horizontalDataBox(destBox));
    CH_assert(destBoxType[SpaceDim-1] == 1);

    // TODO

    if (m_lev0DXi == a_dXi) {
        a_dest.setVal(1.0, a_destComp);
    } else {
        MayDay::Error("Bad a_dXi");
    }
}


// -----------------------------------------------------------------------------
// Reads the bathymetry from file and interpolates onto level 0 grids.
// This also trims the region and smooths the edges.
// -----------------------------------------------------------------------------
void MassBayMap::createBaseMap ()
{
    // This map only works in 3D
    if (SpaceDim != 3) {
        MayDay::Error("MassBayMap only works in 3D");
    }

    const string infile("../geometry/maps/MassBayTopo.hd5");
    const int Nx = 200;           // How many x nodes are there in the topo file?
    const int Ny = 200;           // How many y nodes are there in the topo file?

    const int cropXMin = 40;      // This is the region that we will use for
    const int cropXMax = 180;     // interpolation. It is the data with NaNs
    const int cropYMin = 70;      // cropped out.
    const int cropYMax = 150;

    const Real xmin = -15929.6;   // We will interpolate the topo using the data
    const Real xmax =  50000.0;   // in this region.
    const Real ymin = -28040.2;
    const Real ymax =  30000.0;
    const Real clipDepth = -70.0; // We will clip the ocean at this depth.

    const IntVect hmask = IntVect::Unit - BASISV(SpaceDim-1);
    const ProblemContext* ctx = ProblemContext::getInstance();

    // Open the input HDF5 file.
    pout() << "Opening " << infile << "..." << std::flush;
    // HDF5Handle ifHandle(infile, HDF5Handle::OPEN_RDONLY, "/");
    hid_t fileID = H5Fopen(infile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fileID < 0) {
        ostringstream msg;
        msg << "H5Fopen failed to open " << infile << ". Return value = " << fileID;
        MayDay::Error(msg.str().c_str());
    }
    pout() << "done." << endl;


    // Read X data
    Vector<double> xVect(Nx, quietNAN);
    {
        hid_t xID = H5Dopen(fileID, "/X");
        if (xID < 0) {
            ostringstream msg;
            msg << "H5Dopen failed to open X dataset. Return value = " << xID;
            MayDay::Error(msg.str().c_str());
        }

        herr_t status = H5Dread(xID,                // hid_t dataset_id    IN: Identifier of the dataset read from.
                                H5T_NATIVE_DOUBLE,  // hid_t mem_type_id   IN: Identifier of the memory datatype.
                                H5S_ALL,            // hid_t mem_space_id  IN: Identifier of the memory dataspace.
                                H5S_ALL,            // hid_t file_space_id IN: Identifier of the dataset's dataspace in the file.
                                H5P_DEFAULT,        // hid_t xfer_plist_id     IN: Identifier of a transfer property list for this I/O operation.
                                &xVect[0]);         // void * buf  OUT: Buffer to receive data read from file.
        H5Dclose(xID);
    }

    // Read Y data
    Vector<double> yVect(Ny, quietNAN);
    {
        hid_t yID = H5Dopen(fileID, "/Y");
        if (yID < 0) {
            ostringstream msg;
            msg << "H5Dopen failed to open Y dataset. Return value = " << yID;
            MayDay::Error(msg.str().c_str());
        }

        herr_t status = H5Dread(yID,                // hid_t dataset_id    IN: Identifier of the dataset read from.
                                H5T_NATIVE_DOUBLE,  // hid_t mem_type_id   IN: Identifier of the memory datatype.
                                H5S_ALL,            // hid_t mem_space_id  IN: Identifier of the memory dataspace.
                                H5S_ALL,            // hid_t file_space_id IN: Identifier of the dataset's dataspace in the file.
                                H5P_DEFAULT,        // hid_t xfer_plist_id     IN: Identifier of a transfer property list for this I/O operation.
                                &yVect[0]);         // void * buf  OUT: Buffer to receive data read from file.
        H5Dclose(yID);
    }

    // Read Depth data
    Vector<double> depthVect(Nx*Ny, quietNAN);
    {
        hid_t depthID = H5Dopen(fileID, "/Depth");
        if (depthID < 0) {
            ostringstream msg;
            msg << "H5Dopen failed to open Depth dataset. Return value = " << depthID;
            MayDay::Error(msg.str().c_str());
        }

        herr_t status = H5Dread(depthID,            // hid_t dataset_id    IN: Identifier of the dataset read from.
                                H5T_NATIVE_DOUBLE,  // hid_t mem_type_id   IN: Identifier of the memory datatype.
                                H5S_ALL,            // hid_t mem_space_id  IN: Identifier of the memory dataspace.
                                H5S_ALL,            // hid_t file_space_id IN: Identifier of the dataset's dataspace in the file.
                                H5P_DEFAULT,        // hid_t xfer_plist_id     IN: Identifier of a transfer property list for this I/O operation.
                                &depthVect[0]);     // void * buf  OUT: Buffer to receive data read from file.
        H5Dclose(depthID);
    }

    // We are done reading from file.
    if (!(fileID < 0)) {
        H5Fclose(fileID);
    }

    // Create the level 0 grids
    ProblemDomain domain;
    DisjointBoxLayout grids;
    DataIterator dit;
    {
        // Create the domain
        Box domBox(IntVect::Zero, IntVect(D_DECL(Nx-1,Ny-1,0)), IntVect::Unit);
        domBox.shiftHalf(SpaceDim-1, 1);
        domBox.enclosedCells();
        CH_assert(!domBox.isEmpty());
        CH_assert(domBox.type() == IntVect::Zero);
        domain.define(domBox);

        Vector<Box> vbox(1, domBox);
        Vector<int> vproc(1, 0);
        grids.define(vbox, vproc, domain);

        dit = grids.dataIterator();
    }

    // Package the data into a LevelData<FArrayBox>.
    LevelData<FArrayBox> depthData(grids, 1, IntVect::Zero, FlatNodeDataFactory(BASISV(SpaceDim-1)));
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& depthFAB = depthData[dit];
        BoxIterator bit(depthFAB.box());
        int idx;

        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& nc = bit();
            idx = nc[1] + nc[0]*Ny;
            depthFAB(nc) = depthVect[idx];
        }
    }

    // Create a smoothed mask. This crops the NaNs out.
    LevelData<FArrayBox> mask;
    {
        // Create the LDs
        LevelData<FArrayBox> diffusedMask(grids, 1, hmask);
        LevelData<FArrayBox> sharpMask(grids, 1);

        // Create the mask
        Box cropBox(IntVect(D_DECL(cropXMin,cropYMin,0)),
                    IntVect(D_DECL(cropXMax,cropYMax,0)),
                    IntVect::Zero);

        for (dit.reset(); dit.ok(); ++dit) {
            sharpMask[dit].setVal(1.0);
            sharpMask[dit].setVal(0.0, cropBox, 0);
        }

        // Set initial guess for diffusive solve
        sharpMask.copyTo(diffusedMask);

    //     // Set up a diffusion operator...
    //     BCMethodHolder bc;
    //     RefCountedPtr<BCGhostClass> BCPtr(
    //         new EllipticConstDiriBCGhostClass(RealVect::Unit,
    //                                           RealVect::Unit,
    //                                           hmask,
    //                                           hmask)
    //     ); // extrap order
    //     bc.addBCMethod(BCPtr);

    //     RefCountedPtr<LevelData<FluxBox> > JgupPtr(new LevelData<FluxBox>);
    //     JgupPtr->define(grids1, SpaceDim, hmask);
    //     (*JgupPtr)[dit1].setVal(0.0);
    //     (*JgupPtr)[dit1][0].setVal(1.0, 0);
    //     (*JgupPtr)[dit1][1].setVal(1.0, 1);

    //     MappedAMRPoissonOp op;
    //     op.linearOpDefine(bc,
    //                       RealVect::Unit, // dx
    //                       RealVect::Zero, // dxCrse
    //                       1.0,            // alpha
    //                       -3.0,           // beta (diffusion strength)
    //                       true,           // is diagonal?
    //                       true,           // horizontal op?
    //                       JgupPtr);

    //     BiCGStabSolver<LevelData<FArrayBox> > solver;
    //     solver.define(&op, false);
    //     solver.m_verbosity = 4;
    //     solver.solve(diffusedMask, sharpMask);
    //     {
    //         Box valid = grids1[dit1];

    //         for (int dir = 0; dir < SpaceDim-1; ++dir) {
    //             ExtrapolateFaceNoEV(diffusedMask[dit1], diffusedMask[dit1], valid, dir, Side::Lo, 2);
    //             ExtrapolateFaceNoEV(diffusedMask[dit1], diffusedMask[dit1], valid, dir, Side::Hi, 2);
    //             valid.grow(dir, 1);
    //         }
    //     }

    //     // writeLevelname(&diffusedMask, "diffusedMask.hdf5");
    //     // writeLevelname(&sharpMask, "sharpMask.hdf5");

        mask.define(grids, 1, depthData.ghostVect(), FlatNodeDataFactory(BASISV(SpaceDim-1)));
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& maskFAB = mask[dit];
            maskFAB.shiftHalf(SpaceDim-1, 1);

            pout() << "maskFAB.box() = " << maskFAB.box() << endl;

            FORT_CONVERTFAB(
                CHF_FRA1(maskFAB,0),
                CHF_BOX(maskFAB.box()),
                CHF_CONST_INTVECT(hmask),
                CHF_CONST_FRA1(diffusedMask[dit],0),
                CHF_CONST_INTVECT(IntVect::Zero));

            maskFAB.shiftHalf(SpaceDim-1, -1);
        }
    }

    // Apply the mask to clip the data.
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& depthFAB = depthData[dit];
        const FArrayBox& maskFAB = mask[dit];

        BoxIterator bit(depthFAB.box());
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& nc = bit();
            Real m = maskFAB(nc);
            depthFAB(nc) = Max(depthFAB(nc), clipDepth);
            depthFAB(nc) = clipDepth*m + depthFAB(nc)*(1.0-m);
        }
    }


    // We need the first derivatives at each node. The values come from cubic
    // spline interpolations of the data along all x and y coordinate lines.

    // First the x derivatives...
    LevelData<FArrayBox> dfdx(grids, 1, depthData.ghostVect(), FlatNodeDataFactory(BASISV(SpaceDim-1)));
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& dfdxFAB = dfdx[dit];
        const FArrayBox& depthFAB = depthData[dit];

        Vector<Real> depthSlice(Nx);
        Vector<Real> dfdx(Nx);
        IntVect nc = IntVect::Zero;

        for (int j = 0; j < Ny; ++j) {
            nc[1] = j;

            // Copy slice to local vector
            for (int i = 0; i < Nx; ++i) {
                nc[0] = i;
                depthSlice[i] = depthFAB(nc);
            }

            // // Compute first derivatives
            // CubicSpline cs;
            // cs.solve(depthSlice, xData);
            // cs.interpFirstDeriv(dfdx, xData);

            // Just use a finite difference
            const IntVect& e = BASISV(0);
            nc[0] = 0;
            dfdx[0] = (depthFAB(nc+e) - depthFAB(nc)) / (xVect[1] - xVect[0]);
            for (int i = 1; i < Nx-1; ++i) {
                nc[0] = i;
                dfdx[i] = (depthFAB(nc+e) - depthFAB(nc-e)) / (xVect[i+1] - xVect[i-1]);
            }
            nc[0] = Nx-1;
            dfdx[Nx-1] = (depthFAB(nc) - depthFAB(nc-e)) / (xVect[Nx-1] - xVect[Nx-2]);

            // Copy the data to a permanent holder
            for (int i = 0; i < Nx; ++i) {
                nc[0] = i;
                dfdxFAB(nc) = dfdx[i];
            }
        }
    }

    // Then the y derivatives...
    LevelData<FArrayBox> dfdy(grids, 1, depthData.ghostVect(), FlatNodeDataFactory(BASISV(SpaceDim-1)));
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& dfdyFAB = dfdy[dit];
        const FArrayBox& depthFAB = depthData[dit];

        Vector<Real> depthSlice(Ny);
        Vector<Real> dfdy(Ny);
        IntVect nc = IntVect::Zero;

        for (int i = 0; i < Nx; ++i) {
            nc[0] = i;

            // Copy slice to local vector
            for (int j = 0; j < Ny; ++j) {
                nc[1] = j;
                depthSlice[j] = depthFAB(nc);
            }

            // // Compute first derivatives
            // CubicSpline cs;
            // cs.solve(depthSlice, yData);
            // cs.interpFirstDeriv(dfdy, yData);

            // Just use a finite difference
            const IntVect& e = BASISV(1);
            nc[1] = 0;
            dfdy[0] = (depthFAB(nc+e) - depthFAB(nc)) / (yVect[1] - yVect[0]);
            for (int j = 1; j < Ny-1; ++j) {
                nc[1] = j;
                dfdy[j] = (depthFAB(nc+e) - depthFAB(nc-e)) / (yVect[j+1] - yVect[j-1]);
            }
            nc[1] = Ny-1;
            dfdy[Ny-1] = (depthFAB(nc) - depthFAB(nc-e)) / (yVect[Ny-1] - yVect[Ny-2]);

            // Copy the data to a permanent holder
            for (int j = 0; j < Ny; ++j) {
                nc[1] = j;
                dfdyFAB(nc) = dfdy[j];
            }
        }
    }

    // Where do we *want* data?
    Box interpBox = ctx->domain.domainBox();
    interpBox.surroundingNodes();
    interpBox = flattenBox(interpBox, SpaceDim-1);
    // Box interpBox(IntVect::Zero, IntVect(D_DECL(256,256,0)), IntVect::Unit);
    pout() << "interpBox = " << interpBox << endl;

    FArrayBox xInterp(interpBox, 1);
    FArrayBox yInterp(interpBox, 1);
    {
        BoxIterator bit(interpBox);
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& nc = bit();

            xInterp(nc) = xmin + (xmax - xmin) * Real(nc[0]) / Real(interpBox.size(0));
            yInterp(nc) = ymin + (ymax - ymin) * Real(nc[1]) / Real(interpBox.size(1));
        }
    }

    // // Interpolate onto our interpBox.
    // FArrayBox depthInterp(interpBox, 1);
    // HermiteInterp2D(depthInterp,
    //                 xInterp,
    //                 yInterp,
    //                 interpBox,
    //                 0, // xdir
    //                 1, // ydir
    //                 xData,
    //                 yData,
    //                 depthFAB,
    //                 dfdxFAB,
    //                 dfdyFAB);

    // writeFABname(&depthInterp, "depthInterp.hdf5");

    // // // Extend the data onto a larger FAB.
    // // {
    // //     Box extendedBox = depthInterp.box();
    // //     extendedBox.grow(128*hmask);

    // //     FArrayBox tmpFAB(depthInterp.box(), 1);
    // //     tmpFAB.copy(depthInterp);
    // //     tmpFAB.shift(1, 50); // Changes the location of the ridge center.

    // //     depthInterp.define(extendedBox, 1);
    // //     depthInterp.setVal(clipDepth);
    // //     depthInterp.copy(tmpFAB);
    // //     depthInterp.shift(-extendedBox.smallEnd());
    // // }

}
