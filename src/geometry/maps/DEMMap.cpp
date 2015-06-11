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
// DEMMap is not for public use yet.

#include "DEMMap.H"
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


bool DEMMap::s_fileIsRead = false;
RefCountedPtr<BoxLayoutData<FArrayBox> > DEMMap::s_depthPtr;

// -----------------------------------------------------------------------------
// Construtor
// -----------------------------------------------------------------------------
DEMMap::DEMMap ()
: BathymetricBaseMap()
{
    if (!s_fileIsRead) {
        this->createBaseMap();


        s_fileIsRead = true;
    }
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
DEMMap::~DEMMap ()
{;}


// -----------------------------------------------------------------------------
// Must return the name of the coordinate mapping
// -----------------------------------------------------------------------------
const char* DEMMap::getCoorMapName () const
{
    return "DEMMap";
}


// -----------------------------------------------------------------------------
// Must return whether or not this metric is diagonal
// -----------------------------------------------------------------------------
bool DEMMap::isDiagonal () const
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
void DEMMap::fill_bathymetry (FArrayBox&       a_dest,
                                  const int        a_destComp,
                                  const FArrayBox& a_cartPos,
                                  const RealVect&  a_dXi) const
{
    const Box& destBox = a_dest.box();
    const IntVect destBoxType = destBox.type();
    
    DataIterator dit=(*s_depthPtr).dataIterator();  // this does not work. I need to access the 
                                     // boxLayout which was used to construct this 
    // The holder needs to be flat and nodal in the vertical.
    CH_assert(destBox == horizontalDataBox(destBox));
    CH_assert(destBoxType[SpaceDim-1] == 1);

    // TODO
    // The problem is how to copy the depths. destBox is a box 
    // which is generally contained in the domain. 
    if (m_lev0DXi == a_dXi) {
      for (dit.begin();dit.ok();++dit) {
	// square peg in round hole. Why is a_dest defined as a FArrayBox, when
        // we store the depths in a NodeFArrayBox ?
	FArrayBox& depths=(*s_depthPtr)[dit]; // WHy the ampersand
	
        a_dest.copy(depths);
      }
    } else {
        MayDay::Error("Bad a_dXi");
    }
}


// -----------------------------------------------------------------------------
// Reads the bathymetry from the Digital Elevation Model and interpolates onto level 0 grids.
// This also trims the region and smooths the edges.
// -----------------------------------------------------------------------------
void DEMMap::createBaseMap ()
{
  
    std::string infile;
    int Nx;
    int Ny;
    const ProblemContext* ctx = ProblemContext::getInstance();

    infile = ctx->DemFile;
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
    // Determine size of X data
    {
      hid_t xID = H5Dopen(fileID, "/X");
      if (xID < 0) {
	ostringstream msg;
	msg << "H5Dopen failed to open X dataset. Return value = " << xID;
	MayDay::Error(msg.str().c_str());
      }
    
      Nx=int( H5Sget_simple_extent_npoints(  xID ));				 
      H5Dclose(xID);
    }
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
    // Determine size of Y data
    {hid_t xID = H5Dopen(fileID, "/Y");
      if (xID < 0) {
	//	ostringstream msg;
	//	msg << "H5Dopen failed to open Y dataset. Return value = " << xID;
	//	MayDay::Error(msg.str().c_str());
        // The DEM is one dimensional
	if (SpaceDim == 3) {
	  MayDay::Error("DEMMap is two-dimensional but problem is solved in 3D");
	}
  
	Ny=1;
      }
      else
	{
	  Ny=int( H5Sget_simple_extent_npoints( xID ));				 
	  H5Dclose(xID);
	}
    }
    // Read Y data
    
      Vector<double> yVect(Ny, quietNAN);
      if (Ny>1) {    
            
        hid_t yID = H5Dopen(fileID, "/Y");
        
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
      Box domBox(IntVect::Zero, IntVect(D_DECL(Nx-1,Ny-1,0)), IntVect::Unit); //node centered
        domBox.shiftHalf(SpaceDim-1, 1); // cell centered in vertical
        domBox.enclosedCells(); // pixie dust
        CH_assert(!domBox.isEmpty());
        CH_assert(domBox.type() == IntVect::Zero);
        domain.define(domBox);

        Vector<Box> vbox(1, domBox); // # of boxes and initialization list
        Vector<int> vproc(1, 0); 
        grids.define(vbox, vproc, domain);

        dit = grids.dataIterator();
    }

    if(Ny==1)
      { CubicSpline cs;
	cs.solve(depthVect,xVect);
	Box interpBox = ctx->domain.domainBox();
	interpBox.surroundingNodes();
	interpBox = flattenBox(interpBox, SpaceDim-1);
	// Box interpBox(IntVect::Zero, IntVect(D_DECL(256,256,0)), IntVect::Unit);
	pout() << "interpBox = " << interpBox << endl;
	int nx = ctx->nx[0];
	int nx_offset = ctx->nx_offset[0];
	Real dX = ctx->domainLength[0]/Real(nx);
	Real xmin = (Real(nx_offset))*dX;
	Real xmax = xmin+ctx->domainLength[0];
	Vector<Real> xInterp(nx+1);

	Vector<Real> depthInterp(nx+1);

	  
	for (int i=0; i<nx+1; ++i) {
            
	    
	  xInterp[i] = xmin + dX* Real(i) ;

	}	    
	cs.interp(depthInterp,xInterp);

	// now we chomboize it
	s_depthPtr = RefCountedPtr<BoxLayoutData<FArrayBox> >(new BoxLayoutData<FArrayBox>);
	Vector<Box> vB(numProc(),domain.domainBox());
	Vector<int> vP(numProc());
	for (int p=0; p<numProc(); ++p) {
	  vP[p]=p;}
	BoxLayout bL(vB,vP);
	bL.close();// each process has the same copy
	//	FlatNodeDataFactory DF=FlatNodeDataFactory(BASISV(SpaceDim-1));
	//	s_depthPtr->define(bL,1,(DataFactory<NodeFArrayBox>)DF); // Every processor is good to go
	s_depthPtr->define(bL,1); // Still unallocated the space
	for (DataIterator dit=bL.dataIterator(); dit.ok(); ++dit) {
	  (*s_depthPtr)[dit()].setVal(0.0); 
	  // If I get this right, now we have allocated the space (hopefully, of the right type)
	  for (BoxIterator bit(domain.domainBox()); bit.ok(); ++bit){
	    (*s_depthPtr)[dit()](bit())=depthInterp[(bit()[0]-nx_offset)]; // and now we fill it with the (hopefully) right data
	  }

	}
	
     
	    return;}
    

    // Package the data into a LevelData<FArrayBox>.
    LevelData<FArrayBox> depthData(grids, 1, IntVect::Zero, FlatNodeDataFactory(BASISV(SpaceDim-1))); 
    // cell centered in x,y and flat in z. 
    for (dit.reset(); dit.ok(); ++dit) { //only rank 0 does the work...
        FArrayBox& depthFAB = depthData[dit];
        BoxIterator bit(depthFAB.box());
        int idx;

        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& nc = bit();
            idx = nc[1] + nc[0]*Ny;
            depthFAB(nc) = depthVect[idx];
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
    int nx = ctx->nx[0];
    int nx_offset = ctx->nx_offset[0];
    
    int ny = ctx->nx[1];
    int ny_offset = ctx->nx_offset[1];
    
    Real dX = ctx->domainLength[0]/Real(nx);
    Real dY = ctx->domainLength[1]/Real(ny);
    Real xmin = (Real(nx_offset))*dX;
    Real xmax = xmin+ctx->domainLength[0];
    Real ymin = (Real(ny_offset))*dY;
    Real ymax = ymin+ctx->domainLength[1];
    
    FArrayBox xInterp(interpBox, 1);
    FArrayBox yInterp(interpBox, 1);
    FArrayBox depthInterp;

    
        BoxIterator bit(interpBox);
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& nc = bit();

            xInterp(nc) = xmin + (xmax - xmin) * Real(nc[0]) / Real(interpBox.size(0));
            yInterp(nc) = ymin + (ymax - ymin) * Real(nc[1]) / Real(interpBox.size(1));
 // Interpolate onto our interpBox.
    
    }
    depthInterp.define(interpBox, 1);
    for (dit.reset(); dit.ok(); ++dit) { //only rank 0 does the work...

	
    HermiteInterp2D(depthInterp,
                    xInterp,
                    yInterp,
                    interpBox,
                    0, // xdir
                    1, // ydir
                    xVect,
                    yVect,
                    depthData[dit],
                    dfdx[dit],
                    dfdy[dit]);
    }
     writeFABname(&depthInterp, "depthInterp.hdf5");
    
	// s_depthPtr = RefCountedPtr<BoxLayoutData<NodeFArrayBox> >(new BoxLayoutData<NodeFArrayBox>);
	// Vector<Box> vB(nProc(),domainBox);
	// Vector<int> vP(nProc());
	// for (int p=0; p<nProc(); ++p) {
	//   vP[p]=p;}
	// BoxLayout bL(vB,vP);
	// s_depthPtr->define(bL,1,FlatNodeDataFactory(BASISV(SpaceDim-1))); // Every processor is good to go
	// for (LayoutIterator dit=bL.layoutIterator(); dit.ok(); ++dit) {
	//   for (BoxIterator bit(domainBox); bit.ok(); ++bit){
	//     (*s_depthPtr)[dit](bit())=depthInterp(bit()[0]-nx_offset);
	//   }
	// }      
    
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
