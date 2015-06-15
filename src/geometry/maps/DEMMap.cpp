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
Vector<RefCountedPtr<FArrayBox> >  DEMMap::s_depthPtr; 

// -----------------------------------------------------------------------------
// Construtor
// -----------------------------------------------------------------------------
DEMMap::DEMMap ()
: BathymetricBaseMap()
{
    if (!s_fileIsRead) {
      const ProblemContext* ctx = ProblemContext::getInstance();
      // create vector to contain the DEMs at each level of refinement
             s_depthPtr.resize(ctx->max_level+1);
	     // fills the maps.       
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
    //    static int counter=1;
    // The holder needs to be flat and nodal in the vertical.
    CH_assert(destBox == horizontalDataBox(destBox));
    CH_assert(destBoxType[SpaceDim-1] == 1);

    // we determine what level has called fill_bathymetry
    // and copy the bathymetry from the FArrayBox at the
    // appropriate level.
    const ProblemContext* ctx = ProblemContext::getInstance();
    bool depthIsSet = false;
    
    int levelMax=ctx->max_level;
    RealVect DX=a_dXi;
    for(int level=0; level<levelMax+1; ++level){
      //pout() << "leve "<< level << " m_lev " << m_lev0DXi <<" and DX " << DX << endl; 
      
      if (m_lev0DXi == DX) {// AS: I am a bit uncomfortable here, as we are 
	// comparing Real values. Rounding errors can be a PITA.
	// If there is a problem, will be caught by the MayDay call
	
	if(s_depthPtr[level]->box().contains(destBox)){
	  a_dest.copy(*(s_depthPtr[level]), 0,a_destComp);
	  depthIsSet=true;
	  break;
	}
	else{MayDay::Error("DEMMap:: fill_bathymetry - Need bigger box for cached depth"); 
	}
      } 
	DX*=(ctx->refRatios[level]); // calculates coarsened dx at next level.  
	
      }
    if(depthIsSet==false) {MayDay::Error("DEMMap::fill_bathymetry can't determine level");}
    //const Box& dummyBox= s_depthPtr[0]->box();
    
}
// -----------------------------------------------------------------------------
// Reads the bathymetry from the Digital Elevation Model and interpolates onto level grids depending on the dimensionality of the problem
// -----------------------------------------------------------------------------
void DEMMap::createBaseMap ()
{    

  //  Vector<Real> depthVect,xVect,yVect;

  //Read_DEM_File(depthVect,xVect,yVect);

  const Vector<Real>&     xVect=this->readHDF5Vector("/X");
  const Vector<Real>& depthVect=this->readHDF5Vector("/Depth");
  if(SpaceDim==2)
    {this->Create_Level_DEM_2D(depthVect,xVect);}
  else
    {{const Vector<Real>& yVect=readHDF5Vector("/Y");
	this->Create_Level_DEM_3D(depthVect,xVect,yVect);}
      return;
    }    
}
//---------------------------------------------------------------------------
// fills the 
void DEMMap::Create_Level_DEM_3D(const Vector<Real>& depthVect,
                                 const Vector<Real>& xVect,
				 const Vector<Real>& yVect)
{   
  // get problem context
  const ProblemContext* ctx = ProblemContext::getInstance();
  int Nx=xVect.size();
  int Ny=yVect.size();
  // create a box for the data just read.
  Box domBox(IntVect::Zero, IntVect(D_DECL(Nx,Ny,0)), IntVect::Unit); //node centered
  domBox.shiftHalf(SpaceDim-1, 1); // cell centered in vertical
  domBox.enclosedCells(); // pixie dust
  CH_assert(!domBox.isEmpty());
  CH_assert(domBox.type() == IntVect::Zero);
  
  
  // package the input vectors into FABs
  
  FArrayBox depthFAB(domBox, 1);
  FArrayBox dfdxFAB(domBox,1);
  FArrayBox dfdyFAB(domBox,1);

  // cell centered in x,y and flat in z. 
    
    

  
  for (BoxIterator bit(domBox); bit.ok(); ++bit) {
    const IntVect& nc = bit(); //AS: Why the ampersand? 
    int idx = nc[0] + nc[1]*Nx; // arrays follow Fortran ordering 
    depthFAB(nc) = depthVect[idx];
  }
  {// X derivative
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
  //  pout() << "calculated dfdx " << endl;
  {  // now the y derivative
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
  //  pout() << "calculated dfdy " << endl;



  // Stuff we need on the level=0 grid. We iterate from here on
int nx_offset = ctx->nx_offset[0];
//Real xmin = (Real(nx_offset))*m_lev0DXi[0];
int ny_offset = ctx->nx_offset[1];
//Real ymin = (Real(ny_offset))*m_lev0DXi[1];
int nx = ctx->nx[0];
int ny = ctx->nx[1];
Real dX=m_lev0DXi[0];
Real dY=m_lev0DXi[1];
Box interpBox = ctx->domain.domainBox();

 interpBox.grow(IntVect(D_DECL(5,5,0)));// we enlarge the box by 5, which should accomodate all the necessary 
 // ghost points...
interpBox.surroundingNodes();
interpBox = flattenBox(interpBox, SpaceDim-1);
 FArrayBox xInterp(interpBox,1);
 FArrayBox yInterp(interpBox,1);
 // pout()<< "xmin and ymin " << xmin << " " << ymin<<endl;
for (int level=0; level<ctx->max_level+1; ++level){
   s_depthPtr[level] = RefCountedPtr<FArrayBox> (new FArrayBox);
   s_depthPtr[level]->define(interpBox,1); 
  BoxIterator bit(interpBox);
  for (bit.reset(); bit.ok(); ++bit) {
    const IntVect& nc = bit();
    xInterp(nc) =  nc[0]*dX;
    yInterp(nc) =  nc[1]*dY;
  
  }
  // Interpolate onto our interpBox and store the results in the appropriate container
  HermiteInterp2D(*(s_depthPtr[level]),
     		  xInterp,
     		  yInterp,
     		  interpBox,
     		  0, // xdir
     		  1, // ydir
     		  xVect,
     		  yVect,
     		  depthFAB,
     		  dfdxFAB,
     		  dfdyFAB);
  //    if(level==0) {writeFABname(s_depthPtr[level],"depthPtr.hdf5");
  //  writeFABname(&xInterp,"xInterp.hdf5");
  //  writeFABname(&yInterp,"yInterp.hdf5");

 // move on to the next grid
  nx *= ctx->refRatios[level][0];
  ny *= ctx->refRatios[level][1];
  nx_offset *= ctx->refRatios[level][0]; 
  ny_offset *= ctx->refRatios[level][1]; 
  dX /= Real(ctx->refRatios[level][0]);
  dY /= Real(ctx->refRatios[level][1]);
  interpBox.refine(ctx->refRatios[level]);
  xInterp.resize(interpBox,1);
  yInterp.resize(interpBox,1);

  //  MayDay::Error("Ciao bello");}
 }
 

return;
}

   

 

   

void DEMMap::Create_Level_DEM_2D(const Vector<Real>& depthVect, 
				 const Vector<Real>&xVect)

{  const ProblemContext* ctx = ProblemContext::getInstance();
  // prepare depthVect for interpolation
  CubicSpline cs;
  cs.solve(depthVect,xVect);
  // Stuff we need on the level=0 grid
  int nx_offset = ctx->nx_offset[0];
  Real xmin = (Real(nx_offset))*m_lev0DXi[0];
  int nx = ctx->nx[0];
  Real dX=m_lev0DXi[0];
  Box interpBox = ctx->domain.domainBox();
  interpBox.grow(IntVect(D_DECL(5,0,0)));// we enlarge the box by 5, which should accomodate all the necessary 
  interpBox.surroundingNodes();
  interpBox = flattenBox(interpBox, SpaceDim-1);
  Vector<Real> xInterp(interpBox.size(0));
  Vector<Real> depthInterp(interpBox.size(0));

  for (int level=0; level<ctx->max_level+1; ++level) 
    {// we create the coordinates where we will interpolate


      for (int i=interpBox.smallEnd(0); i<interpBox.bigEnd(0); ++i) {

	xInterp[i-interpBox.smallEnd(0)] =  dX* Real(i) ;
      }       
      
      cs.interp(depthInterp,xInterp);
      
      // now we chomboize the result
      s_depthPtr[level] = RefCountedPtr<FArrayBox> (new FArrayBox);
      s_depthPtr[level]->define(interpBox,1); 
      //pout() << "interpBox level "<< level << " = " << interpBox << endl;
      BoxIterator bit(interpBox);
      for (bit.begin(); bit.ok(); ++bit){
	
	(*(s_depthPtr[level]))(bit())=depthInterp[(bit()[0]-interpBox.smallEnd(0))]; // and now we fill it with the (hopefully) right data
      }
      // now we update the stuff we need for the next level
      nx *= ctx->refRatios[level][0]; 
      nx_offset *= ctx->refRatios[level][0]; 
      dX /= Real(ctx->refRatios[level][0]);
      interpBox.refine(ctx->refRatios[level]);
      xInterp.resize(interpBox.size(0));
      depthInterp.resize(interpBox.size(0));
    }
  return; }

Vector<Real> DEMMap::readHDF5Vector(const char* invector)

{  const ProblemContext* ctx = ProblemContext::getInstance();
  int N; 
  std::string infile;
  Vector<Real> V;
  infile=ctx->DemFile;
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
    // Determine size of data
    {
      hid_t xID = H5Dopen(fileID, invector);
      if (xID < 0) {
	ostringstream msg;
	msg << "H5Dopen failed to open "<< invector << " dataset. Return value = " << xID;
	MayDay::Error(msg.str().c_str());
      }
      hid_t xDSPC=H5Dget_space(xID);
      N=int( H5Sget_simple_extent_npoints( xDSPC ));				 
      H5Sclose(xDSPC);
      H5Dclose(xID);
    }
    // Read  data
    V.resize(N);
    V.assign(quietNAN);
    {
        hid_t xID = H5Dopen(fileID, invector);

        herr_t status = H5Dread(xID,                // hid_t dataset_id    IN: Identifier of the dataset read from.
                                H5T_NATIVE_DOUBLE,  // hid_t mem_type_id   IN: Identifier of the memory datatype.
                                H5S_ALL,            // hid_t mem_space_id  IN: Identifier of the memory dataspace.
                                H5S_ALL,            // hid_t file_space_id IN: Identifier of the dataset's dataspace in the file.
                                H5P_DEFAULT,        // hid_t xfer_plist_id     IN: Identifier of a transfer property list for this I/O operation.
                                &V[0]);         // void * buf  OUT: Buffer to receive data read from file.
        H5Dclose(xID);
	if(status<0) {MayDay::Error("Cannot read data from depth file");} 
    }

    return V;
}


