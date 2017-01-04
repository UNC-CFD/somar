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

// This file contains utilities for printing to the terminal and files.
// Some of these functions have been taken directly from Chombo. Please see
// Chombo's copyright notice.

#ifndef __PRINTING_H__INCLUDED__
#define __PRINTING_H__INCLUDED__

// Standard includes
#include <iostream>
#include <string>
#include <sstream>

// Chombo includes
#include "LevelData.H"
#include "RealVect.H"
#include "LevelGeometry.H"

// Project includes
#include "GNUC_Extensions.H"
class FArrayBox;
class FluxBox;


// -----------------------------------------------------------------------------
// Create a banner
// -----------------------------------------------------------------------------
std::string banner(const std::string& a_str);


// -----------------------------------------------------------------------------
// Print a variable to pout.* with a nice label
// -----------------------------------------------------------------------------
#define POUT(x)         pout() << #x << " = " << x << std::endl;
#define TOUT(n,x)       tout((n)) << #x << " = " << x << std::endl;


// -----------------------------------------------------------------------------
// Print to terminal from one or all procs
// -----------------------------------------------------------------------------
struct null_ostream: std::ostream {
    null_ostream(): std::ios(0), std::ostream(0) {}
};

inline std::ostream& tout(int a_procid = -1) {
#   ifdef CH_MPI
        static null_ostream nullstream;
        if(a_procid == -1) {
            return std::cout << "Proc[" << procID() << "]: ";
        } else if(a_procid <= -2) {
            return std::cout;
        } else if(procID() != a_procid) {
            return nullstream;
        }
#   endif
    return std::cout;
}


// -----------------------------------------------------------------------------
// Color stream manipulators
// NOTE: Empty in release mode!
// -----------------------------------------------------------------------------
namespace color {
#   define USE_COLOR
#   if ((!defined NDEBUG) && (defined USE_COLOR))
        // Return to default color
        static const char __nowarn_unused *none     = "\033[0m";

        // Low intensity colors
        static const char __nowarn_unused *black    = "\033[0;30m";
        static const char __nowarn_unused *blue     = "\033[0;34m";
        static const char __nowarn_unused *green    = "\033[0;32m";
        static const char __nowarn_unused *cyan     = "\033[0;36m";
        static const char __nowarn_unused *red      = "\033[0;31m";
        static const char __nowarn_unused *purple   = "\033[0;35m";
        static const char __nowarn_unused *brown    = "\033[0;33m";
        static const char __nowarn_unused *white    = "\033[0;37m";

        // High intensity colors
        static const char __nowarn_unused *hiblack  = "\033[1;30m";
        static const char __nowarn_unused *hiblue   = "\033[1;34m";
        static const char __nowarn_unused *higreen  = "\033[1;32m";
        static const char __nowarn_unused *hicyan   = "\033[1;36m";
        static const char __nowarn_unused *hired    = "\033[1;31m";
        static const char __nowarn_unused *hipurple = "\033[1;35m";
        static const char __nowarn_unused *hibrown  = "\033[1;33m";
        static const char __nowarn_unused *hiwhite  = "\033[1;37m";

        // Aliases
        static const char __nowarn_unused *yellow   = hibrown;
#   else
        // Return to default color
        static const char __nowarn_unused *none     = "";

        // Low intensity colors
        static const char __nowarn_unused *black    = "";
        static const char __nowarn_unused *blue     = "";
        static const char __nowarn_unused *green    = "";
        static const char __nowarn_unused *cyan     = "";
        static const char __nowarn_unused *red      = "";
        static const char __nowarn_unused *purple   = "";
        static const char __nowarn_unused *brown    = "";
        static const char __nowarn_unused *white    = "";

        // High intensity colors
        static const char __nowarn_unused *hiblack  = "";
        static const char __nowarn_unused *hiblue   = "";
        static const char __nowarn_unused *higreen  = "";
        static const char __nowarn_unused *hicyan   = "";
        static const char __nowarn_unused *hired    = "";
        static const char __nowarn_unused *hipurple = "";
        static const char __nowarn_unused *hibrown  = "";
        static const char __nowarn_unused *hiwhite  = "";

        // Aliases
        static const char __nowarn_unused *yellow   = hibrown;
#   endif
}


// -----------------------------------------------------------------------------
// Converts a time in seconds (double) into the more readable days, hours, mins,
// and seconds (std::string).
// -----------------------------------------------------------------------------
std::string formatTime (const double a_seconds);


// -----------------------------------------------------------------------------
// Print a Vector<> to pout
// -----------------------------------------------------------------------------
#define PrintVector(v)  pout() << #v " = "; _PrintVector( (v) );

template<typename T>
inline void _PrintVector(const Vector<T>& v) {
    pout() << "(";
    for(int i = 0; i < v.size(); ++i) {
        pout() << v[i];
        if(i != v.size()-1) pout() << ", ";
    }
    pout() << ")" << std::endl;
}


// -----------------------------------------------------------------------------
// Sends a RealVect to an ostream.
// -----------------------------------------------------------------------------
inline std::ostream& operator<<(std::ostream& stream, RealVect& rv)
{
    stream << "("
           D_TERM( << rv[0],
                   << ", " << rv[1],
                   << ", " << rv[2] )
           << ")";
    return stream;
}


// -----------------------------------------------------------------------------
// Data to file functions
// -----------------------------------------------------------------------------
// Write a FArrayBox to a text file
void writeTextFile (const FArrayBox&    a_data,
                    const char*         a_filename,
                    const Real&         a_time = -1);

// Write a LevelData to a text file
void writeTextFile (const LevelData<FArrayBox>& a_data,
                    const char*                 a_filename,
                    const Real&                 a_time = -1,
                    const bool                  a_writeGhosts = true);

// Write a FluxBox to a text file
void writeTextFile (const FluxBox&  a_data,
                    const char*     a_filename,
                    const Real&     a_time = -1);

void writeTextFile (const LevelData<FluxBox>&   a_data,
                    const char*                 a_filename,
                    const Real&                 a_time = -1);

void _writeLevelHDF5 (const LevelData<FArrayBox>& a_data,
                      const char*                 a_filename,
                      Real                        a_time = -1.0,
                      bool                        a_oneGhost = false);

void _writeLevelHDF5 (const LevelData<FluxBox>& a_data,
                      const char*               a_filename,
                      Real                      a_time = -1.0,
                      bool                      a_oneGhost = false);

#define writeLevelHDF5(d,t,g) \
do { \
    static int fidx = 0; \
    ostringstream fname; \
    fname << #d << "_" << fidx++ << ".hdf5"; \
    _writeLevelHDF5(d, fname.str().c_str(), (t), g); \
} while(0)


void _writeHDF5 (const Vector<LevelData<FArrayBox>*>& a_data,
                 const LevelGeometry&                 a_levGeo,
                 int                                  a_lmin,
                 int                                  a_lmax,
                 const char*                          a_filename,
                 const Vector<std::string>&           a_names = Vector<std::string>(),
                 const Real                           a_time = -1.0);

void _writeHDF5 (const Vector<LevelData<FluxBox>*>& a_data,
                 const LevelGeometry&               a_levGeo,
                 int                                a_lmin,
                 int                                a_lmax,
                 const char*                        a_filename);

#define writeHDF5CompNames(d,n,lg) \
do { \
    static int fidx = 0; \
    ostringstream fname; \
    fname << #d << "_" << fidx++ << ".hdf5"; \
    _writeHDF5((d), (lg), 0, -1, fname.str().c_str(), (n)); \
} while(0)

#define writeLevelsHDF5(d,lg,min,max) \
do { \
    static int fidx = 0; \
    ostringstream fname; \
    fname << #d << "_" << fidx++ << ".hdf5"; \
    _writeHDF5((d), (lg), (min), (max), fname.str().c_str()); \
} while(0)

#define writeHDF5(d,lg) writeLevelsHDF5(d,(lg),0,-1)

#define writeLevelHDF5CompNames(d,n,lg) \
do { \
    const Vector<LevelData<FArrayBox>*> level_##d(1, (LevelData<FArrayBox>*)&(d)); \
    writeHDF5CompNames(level_##d, (n), (lg)); \
} while(0)


///
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   {\bf Arguments:}\\
 filename  :  file to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_vectNames:  names of variables.\\
 a_domain :  domain at coarsest level.\\
 a_dx     :  grid spacing in each direction at coarsest level.\\
 a_dt     :  time step at coarsest level.\\
 a_time     :  time.\\
 a_vectRatio :  refinement ratio in each direction at all levels
 (ith entry is refinement ratio in each direction between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is blocking.

*/
void
WriteAnisotropicAMRHierarchyHDF5(
    const string& filename,
    const Vector<DisjointBoxLayout>& a_vectGrids,
    const Vector<LevelData<FArrayBox>* > & a_vectData,
    const Vector<string>& a_vectNames,
    const Box& a_domain,
    const RealVect& a_dx,
    const Real& a_dt,
    const Real& a_time,
    const Vector<IntVect>& a_vectRatios,
    const int& a_numLevels);


//
/**
   Writes hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.

   {\bf Arguments:}\\
handle  :  handle to output to.\\
 a_vectGrids : grids at each level.\\
 a_vectData :  data at each level.\\
 a_vectNames:  names of variables.\\
 a_domain :  domain at coarsest level.\\
 a_dx     :  grid spacing in each direction at coarsest level.\\
 a_dt     :  time step at coarsest level.\\
 a_time     :  time.\\
 a_vectRatio :  refinement ratio in each direction at all levels
 (ith entry is refinement ratio between levels i and i + 1).\\
 a_numLevels :  number of levels to output.\\
This is not blocking.

*/
void
WriteAnisotropicAMRHierarchyHDF5(
    HDF5Handle& handle,
    const Vector<DisjointBoxLayout>& a_vectGrids,
    const Vector<LevelData<FArrayBox>* > & a_vectData,
    const Vector<string>& a_vectNames,
    const Box& a_domain,
    const RealVect& a_dx,
    const Real& a_dt,
    const Real& a_time,
    const Vector<IntVect>& a_vectRatios,
    const int& a_numLevels);

///
/**
   Reads hierarchy of levels in ANISOTROPIC HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
filename  :  file to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_vectNames:  names of variables.\\
a_domain :  domain at coarsest level.\\
a_dx     :  grid spacing at coarsest level.\\
a_dt     :  time step at coarsest level.\\
a_time     :  time.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is blocking
*/
int
ReadAnisotropicAMRHierarchyHDF5(const string& filename,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Vector<string>& a_vectNames,
                     Box& a_domain,
                     RealVect& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     Vector<IntVect>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in ANISOTROPIC HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
handle  :  handle to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_vectNames:  names of variables.\\
a_domain :  domain at coarsest level.\\
a_dx     :  grid spacing at coarsest level.\\
a_dt     :  time step at coarsest level.\\
a_time     :  time.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is not blocking
*/
int
ReadAnisotropicAMRHierarchyHDF5(HDF5Handle& handle,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Vector<string>& a_vectNames,
                     Box& a_domain,
                     RealVect& a_dx,
                     Real& a_dt,
                     Real& a_time,
                     Vector<IntVect>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in ANISOTROPIC HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
filename  :  file to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_domain :  domain at coarsest level.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is blocking.
*/
int
ReadAnisotropicAMRHierarchyHDF5(const string& filename,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     Vector<IntVect>& a_refRatio,
                     int& a_numLevels);

///
/**
   Reads hierarchy of levels in HDF5 format.  Only available if the
   preprocessor macro HDF5 is defined at compilation.  Returns 0 on success.

   {\bf Arguments:}\\
handle  :  handle to input from.\\
a_vectGrids : grids at each level.\\
a_vectData :  data at each level.\\
a_domain :  domain at coarsest level.\\
a_vectRatio :  refinement ratio at all levels
(ith entry is refinement ratio between levels i and i + 1).\\
a_numLevels :  number of levels to read.\\

{\bf Returns:} \\
status code with values:\\
0: success \\
-1: bogus number of levels \\
-2: bogus number of components \\
-3: error in readlevel \\
-4: file open failed \\
This is not blocking.
*/
int
ReadAnisotropicAMRHierarchyHDF5(HDF5Handle& handle,
                     Vector<DisjointBoxLayout>& a_vectGrids,
                     Vector<LevelData<FArrayBox>* > & a_vectData,
                     Box& a_domain,
                     Vector<IntVect>& a_refRatio,
                     int& a_numLevels);


#endif //!__PRINTING_H__INCLUDED__

