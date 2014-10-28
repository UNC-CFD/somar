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

// This file includes all of the debugging utilities.

#ifndef __DEBUG_H__INCLUDED__
#define __DEBUG_H__INCLUDED__

#include <string>
#include <sstream>
#include "MayDay.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "LevelData.H"

#ifndef NDEBUG
#   include "DebugOut.H"
#endif


// Nicer than using #ifndef NDEBUG...
#ifdef NDEBUG
#	define DEBUG_RAWCODE(x)
#	define DEBUG_CODE(x)
#else
#	define DEBUG_RAWCODE(x) x
#	define DEBUG_CODE(x) { x }
#endif


// Allows inline assignments that depend on release mode
// For example:
// const char* dbgMsg = IF_DEBUG("Debug mode!", "Release mode!");
#ifndef NDEBUG
#   define IF_DEBUG(d,r) d
#else
#   define IF_DEBUG(d,r) r
#endif

#define DEBUG_BOOL IF_DEBUG(true, false)


// -----------------------------------------------------------------------------
// Puts a checkpoint into pout.*
// -----------------------------------------------------------------------------
#define CHECKPOINT()	\
	DEBUG_RAWCODE(_WriteCheckpoint(__PRETTY_FUNCTION__, __FILE__, __LINE__))
#define CHECKMARK(mark)	\
	DEBUG_RAWCODE(_WriteCheckpoint(__PRETTY_FUNCTION__, __FILE__, __LINE__, mark))
void _WriteCheckpoint(std::string a_func, std::string a_file, unsigned int a_line, std::string a_mark = "");


// -----------------------------------------------------------------------------
// Enables crude use of streams with MayDay::Error
//	NOTE: PrintFunctionStack();\ removed from line just after do {
// -----------------------------------------------------------------------------
#define ERROR(x)							\
	do {									\
		ostringstream ostr;					\
		ostr << x;							\
		MayDay::Error(ostr.str().c_str());	\
	} while(0)


// -----------------------------------------------------------------------------
// Throws an error and reports location
// -----------------------------------------------------------------------------
#define UNDEFINED_FUNCTION()	\
	ERROR(__PRETTY_FUNCTION__ << " is undefined")


// -----------------------------------------------------------------------------
// Put this in all functions that need to be implemented.
// This will not throw and error, just display a message at runtime.
// -----------------------------------------------------------------------------
#define TODO() {static bool w = true; if(w) {_todo(__FILE__, __PRETTY_FUNCTION__, __LINE__); w = false;}}
void _todo(const char* a_filename, const char* a_funcname, const int a_linenumber);


// -----------------------------------------------------------------------------
// Performs an enhanced compatibility check on 2 DisjointBoxLayouts
// -----------------------------------------------------------------------------
#ifndef NDEBUG
#   define DBLCHECK(dbl1, dbl2)                 \
        if (!(dbl1 == dbl2)) {                  \
            ostringstream s1,s2;                \
            s1 << dbl1;                         \
            s2 << dbl2;                         \
            CH_assert(s1.str() == s2.str());    \
        }                                       \
        CH_assert(dbl1.compatible(dbl2))
#else
#   define DBLCHECK(dbl1, dbl2)
#endif


// -----------------------------------------------------------------------------
// Throws an error if a NAN is found in the testBox / valid regions.
// -----------------------------------------------------------------------------
typedef struct {
#ifndef NDEBUG
    bool    problemFound;
    IntVect pos;
    int     comp;
    Box     box;
    Real    val;
    string  msg;
#endif
} naninfo_type;

void checkForNAN (const FArrayBox& a_data, const Box& a_testBox, naninfo_type* a_info = NULL);
void checkForNAN (const FluxBox& a_data, const Box& a_testBox, naninfo_type* a_info = NULL);
void checkForValidNAN (const LevelData<FArrayBox>& a_data);
void checkForValidNAN (const LevelData<FluxBox>& a_data);


// -----------------------------------------------------------------------------
// Initializes data holders to NAN
// -----------------------------------------------------------------------------
#ifndef NDEBUG
    // Set a_levels[a_min:a_max] to NAN.
    void debugInitLevels (Vector<BoxLayoutData<FArrayBox>* >& a_levels,
                          int                                 a_min,
                          int                                 a_max);

    // Set a_level to NAN.
    void debugInitLevel (BoxLayoutData<FArrayBox>& a_level);

    // Set a_level to NAN.
    void debugInitLevel (BoxLayoutData<FluxBox>& a_level);

    // Set FArrayBox to NAN.
    void debugInit (FArrayBox& a_fab);

    // Set FluxBox to NAN.
    void debugInit (FluxBox& a_flub);

#else
#   define debugInitLevels(a,b,c)
#   define debugInitLevel(a)
#   define debugInit(a)
#endif


#endif //!__DEBUG_H__INCLUDED__

