#ifdef CH_LANG_CC
/*******************************************************************************
 *    FILE: Subspace.H
 *    DESC: Utilities for defining and storing data on flattened spaces.
 *    DATE: Sun 12 May 2013 01:07:22 PM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/
#endif


#ifndef __Subspace_H__INCLUDED__
#define __Subspace_H__INCLUDED__
#ifdef CH_LANG_CC

#include "ProblemDomain.H"

// Takes an N-dimensional Box and projects it to a subspace defined by a_mask.
// a_mask should be 0 in flattened directions and 1 in unmodified directions.
// In other words, a_mask flags the tangential directions of the subspace.
Box flattenBox (const Box&     a_box,
                const IntVect& a_mask);

// Takes an N-dimensional Box and projects it to the (N_1)-dimensional surface
// whose normal direction is a_normDir.
Box flattenBox (const Box& a_box,
                const int  a_normDir);

// Returns a box that can define holders for functions that only depend on the
// vertical coordinate.
Box verticalDataBox (const ProblemDomain& a_domain);

// Returns a box that can define holders for functions that only depend on the
// vertical coordinate.
Box verticalDataBox (const Box& a_box);

// Returns a box that can define holders for functions that only depend on the
// horizontal coordinates.
Box horizontalDataBox (const ProblemDomain& a_domain);

// Returns a box that can define holders for functions that only depend on the
// horizontal coordinates.
Box horizontalDataBox (const Box& a_box);

#else // End C++ code, begin Fortran code.

#if CH_SPACEDIM == 2
#  define FLATIX(i)       i##0     , 0
#  define FLATOFFSET(i,o) i##0 o##0, 0
#  define VERTCOMP(i)                i##1
#elif CH_SPACEDIM == 3
#  define FLATIX(i)       i##0     , i##1     , 0
#  define FLATOFFSET(i,o) i##0 o##0, i##1 o##1, 0
#  define VERTCOMP(i)                           i##2
#else
#  error Bad SpaceDim
#endif

#endif // End of fortran code
#endif //!__Subspace_H__INCLUDED__
