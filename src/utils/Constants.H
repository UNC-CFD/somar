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
#ifndef __CONSTANTS_H__INCLUDED__
#define __CONSTANTS_H__INCLUDED__

#include "REAL.H"
#include "CONSTANTS.H"
#include "Box.H"


// Mathematical constants
static const Real E		= 2.718281828459045;
static const Real SQRT2	= 1.414213562373095;
static const Real SQRT3	= 1.732050807568877;


// NaN
#define _USE_NAN_

#ifndef _USE_NAN_
#   define quietNAN 4.321e300
#   define loudNAN  4.321e300
#else // !_USE_NAN_
#   include <limits>
#   define quietNAN std::numeric_limits<Real>::quiet_NaN()
#   define loudNAN  std::numeric_limits<Real>::signaling_NaN()
#endif // !_USE_NAN_


// The canonical empty box.
const Box EMPTYBOX;


#endif //!__CONSTANTS_H__INCLUDED__

