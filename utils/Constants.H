/*******************************************************************************
 *	FILE: Constants.hpp
 *	DESC: Some useful definitions.
 *	DATE: Thu 18 Aug 2011 12:31:14 AM EDT
 *	MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __CONSTANTS_HPP__INCLUDED__
#define __CONSTANTS_HPP__INCLUDED__

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


#endif //!__CONSTANTS_HPP__INCLUDED__

