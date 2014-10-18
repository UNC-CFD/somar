/*******************************************************************************
 *    FILE: HermiteInterp.H
 *    DESC: Simple Hermite polynomial interpolation utilities.
 *    DATE: Mon 26 May 2014 06:17:58 PM EDT  (UTC-0400)
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef HermiteInterp_H__INCLUDED__
#define HermiteInterp_H__INCLUDED__

#include "FArrayBox.H"


void HermiteInterp2D (FArrayBox&          a_fInterp,
                      const FArrayBox&    a_xInterp,
                      const FArrayBox&    a_yInterp,
                      const Box&          a_interpBox,
                      const int           a_xdir,
                      const int           a_ydir,
                      const Vector<Real>& a_x,
                      const Vector<Real>& a_y,
                      const FArrayBox&    a_f,
                      const FArrayBox&    a_dfdx,
                      const FArrayBox&    a_dfdy);


#endif //!HermiteInterp_H__INCLUDED__
