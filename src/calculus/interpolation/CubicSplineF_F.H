#ifndef _CUBICSPLINEF_F_H_
#define _CUBICSPLINEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCUBICSPLINE_SOLVESECONDDERIV 
#define GUARDCUBICSPLINE_SOLVESECONDDERIV 
// Prototype for Fortran procedure CubicSpline_SolveSecondDeriv ...
//
void FORTRAN_NAME( CUBICSPLINE_SOLVESECONDDERIV ,cubicspline_solvesecondderiv )(
      CHFp_VR(d2f)
      ,CHFp_CONST_VR(x)
      ,CHFp_CONST_VR(f)
      ,CHFp_CONST_REAL(lofbc)
      ,CHFp_CONST_REAL(hifbc)
      ,CHFp_VR(u) );

#define FORT_CUBICSPLINE_SOLVESECONDDERIV FORTRAN_NAME( inlineCUBICSPLINE_SOLVESECONDDERIV, inlineCUBICSPLINE_SOLVESECONDDERIV)
#define FORTNT_CUBICSPLINE_SOLVESECONDDERIV FORTRAN_NAME( CUBICSPLINE_SOLVESECONDDERIV, cubicspline_solvesecondderiv)

inline void FORTRAN_NAME(inlineCUBICSPLINE_SOLVESECONDDERIV, inlineCUBICSPLINE_SOLVESECONDDERIV)(
      CHFp_VR(d2f)
      ,CHFp_CONST_VR(x)
      ,CHFp_CONST_VR(f)
      ,CHFp_CONST_REAL(lofbc)
      ,CHFp_CONST_REAL(hifbc)
      ,CHFp_VR(u) )
{
 CH_TIMELEAF("FORT_CUBICSPLINE_SOLVESECONDDERIV");
 FORTRAN_NAME( CUBICSPLINE_SOLVESECONDDERIV ,cubicspline_solvesecondderiv )(
      CHFt_VR(d2f)
      ,CHFt_CONST_VR(x)
      ,CHFt_CONST_VR(f)
      ,CHFt_CONST_REAL(lofbc)
      ,CHFt_CONST_REAL(hifbc)
      ,CHFt_VR(u) );
}
#endif  // GUARDCUBICSPLINE_SOLVESECONDDERIV 

}

#endif
