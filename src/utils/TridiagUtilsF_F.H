#ifndef _TRIDIAGUTILSF_F_H_
#define _TRIDIAGUTILSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDTRIDIAGPOISSONNN1D 
#define GUARDTRIDIAGPOISSONNN1D 
// Prototype for Fortran procedure TriDiagPoissonNN1D ...
//
void FORTRAN_NAME( TRIDIAGPOISSONNN1D ,tridiagpoissonnn1d )(
      CHFp_VR(x)
      ,CHFp_CONST_VR(a)
      ,CHFp_CONST_VR(b)
      ,CHFp_CONST_VR(c)
      ,CHFp_VR(gam) );

#define FORT_TRIDIAGPOISSONNN1D FORTRAN_NAME( inlineTRIDIAGPOISSONNN1D, inlineTRIDIAGPOISSONNN1D)
#define FORTNT_TRIDIAGPOISSONNN1D FORTRAN_NAME( TRIDIAGPOISSONNN1D, tridiagpoissonnn1d)

inline void FORTRAN_NAME(inlineTRIDIAGPOISSONNN1D, inlineTRIDIAGPOISSONNN1D)(
      CHFp_VR(x)
      ,CHFp_CONST_VR(a)
      ,CHFp_CONST_VR(b)
      ,CHFp_CONST_VR(c)
      ,CHFp_VR(gam) )
{
 CH_TIMELEAF("FORT_TRIDIAGPOISSONNN1D");
 FORTRAN_NAME( TRIDIAGPOISSONNN1D ,tridiagpoissonnn1d )(
      CHFt_VR(x)
      ,CHFt_CONST_VR(a)
      ,CHFt_CONST_VR(b)
      ,CHFt_CONST_VR(c)
      ,CHFt_VR(gam) );
}
#endif  // GUARDTRIDIAGPOISSONNN1D 

#ifndef GUARDTRIDIAGPOISSONNN1DFAB 
#define GUARDTRIDIAGPOISSONNN1DFAB 
// Prototype for Fortran procedure TriDiagPoissonNN1DFAB ...
//
void FORTRAN_NAME( TRIDIAGPOISSONNN1DFAB ,tridiagpoissonnn1dfab )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA1(sigma)
      ,CHFp_BOX(bottomBox)
      ,CHFp_CONST_INT(Nx)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir) );

#define FORT_TRIDIAGPOISSONNN1DFAB FORTRAN_NAME( inlineTRIDIAGPOISSONNN1DFAB, inlineTRIDIAGPOISSONNN1DFAB)
#define FORTNT_TRIDIAGPOISSONNN1DFAB FORTRAN_NAME( TRIDIAGPOISSONNN1DFAB, tridiagpoissonnn1dfab)

inline void FORTRAN_NAME(inlineTRIDIAGPOISSONNN1DFAB, inlineTRIDIAGPOISSONNN1DFAB)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA1(sigma)
      ,CHFp_BOX(bottomBox)
      ,CHFp_CONST_INT(Nx)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_TRIDIAGPOISSONNN1DFAB");
 FORTRAN_NAME( TRIDIAGPOISSONNN1DFAB ,tridiagpoissonnn1dfab )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA1(sigma)
      ,CHFt_BOX(bottomBox)
      ,CHFt_CONST_INT(Nx)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDTRIDIAGPOISSONNN1DFAB 

}

#endif
