#ifndef _GEOSOURCEINTERFACEF_F_H_
#define _GEOSOURCEINTERFACEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDDEFAULT_FILL_J_2D 
#define GUARDDEFAULT_FILL_J_2D 
// Prototype for Fortran procedure DEFAULT_FILL_J_2D ...
//
void FORTRAN_NAME( DEFAULT_FILL_J_2D ,default_fill_j_2d )(
      CHFp_FRA1(J)
      ,CHFp_CONST_FRA(x0)
      ,CHFp_CONST_FRA(x1)
      ,CHFp_BOX(destBox) );

#define FORT_DEFAULT_FILL_J_2D FORTRAN_NAME( inlineDEFAULT_FILL_J_2D, inlineDEFAULT_FILL_J_2D)
#define FORTNT_DEFAULT_FILL_J_2D FORTRAN_NAME( DEFAULT_FILL_J_2D, default_fill_j_2d)

inline void FORTRAN_NAME(inlineDEFAULT_FILL_J_2D, inlineDEFAULT_FILL_J_2D)(
      CHFp_FRA1(J)
      ,CHFp_CONST_FRA(x0)
      ,CHFp_CONST_FRA(x1)
      ,CHFp_BOX(destBox) )
{
 CH_TIMELEAF("FORT_DEFAULT_FILL_J_2D");
 FORTRAN_NAME( DEFAULT_FILL_J_2D ,default_fill_j_2d )(
      CHFt_FRA1(J)
      ,CHFt_CONST_FRA(x0)
      ,CHFt_CONST_FRA(x1)
      ,CHFt_BOX(destBox) );
}
#endif  // GUARDDEFAULT_FILL_J_2D 

#ifndef GUARDDEFAULT_FILL_J_3D 
#define GUARDDEFAULT_FILL_J_3D 
// Prototype for Fortran procedure DEFAULT_FILL_J_3D ...
//
void FORTRAN_NAME( DEFAULT_FILL_J_3D ,default_fill_j_3d )(
      CHFp_FRA1(J)
      ,CHFp_CONST_FRA(x0)
      ,CHFp_CONST_FRA(x1)
      ,CHFp_CONST_FRA(x2)
      ,CHFp_BOX(destBox) );

#define FORT_DEFAULT_FILL_J_3D FORTRAN_NAME( inlineDEFAULT_FILL_J_3D, inlineDEFAULT_FILL_J_3D)
#define FORTNT_DEFAULT_FILL_J_3D FORTRAN_NAME( DEFAULT_FILL_J_3D, default_fill_j_3d)

inline void FORTRAN_NAME(inlineDEFAULT_FILL_J_3D, inlineDEFAULT_FILL_J_3D)(
      CHFp_FRA1(J)
      ,CHFp_CONST_FRA(x0)
      ,CHFp_CONST_FRA(x1)
      ,CHFp_CONST_FRA(x2)
      ,CHFp_BOX(destBox) )
{
 CH_TIMELEAF("FORT_DEFAULT_FILL_J_3D");
 FORTRAN_NAME( DEFAULT_FILL_J_3D ,default_fill_j_3d )(
      CHFt_FRA1(J)
      ,CHFt_CONST_FRA(x0)
      ,CHFt_CONST_FRA(x1)
      ,CHFt_CONST_FRA(x2)
      ,CHFt_BOX(destBox) );
}
#endif  // GUARDDEFAULT_FILL_J_3D 

#ifndef GUARDDEFAULT_FILL_GAMMA 
#define GUARDDEFAULT_FILL_GAMMA 
// Prototype for Fortran procedure DEFAULT_FILL_GAMMA ...
//
void FORTRAN_NAME( DEFAULT_FILL_GAMMA ,default_fill_gamma )(
      CHFp_FRA1(gamma)
      ,CHFp_CONST_FRA(gup)
      ,CHFp_CONST_FRA(gdn)
      ,CHFp_CONST_INT(a)
      ,CHFp_CONST_INT(b)
      ,CHFp_CONST_INT(c)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REAL(multiplier)
      ,CHFp_BOX(destBox) );

#define FORT_DEFAULT_FILL_GAMMA FORTRAN_NAME( inlineDEFAULT_FILL_GAMMA, inlineDEFAULT_FILL_GAMMA)
#define FORTNT_DEFAULT_FILL_GAMMA FORTRAN_NAME( DEFAULT_FILL_GAMMA, default_fill_gamma)

inline void FORTRAN_NAME(inlineDEFAULT_FILL_GAMMA, inlineDEFAULT_FILL_GAMMA)(
      CHFp_FRA1(gamma)
      ,CHFp_CONST_FRA(gup)
      ,CHFp_CONST_FRA(gdn)
      ,CHFp_CONST_INT(a)
      ,CHFp_CONST_INT(b)
      ,CHFp_CONST_INT(c)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REAL(multiplier)
      ,CHFp_BOX(destBox) )
{
 CH_TIMELEAF("FORT_DEFAULT_FILL_GAMMA");
 FORTRAN_NAME( DEFAULT_FILL_GAMMA ,default_fill_gamma )(
      CHFt_FRA1(gamma)
      ,CHFt_CONST_FRA(gup)
      ,CHFt_CONST_FRA(gdn)
      ,CHFt_CONST_INT(a)
      ,CHFt_CONST_INT(b)
      ,CHFt_CONST_INT(c)
      ,CHFt_CONST_REALVECT(dXi)
      ,CHFt_CONST_REAL(multiplier)
      ,CHFt_BOX(destBox) );
}
#endif  // GUARDDEFAULT_FILL_GAMMA 

#ifndef GUARDSIMPLECCDERIV 
#define GUARDSIMPLECCDERIV 
// Prototype for Fortran procedure SIMPLECCDERIV ...
//
void FORTRAN_NAME( SIMPLECCDERIV ,simpleccderiv )(
      CHFp_FRA1(dest)
      ,CHFp_CONST_FRA1(src)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(destBox) );

#define FORT_SIMPLECCDERIV FORTRAN_NAME( inlineSIMPLECCDERIV, inlineSIMPLECCDERIV)
#define FORTNT_SIMPLECCDERIV FORTRAN_NAME( SIMPLECCDERIV, simpleccderiv)

inline void FORTRAN_NAME(inlineSIMPLECCDERIV, inlineSIMPLECCDERIV)(
      CHFp_FRA1(dest)
      ,CHFp_CONST_FRA1(src)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(destBox) )
{
 CH_TIMELEAF("FORT_SIMPLECCDERIV");
 FORTRAN_NAME( SIMPLECCDERIV ,simpleccderiv )(
      CHFt_FRA1(dest)
      ,CHFt_CONST_FRA1(src)
      ,CHFt_CONST_REAL(scale)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(destBox) );
}
#endif  // GUARDSIMPLECCDERIV 

#ifndef GUARDSIMPLEFCDERIV 
#define GUARDSIMPLEFCDERIV 
// Prototype for Fortran procedure SIMPLEFCDERIV ...
//
void FORTRAN_NAME( SIMPLEFCDERIV ,simplefcderiv )(
      CHFp_FRA1(dest)
      ,CHFp_CONST_FRA1(src)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(destBox) );

#define FORT_SIMPLEFCDERIV FORTRAN_NAME( inlineSIMPLEFCDERIV, inlineSIMPLEFCDERIV)
#define FORTNT_SIMPLEFCDERIV FORTRAN_NAME( SIMPLEFCDERIV, simplefcderiv)

inline void FORTRAN_NAME(inlineSIMPLEFCDERIV, inlineSIMPLEFCDERIV)(
      CHFp_FRA1(dest)
      ,CHFp_CONST_FRA1(src)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(destBox) )
{
 CH_TIMELEAF("FORT_SIMPLEFCDERIV");
 FORTRAN_NAME( SIMPLEFCDERIV ,simplefcderiv )(
      CHFt_FRA1(dest)
      ,CHFt_CONST_FRA1(src)
      ,CHFt_CONST_REAL(scale)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(destBox) );
}
#endif  // GUARDSIMPLEFCDERIV 

#ifndef GUARDADDPROD2 
#define GUARDADDPROD2 
// Prototype for Fortran procedure ADDPROD2 ...
//
void FORTRAN_NAME( ADDPROD2 ,addprod2 )(
      CHFp_FRA1(dest)
      ,CHFp_CONST_FRA1(p0)
      ,CHFp_CONST_FRA1(p1)
      ,CHFp_BOX(destBox) );

#define FORT_ADDPROD2 FORTRAN_NAME( inlineADDPROD2, inlineADDPROD2)
#define FORTNT_ADDPROD2 FORTRAN_NAME( ADDPROD2, addprod2)

inline void FORTRAN_NAME(inlineADDPROD2, inlineADDPROD2)(
      CHFp_FRA1(dest)
      ,CHFp_CONST_FRA1(p0)
      ,CHFp_CONST_FRA1(p1)
      ,CHFp_BOX(destBox) )
{
 CH_TIMELEAF("FORT_ADDPROD2");
 FORTRAN_NAME( ADDPROD2 ,addprod2 )(
      CHFt_FRA1(dest)
      ,CHFt_CONST_FRA1(p0)
      ,CHFt_CONST_FRA1(p1)
      ,CHFt_BOX(destBox) );
}
#endif  // GUARDADDPROD2 

#ifndef GUARDSUBPROD2 
#define GUARDSUBPROD2 
// Prototype for Fortran procedure SUBPROD2 ...
//
void FORTRAN_NAME( SUBPROD2 ,subprod2 )(
      CHFp_FRA1(dest)
      ,CHFp_CONST_FRA1(p0)
      ,CHFp_CONST_FRA1(p1)
      ,CHFp_BOX(destBox) );

#define FORT_SUBPROD2 FORTRAN_NAME( inlineSUBPROD2, inlineSUBPROD2)
#define FORTNT_SUBPROD2 FORTRAN_NAME( SUBPROD2, subprod2)

inline void FORTRAN_NAME(inlineSUBPROD2, inlineSUBPROD2)(
      CHFp_FRA1(dest)
      ,CHFp_CONST_FRA1(p0)
      ,CHFp_CONST_FRA1(p1)
      ,CHFp_BOX(destBox) )
{
 CH_TIMELEAF("FORT_SUBPROD2");
 FORTRAN_NAME( SUBPROD2 ,subprod2 )(
      CHFt_FRA1(dest)
      ,CHFt_CONST_FRA1(p0)
      ,CHFt_CONST_FRA1(p1)
      ,CHFt_BOX(destBox) );
}
#endif  // GUARDSUBPROD2 

}

#endif
