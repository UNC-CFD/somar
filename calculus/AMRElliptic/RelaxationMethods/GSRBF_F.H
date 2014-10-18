#ifndef _GSRBF_F_H_
#define _GSRBF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGSRBITER1D 
#define GUARDGSRBITER1D 
// Prototype for Fortran procedure GSRBITER1D ...
//
void FORTRAN_NAME( GSRBITER1D ,gsrbiter1d )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(JgDir)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_FRA1(lapDiag)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(lineDir)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBITER1D FORTRAN_NAME( inlineGSRBITER1D, inlineGSRBITER1D)
#define FORTNT_GSRBITER1D FORTRAN_NAME( GSRBITER1D, gsrbiter1d)

inline void FORTRAN_NAME(inlineGSRBITER1D, inlineGSRBITER1D)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(JgDir)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_FRA1(lapDiag)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(lineDir)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBITER1D");
 FORTRAN_NAME( GSRBITER1D ,gsrbiter1d )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA(JgDir)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_CONST_FRA1(lapDiag)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(lineDir)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBITER1D 

#ifndef GUARDGSRBITER2D 
#define GUARDGSRBITER2D 
// Prototype for Fortran procedure GSRBITER2D ...
//
void FORTRAN_NAME( GSRBITER2D ,gsrbiter2d )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(extrap)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_FRA1(lapDiag)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBITER2D FORTRAN_NAME( inlineGSRBITER2D, inlineGSRBITER2D)
#define FORTNT_GSRBITER2D FORTRAN_NAME( GSRBITER2D, gsrbiter2d)

inline void FORTRAN_NAME(inlineGSRBITER2D, inlineGSRBITER2D)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(extrap)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_FRA1(lapDiag)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBITER2D");
 FORTRAN_NAME( GSRBITER2D ,gsrbiter2d )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(extrap)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA(Jg0)
      ,CHFt_CONST_FRA(Jg1)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_CONST_FRA1(lapDiag)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBITER2D 

#ifndef GUARDGSRBITER3D 
#define GUARDGSRBITER3D 
// Prototype for Fortran procedure GSRBITER3D ...
//
void FORTRAN_NAME( GSRBITER3D ,gsrbiter3d )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(extrap)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA(Jg2)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_FRA1(lapDiag)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBITER3D FORTRAN_NAME( inlineGSRBITER3D, inlineGSRBITER3D)
#define FORTNT_GSRBITER3D FORTRAN_NAME( GSRBITER3D, gsrbiter3d)

inline void FORTRAN_NAME(inlineGSRBITER3D, inlineGSRBITER3D)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(extrap)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA(Jg2)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_FRA1(lapDiag)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBITER3D");
 FORTRAN_NAME( GSRBITER3D ,gsrbiter3d )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(extrap)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA(Jg0)
      ,CHFt_CONST_FRA(Jg1)
      ,CHFt_CONST_FRA(Jg2)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_CONST_FRA1(lapDiag)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBITER3D 

#ifndef GUARDGSRBITER2DORTHO 
#define GUARDGSRBITER2DORTHO 
// Prototype for Fortran procedure GSRBITER2DORTHO ...
//
void FORTRAN_NAME( GSRBITER2DORTHO ,gsrbiter2dortho )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA1(Jgxx)
      ,CHFp_CONST_FRA1(Jgyy)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_FRA1(lapDiag)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBITER2DORTHO FORTRAN_NAME( inlineGSRBITER2DORTHO, inlineGSRBITER2DORTHO)
#define FORTNT_GSRBITER2DORTHO FORTRAN_NAME( GSRBITER2DORTHO, gsrbiter2dortho)

inline void FORTRAN_NAME(inlineGSRBITER2DORTHO, inlineGSRBITER2DORTHO)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA1(Jgxx)
      ,CHFp_CONST_FRA1(Jgyy)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_FRA1(lapDiag)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBITER2DORTHO");
 FORTRAN_NAME( GSRBITER2DORTHO ,gsrbiter2dortho )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA1(Jgxx)
      ,CHFt_CONST_FRA1(Jgyy)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_CONST_FRA1(lapDiag)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBITER2DORTHO 

#ifndef GUARDGSRBITER3DORTHO 
#define GUARDGSRBITER3DORTHO 
// Prototype for Fortran procedure GSRBITER3DORTHO ...
//
void FORTRAN_NAME( GSRBITER3DORTHO ,gsrbiter3dortho )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA1(Jgxx)
      ,CHFp_CONST_FRA1(Jgyy)
      ,CHFp_CONST_FRA1(Jgzz)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_FRA1(lapDiag)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBITER3DORTHO FORTRAN_NAME( inlineGSRBITER3DORTHO, inlineGSRBITER3DORTHO)
#define FORTNT_GSRBITER3DORTHO FORTRAN_NAME( GSRBITER3DORTHO, gsrbiter3dortho)

inline void FORTRAN_NAME(inlineGSRBITER3DORTHO, inlineGSRBITER3DORTHO)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA1(Jgxx)
      ,CHFp_CONST_FRA1(Jgyy)
      ,CHFp_CONST_FRA1(Jgzz)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_FRA1(lapDiag)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBITER3DORTHO");
 FORTRAN_NAME( GSRBITER3DORTHO ,gsrbiter3dortho )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA1(Jgxx)
      ,CHFt_CONST_FRA1(Jgyy)
      ,CHFt_CONST_FRA1(Jgzz)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_CONST_FRA1(lapDiag)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBITER3DORTHO 

#ifndef GUARDGSRBBOUNDARYITER1D 
#define GUARDGSRBBOUNDARYITER1D 
// Prototype for Fortran procedure GSRBBOUNDARYITER1D ...
//
void FORTRAN_NAME( GSRBBOUNDARYITER1D ,gsrbboundaryiter1d )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(loStencil)
      ,CHFp_CONST_INT(hiStencil)
      ,CHFp_CONST_INT(lineDir)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBBOUNDARYITER1D FORTRAN_NAME( inlineGSRBBOUNDARYITER1D, inlineGSRBBOUNDARYITER1D)
#define FORTNT_GSRBBOUNDARYITER1D FORTRAN_NAME( GSRBBOUNDARYITER1D, gsrbboundaryiter1d)

inline void FORTRAN_NAME(inlineGSRBBOUNDARYITER1D, inlineGSRBBOUNDARYITER1D)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(loStencil)
      ,CHFp_CONST_INT(hiStencil)
      ,CHFp_CONST_INT(lineDir)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBBOUNDARYITER1D");
 FORTRAN_NAME( GSRBBOUNDARYITER1D ,gsrbboundaryiter1d )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA(Jg0)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(loStencil)
      ,CHFt_CONST_INT(hiStencil)
      ,CHFt_CONST_INT(lineDir)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBBOUNDARYITER1D 

#ifndef GUARDGSRBBOUNDARYITER2D 
#define GUARDGSRBBOUNDARYITER2D 
// Prototype for Fortran procedure GSRBBOUNDARYITER2D ...
//
void FORTRAN_NAME( GSRBBOUNDARYITER2D ,gsrbboundaryiter2d )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(extrap)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(loXStencil)
      ,CHFp_CONST_INT(hiXStencil)
      ,CHFp_CONST_INT(loYStencil)
      ,CHFp_CONST_INT(hiYStencil)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBBOUNDARYITER2D FORTRAN_NAME( inlineGSRBBOUNDARYITER2D, inlineGSRBBOUNDARYITER2D)
#define FORTNT_GSRBBOUNDARYITER2D FORTRAN_NAME( GSRBBOUNDARYITER2D, gsrbboundaryiter2d)

inline void FORTRAN_NAME(inlineGSRBBOUNDARYITER2D, inlineGSRBBOUNDARYITER2D)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(extrap)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(loXStencil)
      ,CHFp_CONST_INT(hiXStencil)
      ,CHFp_CONST_INT(loYStencil)
      ,CHFp_CONST_INT(hiYStencil)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBBOUNDARYITER2D");
 FORTRAN_NAME( GSRBBOUNDARYITER2D ,gsrbboundaryiter2d )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(extrap)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA(Jg0)
      ,CHFt_CONST_FRA(Jg1)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(loXStencil)
      ,CHFt_CONST_INT(hiXStencil)
      ,CHFt_CONST_INT(loYStencil)
      ,CHFt_CONST_INT(hiYStencil)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBBOUNDARYITER2D 

#ifndef GUARDGSRBBOUNDARYITER3D 
#define GUARDGSRBBOUNDARYITER3D 
// Prototype for Fortran procedure GSRBBOUNDARYITER3D ...
//
void FORTRAN_NAME( GSRBBOUNDARYITER3D ,gsrbboundaryiter3d )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(extrap)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA(Jg2)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(loXStencil)
      ,CHFp_CONST_INT(hiXStencil)
      ,CHFp_CONST_INT(loYStencil)
      ,CHFp_CONST_INT(hiYStencil)
      ,CHFp_CONST_INT(loZStencil)
      ,CHFp_CONST_INT(hiZStencil)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBBOUNDARYITER3D FORTRAN_NAME( inlineGSRBBOUNDARYITER3D, inlineGSRBBOUNDARYITER3D)
#define FORTNT_GSRBBOUNDARYITER3D FORTRAN_NAME( GSRBBOUNDARYITER3D, gsrbboundaryiter3d)

inline void FORTRAN_NAME(inlineGSRBBOUNDARYITER3D, inlineGSRBBOUNDARYITER3D)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(extrap)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA(Jg2)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(loXStencil)
      ,CHFp_CONST_INT(hiXStencil)
      ,CHFp_CONST_INT(loYStencil)
      ,CHFp_CONST_INT(hiYStencil)
      ,CHFp_CONST_INT(loZStencil)
      ,CHFp_CONST_INT(hiZStencil)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBBOUNDARYITER3D");
 FORTRAN_NAME( GSRBBOUNDARYITER3D ,gsrbboundaryiter3d )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(extrap)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA(Jg0)
      ,CHFt_CONST_FRA(Jg1)
      ,CHFt_CONST_FRA(Jg2)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(loXStencil)
      ,CHFt_CONST_INT(hiXStencil)
      ,CHFt_CONST_INT(loYStencil)
      ,CHFt_CONST_INT(hiYStencil)
      ,CHFt_CONST_INT(loZStencil)
      ,CHFt_CONST_INT(hiZStencil)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBBOUNDARYITER3D 

#ifndef GUARDGSRBBOUNDARYITER2DORTHO 
#define GUARDGSRBBOUNDARYITER2DORTHO 
// Prototype for Fortran procedure GSRBBOUNDARYITER2DORTHO ...
//
void FORTRAN_NAME( GSRBBOUNDARYITER2DORTHO ,gsrbboundaryiter2dortho )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA1(Jgxx)
      ,CHFp_CONST_FRA1(Jgyy)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(loXStencil)
      ,CHFp_CONST_INT(hiXStencil)
      ,CHFp_CONST_INT(loYStencil)
      ,CHFp_CONST_INT(hiYStencil)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBBOUNDARYITER2DORTHO FORTRAN_NAME( inlineGSRBBOUNDARYITER2DORTHO, inlineGSRBBOUNDARYITER2DORTHO)
#define FORTNT_GSRBBOUNDARYITER2DORTHO FORTRAN_NAME( GSRBBOUNDARYITER2DORTHO, gsrbboundaryiter2dortho)

inline void FORTRAN_NAME(inlineGSRBBOUNDARYITER2DORTHO, inlineGSRBBOUNDARYITER2DORTHO)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA1(Jgxx)
      ,CHFp_CONST_FRA1(Jgyy)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(loXStencil)
      ,CHFp_CONST_INT(hiXStencil)
      ,CHFp_CONST_INT(loYStencil)
      ,CHFp_CONST_INT(hiYStencil)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBBOUNDARYITER2DORTHO");
 FORTRAN_NAME( GSRBBOUNDARYITER2DORTHO ,gsrbboundaryiter2dortho )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA1(Jgxx)
      ,CHFt_CONST_FRA1(Jgyy)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(loXStencil)
      ,CHFt_CONST_INT(hiXStencil)
      ,CHFt_CONST_INT(loYStencil)
      ,CHFt_CONST_INT(hiYStencil)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBBOUNDARYITER2DORTHO 

#ifndef GUARDGSRBBOUNDARYITER3DORTHO 
#define GUARDGSRBBOUNDARYITER3DORTHO 
// Prototype for Fortran procedure GSRBBOUNDARYITER3DORTHO ...
//
void FORTRAN_NAME( GSRBBOUNDARYITER3DORTHO ,gsrbboundaryiter3dortho )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA1(Jgxx)
      ,CHFp_CONST_FRA1(Jgyy)
      ,CHFp_CONST_FRA1(Jgzz)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(loXStencil)
      ,CHFp_CONST_INT(hiXStencil)
      ,CHFp_CONST_INT(loYStencil)
      ,CHFp_CONST_INT(hiYStencil)
      ,CHFp_CONST_INT(loZStencil)
      ,CHFp_CONST_INT(hiZStencil)
      ,CHFp_CONST_INT(redBlack) );

#define FORT_GSRBBOUNDARYITER3DORTHO FORTRAN_NAME( inlineGSRBBOUNDARYITER3DORTHO, inlineGSRBBOUNDARYITER3DORTHO)
#define FORTNT_GSRBBOUNDARYITER3DORTHO FORTRAN_NAME( GSRBBOUNDARYITER3DORTHO, gsrbboundaryiter3dortho)

inline void FORTRAN_NAME(inlineGSRBBOUNDARYITER3DORTHO, inlineGSRBBOUNDARYITER3DORTHO)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA1(Jgxx)
      ,CHFp_CONST_FRA1(Jgyy)
      ,CHFp_CONST_FRA1(Jgzz)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(loXStencil)
      ,CHFp_CONST_INT(hiXStencil)
      ,CHFp_CONST_INT(loYStencil)
      ,CHFp_CONST_INT(hiYStencil)
      ,CHFp_CONST_INT(loZStencil)
      ,CHFp_CONST_INT(hiZStencil)
      ,CHFp_CONST_INT(redBlack) )
{
 CH_TIMELEAF("FORT_GSRBBOUNDARYITER3DORTHO");
 FORTRAN_NAME( GSRBBOUNDARYITER3DORTHO ,gsrbboundaryiter3dortho )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA1(Jgxx)
      ,CHFt_CONST_FRA1(Jgyy)
      ,CHFt_CONST_FRA1(Jgzz)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(loXStencil)
      ,CHFt_CONST_INT(hiXStencil)
      ,CHFt_CONST_INT(loYStencil)
      ,CHFt_CONST_INT(hiYStencil)
      ,CHFt_CONST_INT(loZStencil)
      ,CHFt_CONST_INT(hiZStencil)
      ,CHFt_CONST_INT(redBlack) );
}
#endif  // GUARDGSRBBOUNDARYITER3DORTHO 

#ifndef GUARDLINEGSRBITER2D 
#define GUARDLINEGSRBITER2D 
// Prototype for Fortran procedure LineGSRBIter2D ...
//
void FORTRAN_NAME( LINEGSRBITER2D ,linegsrbiter2d )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(dzCrse)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack)
      ,CHFp_FRA1(DU)
      ,CHFp_FRA1(D)
      ,CHFp_FRA1(DL)
      ,CHFp_FRA1(B)
      ,CHFp_CONST_INT(loXBCType)
      ,CHFp_CONST_INT(hiXBCType)
      ,CHFp_CONST_INT(loYBCType)
      ,CHFp_CONST_INT(hiYBCType) );

#define FORT_LINEGSRBITER2D FORTRAN_NAME( inlineLINEGSRBITER2D, inlineLINEGSRBITER2D)
#define FORTNT_LINEGSRBITER2D FORTRAN_NAME( LINEGSRBITER2D, linegsrbiter2d)

inline void FORTRAN_NAME(inlineLINEGSRBITER2D, inlineLINEGSRBITER2D)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(dzCrse)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack)
      ,CHFp_FRA1(DU)
      ,CHFp_FRA1(D)
      ,CHFp_FRA1(DL)
      ,CHFp_FRA1(B)
      ,CHFp_CONST_INT(loXBCType)
      ,CHFp_CONST_INT(hiXBCType)
      ,CHFp_CONST_INT(loYBCType)
      ,CHFp_CONST_INT(hiYBCType) )
{
 CH_TIMELEAF("FORT_LINEGSRBITER2D");
 FORTRAN_NAME( LINEGSRBITER2D ,linegsrbiter2d )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(extrap)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_FRA(Jg0)
      ,CHFt_CONST_FRA(Jg1)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(dzCrse)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(redBlack)
      ,CHFt_FRA1(DU)
      ,CHFt_FRA1(D)
      ,CHFt_FRA1(DL)
      ,CHFt_FRA1(B)
      ,CHFt_CONST_INT(loXBCType)
      ,CHFt_CONST_INT(hiXBCType)
      ,CHFt_CONST_INT(loYBCType)
      ,CHFt_CONST_INT(hiYBCType) );
}
#endif  // GUARDLINEGSRBITER2D 

#ifndef GUARDLINEGSRBITER3D 
#define GUARDLINEGSRBITER3D 
// Prototype for Fortran procedure LineGSRBIter3D ...
//
void FORTRAN_NAME( LINEGSRBITER3D ,linegsrbiter3d )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA(Jg2)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(dzCrse)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack)
      ,CHFp_FRA1(DU)
      ,CHFp_FRA1(D)
      ,CHFp_FRA1(DL)
      ,CHFp_FRA1(B)
      ,CHFp_CONST_INT(loXBCType)
      ,CHFp_CONST_INT(hiXBCType)
      ,CHFp_CONST_INT(loYBCType)
      ,CHFp_CONST_INT(hiYBCType)
      ,CHFp_CONST_INT(loZBCType)
      ,CHFp_CONST_INT(hiZBCType) );

#define FORT_LINEGSRBITER3D FORTRAN_NAME( inlineLINEGSRBITER3D, inlineLINEGSRBITER3D)
#define FORTNT_LINEGSRBITER3D FORTRAN_NAME( LINEGSRBITER3D, linegsrbiter3d)

inline void FORTRAN_NAME(inlineLINEGSRBITER3D, inlineLINEGSRBITER3D)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA(Jg2)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(dzCrse)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_INT(redBlack)
      ,CHFp_FRA1(DU)
      ,CHFp_FRA1(D)
      ,CHFp_FRA1(DL)
      ,CHFp_FRA1(B)
      ,CHFp_CONST_INT(loXBCType)
      ,CHFp_CONST_INT(hiXBCType)
      ,CHFp_CONST_INT(loYBCType)
      ,CHFp_CONST_INT(hiYBCType)
      ,CHFp_CONST_INT(loZBCType)
      ,CHFp_CONST_INT(hiZBCType) )
{
 CH_TIMELEAF("FORT_LINEGSRBITER3D");
 FORTRAN_NAME( LINEGSRBITER3D ,linegsrbiter3d )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(extrap)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_FRA(Jg0)
      ,CHFt_CONST_FRA(Jg1)
      ,CHFt_CONST_FRA(Jg2)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(dzCrse)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_INT(redBlack)
      ,CHFt_FRA1(DU)
      ,CHFt_FRA1(D)
      ,CHFt_FRA1(DL)
      ,CHFt_FRA1(B)
      ,CHFt_CONST_INT(loXBCType)
      ,CHFt_CONST_INT(hiXBCType)
      ,CHFt_CONST_INT(loYBCType)
      ,CHFt_CONST_INT(hiYBCType)
      ,CHFt_CONST_INT(loZBCType)
      ,CHFt_CONST_INT(hiZBCType) );
}
#endif  // GUARDLINEGSRBITER3D 

}

#endif
