#ifndef _BATHYMETRICBASEMAPF_F_H_
#define _BATHYMETRICBASEMAPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDFILL_BATHYHORIZMAP 
#define GUARDFILL_BATHYHORIZMAP 
// Prototype for Fortran procedure FILL_BATHYHORIZMAP ...
//
void FORTRAN_NAME( FILL_BATHYHORIZMAP ,fill_bathyhorizmap )(
      CHFp_FRA1(xFAB)
      ,CHFp_BOX(xBox)
      ,CHFp_CONST_INT(xBoxDirType)
      ,CHFp_CONST_INT(xDir)
      ,CHFp_CONST_REAL(L)
      ,CHFp_CONST_REAL(dXi) );

#define FORT_FILL_BATHYHORIZMAP FORTRAN_NAME( inlineFILL_BATHYHORIZMAP, inlineFILL_BATHYHORIZMAP)
#define FORTNT_FILL_BATHYHORIZMAP FORTRAN_NAME( FILL_BATHYHORIZMAP, fill_bathyhorizmap)

inline void FORTRAN_NAME(inlineFILL_BATHYHORIZMAP, inlineFILL_BATHYHORIZMAP)(
      CHFp_FRA1(xFAB)
      ,CHFp_BOX(xBox)
      ,CHFp_CONST_INT(xBoxDirType)
      ,CHFp_CONST_INT(xDir)
      ,CHFp_CONST_REAL(L)
      ,CHFp_CONST_REAL(dXi) )
{
 CH_TIMELEAF("FORT_FILL_BATHYHORIZMAP");
 FORTRAN_NAME( FILL_BATHYHORIZMAP ,fill_bathyhorizmap )(
      CHFt_FRA1(xFAB)
      ,CHFt_BOX(xBox)
      ,CHFt_CONST_INT(xBoxDirType)
      ,CHFt_CONST_INT(xDir)
      ,CHFt_CONST_REAL(L)
      ,CHFt_CONST_REAL(dXi) );
}
#endif  // GUARDFILL_BATHYHORIZMAP 

#ifndef GUARDFILL_BATHYVERTMAP 
#define GUARDFILL_BATHYVERTMAP 
// Prototype for Fortran procedure FILL_BATHYVERTMAP ...
//
void FORTRAN_NAME( FILL_BATHYVERTMAP ,fill_bathyvertmap )(
      CHFp_FRA1(zFAB)
      ,CHFp_BOX(zBox)
      ,CHFp_CONST_INTVECT(zBoxType)
      ,CHFp_CONST_FRA1(depth)
      ,CHFp_CONST_REALVECT(L)
      ,CHFp_CONST_REALVECT(dXi) );

#define FORT_FILL_BATHYVERTMAP FORTRAN_NAME( inlineFILL_BATHYVERTMAP, inlineFILL_BATHYVERTMAP)
#define FORTNT_FILL_BATHYVERTMAP FORTRAN_NAME( FILL_BATHYVERTMAP, fill_bathyvertmap)

inline void FORTRAN_NAME(inlineFILL_BATHYVERTMAP, inlineFILL_BATHYVERTMAP)(
      CHFp_FRA1(zFAB)
      ,CHFp_BOX(zBox)
      ,CHFp_CONST_INTVECT(zBoxType)
      ,CHFp_CONST_FRA1(depth)
      ,CHFp_CONST_REALVECT(L)
      ,CHFp_CONST_REALVECT(dXi) )
{
 CH_TIMELEAF("FORT_FILL_BATHYVERTMAP");
 FORTRAN_NAME( FILL_BATHYVERTMAP ,fill_bathyvertmap )(
      CHFt_FRA1(zFAB)
      ,CHFt_BOX(zBox)
      ,CHFt_CONST_INTVECT(zBoxType)
      ,CHFt_CONST_FRA1(depth)
      ,CHFt_CONST_REALVECT(L)
      ,CHFt_CONST_REALVECT(dXi) );
}
#endif  // GUARDFILL_BATHYVERTMAP 

#ifndef GUARDFILL_BATHYDXDXI 
#define GUARDFILL_BATHYDXDXI 
// Prototype for Fortran procedure FILL_BATHYDXDXI ...
//
void FORTRAN_NAME( FILL_BATHYDXDXI ,fill_bathydxdxi )(
      CHFp_FRA1(dxdXi)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType)
      ,CHFp_CONST_INT(nu)
      ,CHFp_CONST_REALVECT(L)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REAL(scale) );

#define FORT_FILL_BATHYDXDXI FORTRAN_NAME( inlineFILL_BATHYDXDXI, inlineFILL_BATHYDXDXI)
#define FORTNT_FILL_BATHYDXDXI FORTRAN_NAME( FILL_BATHYDXDXI, fill_bathydxdxi)

inline void FORTRAN_NAME(inlineFILL_BATHYDXDXI, inlineFILL_BATHYDXDXI)(
      CHFp_FRA1(dxdXi)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType)
      ,CHFp_CONST_INT(nu)
      ,CHFp_CONST_REALVECT(L)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REAL(scale) )
{
 CH_TIMELEAF("FORT_FILL_BATHYDXDXI");
 FORTRAN_NAME( FILL_BATHYDXDXI ,fill_bathydxdxi )(
      CHFt_FRA1(dxdXi)
      ,CHFt_BOX(destBox)
      ,CHFt_CONST_INTVECT(destBoxType)
      ,CHFt_CONST_INT(nu)
      ,CHFt_CONST_REALVECT(L)
      ,CHFt_CONST_REALVECT(dXi)
      ,CHFt_CONST_REAL(scale) );
}
#endif  // GUARDFILL_BATHYDXDXI 

#ifndef GUARDFILL_BATHYDZDXI 
#define GUARDFILL_BATHYDZDXI 
// Prototype for Fortran procedure FILL_BATHYDZDXI ...
//
void FORTRAN_NAME( FILL_BATHYDZDXI ,fill_bathydzdxi )(
      CHFp_FRA1(dxdXi)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INT(nu)
      ,CHFp_CONST_FRA1(depth)
      ,CHFp_CONST_REALVECT(L)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REAL(scale) );

#define FORT_FILL_BATHYDZDXI FORTRAN_NAME( inlineFILL_BATHYDZDXI, inlineFILL_BATHYDZDXI)
#define FORTNT_FILL_BATHYDZDXI FORTRAN_NAME( FILL_BATHYDZDXI, fill_bathydzdxi)

inline void FORTRAN_NAME(inlineFILL_BATHYDZDXI, inlineFILL_BATHYDZDXI)(
      CHFp_FRA1(dxdXi)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INT(nu)
      ,CHFp_CONST_FRA1(depth)
      ,CHFp_CONST_REALVECT(L)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REAL(scale) )
{
 CH_TIMELEAF("FORT_FILL_BATHYDZDXI");
 FORTRAN_NAME( FILL_BATHYDZDXI ,fill_bathydzdxi )(
      CHFt_FRA1(dxdXi)
      ,CHFt_BOX(destBox)
      ,CHFt_CONST_INT(nu)
      ,CHFt_CONST_FRA1(depth)
      ,CHFt_CONST_REALVECT(L)
      ,CHFt_CONST_REALVECT(dXi)
      ,CHFt_CONST_REAL(scale) );
}
#endif  // GUARDFILL_BATHYDZDXI 

#ifndef GUARDFILL_BATHYDZDZETA 
#define GUARDFILL_BATHYDZDZETA 
// Prototype for Fortran procedure FILL_BATHYDZDZETA ...
//
void FORTRAN_NAME( FILL_BATHYDZDZETA ,fill_bathydzdzeta )(
      CHFp_FRA1(dxdZeta)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_FRA1(depth)
      ,CHFp_CONST_REALVECT(L)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REAL(scale) );

#define FORT_FILL_BATHYDZDZETA FORTRAN_NAME( inlineFILL_BATHYDZDZETA, inlineFILL_BATHYDZDZETA)
#define FORTNT_FILL_BATHYDZDZETA FORTRAN_NAME( FILL_BATHYDZDZETA, fill_bathydzdzeta)

inline void FORTRAN_NAME(inlineFILL_BATHYDZDZETA, inlineFILL_BATHYDZDZETA)(
      CHFp_FRA1(dxdZeta)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_FRA1(depth)
      ,CHFp_CONST_REALVECT(L)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REAL(scale) )
{
 CH_TIMELEAF("FORT_FILL_BATHYDZDZETA");
 FORTRAN_NAME( FILL_BATHYDZDZETA ,fill_bathydzdzeta )(
      CHFt_FRA1(dxdZeta)
      ,CHFt_BOX(destBox)
      ,CHFt_CONST_FRA1(depth)
      ,CHFt_CONST_REALVECT(L)
      ,CHFt_CONST_REALVECT(dXi)
      ,CHFt_CONST_REAL(scale) );
}
#endif  // GUARDFILL_BATHYDZDZETA 

}

#endif
