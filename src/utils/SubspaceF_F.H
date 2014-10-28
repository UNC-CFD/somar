#ifndef _SUBSPACEF_F_H_
#define _SUBSPACEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDUNMAPPEDVERTINTEGRAL 
#define GUARDUNMAPPEDVERTINTEGRAL 
// Prototype for Fortran procedure UnmappedVertIntegral ...
//
void FORTRAN_NAME( UNMAPPEDVERTINTEGRAL ,unmappedvertintegral )(
      CHFp_FRA1(levelSum)
      ,CHFp_CONST_FRA1(func)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dz) );

#define FORT_UNMAPPEDVERTINTEGRAL FORTRAN_NAME( inlineUNMAPPEDVERTINTEGRAL, inlineUNMAPPEDVERTINTEGRAL)
#define FORTNT_UNMAPPEDVERTINTEGRAL FORTRAN_NAME( UNMAPPEDVERTINTEGRAL, unmappedvertintegral)

inline void FORTRAN_NAME(inlineUNMAPPEDVERTINTEGRAL, inlineUNMAPPEDVERTINTEGRAL)(
      CHFp_FRA1(levelSum)
      ,CHFp_CONST_FRA1(func)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dz) )
{
 CH_TIMELEAF("FORT_UNMAPPEDVERTINTEGRAL");
 FORTRAN_NAME( UNMAPPEDVERTINTEGRAL ,unmappedvertintegral )(
      CHFt_FRA1(levelSum)
      ,CHFt_CONST_FRA1(func)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dz) );
}
#endif  // GUARDUNMAPPEDVERTINTEGRAL 

#ifndef GUARDADDEXTRUSION 
#define GUARDADDEXTRUSION 
// Prototype for Fortran procedure AddExtrusion ...
//
void FORTRAN_NAME( ADDEXTRUSION ,addextrusion )(
      CHFp_FRA1(dest)
      ,CHFp_CONST_FRA1(src)
      ,CHFp_BOX(region) );

#define FORT_ADDEXTRUSION FORTRAN_NAME( inlineADDEXTRUSION, inlineADDEXTRUSION)
#define FORTNT_ADDEXTRUSION FORTRAN_NAME( ADDEXTRUSION, addextrusion)

inline void FORTRAN_NAME(inlineADDEXTRUSION, inlineADDEXTRUSION)(
      CHFp_FRA1(dest)
      ,CHFp_CONST_FRA1(src)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_ADDEXTRUSION");
 FORTRAN_NAME( ADDEXTRUSION ,addextrusion )(
      CHFt_FRA1(dest)
      ,CHFt_CONST_FRA1(src)
      ,CHFt_BOX(region) );
}
#endif  // GUARDADDEXTRUSION 

}

#endif
