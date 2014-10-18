#ifndef _SOLITARYWAVEBCUTILF_F_H_
#define _SOLITARYWAVEBCUTILF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSOLITARYWAVE_SETVELIC 
#define GUARDSOLITARYWAVE_SETVELIC 
// Prototype for Fortran procedure SOLITARYWAVE_SETVELIC ...
//
void FORTRAN_NAME( SOLITARYWAVE_SETVELIC ,solitarywave_setvelic )(
      CHFp_FRA(vel)
      ,CHFp_CONST_INT(velDir)
      ,CHFp_CONST_FRA1(CCA)
      ,CHFp_CONST_FRA1(FCA)
      ,CHFp_CONST_FRA1(CCphi)
      ,CHFp_BOX(valid)
      ,CHFp_BOX(horizCCBox)
      ,CHFp_BOX(vertCCBox)
      ,CHFp_CONST_REALVECT(physDx) );

#define FORT_SOLITARYWAVE_SETVELIC FORTRAN_NAME( inlineSOLITARYWAVE_SETVELIC, inlineSOLITARYWAVE_SETVELIC)
#define FORTNT_SOLITARYWAVE_SETVELIC FORTRAN_NAME( SOLITARYWAVE_SETVELIC, solitarywave_setvelic)

inline void FORTRAN_NAME(inlineSOLITARYWAVE_SETVELIC, inlineSOLITARYWAVE_SETVELIC)(
      CHFp_FRA(vel)
      ,CHFp_CONST_INT(velDir)
      ,CHFp_CONST_FRA1(CCA)
      ,CHFp_CONST_FRA1(FCA)
      ,CHFp_CONST_FRA1(CCphi)
      ,CHFp_BOX(valid)
      ,CHFp_BOX(horizCCBox)
      ,CHFp_BOX(vertCCBox)
      ,CHFp_CONST_REALVECT(physDx) )
{
 CH_TIMELEAF("FORT_SOLITARYWAVE_SETVELIC");
 FORTRAN_NAME( SOLITARYWAVE_SETVELIC ,solitarywave_setvelic )(
      CHFt_FRA(vel)
      ,CHFt_CONST_INT(velDir)
      ,CHFt_CONST_FRA1(CCA)
      ,CHFt_CONST_FRA1(FCA)
      ,CHFt_CONST_FRA1(CCphi)
      ,CHFt_BOX(valid)
      ,CHFt_BOX(horizCCBox)
      ,CHFt_BOX(vertCCBox)
      ,CHFt_CONST_REALVECT(physDx) );
}
#endif  // GUARDSOLITARYWAVE_SETVELIC 

#ifndef GUARDSOLITARYWAVE_SETSCALARIC 
#define GUARDSOLITARYWAVE_SETSCALARIC 
// Prototype for Fortran procedure SOLITARYWAVE_SETSCALARIC ...
//
void FORTRAN_NAME( SOLITARYWAVE_SETSCALARIC ,solitarywave_setscalaric )(
      CHFp_FRA1(CCB)
      ,CHFp_CONST_FRA1(CCBbg)
      ,CHFp_CONST_FRA1(CCA)
      ,CHFp_CONST_FRA1(CCphi)
      ,CHFp_CONST_FRA1(CCNsq)
      ,CHFp_CONST_REAL(c0)
      ,CHFp_BOX(valid)
      ,CHFp_BOX(horizCCBox)
      ,CHFp_BOX(vertCCBox) );

#define FORT_SOLITARYWAVE_SETSCALARIC FORTRAN_NAME( inlineSOLITARYWAVE_SETSCALARIC, inlineSOLITARYWAVE_SETSCALARIC)
#define FORTNT_SOLITARYWAVE_SETSCALARIC FORTRAN_NAME( SOLITARYWAVE_SETSCALARIC, solitarywave_setscalaric)

inline void FORTRAN_NAME(inlineSOLITARYWAVE_SETSCALARIC, inlineSOLITARYWAVE_SETSCALARIC)(
      CHFp_FRA1(CCB)
      ,CHFp_CONST_FRA1(CCBbg)
      ,CHFp_CONST_FRA1(CCA)
      ,CHFp_CONST_FRA1(CCphi)
      ,CHFp_CONST_FRA1(CCNsq)
      ,CHFp_CONST_REAL(c0)
      ,CHFp_BOX(valid)
      ,CHFp_BOX(horizCCBox)
      ,CHFp_BOX(vertCCBox) )
{
 CH_TIMELEAF("FORT_SOLITARYWAVE_SETSCALARIC");
 FORTRAN_NAME( SOLITARYWAVE_SETSCALARIC ,solitarywave_setscalaric )(
      CHFt_FRA1(CCB)
      ,CHFt_CONST_FRA1(CCBbg)
      ,CHFt_CONST_FRA1(CCA)
      ,CHFt_CONST_FRA1(CCphi)
      ,CHFt_CONST_FRA1(CCNsq)
      ,CHFt_CONST_REAL(c0)
      ,CHFt_BOX(valid)
      ,CHFt_BOX(horizCCBox)
      ,CHFt_BOX(vertCCBox) );
}
#endif  // GUARDSOLITARYWAVE_SETSCALARIC 

}

#endif
