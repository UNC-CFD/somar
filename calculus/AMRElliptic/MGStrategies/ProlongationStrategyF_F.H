#ifndef _PROLONGATIONSTRATEGYF_F_H_
#define _PROLONGATIONSTRATEGYF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCONSTINTERPPS 
#define GUARDCONSTINTERPPS 
// Prototype for Fortran procedure ConstInterpPS ...
//
void FORTRAN_NAME( CONSTINTERPPS ,constinterpps )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(fineRegion)
      ,CHFp_CONST_INTVECT(m) );

#define FORT_CONSTINTERPPS FORTRAN_NAME( inlineCONSTINTERPPS, inlineCONSTINTERPPS)
#define FORTNT_CONSTINTERPPS FORTRAN_NAME( CONSTINTERPPS, constinterpps)

inline void FORTRAN_NAME(inlineCONSTINTERPPS, inlineCONSTINTERPPS)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(fineRegion)
      ,CHFp_CONST_INTVECT(m) )
{
 CH_TIMELEAF("FORT_CONSTINTERPPS");
 FORTRAN_NAME( CONSTINTERPPS ,constinterpps )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_BOX(fineRegion)
      ,CHFt_CONST_INTVECT(m) );
}
#endif  // GUARDCONSTINTERPPS 

#ifndef GUARDCONSTINTERPWITHAVGPS 
#define GUARDCONSTINTERPWITHAVGPS 
// Prototype for Fortran procedure ConstInterpWithAvgPS ...
//
void FORTRAN_NAME( CONSTINTERPWITHAVGPS ,constinterpwithavgps )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(fineRegion)
      ,CHFp_CONST_INTVECT(m)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_REAL(dxProduct)
      ,CHFp_REAL(vol)
      ,CHFp_REAL(sum) );

#define FORT_CONSTINTERPWITHAVGPS FORTRAN_NAME( inlineCONSTINTERPWITHAVGPS, inlineCONSTINTERPWITHAVGPS)
#define FORTNT_CONSTINTERPWITHAVGPS FORTRAN_NAME( CONSTINTERPWITHAVGPS, constinterpwithavgps)

inline void FORTRAN_NAME(inlineCONSTINTERPWITHAVGPS, inlineCONSTINTERPWITHAVGPS)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(fineRegion)
      ,CHFp_CONST_INTVECT(m)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_REAL(dxProduct)
      ,CHFp_REAL(vol)
      ,CHFp_REAL(sum) )
{
 CH_TIMELEAF("FORT_CONSTINTERPWITHAVGPS");
 FORTRAN_NAME( CONSTINTERPWITHAVGPS ,constinterpwithavgps )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_BOX(fineRegion)
      ,CHFt_CONST_INTVECT(m)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_CONST_REAL(dxProduct)
      ,CHFt_REAL(vol)
      ,CHFt_REAL(sum) );
}
#endif  // GUARDCONSTINTERPWITHAVGPS 

}

#endif
