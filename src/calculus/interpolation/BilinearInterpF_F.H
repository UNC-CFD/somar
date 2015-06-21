#ifndef _BILINEARINTERPF_F_H_
#define _BILINEARINTERPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDBILINEARINTERP2DF 
#define GUARDBILINEARINTERP2DF 
// Prototype for Fortran procedure BilinearInterp2DF ...
//
void FORTRAN_NAME( BILINEARINTERP2DF ,bilinearinterp2df )(
      CHFp_FRA(fInterp)
      ,CHFp_CONST_FRA1(xInterp)
      ,CHFp_CONST_FRA1(yInterp)
      ,CHFp_BOX(interpBox)
      ,CHFp_CONST_INT(xdir)
      ,CHFp_CONST_INT(ydir)
      ,CHFp_CONST_VR(x)
      ,CHFp_CONST_VR(y)
      ,CHFp_CONST_FRA(f));

#define FORT_BILINEARINTERP2DF FORTRAN_NAME( inlineBILINEARINTERP2DF, inlineBILINEARINTERP2DF)
#define FORTNT_BILINEARINTERP2DF FORTRAN_NAME( BILINEARINTERP2DF, bilinearinterp2df)

inline void FORTRAN_NAME(inlineBILINEARINTERP2DF, inlineBILINEARINTERP2DF)(
      CHFp_FRA(fInterp)
      ,CHFp_CONST_FRA1(xInterp)
      ,CHFp_CONST_FRA1(yInterp)
      ,CHFp_BOX(interpBox)
      ,CHFp_CONST_INT(xdir)
      ,CHFp_CONST_INT(ydir)
      ,CHFp_CONST_VR(x)
      ,CHFp_CONST_VR(y)
      ,CHFp_CONST_FRA(f))
{
 CH_TIMELEAF("FORT_BILINEARINTERP2DF");
 FORTRAN_NAME( BILINEARINTERP2DF ,bilinearinterp2df )(
      CHFt_FRA(fInterp)
      ,CHFt_CONST_FRA1(xInterp)
      ,CHFt_CONST_FRA1(yInterp)
      ,CHFt_BOX(interpBox)
      ,CHFt_CONST_INT(xdir)
      ,CHFt_CONST_INT(ydir)
      ,CHFt_CONST_VR(x)
      ,CHFt_CONST_VR(y)
      ,CHFt_CONST_FRA(f) );
}
#endif  // GUARDBILINEARINTERP2DF 

}

#endif
