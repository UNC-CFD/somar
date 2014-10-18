#ifndef _HERMITEINTERPF_F_H_
#define _HERMITEINTERPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDHERMITEINTERP2DF 
#define GUARDHERMITEINTERP2DF 
// Prototype for Fortran procedure HermiteInterp2DF ...
//
void FORTRAN_NAME( HERMITEINTERP2DF ,hermiteinterp2df )(
      CHFp_FRA(fInterp)
      ,CHFp_CONST_FRA1(xInterp)
      ,CHFp_CONST_FRA1(yInterp)
      ,CHFp_BOX(interpBox)
      ,CHFp_CONST_INT(xdir)
      ,CHFp_CONST_INT(ydir)
      ,CHFp_CONST_VR(x)
      ,CHFp_CONST_VR(y)
      ,CHFp_CONST_FRA(f)
      ,CHFp_CONST_FRA(dfdx)
      ,CHFp_CONST_FRA(dfdy) );

#define FORT_HERMITEINTERP2DF FORTRAN_NAME( inlineHERMITEINTERP2DF, inlineHERMITEINTERP2DF)
#define FORTNT_HERMITEINTERP2DF FORTRAN_NAME( HERMITEINTERP2DF, hermiteinterp2df)

inline void FORTRAN_NAME(inlineHERMITEINTERP2DF, inlineHERMITEINTERP2DF)(
      CHFp_FRA(fInterp)
      ,CHFp_CONST_FRA1(xInterp)
      ,CHFp_CONST_FRA1(yInterp)
      ,CHFp_BOX(interpBox)
      ,CHFp_CONST_INT(xdir)
      ,CHFp_CONST_INT(ydir)
      ,CHFp_CONST_VR(x)
      ,CHFp_CONST_VR(y)
      ,CHFp_CONST_FRA(f)
      ,CHFp_CONST_FRA(dfdx)
      ,CHFp_CONST_FRA(dfdy) )
{
 CH_TIMELEAF("FORT_HERMITEINTERP2DF");
 FORTRAN_NAME( HERMITEINTERP2DF ,hermiteinterp2df )(
      CHFt_FRA(fInterp)
      ,CHFt_CONST_FRA1(xInterp)
      ,CHFt_CONST_FRA1(yInterp)
      ,CHFt_BOX(interpBox)
      ,CHFt_CONST_INT(xdir)
      ,CHFt_CONST_INT(ydir)
      ,CHFt_CONST_VR(x)
      ,CHFt_CONST_VR(y)
      ,CHFt_CONST_FRA(f)
      ,CHFt_CONST_FRA(dfdx)
      ,CHFt_CONST_FRA(dfdy) );
}
#endif  // GUARDHERMITEINTERP2DF 

}

#endif
