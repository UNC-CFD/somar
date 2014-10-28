#ifndef _NODEINTERPF_F_H_
#define _NODEINTERPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDLINEARNODEINTERP 
#define GUARDLINEARNODEINTERP 
// Prototype for Fortran procedure LinearNodeInterp ...
//
void FORTRAN_NAME( LINEARNODEINTERP ,linearnodeinterp )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(crse)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INTVECT(ref)
      ,CHFp_BOX(refBox) );

#define FORT_LINEARNODEINTERP FORTRAN_NAME( inlineLINEARNODEINTERP, inlineLINEARNODEINTERP)
#define FORTNT_LINEARNODEINTERP FORTRAN_NAME( LINEARNODEINTERP, linearnodeinterp)

inline void FORTRAN_NAME(inlineLINEARNODEINTERP, inlineLINEARNODEINTERP)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(crse)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INTVECT(ref)
      ,CHFp_BOX(refBox) )
{
 CH_TIMELEAF("FORT_LINEARNODEINTERP");
 FORTRAN_NAME( LINEARNODEINTERP ,linearnodeinterp )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(crse)
      ,CHFt_BOX(crseBox)
      ,CHFt_CONST_INTVECT(ref)
      ,CHFt_BOX(refBox) );
}
#endif  // GUARDLINEARNODEINTERP 

}

#endif
