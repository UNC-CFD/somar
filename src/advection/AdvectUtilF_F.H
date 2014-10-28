#ifndef _ADVECTUTILF_F_H_
#define _ADVECTUTILF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDRIEMANNSOLVER 
#define GUARDRIEMANNSOLVER 
// Prototype for Fortran procedure RiemannSolver ...
//
void FORTRAN_NAME( RIEMANNSOLVER ,riemannsolver )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_FRA1(advVel)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_RIEMANNSOLVER FORTRAN_NAME( inlineRIEMANNSOLVER, inlineRIEMANNSOLVER)
#define FORTNT_RIEMANNSOLVER FORTRAN_NAME( RIEMANNSOLVER, riemannsolver)

inline void FORTRAN_NAME(inlineRIEMANNSOLVER, inlineRIEMANNSOLVER)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_FRA1(advVel)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RIEMANNSOLVER");
 FORTRAN_NAME( RIEMANNSOLVER ,riemannsolver )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_FRA1(advVel)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRIEMANNSOLVER 

#ifndef GUARDAVGSTATES 
#define GUARDAVGSTATES 
// Prototype for Fortran procedure AvgStates ...
//
void FORTRAN_NAME( AVGSTATES ,avgstates )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_AVGSTATES FORTRAN_NAME( inlineAVGSTATES, inlineAVGSTATES)
#define FORTNT_AVGSTATES FORTRAN_NAME( AVGSTATES, avgstates)

inline void FORTRAN_NAME(inlineAVGSTATES, inlineAVGSTATES)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_AVGSTATES");
 FORTRAN_NAME( AVGSTATES ,avgstates )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDAVGSTATES 

#ifndef GUARDQUASILINEARUPDATEF 
#define GUARDQUASILINEARUPDATEF 
// Prototype for Fortran procedure QuasilinearUpdateF ...
//
void FORTRAN_NAME( QUASILINEARUPDATEF ,quasilinearupdatef )(
      CHFp_FRA(AdWdx)
      ,CHFp_CONST_FRA(WHalf)
      ,CHFp_CONST_FRA1(cellVel)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_QUASILINEARUPDATEF FORTRAN_NAME( inlineQUASILINEARUPDATEF, inlineQUASILINEARUPDATEF)
#define FORTNT_QUASILINEARUPDATEF FORTRAN_NAME( QUASILINEARUPDATEF, quasilinearupdatef)

inline void FORTRAN_NAME(inlineQUASILINEARUPDATEF, inlineQUASILINEARUPDATEF)(
      CHFp_FRA(AdWdx)
      ,CHFp_CONST_FRA(WHalf)
      ,CHFp_CONST_FRA1(cellVel)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_QUASILINEARUPDATEF");
 FORTRAN_NAME( QUASILINEARUPDATEF ,quasilinearupdatef )(
      CHFt_FRA(AdWdx)
      ,CHFt_CONST_FRA(WHalf)
      ,CHFt_CONST_FRA1(cellVel)
      ,CHFt_CONST_REAL(scale)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDQUASILINEARUPDATEF 

}

#endif
