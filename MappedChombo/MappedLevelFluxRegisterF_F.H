#ifndef _MAPPEDLEVELFLUXREGISTERF_F_H_
#define _MAPPEDLEVELFLUXREGISTERF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMAPPEDINCREMENTFINE 
#define GUARDMAPPEDINCREMENTFINE 
// Prototype for Fortran procedure MAPPEDINCREMENTFINE ...
//
void FORTRAN_NAME( MAPPEDINCREMENTFINE ,mappedincrementfine )(
      CHFp_CONST_FRA(fine)
      ,CHFp_FRA(cFine)
      ,CHFp_BOX(fineBox)
      ,CHFp_CONST_INTVECT(nRef)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(srcStart)
      ,CHFp_CONST_INT(destStart)
      ,CHFp_CONST_INT(ncomp) );

#define FORT_MAPPEDINCREMENTFINE FORTRAN_NAME( inlineMAPPEDINCREMENTFINE, inlineMAPPEDINCREMENTFINE)
#define FORTNT_MAPPEDINCREMENTFINE FORTRAN_NAME( MAPPEDINCREMENTFINE, mappedincrementfine)

inline void FORTRAN_NAME(inlineMAPPEDINCREMENTFINE, inlineMAPPEDINCREMENTFINE)(
      CHFp_CONST_FRA(fine)
      ,CHFp_FRA(cFine)
      ,CHFp_BOX(fineBox)
      ,CHFp_CONST_INTVECT(nRef)
      ,CHFp_CONST_REAL(scale)
      ,CHFp_CONST_INT(srcStart)
      ,CHFp_CONST_INT(destStart)
      ,CHFp_CONST_INT(ncomp) )
{
 CH_TIMELEAF("FORT_MAPPEDINCREMENTFINE");
 FORTRAN_NAME( MAPPEDINCREMENTFINE ,mappedincrementfine )(
      CHFt_CONST_FRA(fine)
      ,CHFt_FRA(cFine)
      ,CHFt_BOX(fineBox)
      ,CHFt_CONST_INTVECT(nRef)
      ,CHFt_CONST_REAL(scale)
      ,CHFt_CONST_INT(srcStart)
      ,CHFt_CONST_INT(destStart)
      ,CHFt_CONST_INT(ncomp) );
}
#endif  // GUARDMAPPEDINCREMENTFINE 

}

#endif
