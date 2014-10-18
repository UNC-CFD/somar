#ifndef _STRESSMETRICF_F_H_
#define _STRESSMETRICF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDAPPLYVISCOUSFACTOR 
#define GUARDAPPLYVISCOUSFACTOR 
// Prototype for Fortran procedure ApplyViscousFactor ...
//
void FORTRAN_NAME( APPLYVISCOUSFACTOR ,applyviscousfactor )(
      CHFp_FRA1(stress)
      ,CHFp_CONST_FRA(physCoor)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REAL(dynVisc)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType) );

#define FORT_APPLYVISCOUSFACTOR FORTRAN_NAME( inlineAPPLYVISCOUSFACTOR, inlineAPPLYVISCOUSFACTOR)
#define FORTNT_APPLYVISCOUSFACTOR FORTRAN_NAME( APPLYVISCOUSFACTOR, applyviscousfactor)

inline void FORTRAN_NAME(inlineAPPLYVISCOUSFACTOR, inlineAPPLYVISCOUSFACTOR)(
      CHFp_FRA1(stress)
      ,CHFp_CONST_FRA(physCoor)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REAL(dynVisc)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType) )
{
 CH_TIMELEAF("FORT_APPLYVISCOUSFACTOR");
 FORTRAN_NAME( APPLYVISCOUSFACTOR ,applyviscousfactor )(
      CHFt_FRA1(stress)
      ,CHFt_CONST_FRA(physCoor)
      ,CHFt_CONST_REALVECT(dXi)
      ,CHFt_CONST_REAL(dynVisc)
      ,CHFt_BOX(destBox)
      ,CHFt_CONST_INTVECT(destBoxType) );
}
#endif  // GUARDAPPLYVISCOUSFACTOR 

}

#endif
