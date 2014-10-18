#ifndef _TWISTEDMAPF_F_H_
#define _TWISTEDMAPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDTWISTED0_FILL_PHYSCOOR 
#define GUARDTWISTED0_FILL_PHYSCOOR 
// Prototype for Fortran procedure TWISTED0_FILL_PHYSCOOR ...
//
void FORTRAN_NAME( TWISTED0_FILL_PHYSCOOR ,twisted0_fill_physcoor )(
      CHFp_FRA1(dest)
      ,CHFp_CONST_INT(mu)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REAL(pert)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType) );

#define FORT_TWISTED0_FILL_PHYSCOOR FORTRAN_NAME( inlineTWISTED0_FILL_PHYSCOOR, inlineTWISTED0_FILL_PHYSCOOR)
#define FORTNT_TWISTED0_FILL_PHYSCOOR FORTRAN_NAME( TWISTED0_FILL_PHYSCOOR, twisted0_fill_physcoor)

inline void FORTRAN_NAME(inlineTWISTED0_FILL_PHYSCOOR, inlineTWISTED0_FILL_PHYSCOOR)(
      CHFp_FRA1(dest)
      ,CHFp_CONST_INT(mu)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REAL(pert)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType) )
{
 CH_TIMELEAF("FORT_TWISTED0_FILL_PHYSCOOR");
 FORTRAN_NAME( TWISTED0_FILL_PHYSCOOR ,twisted0_fill_physcoor )(
      CHFt_FRA1(dest)
      ,CHFt_CONST_INT(mu)
      ,CHFt_CONST_REALVECT(dXi)
      ,CHFt_CONST_REAL(pert)
      ,CHFt_BOX(destBox)
      ,CHFt_CONST_INTVECT(destBoxType) );
}
#endif  // GUARDTWISTED0_FILL_PHYSCOOR 

#ifndef GUARDTWISTED0_FILL_PHYSCOOR_ALL_COMPS 
#define GUARDTWISTED0_FILL_PHYSCOOR_ALL_COMPS 
// Prototype for Fortran procedure TWISTED0_FILL_PHYSCOOR_ALL_COMPS ...
//
void FORTRAN_NAME( TWISTED0_FILL_PHYSCOOR_ALL_COMPS ,twisted0_fill_physcoor_all_comps )(
      CHFp_FRA(dest)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REALVECT(pert)
      ,CHFp_CONST_REALVECT(mult)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType) );

#define FORT_TWISTED0_FILL_PHYSCOOR_ALL_COMPS FORTRAN_NAME( inlineTWISTED0_FILL_PHYSCOOR_ALL_COMPS, inlineTWISTED0_FILL_PHYSCOOR_ALL_COMPS)
#define FORTNT_TWISTED0_FILL_PHYSCOOR_ALL_COMPS FORTRAN_NAME( TWISTED0_FILL_PHYSCOOR_ALL_COMPS, twisted0_fill_physcoor_all_comps)

inline void FORTRAN_NAME(inlineTWISTED0_FILL_PHYSCOOR_ALL_COMPS, inlineTWISTED0_FILL_PHYSCOOR_ALL_COMPS)(
      CHFp_FRA(dest)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REALVECT(pert)
      ,CHFp_CONST_REALVECT(mult)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType) )
{
 CH_TIMELEAF("FORT_TWISTED0_FILL_PHYSCOOR_ALL_COMPS");
 FORTRAN_NAME( TWISTED0_FILL_PHYSCOOR_ALL_COMPS ,twisted0_fill_physcoor_all_comps )(
      CHFt_FRA(dest)
      ,CHFt_CONST_REALVECT(dXi)
      ,CHFt_CONST_REALVECT(pert)
      ,CHFt_CONST_REALVECT(mult)
      ,CHFt_BOX(destBox)
      ,CHFt_CONST_INTVECT(destBoxType) );
}
#endif  // GUARDTWISTED0_FILL_PHYSCOOR_ALL_COMPS 

#ifndef GUARDTWISTED0_FILL_DXDXI 
#define GUARDTWISTED0_FILL_DXDXI 
// Prototype for Fortran procedure TWISTED0_FILL_DXDXI ...
//
void FORTRAN_NAME( TWISTED0_FILL_DXDXI ,twisted0_fill_dxdxi )(
      CHFp_FRA1(dest)
      ,CHFp_CONST_INT(mu)
      ,CHFp_CONST_INT(nu)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REAL(pert)
      ,CHFp_CONST_REAL(mult)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType) );

#define FORT_TWISTED0_FILL_DXDXI FORTRAN_NAME( inlineTWISTED0_FILL_DXDXI, inlineTWISTED0_FILL_DXDXI)
#define FORTNT_TWISTED0_FILL_DXDXI FORTRAN_NAME( TWISTED0_FILL_DXDXI, twisted0_fill_dxdxi)

inline void FORTRAN_NAME(inlineTWISTED0_FILL_DXDXI, inlineTWISTED0_FILL_DXDXI)(
      CHFp_FRA1(dest)
      ,CHFp_CONST_INT(mu)
      ,CHFp_CONST_INT(nu)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REAL(pert)
      ,CHFp_CONST_REAL(mult)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType) )
{
 CH_TIMELEAF("FORT_TWISTED0_FILL_DXDXI");
 FORTRAN_NAME( TWISTED0_FILL_DXDXI ,twisted0_fill_dxdxi )(
      CHFt_FRA1(dest)
      ,CHFt_CONST_INT(mu)
      ,CHFt_CONST_INT(nu)
      ,CHFt_CONST_REALVECT(dXi)
      ,CHFt_CONST_REAL(pert)
      ,CHFt_CONST_REAL(mult)
      ,CHFt_BOX(destBox)
      ,CHFt_CONST_INTVECT(destBoxType) );
}
#endif  // GUARDTWISTED0_FILL_DXDXI 

#ifndef GUARDTWISTED0_FILL_J 
#define GUARDTWISTED0_FILL_J 
// Prototype for Fortran procedure TWISTED0_FILL_J ...
//
void FORTRAN_NAME( TWISTED0_FILL_J ,twisted0_fill_j )(
      CHFp_FRA1(dest)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REALVECT(pert)
      ,CHFp_CONST_REAL(mult)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType) );

#define FORT_TWISTED0_FILL_J FORTRAN_NAME( inlineTWISTED0_FILL_J, inlineTWISTED0_FILL_J)
#define FORTNT_TWISTED0_FILL_J FORTRAN_NAME( TWISTED0_FILL_J, twisted0_fill_j)

inline void FORTRAN_NAME(inlineTWISTED0_FILL_J, inlineTWISTED0_FILL_J)(
      CHFp_FRA1(dest)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REALVECT(pert)
      ,CHFp_CONST_REAL(mult)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType) )
{
 CH_TIMELEAF("FORT_TWISTED0_FILL_J");
 FORTRAN_NAME( TWISTED0_FILL_J ,twisted0_fill_j )(
      CHFt_FRA1(dest)
      ,CHFt_CONST_REALVECT(dXi)
      ,CHFt_CONST_REALVECT(pert)
      ,CHFt_CONST_REAL(mult)
      ,CHFt_BOX(destBox)
      ,CHFt_CONST_INTVECT(destBoxType) );
}
#endif  // GUARDTWISTED0_FILL_J 

#ifndef GUARDTWISTED1_FILL_PHYSCOOR 
#define GUARDTWISTED1_FILL_PHYSCOOR 
// Prototype for Fortran procedure TWISTED1_FILL_PHYSCOOR ...
//
void FORTRAN_NAME( TWISTED1_FILL_PHYSCOOR ,twisted1_fill_physcoor )(
      CHFp_FRA1(dest)
      ,CHFp_CONST_INT(mu)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REALVECT(domLen)
      ,CHFp_CONST_REAL(pert)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType) );

#define FORT_TWISTED1_FILL_PHYSCOOR FORTRAN_NAME( inlineTWISTED1_FILL_PHYSCOOR, inlineTWISTED1_FILL_PHYSCOOR)
#define FORTNT_TWISTED1_FILL_PHYSCOOR FORTRAN_NAME( TWISTED1_FILL_PHYSCOOR, twisted1_fill_physcoor)

inline void FORTRAN_NAME(inlineTWISTED1_FILL_PHYSCOOR, inlineTWISTED1_FILL_PHYSCOOR)(
      CHFp_FRA1(dest)
      ,CHFp_CONST_INT(mu)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_CONST_REALVECT(domLen)
      ,CHFp_CONST_REAL(pert)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType) )
{
 CH_TIMELEAF("FORT_TWISTED1_FILL_PHYSCOOR");
 FORTRAN_NAME( TWISTED1_FILL_PHYSCOOR ,twisted1_fill_physcoor )(
      CHFt_FRA1(dest)
      ,CHFt_CONST_INT(mu)
      ,CHFt_CONST_REALVECT(dXi)
      ,CHFt_CONST_REALVECT(domLen)
      ,CHFt_CONST_REAL(pert)
      ,CHFt_BOX(destBox)
      ,CHFt_CONST_INTVECT(destBoxType) );
}
#endif  // GUARDTWISTED1_FILL_PHYSCOOR 

}

#endif
