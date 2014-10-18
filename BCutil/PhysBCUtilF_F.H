#ifndef _PHYSBCUTILF_F_H_
#define _PHYSBCUTILF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTE_CCNSQ 
#define GUARDCOMPUTE_CCNSQ 
// Prototype for Fortran procedure COMPUTE_CCNSQ ...
//
void FORTRAN_NAME( COMPUTE_CCNSQ ,compute_ccnsq )(
      CHFp_FRA1(Nsq)
      ,CHFp_CONST_FRA1(buoyancy)
      ,CHFp_CONST_FRA(dXidz)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_BOX(destBox) );

#define FORT_COMPUTE_CCNSQ FORTRAN_NAME( inlineCOMPUTE_CCNSQ, inlineCOMPUTE_CCNSQ)
#define FORTNT_COMPUTE_CCNSQ FORTRAN_NAME( COMPUTE_CCNSQ, compute_ccnsq)

inline void FORTRAN_NAME(inlineCOMPUTE_CCNSQ, inlineCOMPUTE_CCNSQ)(
      CHFp_FRA1(Nsq)
      ,CHFp_CONST_FRA1(buoyancy)
      ,CHFp_CONST_FRA(dXidz)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_BOX(destBox) )
{
 CH_TIMELEAF("FORT_COMPUTE_CCNSQ");
 FORTRAN_NAME( COMPUTE_CCNSQ ,compute_ccnsq )(
      CHFt_FRA1(Nsq)
      ,CHFt_CONST_FRA1(buoyancy)
      ,CHFt_CONST_FRA(dXidz)
      ,CHFt_CONST_REALVECT(dXi)
      ,CHFt_BOX(destBox) );
}
#endif  // GUARDCOMPUTE_CCNSQ 

}

#endif
