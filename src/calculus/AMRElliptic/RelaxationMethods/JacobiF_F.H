#ifndef _JACOBIF_F_H_
#define _JACOBIF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDJACOBIITER 
#define GUARDJACOBIITER 
// Prototype for Fortran procedure JACOBIITER ...
//
void FORTRAN_NAME( JACOBIITER ,jacobiiter )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(res)
      ,CHFp_CONST_FRA1(lapDiag)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) );

#define FORT_JACOBIITER FORTRAN_NAME( inlineJACOBIITER, inlineJACOBIITER)
#define FORTNT_JACOBIITER FORTRAN_NAME( JACOBIITER, jacobiiter)

inline void FORTRAN_NAME(inlineJACOBIITER, inlineJACOBIITER)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(res)
      ,CHFp_CONST_FRA1(lapDiag)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) )
{
 CH_TIMELEAF("FORT_JACOBIITER");
 FORTRAN_NAME( JACOBIITER ,jacobiiter )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(res)
      ,CHFt_CONST_FRA1(lapDiag)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta) );
}
#endif  // GUARDJACOBIITER 

}

#endif
