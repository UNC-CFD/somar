#ifndef _MAPPEDAMRPOISSONOPORTHOF_F_H_
#define _MAPPEDAMRPOISSONOPORTHOF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMAPPEDGETFLUXORTHO 
#define GUARDMAPPEDGETFLUXORTHO 
// Prototype for Fortran procedure MAPPEDGETFLUXORTHO ...
//
void FORTRAN_NAME( MAPPEDGETFLUXORTHO ,mappedgetfluxortho )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(Jgaa)
      ,CHFp_BOX(FCRegion)
      ,CHFp_CONST_REAL(beta_dx)
      ,CHFp_CONST_INT(adir) );

#define FORT_MAPPEDGETFLUXORTHO FORTRAN_NAME( inlineMAPPEDGETFLUXORTHO, inlineMAPPEDGETFLUXORTHO)
#define FORTNT_MAPPEDGETFLUXORTHO FORTRAN_NAME( MAPPEDGETFLUXORTHO, mappedgetfluxortho)

inline void FORTRAN_NAME(inlineMAPPEDGETFLUXORTHO, inlineMAPPEDGETFLUXORTHO)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(Jgaa)
      ,CHFp_BOX(FCRegion)
      ,CHFp_CONST_REAL(beta_dx)
      ,CHFp_CONST_INT(adir) )
{
 CH_TIMELEAF("FORT_MAPPEDGETFLUXORTHO");
 FORTRAN_NAME( MAPPEDGETFLUXORTHO ,mappedgetfluxortho )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA1(Jgaa)
      ,CHFt_BOX(FCRegion)
      ,CHFt_CONST_REAL(beta_dx)
      ,CHFt_CONST_INT(adir) );
}
#endif  // GUARDMAPPEDGETFLUXORTHO 

}

#endif
