#ifndef _MAPPEDQUADCFINTERPF_F_H_
#define _MAPPEDQUADCFINTERPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMAPPEDQUADINTERP 
#define GUARDMAPPEDQUADINTERP 
// Prototype for Fortran procedure mappedquadinterp ...
//
void FORTRAN_NAME( MAPPEDQUADINTERP ,mappedquadinterp )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(phistar)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(ihilo)
      ,CHFp_CONST_REAL(h)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(scomp)
      ,CHFp_CONST_INT(ecomp)
      ,CHFp_CONST_INT(nref) );

#define FORT_MAPPEDQUADINTERP FORTRAN_NAME( inlineMAPPEDQUADINTERP, inlineMAPPEDQUADINTERP)
#define FORTNT_MAPPEDQUADINTERP FORTRAN_NAME( MAPPEDQUADINTERP, mappedquadinterp)

inline void FORTRAN_NAME(inlineMAPPEDQUADINTERP, inlineMAPPEDQUADINTERP)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(phistar)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(ihilo)
      ,CHFp_CONST_REAL(h)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(scomp)
      ,CHFp_CONST_INT(ecomp)
      ,CHFp_CONST_INT(nref) )
{
 CH_TIMELEAF("FORT_MAPPEDQUADINTERP");
 FORTRAN_NAME( MAPPEDQUADINTERP ,mappedquadinterp )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(phistar)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(ihilo)
      ,CHFt_CONST_REAL(h)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(scomp)
      ,CHFt_CONST_INT(ecomp)
      ,CHFt_CONST_INT(nref) );
}
#endif  // GUARDMAPPEDQUADINTERP 

#ifndef GUARDMAPPEDPHISTAR 
#define GUARDMAPPEDPHISTAR 
// Prototype for Fortran procedure MAPPEDPHISTAR ...
//
void FORTRAN_NAME( MAPPEDPHISTAR ,mappedphistar )(
      CHFp_FRA(fPhiStar)
      ,CHFp_BOX(region)
      ,CHFp_CONST_FRA(phic)
      ,CHFp_FRA(coarslope)
      ,CHFp_FRA(coarcurva)
      ,CHFp_FRA(coarmixed)
      ,CHFp_CONST_REALVECT(dxf)
      ,CHFp_CONST_INT(ivar)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(sign)
      ,CHFp_CONST_INTVECT(nRef) );

#define FORT_MAPPEDPHISTAR FORTRAN_NAME( inlineMAPPEDPHISTAR, inlineMAPPEDPHISTAR)
#define FORTNT_MAPPEDPHISTAR FORTRAN_NAME( MAPPEDPHISTAR, mappedphistar)

inline void FORTRAN_NAME(inlineMAPPEDPHISTAR, inlineMAPPEDPHISTAR)(
      CHFp_FRA(fPhiStar)
      ,CHFp_BOX(region)
      ,CHFp_CONST_FRA(phic)
      ,CHFp_FRA(coarslope)
      ,CHFp_FRA(coarcurva)
      ,CHFp_FRA(coarmixed)
      ,CHFp_CONST_REALVECT(dxf)
      ,CHFp_CONST_INT(ivar)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(sign)
      ,CHFp_CONST_INTVECT(nRef) )
{
 CH_TIMELEAF("FORT_MAPPEDPHISTAR");
 FORTRAN_NAME( MAPPEDPHISTAR ,mappedphistar )(
      CHFt_FRA(fPhiStar)
      ,CHFt_BOX(region)
      ,CHFt_CONST_FRA(phic)
      ,CHFt_FRA(coarslope)
      ,CHFt_FRA(coarcurva)
      ,CHFt_FRA(coarmixed)
      ,CHFt_CONST_REALVECT(dxf)
      ,CHFt_CONST_INT(ivar)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(sign)
      ,CHFt_CONST_INTVECT(nRef) );
}
#endif  // GUARDMAPPEDPHISTAR 

}

#endif
