#ifndef _MAPPEDAMRPOISSONOPF_F_H_
#define _MAPPEDAMRPOISSONOPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSUBTRACTOP 
#define GUARDSUBTRACTOP 
// Prototype for Fortran procedure SUBTRACTOP ...
//
void FORTRAN_NAME( SUBTRACTOP ,subtractop )(
      CHFp_FRA(subres)
      ,CHFp_CONST_FRA(a1)
      ,CHFp_CONST_FRA(a2)
      ,CHFp_BOX(region) );

#define FORT_SUBTRACTOP FORTRAN_NAME( inlineSUBTRACTOP, inlineSUBTRACTOP)
#define FORTNT_SUBTRACTOP FORTRAN_NAME( SUBTRACTOP, subtractop)

inline void FORTRAN_NAME(inlineSUBTRACTOP, inlineSUBTRACTOP)(
      CHFp_FRA(subres)
      ,CHFp_CONST_FRA(a1)
      ,CHFp_CONST_FRA(a2)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_SUBTRACTOP");
 FORTRAN_NAME( SUBTRACTOP ,subtractop )(
      CHFt_FRA(subres)
      ,CHFt_CONST_FRA(a1)
      ,CHFt_CONST_FRA(a2)
      ,CHFt_BOX(region) );
}
#endif  // GUARDSUBTRACTOP 

#ifndef GUARDAXBYIP 
#define GUARDAXBYIP 
// Prototype for Fortran procedure AXBYIP ...
//
void FORTRAN_NAME( AXBYIP ,axbyip )(
      CHFp_FRA(lhs)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(region) );

#define FORT_AXBYIP FORTRAN_NAME( inlineAXBYIP, inlineAXBYIP)
#define FORTNT_AXBYIP FORTRAN_NAME( AXBYIP, axbyip)

inline void FORTRAN_NAME(inlineAXBYIP, inlineAXBYIP)(
      CHFp_FRA(lhs)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_AXBYIP");
 FORTRAN_NAME( AXBYIP ,axbyip )(
      CHFt_FRA(lhs)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_BOX(region) );
}
#endif  // GUARDAXBYIP 

#ifndef GUARDFILLMAPPEDLAPDIAG1D 
#define GUARDFILLMAPPEDLAPDIAG1D 
// Prototype for Fortran procedure FILLMAPPEDLAPDIAG1D ...
//
void FORTRAN_NAME( FILLMAPPEDLAPDIAG1D ,fillmappedlapdiag1d )(
      CHFp_FRA1(lapDiag)
      ,CHFp_CONST_FRA1(Jgxx)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_FILLMAPPEDLAPDIAG1D FORTRAN_NAME( inlineFILLMAPPEDLAPDIAG1D, inlineFILLMAPPEDLAPDIAG1D)
#define FORTNT_FILLMAPPEDLAPDIAG1D FORTRAN_NAME( FILLMAPPEDLAPDIAG1D, fillmappedlapdiag1d)

inline void FORTRAN_NAME(inlineFILLMAPPEDLAPDIAG1D, inlineFILLMAPPEDLAPDIAG1D)(
      CHFp_FRA1(lapDiag)
      ,CHFp_CONST_FRA1(Jgxx)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_FILLMAPPEDLAPDIAG1D");
 FORTRAN_NAME( FILLMAPPEDLAPDIAG1D ,fillmappedlapdiag1d )(
      CHFt_FRA1(lapDiag)
      ,CHFt_CONST_FRA1(Jgxx)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDFILLMAPPEDLAPDIAG1D 

#ifndef GUARDFILLMAPPEDLAPDIAG2D 
#define GUARDFILLMAPPEDLAPDIAG2D 
// Prototype for Fortran procedure FILLMAPPEDLAPDIAG2D ...
//
void FORTRAN_NAME( FILLMAPPEDLAPDIAG2D ,fillmappedlapdiag2d )(
      CHFp_FRA1(lapDiag)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx) );

#define FORT_FILLMAPPEDLAPDIAG2D FORTRAN_NAME( inlineFILLMAPPEDLAPDIAG2D, inlineFILLMAPPEDLAPDIAG2D)
#define FORTNT_FILLMAPPEDLAPDIAG2D FORTRAN_NAME( FILLMAPPEDLAPDIAG2D, fillmappedlapdiag2d)

inline void FORTRAN_NAME(inlineFILLMAPPEDLAPDIAG2D, inlineFILLMAPPEDLAPDIAG2D)(
      CHFp_FRA1(lapDiag)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx) )
{
 CH_TIMELEAF("FORT_FILLMAPPEDLAPDIAG2D");
 FORTRAN_NAME( FILLMAPPEDLAPDIAG2D ,fillmappedlapdiag2d )(
      CHFt_FRA1(lapDiag)
      ,CHFt_CONST_FRA(Jg0)
      ,CHFt_CONST_FRA(Jg1)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx) );
}
#endif  // GUARDFILLMAPPEDLAPDIAG2D 

#ifndef GUARDFILLMAPPEDLAPDIAG3D 
#define GUARDFILLMAPPEDLAPDIAG3D 
// Prototype for Fortran procedure FILLMAPPEDLAPDIAG3D ...
//
void FORTRAN_NAME( FILLMAPPEDLAPDIAG3D ,fillmappedlapdiag3d )(
      CHFp_FRA1(lapDiag)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA(Jg2)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx) );

#define FORT_FILLMAPPEDLAPDIAG3D FORTRAN_NAME( inlineFILLMAPPEDLAPDIAG3D, inlineFILLMAPPEDLAPDIAG3D)
#define FORTNT_FILLMAPPEDLAPDIAG3D FORTRAN_NAME( FILLMAPPEDLAPDIAG3D, fillmappedlapdiag3d)

inline void FORTRAN_NAME(inlineFILLMAPPEDLAPDIAG3D, inlineFILLMAPPEDLAPDIAG3D)(
      CHFp_FRA1(lapDiag)
      ,CHFp_CONST_FRA(Jg0)
      ,CHFp_CONST_FRA(Jg1)
      ,CHFp_CONST_FRA(Jg2)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx) )
{
 CH_TIMELEAF("FORT_FILLMAPPEDLAPDIAG3D");
 FORTRAN_NAME( FILLMAPPEDLAPDIAG3D ,fillmappedlapdiag3d )(
      CHFt_FRA1(lapDiag)
      ,CHFt_CONST_FRA(Jg0)
      ,CHFt_CONST_FRA(Jg1)
      ,CHFt_CONST_FRA(Jg2)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx) );
}
#endif  // GUARDFILLMAPPEDLAPDIAG3D 

#ifndef GUARDDIAGPRECOND 
#define GUARDDIAGPRECOND 
// Prototype for Fortran procedure DIAGPRECOND ...
//
void FORTRAN_NAME( DIAGPRECOND ,diagprecond )(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA1(lapDiag)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) );

#define FORT_DIAGPRECOND FORTRAN_NAME( inlineDIAGPRECOND, inlineDIAGPRECOND)
#define FORTNT_DIAGPRECOND FORTRAN_NAME( DIAGPRECOND, diagprecond)

inline void FORTRAN_NAME(inlineDIAGPRECOND, inlineDIAGPRECOND)(
      CHFp_FRA(phi)
      ,CHFp_CONST_FRA(rhs)
      ,CHFp_CONST_FRA1(lapDiag)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta) )
{
 CH_TIMELEAF("FORT_DIAGPRECOND");
 FORTRAN_NAME( DIAGPRECOND ,diagprecond )(
      CHFt_FRA(phi)
      ,CHFt_CONST_FRA(rhs)
      ,CHFt_CONST_FRA1(lapDiag)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta) );
}
#endif  // GUARDDIAGPRECOND 

#ifndef GUARDMAPPEDGETFLUX 
#define GUARDMAPPEDGETFLUX 
// Prototype for Fortran procedure MAPPEDGETFLUX ...
//
void FORTRAN_NAME( MAPPEDGETFLUX ,mappedgetflux )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(extrap)
      ,CHFp_CONST_FRA(Jga)
      ,CHFp_BOX(FCRegion)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(adir) );

#define FORT_MAPPEDGETFLUX FORTRAN_NAME( inlineMAPPEDGETFLUX, inlineMAPPEDGETFLUX)
#define FORTNT_MAPPEDGETFLUX FORTRAN_NAME( MAPPEDGETFLUX, mappedgetflux)

inline void FORTRAN_NAME(inlineMAPPEDGETFLUX, inlineMAPPEDGETFLUX)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(extrap)
      ,CHFp_CONST_FRA(Jga)
      ,CHFp_BOX(FCRegion)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(adir) )
{
 CH_TIMELEAF("FORT_MAPPEDGETFLUX");
 FORTRAN_NAME( MAPPEDGETFLUX ,mappedgetflux )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA(extrap)
      ,CHFt_CONST_FRA(Jga)
      ,CHFt_BOX(FCRegion)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(adir) );
}
#endif  // GUARDMAPPEDGETFLUX 

#ifndef GUARDMAPPEDGETFLATFLUX2D 
#define GUARDMAPPEDGETFLATFLUX2D 
// Prototype for Fortran procedure MAPPEDGETFLATFLUX2D ...
//
void FORTRAN_NAME( MAPPEDGETFLATFLUX2D ,mappedgetflatflux2d )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(Jgxx)
      ,CHFp_BOX(FCRegion)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx) );

#define FORT_MAPPEDGETFLATFLUX2D FORTRAN_NAME( inlineMAPPEDGETFLATFLUX2D, inlineMAPPEDGETFLATFLUX2D)
#define FORTNT_MAPPEDGETFLATFLUX2D FORTRAN_NAME( MAPPEDGETFLATFLUX2D, mappedgetflatflux2d)

inline void FORTRAN_NAME(inlineMAPPEDGETFLATFLUX2D, inlineMAPPEDGETFLATFLUX2D)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA1(Jgxx)
      ,CHFp_BOX(FCRegion)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_MAPPEDGETFLATFLUX2D");
 FORTRAN_NAME( MAPPEDGETFLATFLUX2D ,mappedgetflatflux2d )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA1(Jgxx)
      ,CHFt_BOX(FCRegion)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDMAPPEDGETFLATFLUX2D 

#ifndef GUARDMAPPEDGETFLATFLUX3D 
#define GUARDMAPPEDGETFLATFLUX3D 
// Prototype for Fortran procedure MAPPEDGETFLATFLUX3D ...
//
void FORTRAN_NAME( MAPPEDGETFLATFLUX3D ,mappedgetflatflux3d )(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(extrap)
      ,CHFp_CONST_FRA(Jga)
      ,CHFp_BOX(FCRegion)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(adir) );

#define FORT_MAPPEDGETFLATFLUX3D FORTRAN_NAME( inlineMAPPEDGETFLATFLUX3D, inlineMAPPEDGETFLATFLUX3D)
#define FORTNT_MAPPEDGETFLATFLUX3D FORTRAN_NAME( MAPPEDGETFLATFLUX3D, mappedgetflatflux3d)

inline void FORTRAN_NAME(inlineMAPPEDGETFLATFLUX3D, inlineMAPPEDGETFLATFLUX3D)(
      CHFp_FRA(flux)
      ,CHFp_CONST_FRA(phi)
      ,CHFp_CONST_FRA(extrap)
      ,CHFp_CONST_FRA(Jga)
      ,CHFp_BOX(FCRegion)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(adir) )
{
 CH_TIMELEAF("FORT_MAPPEDGETFLATFLUX3D");
 FORTRAN_NAME( MAPPEDGETFLATFLUX3D ,mappedgetflatflux3d )(
      CHFt_FRA(flux)
      ,CHFt_CONST_FRA(phi)
      ,CHFt_CONST_FRA(extrap)
      ,CHFt_CONST_FRA(Jga)
      ,CHFt_BOX(FCRegion)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(adir) );
}
#endif  // GUARDMAPPEDGETFLATFLUX3D 

}

#endif
