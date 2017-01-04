#ifndef _STRATUTILSF_F_H_
#define _STRATUTILSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTEMAXBVFREQ2D 
#define GUARDCOMPUTEMAXBVFREQ2D 
// Prototype for Fortran procedure ComputeMaxBVFreq2D ...
//
void FORTRAN_NAME( COMPUTEMAXBVFREQ2D ,computemaxbvfreq2d )(
      CHFp_REAL(maxN)
      ,CHFp_CONST_FRA1(buoyancy0)
      ,CHFp_CONST_FRA1(buoyancy1)
      ,CHFp_CONST_FRA(dXidz)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(p) );

#define FORT_COMPUTEMAXBVFREQ2D FORTRAN_NAME( inlineCOMPUTEMAXBVFREQ2D, inlineCOMPUTEMAXBVFREQ2D)
#define FORTNT_COMPUTEMAXBVFREQ2D FORTRAN_NAME( COMPUTEMAXBVFREQ2D, computemaxbvfreq2d)

inline void FORTRAN_NAME(inlineCOMPUTEMAXBVFREQ2D, inlineCOMPUTEMAXBVFREQ2D)(
      CHFp_REAL(maxN)
      ,CHFp_CONST_FRA1(buoyancy0)
      ,CHFp_CONST_FRA1(buoyancy1)
      ,CHFp_CONST_FRA(dXidz)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(p) )
{
 CH_TIMELEAF("FORT_COMPUTEMAXBVFREQ2D");
 FORTRAN_NAME( COMPUTEMAXBVFREQ2D ,computemaxbvfreq2d )(
      CHFt_REAL(maxN)
      ,CHFt_CONST_FRA1(buoyancy0)
      ,CHFt_CONST_FRA1(buoyancy1)
      ,CHFt_CONST_FRA(dXidz)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(p) );
}
#endif  // GUARDCOMPUTEMAXBVFREQ2D 

#ifndef GUARDCOMPUTEMAXBVFREQ3D 
#define GUARDCOMPUTEMAXBVFREQ3D 
// Prototype for Fortran procedure ComputeMaxBVFreq3D ...
//
void FORTRAN_NAME( COMPUTEMAXBVFREQ3D ,computemaxbvfreq3d )(
      CHFp_REAL(maxN)
      ,CHFp_CONST_FRA1(buoyancy0)
      ,CHFp_CONST_FRA1(buoyancy1)
      ,CHFp_CONST_FRA1(buoyancy2)
      ,CHFp_CONST_FRA(dXidz)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(p) );

#define FORT_COMPUTEMAXBVFREQ3D FORTRAN_NAME( inlineCOMPUTEMAXBVFREQ3D, inlineCOMPUTEMAXBVFREQ3D)
#define FORTNT_COMPUTEMAXBVFREQ3D FORTRAN_NAME( COMPUTEMAXBVFREQ3D, computemaxbvfreq3d)

inline void FORTRAN_NAME(inlineCOMPUTEMAXBVFREQ3D, inlineCOMPUTEMAXBVFREQ3D)(
      CHFp_REAL(maxN)
      ,CHFp_CONST_FRA1(buoyancy0)
      ,CHFp_CONST_FRA1(buoyancy1)
      ,CHFp_CONST_FRA1(buoyancy2)
      ,CHFp_CONST_FRA(dXidz)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(p) )
{
 CH_TIMELEAF("FORT_COMPUTEMAXBVFREQ3D");
 FORTRAN_NAME( COMPUTEMAXBVFREQ3D ,computemaxbvfreq3d )(
      CHFt_REAL(maxN)
      ,CHFt_CONST_FRA1(buoyancy0)
      ,CHFt_CONST_FRA1(buoyancy1)
      ,CHFt_CONST_FRA1(buoyancy2)
      ,CHFt_CONST_FRA(dXidz)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(p) );
}
#endif  // GUARDCOMPUTEMAXBVFREQ3D 

#ifndef GUARDCOMPUTEBVFREQ2D 
#define GUARDCOMPUTEBVFREQ2D 
// Prototype for Fortran procedure ComputeBVFreq2D ...
//
void FORTRAN_NAME( COMPUTEBVFREQ2D ,computebvfreq2d )(
      CHFp_FRA1(BVFreq)
      ,CHFp_CONST_FRA1(buoyancy0)
      ,CHFp_CONST_FRA1(buoyancy1)
      ,CHFp_CONST_FRA1(J01)
      ,CHFp_CONST_FRA1(J11)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(p)
      ,CHFp_CONST_INT(doAbs) );

#define FORT_COMPUTEBVFREQ2D FORTRAN_NAME( inlineCOMPUTEBVFREQ2D, inlineCOMPUTEBVFREQ2D)
#define FORTNT_COMPUTEBVFREQ2D FORTRAN_NAME( COMPUTEBVFREQ2D, computebvfreq2d)

inline void FORTRAN_NAME(inlineCOMPUTEBVFREQ2D, inlineCOMPUTEBVFREQ2D)(
      CHFp_FRA1(BVFreq)
      ,CHFp_CONST_FRA1(buoyancy0)
      ,CHFp_CONST_FRA1(buoyancy1)
      ,CHFp_CONST_FRA1(J01)
      ,CHFp_CONST_FRA1(J11)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(p)
      ,CHFp_CONST_INT(doAbs) )
{
 CH_TIMELEAF("FORT_COMPUTEBVFREQ2D");
 FORTRAN_NAME( COMPUTEBVFREQ2D ,computebvfreq2d )(
      CHFt_FRA1(BVFreq)
      ,CHFt_CONST_FRA1(buoyancy0)
      ,CHFt_CONST_FRA1(buoyancy1)
      ,CHFt_CONST_FRA1(J01)
      ,CHFt_CONST_FRA1(J11)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(p)
      ,CHFt_CONST_INT(doAbs) );
}
#endif  // GUARDCOMPUTEBVFREQ2D 

#ifndef GUARDCOMPUTEBVFREQ3D 
#define GUARDCOMPUTEBVFREQ3D 
// Prototype for Fortran procedure ComputeBVFreq3D ...
//
void FORTRAN_NAME( COMPUTEBVFREQ3D ,computebvfreq3d )(
      CHFp_FRA1(BVFreq)
      ,CHFp_CONST_FRA1(buoyancy0)
      ,CHFp_CONST_FRA1(buoyancy1)
      ,CHFp_CONST_FRA1(buoyancy2)
      ,CHFp_CONST_FRA1(J02)
      ,CHFp_CONST_FRA1(J12)
      ,CHFp_CONST_FRA1(J22)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(p)
      ,CHFp_CONST_INT(doAbs) );

#define FORT_COMPUTEBVFREQ3D FORTRAN_NAME( inlineCOMPUTEBVFREQ3D, inlineCOMPUTEBVFREQ3D)
#define FORTNT_COMPUTEBVFREQ3D FORTRAN_NAME( COMPUTEBVFREQ3D, computebvfreq3d)

inline void FORTRAN_NAME(inlineCOMPUTEBVFREQ3D, inlineCOMPUTEBVFREQ3D)(
      CHFp_FRA1(BVFreq)
      ,CHFp_CONST_FRA1(buoyancy0)
      ,CHFp_CONST_FRA1(buoyancy1)
      ,CHFp_CONST_FRA1(buoyancy2)
      ,CHFp_CONST_FRA1(J02)
      ,CHFp_CONST_FRA1(J12)
      ,CHFp_CONST_FRA1(J22)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(p)
      ,CHFp_CONST_INT(doAbs) )
{
 CH_TIMELEAF("FORT_COMPUTEBVFREQ3D");
 FORTRAN_NAME( COMPUTEBVFREQ3D ,computebvfreq3d )(
      CHFt_FRA1(BVFreq)
      ,CHFt_CONST_FRA1(buoyancy0)
      ,CHFt_CONST_FRA1(buoyancy1)
      ,CHFt_CONST_FRA1(buoyancy2)
      ,CHFt_CONST_FRA1(J02)
      ,CHFt_CONST_FRA1(J12)
      ,CHFt_CONST_FRA1(J22)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(p)
      ,CHFt_CONST_INT(doAbs) );
}
#endif  // GUARDCOMPUTEBVFREQ3D 

#ifndef GUARDSOLVEVERTEIGENPROBLEM 
#define GUARDSOLVEVERTEIGENPROBLEM 
// Prototype for Fortran procedure SOLVEVERTEIGENPROBLEM ...
//
void FORTRAN_NAME( SOLVEVERTEIGENPROBLEM ,solveverteigenproblem )(
      CHFp_REAL(c0)
      ,CHFp_FRA1(phi0)
      ,CHFp_CONST_FRA1(Nsq)
      ,CHFp_CONST_FRA1(jacFC)
      ,CHFp_CONST_FRA1(jacCC)
      ,CHFp_CONST_INTVECT(loIdx)
      ,CHFp_CONST_INT(N)
      ,CHFp_CONST_REAL(dz)
      ,CHFp_CONST_REAL(abstol)
      ,CHFp_FIA1(IFAIL)
      ,CHFp_FIA1(IWORK)
      ,CHFp_FRA1(AB)
      ,CHFp_FRA1(BB)
      ,CHFp_FRA1(Q)
      ,CHFp_FRA1(W)
      ,CHFp_FRA1(WORK)
      ,CHFp_FRA1(Z) );

#define FORT_SOLVEVERTEIGENPROBLEM FORTRAN_NAME( inlineSOLVEVERTEIGENPROBLEM, inlineSOLVEVERTEIGENPROBLEM)
#define FORTNT_SOLVEVERTEIGENPROBLEM FORTRAN_NAME( SOLVEVERTEIGENPROBLEM, solveverteigenproblem)

inline void FORTRAN_NAME(inlineSOLVEVERTEIGENPROBLEM, inlineSOLVEVERTEIGENPROBLEM)(
      CHFp_REAL(c0)
      ,CHFp_FRA1(phi0)
      ,CHFp_CONST_FRA1(Nsq)
      ,CHFp_CONST_FRA1(jacFC)
      ,CHFp_CONST_FRA1(jacCC)
      ,CHFp_CONST_INTVECT(loIdx)
      ,CHFp_CONST_INT(N)
      ,CHFp_CONST_REAL(dz)
      ,CHFp_CONST_REAL(abstol)
      ,CHFp_FIA1(IFAIL)
      ,CHFp_FIA1(IWORK)
      ,CHFp_FRA1(AB)
      ,CHFp_FRA1(BB)
      ,CHFp_FRA1(Q)
      ,CHFp_FRA1(W)
      ,CHFp_FRA1(WORK)
      ,CHFp_FRA1(Z) )
{
 CH_TIMELEAF("FORT_SOLVEVERTEIGENPROBLEM");
 FORTRAN_NAME( SOLVEVERTEIGENPROBLEM ,solveverteigenproblem )(
      CHFt_REAL(c0)
      ,CHFt_FRA1(phi0)
      ,CHFt_CONST_FRA1(Nsq)
      ,CHFt_CONST_FRA1(jacFC)
      ,CHFt_CONST_FRA1(jacCC)
      ,CHFt_CONST_INTVECT(loIdx)
      ,CHFt_CONST_INT(N)
      ,CHFt_CONST_REAL(dz)
      ,CHFt_CONST_REAL(abstol)
      ,CHFt_FIA1(IFAIL)
      ,CHFt_FIA1(IWORK)
      ,CHFt_FRA1(AB)
      ,CHFt_FRA1(BB)
      ,CHFt_FRA1(Q)
      ,CHFt_FRA1(W)
      ,CHFt_FRA1(WORK)
      ,CHFt_FRA1(Z) );
}
#endif  // GUARDSOLVEVERTEIGENPROBLEM 

#ifndef GUARDPROJECTPHASESPEED 
#define GUARDPROJECTPHASESPEED 
// Prototype for Fortran procedure ProjectPhaseSpeed ...
//
void FORTRAN_NAME( PROJECTPHASESPEED ,projectphasespeed )(
      CHFp_FRA(c0i)
      ,CHFp_CONST_FRA1(c0)
      ,CHFp_CONST_FRA(dXidx)
      ,CHFp_BOX(region) );

#define FORT_PROJECTPHASESPEED FORTRAN_NAME( inlinePROJECTPHASESPEED, inlinePROJECTPHASESPEED)
#define FORTNT_PROJECTPHASESPEED FORTRAN_NAME( PROJECTPHASESPEED, projectphasespeed)

inline void FORTRAN_NAME(inlinePROJECTPHASESPEED, inlinePROJECTPHASESPEED)(
      CHFp_FRA(c0i)
      ,CHFp_CONST_FRA1(c0)
      ,CHFp_CONST_FRA(dXidx)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_PROJECTPHASESPEED");
 FORTRAN_NAME( PROJECTPHASESPEED ,projectphasespeed )(
      CHFt_FRA(c0i)
      ,CHFt_CONST_FRA1(c0)
      ,CHFt_CONST_FRA(dXidx)
      ,CHFt_BOX(region) );
}
#endif  // GUARDPROJECTPHASESPEED 

#ifndef GUARDCOMPUTEBVFREQCARTESIAN 
#define GUARDCOMPUTEBVFREQCARTESIAN 
// Prototype for Fortran procedure ComputeBVFreqCartesian ...
//
void FORTRAN_NAME( COMPUTEBVFREQCARTESIAN ,computebvfreqcartesian )(
      CHFp_FRA1(freq)
      ,CHFp_CONST_FRA1(bouyancy)
      ,CHFp_CONST_REAL(dz)
      ,CHFp_BOX(NBox)
      ,CHFp_CONST_INTVECT(NBoxType)
      ,CHFp_CONST_INTVECT(BBoxType)
      ,CHFp_CONST_REAL(p) );

#define FORT_COMPUTEBVFREQCARTESIAN FORTRAN_NAME( inlineCOMPUTEBVFREQCARTESIAN, inlineCOMPUTEBVFREQCARTESIAN)
#define FORTNT_COMPUTEBVFREQCARTESIAN FORTRAN_NAME( COMPUTEBVFREQCARTESIAN, computebvfreqcartesian)

inline void FORTRAN_NAME(inlineCOMPUTEBVFREQCARTESIAN, inlineCOMPUTEBVFREQCARTESIAN)(
      CHFp_FRA1(freq)
      ,CHFp_CONST_FRA1(bouyancy)
      ,CHFp_CONST_REAL(dz)
      ,CHFp_BOX(NBox)
      ,CHFp_CONST_INTVECT(NBoxType)
      ,CHFp_CONST_INTVECT(BBoxType)
      ,CHFp_CONST_REAL(p) )
{
 CH_TIMELEAF("FORT_COMPUTEBVFREQCARTESIAN");
 FORTRAN_NAME( COMPUTEBVFREQCARTESIAN ,computebvfreqcartesian )(
      CHFt_FRA1(freq)
      ,CHFt_CONST_FRA1(bouyancy)
      ,CHFt_CONST_REAL(dz)
      ,CHFt_BOX(NBox)
      ,CHFt_CONST_INTVECT(NBoxType)
      ,CHFt_CONST_INTVECT(BBoxType)
      ,CHFt_CONST_REAL(p) );
}
#endif  // GUARDCOMPUTEBVFREQCARTESIAN 

#ifndef GUARDSOLVEVERTEIGENPROBLEMCARTESIAN 
#define GUARDSOLVEVERTEIGENPROBLEMCARTESIAN 
// Prototype for Fortran procedure SOLVEVERTEIGENPROBLEMCARTESIAN ...
//
void FORTRAN_NAME( SOLVEVERTEIGENPROBLEMCARTESIAN ,solveverteigenproblemcartesian )(
      CHFp_REAL(c0)
      ,CHFp_FRA1(phi0)
      ,CHFp_CONST_FRA1(Nsq)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(N)
      ,CHFp_CONST_REAL(dz)
      ,CHFp_CONST_REAL(abstol)
      ,CHFp_FIA1(IFAIL)
      ,CHFp_FIA1(IWORK)
      ,CHFp_FRA1(AB)
      ,CHFp_FRA1(BB)
      ,CHFp_FRA1(Q)
      ,CHFp_FRA1(W)
      ,CHFp_FRA1(WORK)
      ,CHFp_FRA1(Z) );

#define FORT_SOLVEVERTEIGENPROBLEMCARTESIAN FORTRAN_NAME( inlineSOLVEVERTEIGENPROBLEMCARTESIAN, inlineSOLVEVERTEIGENPROBLEMCARTESIAN)
#define FORTNT_SOLVEVERTEIGENPROBLEMCARTESIAN FORTRAN_NAME( SOLVEVERTEIGENPROBLEMCARTESIAN, solveverteigenproblemcartesian)

inline void FORTRAN_NAME(inlineSOLVEVERTEIGENPROBLEMCARTESIAN, inlineSOLVEVERTEIGENPROBLEMCARTESIAN)(
      CHFp_REAL(c0)
      ,CHFp_FRA1(phi0)
      ,CHFp_CONST_FRA1(Nsq)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(N)
      ,CHFp_CONST_REAL(dz)
      ,CHFp_CONST_REAL(abstol)
      ,CHFp_FIA1(IFAIL)
      ,CHFp_FIA1(IWORK)
      ,CHFp_FRA1(AB)
      ,CHFp_FRA1(BB)
      ,CHFp_FRA1(Q)
      ,CHFp_FRA1(W)
      ,CHFp_FRA1(WORK)
      ,CHFp_FRA1(Z) )
{
 CH_TIMELEAF("FORT_SOLVEVERTEIGENPROBLEMCARTESIAN");
 FORTRAN_NAME( SOLVEVERTEIGENPROBLEMCARTESIAN ,solveverteigenproblemcartesian )(
      CHFt_REAL(c0)
      ,CHFt_FRA1(phi0)
      ,CHFt_CONST_FRA1(Nsq)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(N)
      ,CHFt_CONST_REAL(dz)
      ,CHFt_CONST_REAL(abstol)
      ,CHFt_FIA1(IFAIL)
      ,CHFt_FIA1(IWORK)
      ,CHFt_FRA1(AB)
      ,CHFt_FRA1(BB)
      ,CHFt_FRA1(Q)
      ,CHFt_FRA1(W)
      ,CHFt_FRA1(WORK)
      ,CHFt_FRA1(Z) );
}
#endif  // GUARDSOLVEVERTEIGENPROBLEMCARTESIAN 

#ifndef GUARDCOMPUTEHORIZCOEFFSCARTESIAN 
#define GUARDCOMPUTEHORIZCOEFFSCARTESIAN 
// Prototype for Fortran procedure COMPUTEHORIZCOEFFSCARTESIAN ...
//
void FORTRAN_NAME( COMPUTEHORIZCOEFFSCARTESIAN ,computehorizcoeffscartesian )(
      CHFp_REAL(alpha)
      ,CHFp_REAL(beta)
      ,CHFp_REAL(gamma)
      ,CHFp_CONST_FRA1(CCphi)
      ,CHFp_CONST_FRA1(CCNsq)
      ,CHFp_CONST_FRA1(FCNsq)
      ,CHFp_CONST_REAL(c0)
      ,CHFp_BOX(vertBox)
      ,CHFp_CONST_REAL(dz) );

#define FORT_COMPUTEHORIZCOEFFSCARTESIAN FORTRAN_NAME( inlineCOMPUTEHORIZCOEFFSCARTESIAN, inlineCOMPUTEHORIZCOEFFSCARTESIAN)
#define FORTNT_COMPUTEHORIZCOEFFSCARTESIAN FORTRAN_NAME( COMPUTEHORIZCOEFFSCARTESIAN, computehorizcoeffscartesian)

inline void FORTRAN_NAME(inlineCOMPUTEHORIZCOEFFSCARTESIAN, inlineCOMPUTEHORIZCOEFFSCARTESIAN)(
      CHFp_REAL(alpha)
      ,CHFp_REAL(beta)
      ,CHFp_REAL(gamma)
      ,CHFp_CONST_FRA1(CCphi)
      ,CHFp_CONST_FRA1(CCNsq)
      ,CHFp_CONST_FRA1(FCNsq)
      ,CHFp_CONST_REAL(c0)
      ,CHFp_BOX(vertBox)
      ,CHFp_CONST_REAL(dz) )
{
 CH_TIMELEAF("FORT_COMPUTEHORIZCOEFFSCARTESIAN");
 FORTRAN_NAME( COMPUTEHORIZCOEFFSCARTESIAN ,computehorizcoeffscartesian )(
      CHFt_REAL(alpha)
      ,CHFt_REAL(beta)
      ,CHFt_REAL(gamma)
      ,CHFt_CONST_FRA1(CCphi)
      ,CHFt_CONST_FRA1(CCNsq)
      ,CHFt_CONST_FRA1(FCNsq)
      ,CHFt_CONST_REAL(c0)
      ,CHFt_BOX(vertBox)
      ,CHFt_CONST_REAL(dz) );
}
#endif  // GUARDCOMPUTEHORIZCOEFFSCARTESIAN 

#ifndef GUARDFILLHORIZSTRUCTURE 
#define GUARDFILLHORIZSTRUCTURE 
// Prototype for Fortran procedure FILLHORIZSTRUCTURE ...
//
void FORTRAN_NAME( FILLHORIZSTRUCTURE ,fillhorizstructure )(
      CHFp_FRA1(A0)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(gamma)
      ,CHFp_CONST_REAL(c0)
      ,CHFp_CONST_REAL(xcenter)
      ,CHFp_CONST_REAL(amp)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INTVECT(boxType)
      ,CHFp_CONST_REAL(physDx) );

#define FORT_FILLHORIZSTRUCTURE FORTRAN_NAME( inlineFILLHORIZSTRUCTURE, inlineFILLHORIZSTRUCTURE)
#define FORTNT_FILLHORIZSTRUCTURE FORTRAN_NAME( FILLHORIZSTRUCTURE, fillhorizstructure)

inline void FORTRAN_NAME(inlineFILLHORIZSTRUCTURE, inlineFILLHORIZSTRUCTURE)(
      CHFp_FRA1(A0)
      ,CHFp_CONST_REAL(alpha)
      ,CHFp_CONST_REAL(beta)
      ,CHFp_CONST_REAL(gamma)
      ,CHFp_CONST_REAL(c0)
      ,CHFp_CONST_REAL(xcenter)
      ,CHFp_CONST_REAL(amp)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INTVECT(boxType)
      ,CHFp_CONST_REAL(physDx) )
{
 CH_TIMELEAF("FORT_FILLHORIZSTRUCTURE");
 FORTRAN_NAME( FILLHORIZSTRUCTURE ,fillhorizstructure )(
      CHFt_FRA1(A0)
      ,CHFt_CONST_REAL(alpha)
      ,CHFt_CONST_REAL(beta)
      ,CHFt_CONST_REAL(gamma)
      ,CHFt_CONST_REAL(c0)
      ,CHFt_CONST_REAL(xcenter)
      ,CHFt_CONST_REAL(amp)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INTVECT(boxType)
      ,CHFt_CONST_REAL(physDx) );
}
#endif  // GUARDFILLHORIZSTRUCTURE 

#ifndef GUARDCOMPUTESTRAINRATESQ 
#define GUARDCOMPUTESTRAINRATESQ 
// Prototype for Fortran procedure ComputeStrainRateSq ...
//
void FORTRAN_NAME( COMPUTESTRAINRATESQ ,computestrainratesq )(
      CHFp_FRA1(Ssq)
      ,CHFp_CONST_FRA(cartVel)
      ,CHFp_CONST_FRA(dXidx)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region) );

#define FORT_COMPUTESTRAINRATESQ FORTRAN_NAME( inlineCOMPUTESTRAINRATESQ, inlineCOMPUTESTRAINRATESQ)
#define FORTNT_COMPUTESTRAINRATESQ FORTRAN_NAME( COMPUTESTRAINRATESQ, computestrainratesq)

inline void FORTRAN_NAME(inlineCOMPUTESTRAINRATESQ, inlineCOMPUTESTRAINRATESQ)(
      CHFp_FRA1(Ssq)
      ,CHFp_CONST_FRA(cartVel)
      ,CHFp_CONST_FRA(dXidx)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_COMPUTESTRAINRATESQ");
 FORTRAN_NAME( COMPUTESTRAINRATESQ ,computestrainratesq )(
      CHFt_FRA1(Ssq)
      ,CHFt_CONST_FRA(cartVel)
      ,CHFt_CONST_FRA(dXidx)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(region) );
}
#endif  // GUARDCOMPUTESTRAINRATESQ 

#ifndef GUARDCOMPUTEHORIZSTRAINRATESQ 
#define GUARDCOMPUTEHORIZSTRAINRATESQ 
// Prototype for Fortran procedure ComputeHorizStrainRateSq ...
//
void FORTRAN_NAME( COMPUTEHORIZSTRAINRATESQ ,computehorizstrainratesq )(
      CHFp_FRA1(Ssq)
      ,CHFp_CONST_FRA(cartVel)
      ,CHFp_CONST_FRA(dXidx)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region) );

#define FORT_COMPUTEHORIZSTRAINRATESQ FORTRAN_NAME( inlineCOMPUTEHORIZSTRAINRATESQ, inlineCOMPUTEHORIZSTRAINRATESQ)
#define FORTNT_COMPUTEHORIZSTRAINRATESQ FORTRAN_NAME( COMPUTEHORIZSTRAINRATESQ, computehorizstrainratesq)

inline void FORTRAN_NAME(inlineCOMPUTEHORIZSTRAINRATESQ, inlineCOMPUTEHORIZSTRAINRATESQ)(
      CHFp_FRA1(Ssq)
      ,CHFp_CONST_FRA(cartVel)
      ,CHFp_CONST_FRA(dXidx)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_COMPUTEHORIZSTRAINRATESQ");
 FORTRAN_NAME( COMPUTEHORIZSTRAINRATESQ ,computehorizstrainratesq )(
      CHFt_FRA1(Ssq)
      ,CHFt_CONST_FRA(cartVel)
      ,CHFt_CONST_FRA(dXidx)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(region) );
}
#endif  // GUARDCOMPUTEHORIZSTRAINRATESQ 

#ifndef GUARDCOMPUTEMASKEDRIG 
#define GUARDCOMPUTEMASKEDRIG 
// Prototype for Fortran procedure ComputeMaskedRig ...
//
void FORTRAN_NAME( COMPUTEMASKEDRIG ,computemaskedrig )(
      CHFp_FRA1(Rig)
      ,CHFp_CONST_FRA1(Nsq)
      ,CHFp_CONST_FRA1(Ssq)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(tol)
      ,CHFp_CONST_REAL(maskVal) );

#define FORT_COMPUTEMASKEDRIG FORTRAN_NAME( inlineCOMPUTEMASKEDRIG, inlineCOMPUTEMASKEDRIG)
#define FORTNT_COMPUTEMASKEDRIG FORTRAN_NAME( COMPUTEMASKEDRIG, computemaskedrig)

inline void FORTRAN_NAME(inlineCOMPUTEMASKEDRIG, inlineCOMPUTEMASKEDRIG)(
      CHFp_FRA1(Rig)
      ,CHFp_CONST_FRA1(Nsq)
      ,CHFp_CONST_FRA1(Ssq)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(tol)
      ,CHFp_CONST_REAL(maskVal) )
{
 CH_TIMELEAF("FORT_COMPUTEMASKEDRIG");
 FORTRAN_NAME( COMPUTEMASKEDRIG ,computemaskedrig )(
      CHFt_FRA1(Rig)
      ,CHFt_CONST_FRA1(Nsq)
      ,CHFt_CONST_FRA1(Ssq)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(tol)
      ,CHFt_CONST_REAL(maskVal) );
}
#endif  // GUARDCOMPUTEMASKEDRIG 

}

#endif
