#ifndef _AMRNSF_F_H_
#define _AMRNSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTEBVUPDATEVEC3D 
#define GUARDCOMPUTEBVUPDATEVEC3D 
// Prototype for Fortran procedure ComputeBVUpdateVec3D ...
//
void FORTRAN_NAME( COMPUTEBVUPDATEVEC3D ,computebvupdatevec3d )(
      CHFp_FRA(updateVec)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_FRA1(buoyancy)
      ,CHFp_CONST_FRA1(Nsq)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(theta)
      ,CHFp_CONST_REAL(coriolisF)
      ,CHFp_CONST_REAL(gravityCoef)
      ,CHFp_BOX(region) );

#define FORT_COMPUTEBVUPDATEVEC3D FORTRAN_NAME( inlineCOMPUTEBVUPDATEVEC3D, inlineCOMPUTEBVUPDATEVEC3D)
#define FORTNT_COMPUTEBVUPDATEVEC3D FORTRAN_NAME( COMPUTEBVUPDATEVEC3D, computebvupdatevec3d)

inline void FORTRAN_NAME(inlineCOMPUTEBVUPDATEVEC3D, inlineCOMPUTEBVUPDATEVEC3D)(
      CHFp_FRA(updateVec)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_FRA1(buoyancy)
      ,CHFp_CONST_FRA1(Nsq)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(theta)
      ,CHFp_CONST_REAL(coriolisF)
      ,CHFp_CONST_REAL(gravityCoef)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_COMPUTEBVUPDATEVEC3D");
 FORTRAN_NAME( COMPUTEBVUPDATEVEC3D ,computebvupdatevec3d )(
      CHFt_FRA(updateVec)
      ,CHFt_CONST_FRA(vel)
      ,CHFt_CONST_FRA1(buoyancy)
      ,CHFt_CONST_FRA1(Nsq)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_REAL(theta)
      ,CHFt_CONST_REAL(coriolisF)
      ,CHFt_CONST_REAL(gravityCoef)
      ,CHFt_BOX(region) );
}
#endif  // GUARDCOMPUTEBVUPDATEVEC3D 

#ifndef GUARDBVUPDATEBUOYANCY 
#define GUARDBVUPDATEBUOYANCY 
// Prototype for Fortran procedure BVUpdateBuoyancy ...
//
void FORTRAN_NAME( BVUPDATEBUOYANCY ,bvupdatebuoyancy )(
      CHFp_FRA1(buoyancy)
      ,CHFp_CONST_FRA(thetaVel)
      ,CHFp_CONST_FRA(dzdXi)
      ,CHFp_CONST_FRA1(Nsq)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(region) );

#define FORT_BVUPDATEBUOYANCY FORTRAN_NAME( inlineBVUPDATEBUOYANCY, inlineBVUPDATEBUOYANCY)
#define FORTNT_BVUPDATEBUOYANCY FORTRAN_NAME( BVUPDATEBUOYANCY, bvupdatebuoyancy)

inline void FORTRAN_NAME(inlineBVUPDATEBUOYANCY, inlineBVUPDATEBUOYANCY)(
      CHFp_FRA1(buoyancy)
      ,CHFp_CONST_FRA(thetaVel)
      ,CHFp_CONST_FRA(dzdXi)
      ,CHFp_CONST_FRA1(Nsq)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_BVUPDATEBUOYANCY");
 FORTRAN_NAME( BVUPDATEBUOYANCY ,bvupdatebuoyancy )(
      CHFt_FRA1(buoyancy)
      ,CHFt_CONST_FRA(thetaVel)
      ,CHFt_CONST_FRA(dzdXi)
      ,CHFt_CONST_FRA1(Nsq)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_BOX(region) );
}
#endif  // GUARDBVUPDATEBUOYANCY 

#ifndef GUARDCOMPUTENSQANDDINV2D 
#define GUARDCOMPUTENSQANDDINV2D 
// Prototype for Fortran procedure ComputeNsqAndDinv2D ...
//
void FORTRAN_NAME( COMPUTENSQANDDINV2D ,computensqanddinv2d )(
      CHFp_FRA1(Nsq)
      ,CHFp_FRA1(Dinv)
      ,CHFp_CONST_FRA1(bbar0)
      ,CHFp_CONST_FRA1(bbar1)
      ,CHFp_CONST_FRA(dXidz)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(theta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region) );

#define FORT_COMPUTENSQANDDINV2D FORTRAN_NAME( inlineCOMPUTENSQANDDINV2D, inlineCOMPUTENSQANDDINV2D)
#define FORTNT_COMPUTENSQANDDINV2D FORTRAN_NAME( COMPUTENSQANDDINV2D, computensqanddinv2d)

inline void FORTRAN_NAME(inlineCOMPUTENSQANDDINV2D, inlineCOMPUTENSQANDDINV2D)(
      CHFp_FRA1(Nsq)
      ,CHFp_FRA1(Dinv)
      ,CHFp_CONST_FRA1(bbar0)
      ,CHFp_CONST_FRA1(bbar1)
      ,CHFp_CONST_FRA(dXidz)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(theta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_COMPUTENSQANDDINV2D");
 FORTRAN_NAME( COMPUTENSQANDDINV2D ,computensqanddinv2d )(
      CHFt_FRA1(Nsq)
      ,CHFt_FRA1(Dinv)
      ,CHFt_CONST_FRA1(bbar0)
      ,CHFt_CONST_FRA1(bbar1)
      ,CHFt_CONST_FRA(dXidz)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_REAL(theta)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(region) );
}
#endif  // GUARDCOMPUTENSQANDDINV2D 

#ifndef GUARDCOMPUTENSQANDDINV3D 
#define GUARDCOMPUTENSQANDDINV3D 
// Prototype for Fortran procedure ComputeNsqAndDinv3D ...
//
void FORTRAN_NAME( COMPUTENSQANDDINV3D ,computensqanddinv3d )(
      CHFp_FRA1(Nsq)
      ,CHFp_FRA1(Dinv)
      ,CHFp_CONST_FRA1(bbar0)
      ,CHFp_CONST_FRA1(bbar1)
      ,CHFp_CONST_FRA1(bbar2)
      ,CHFp_CONST_FRA(dXidz)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(theta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region) );

#define FORT_COMPUTENSQANDDINV3D FORTRAN_NAME( inlineCOMPUTENSQANDDINV3D, inlineCOMPUTENSQANDDINV3D)
#define FORTNT_COMPUTENSQANDDINV3D FORTRAN_NAME( COMPUTENSQANDDINV3D, computensqanddinv3d)

inline void FORTRAN_NAME(inlineCOMPUTENSQANDDINV3D, inlineCOMPUTENSQANDDINV3D)(
      CHFp_FRA1(Nsq)
      ,CHFp_FRA1(Dinv)
      ,CHFp_CONST_FRA1(bbar0)
      ,CHFp_CONST_FRA1(bbar1)
      ,CHFp_CONST_FRA1(bbar2)
      ,CHFp_CONST_FRA(dXidz)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(theta)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_COMPUTENSQANDDINV3D");
 FORTRAN_NAME( COMPUTENSQANDDINV3D ,computensqanddinv3d )(
      CHFt_FRA1(Nsq)
      ,CHFt_FRA1(Dinv)
      ,CHFt_CONST_FRA1(bbar0)
      ,CHFt_CONST_FRA1(bbar1)
      ,CHFt_CONST_FRA1(bbar2)
      ,CHFt_CONST_FRA(dXidz)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_REAL(theta)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_BOX(region) );
}
#endif  // GUARDCOMPUTENSQANDDINV3D 

#ifndef GUARDCOMPUTEMINBVDT 
#define GUARDCOMPUTEMINBVDT 
// Prototype for Fortran procedure ComputeMinBVDt ...
//
void FORTRAN_NAME( COMPUTEMINBVDT ,computeminbvdt )(
      CHFp_REAL(minDt)
      ,CHFp_CONST_FRA(c0i)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_BOX(region) );

#define FORT_COMPUTEMINBVDT FORTRAN_NAME( inlineCOMPUTEMINBVDT, inlineCOMPUTEMINBVDT)
#define FORTNT_COMPUTEMINBVDT FORTRAN_NAME( COMPUTEMINBVDT, computeminbvdt)

inline void FORTRAN_NAME(inlineCOMPUTEMINBVDT, inlineCOMPUTEMINBVDT)(
      CHFp_REAL(minDt)
      ,CHFp_CONST_FRA(c0i)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_COMPUTEMINBVDT");
 FORTRAN_NAME( COMPUTEMINBVDT ,computeminbvdt )(
      CHFt_REAL(minDt)
      ,CHFt_CONST_FRA(c0i)
      ,CHFt_CONST_FRA(vel)
      ,CHFt_CONST_REALVECT(dXi)
      ,CHFt_BOX(region) );
}
#endif  // GUARDCOMPUTEMINBVDT 

#ifndef GUARDWEIGHTEDAVG 
#define GUARDWEIGHTEDAVG 
// Prototype for Fortran procedure WeightedAvg ...
//
void FORTRAN_NAME( WEIGHTEDAVG ,weightedavg )(
      CHFp_FRA(avgState)
      ,CHFp_CONST_FRA(state1)
      ,CHFp_CONST_FRA(state2)
      ,CHFp_CONST_REAL(theta)
      ,CHFp_BOX(region) );

#define FORT_WEIGHTEDAVG FORTRAN_NAME( inlineWEIGHTEDAVG, inlineWEIGHTEDAVG)
#define FORTNT_WEIGHTEDAVG FORTRAN_NAME( WEIGHTEDAVG, weightedavg)

inline void FORTRAN_NAME(inlineWEIGHTEDAVG, inlineWEIGHTEDAVG)(
      CHFp_FRA(avgState)
      ,CHFp_CONST_FRA(state1)
      ,CHFp_CONST_FRA(state2)
      ,CHFp_CONST_REAL(theta)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_WEIGHTEDAVG");
 FORTRAN_NAME( WEIGHTEDAVG ,weightedavg )(
      CHFt_FRA(avgState)
      ,CHFt_CONST_FRA(state1)
      ,CHFt_CONST_FRA(state2)
      ,CHFt_CONST_REAL(theta)
      ,CHFt_BOX(region) );
}
#endif  // GUARDWEIGHTEDAVG 

#ifndef GUARDTWOFORMMAG2D 
#define GUARDTWOFORMMAG2D 
// Prototype for Fortran procedure TWOFORMMAG2D ...
//
void FORTRAN_NAME( TWOFORMMAG2D ,twoformmag2d )(
      CHFp_FRA1(mag)
      ,CHFp_CONST_FRA1(twoform)
      ,CHFp_CONST_FRA1(gdnzz)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region) );

#define FORT_TWOFORMMAG2D FORTRAN_NAME( inlineTWOFORMMAG2D, inlineTWOFORMMAG2D)
#define FORTNT_TWOFORMMAG2D FORTRAN_NAME( TWOFORMMAG2D, twoformmag2d)

inline void FORTRAN_NAME(inlineTWOFORMMAG2D, inlineTWOFORMMAG2D)(
      CHFp_FRA1(mag)
      ,CHFp_CONST_FRA1(twoform)
      ,CHFp_CONST_FRA1(gdnzz)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_TWOFORMMAG2D");
 FORTRAN_NAME( TWOFORMMAG2D ,twoformmag2d )(
      CHFt_FRA1(mag)
      ,CHFt_CONST_FRA1(twoform)
      ,CHFt_CONST_FRA1(gdnzz)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region) );
}
#endif  // GUARDTWOFORMMAG2D 

#ifndef GUARDPOWFAB 
#define GUARDPOWFAB 
// Prototype for Fortran procedure POWFAB ...
//
void FORTRAN_NAME( POWFAB ,powfab )(
      CHFp_FRA(state)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(power) );

#define FORT_POWFAB FORTRAN_NAME( inlinePOWFAB, inlinePOWFAB)
#define FORTNT_POWFAB FORTRAN_NAME( POWFAB, powfab)

inline void FORTRAN_NAME(inlinePOWFAB, inlinePOWFAB)(
      CHFp_FRA(state)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(power) )
{
 CH_TIMELEAF("FORT_POWFAB");
 FORTRAN_NAME( POWFAB ,powfab )(
      CHFt_FRA(state)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(power) );
}
#endif  // GUARDPOWFAB 

#ifndef GUARDCOMPUTEENERGY 
#define GUARDCOMPUTEENERGY 
// Prototype for Fortran procedure COMPUTEENERGY ...
//
void FORTRAN_NAME( COMPUTEENERGY ,computeenergy )(
      CHFp_FRA1(energy)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_FRA1(density)
      ,CHFp_BOX(region)
      ,CHFp_CONST_FRA(gdn)
      ,CHFp_CONST_REALVECT(physDx) );

#define FORT_COMPUTEENERGY FORTRAN_NAME( inlineCOMPUTEENERGY, inlineCOMPUTEENERGY)
#define FORTNT_COMPUTEENERGY FORTRAN_NAME( COMPUTEENERGY, computeenergy)

inline void FORTRAN_NAME(inlineCOMPUTEENERGY, inlineCOMPUTEENERGY)(
      CHFp_FRA1(energy)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_CONST_FRA1(density)
      ,CHFp_BOX(region)
      ,CHFp_CONST_FRA(gdn)
      ,CHFp_CONST_REALVECT(physDx) )
{
 CH_TIMELEAF("FORT_COMPUTEENERGY");
 FORTRAN_NAME( COMPUTEENERGY ,computeenergy )(
      CHFt_FRA1(energy)
      ,CHFt_CONST_FRA(vel)
      ,CHFt_CONST_FRA1(density)
      ,CHFt_BOX(region)
      ,CHFt_CONST_FRA(gdn)
      ,CHFt_CONST_REALVECT(physDx) );
}
#endif  // GUARDCOMPUTEENERGY 

#ifndef GUARDCOMPUTEKINETICENERGY 
#define GUARDCOMPUTEKINETICENERGY 
// Prototype for Fortran procedure COMPUTEKINETICENERGY ...
//
void FORTRAN_NAME( COMPUTEKINETICENERGY ,computekineticenergy )(
      CHFp_FRA1(energy)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_BOX(region)
      ,CHFp_CONST_FRA(gdn)
      ,CHFp_CONST_REALVECT(physDx) );

#define FORT_COMPUTEKINETICENERGY FORTRAN_NAME( inlineCOMPUTEKINETICENERGY, inlineCOMPUTEKINETICENERGY)
#define FORTNT_COMPUTEKINETICENERGY FORTRAN_NAME( COMPUTEKINETICENERGY, computekineticenergy)

inline void FORTRAN_NAME(inlineCOMPUTEKINETICENERGY, inlineCOMPUTEKINETICENERGY)(
      CHFp_FRA1(energy)
      ,CHFp_CONST_FRA(vel)
      ,CHFp_BOX(region)
      ,CHFp_CONST_FRA(gdn)
      ,CHFp_CONST_REALVECT(physDx) )
{
 CH_TIMELEAF("FORT_COMPUTEKINETICENERGY");
 FORTRAN_NAME( COMPUTEKINETICENERGY ,computekineticenergy )(
      CHFt_FRA1(energy)
      ,CHFt_CONST_FRA(vel)
      ,CHFt_BOX(region)
      ,CHFt_CONST_FRA(gdn)
      ,CHFt_CONST_REALVECT(physDx) );
}
#endif  // GUARDCOMPUTEKINETICENERGY 

}

#endif
