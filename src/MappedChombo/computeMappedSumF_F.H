#ifndef _COMPUTEMAPPEDSUMF_F_H_
#define _COMPUTEMAPPEDSUMF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCOMPUTEMAPPEDSUM 
#define GUARDCOMPUTEMAPPEDSUM 
// Prototype for Fortran procedure COMPUTEMAPPEDSUM ...
//
void FORTRAN_NAME( COMPUTEMAPPEDSUM ,computemappedsum )(
      CHFp_REAL(sum)
      ,CHFp_REAL(vol)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(J)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_BOX(stateBox)
      ,CHFp_CONST_REAL(dxScale) );

#define FORT_COMPUTEMAPPEDSUM FORTRAN_NAME( inlineCOMPUTEMAPPEDSUM, inlineCOMPUTEMAPPEDSUM)
#define FORTNT_COMPUTEMAPPEDSUM FORTRAN_NAME( COMPUTEMAPPEDSUM, computemappedsum)

inline void FORTRAN_NAME(inlineCOMPUTEMAPPEDSUM, inlineCOMPUTEMAPPEDSUM)(
      CHFp_REAL(sum)
      ,CHFp_REAL(vol)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(J)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_BOX(stateBox)
      ,CHFp_CONST_REAL(dxScale) )
{
 CH_TIMELEAF("FORT_COMPUTEMAPPEDSUM");
 FORTRAN_NAME( COMPUTEMAPPEDSUM ,computemappedsum )(
      CHFt_REAL(sum)
      ,CHFt_REAL(vol)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(J)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_BOX(stateBox)
      ,CHFt_CONST_REAL(dxScale) );
}
#endif  // GUARDCOMPUTEMAPPEDSUM 

#ifndef GUARDCOMPUTEMAPPEDSUMINV 
#define GUARDCOMPUTEMAPPEDSUMINV 
// Prototype for Fortran procedure COMPUTEMAPPEDSUMINV ...
//
void FORTRAN_NAME( COMPUTEMAPPEDSUMINV ,computemappedsuminv )(
      CHFp_REAL(sum)
      ,CHFp_REAL(vol)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_BOX(stateBox)
      ,CHFp_CONST_REAL(dxScale) );

#define FORT_COMPUTEMAPPEDSUMINV FORTRAN_NAME( inlineCOMPUTEMAPPEDSUMINV, inlineCOMPUTEMAPPEDSUMINV)
#define FORTNT_COMPUTEMAPPEDSUMINV FORTRAN_NAME( COMPUTEMAPPEDSUMINV, computemappedsuminv)

inline void FORTRAN_NAME(inlineCOMPUTEMAPPEDSUMINV, inlineCOMPUTEMAPPEDSUMINV)(
      CHFp_REAL(sum)
      ,CHFp_REAL(vol)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_BOX(stateBox)
      ,CHFp_CONST_REAL(dxScale) )
{
 CH_TIMELEAF("FORT_COMPUTEMAPPEDSUMINV");
 FORTRAN_NAME( COMPUTEMAPPEDSUMINV ,computemappedsuminv )(
      CHFt_REAL(sum)
      ,CHFt_REAL(vol)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_BOX(stateBox)
      ,CHFt_CONST_REAL(dxScale) );
}
#endif  // GUARDCOMPUTEMAPPEDSUMINV 

#ifndef GUARDCOMPUTEUNMAPPEDSUM 
#define GUARDCOMPUTEUNMAPPEDSUM 
// Prototype for Fortran procedure COMPUTEUNMAPPEDSUM ...
//
void FORTRAN_NAME( COMPUTEUNMAPPEDSUM ,computeunmappedsum )(
      CHFp_REAL(sum)
      ,CHFp_REAL(vol)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(J)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_BOX(stateBox)
      ,CHFp_CONST_REAL(dxScale) );

#define FORT_COMPUTEUNMAPPEDSUM FORTRAN_NAME( inlineCOMPUTEUNMAPPEDSUM, inlineCOMPUTEUNMAPPEDSUM)
#define FORTNT_COMPUTEUNMAPPEDSUM FORTRAN_NAME( COMPUTEUNMAPPEDSUM, computeunmappedsum)

inline void FORTRAN_NAME(inlineCOMPUTEUNMAPPEDSUM, inlineCOMPUTEUNMAPPEDSUM)(
      CHFp_REAL(sum)
      ,CHFp_REAL(vol)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(J)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_BOX(stateBox)
      ,CHFp_CONST_REAL(dxScale) )
{
 CH_TIMELEAF("FORT_COMPUTEUNMAPPEDSUM");
 FORTRAN_NAME( COMPUTEUNMAPPEDSUM ,computeunmappedsum )(
      CHFt_REAL(sum)
      ,CHFt_REAL(vol)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(J)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_BOX(stateBox)
      ,CHFt_CONST_REAL(dxScale) );
}
#endif  // GUARDCOMPUTEUNMAPPEDSUM 

#ifndef GUARDCOMPUTEUNMAPPEDSUMINV 
#define GUARDCOMPUTEUNMAPPEDSUMINV 
// Prototype for Fortran procedure COMPUTEUNMAPPEDSUMINV ...
//
void FORTRAN_NAME( COMPUTEUNMAPPEDSUMINV ,computeunmappedsuminv )(
      CHFp_REAL(sum)
      ,CHFp_REAL(vol)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_BOX(stateBox)
      ,CHFp_CONST_REAL(dxScale) );

#define FORT_COMPUTEUNMAPPEDSUMINV FORTRAN_NAME( inlineCOMPUTEUNMAPPEDSUMINV, inlineCOMPUTEUNMAPPEDSUMINV)
#define FORTNT_COMPUTEUNMAPPEDSUMINV FORTRAN_NAME( COMPUTEUNMAPPEDSUMINV, computeunmappedsuminv)

inline void FORTRAN_NAME(inlineCOMPUTEUNMAPPEDSUMINV, inlineCOMPUTEUNMAPPEDSUMINV)(
      CHFp_REAL(sum)
      ,CHFp_REAL(vol)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_BOX(stateBox)
      ,CHFp_CONST_REAL(dxScale) )
{
 CH_TIMELEAF("FORT_COMPUTEUNMAPPEDSUMINV");
 FORTRAN_NAME( COMPUTEUNMAPPEDSUMINV ,computeunmappedsuminv )(
      CHFt_REAL(sum)
      ,CHFt_REAL(vol)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_BOX(stateBox)
      ,CHFt_CONST_REAL(dxScale) );
}
#endif  // GUARDCOMPUTEUNMAPPEDSUMINV 

#ifndef GUARDCOMPUTEMAPPEDFLUXSUM 
#define GUARDCOMPUTEMAPPEDFLUXSUM 
// Prototype for Fortran procedure COMPUTEMAPPEDFLUXSUM ...
//
void FORTRAN_NAME( COMPUTEMAPPEDFLUXSUM ,computemappedfluxsum )(
      CHFp_REAL(sum)
      ,CHFp_REAL(area)
      ,CHFp_CONST_FRA1(flux)
      ,CHFp_CONST_FRA1(J)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_BOX(stateBox)
      ,CHFp_CONST_REAL(areaScale) );

#define FORT_COMPUTEMAPPEDFLUXSUM FORTRAN_NAME( inlineCOMPUTEMAPPEDFLUXSUM, inlineCOMPUTEMAPPEDFLUXSUM)
#define FORTNT_COMPUTEMAPPEDFLUXSUM FORTRAN_NAME( COMPUTEMAPPEDFLUXSUM, computemappedfluxsum)

inline void FORTRAN_NAME(inlineCOMPUTEMAPPEDFLUXSUM, inlineCOMPUTEMAPPEDFLUXSUM)(
      CHFp_REAL(sum)
      ,CHFp_REAL(area)
      ,CHFp_CONST_FRA1(flux)
      ,CHFp_CONST_FRA1(J)
      ,CHFp_CONST_FIA1(mask)
      ,CHFp_BOX(stateBox)
      ,CHFp_CONST_REAL(areaScale) )
{
 CH_TIMELEAF("FORT_COMPUTEMAPPEDFLUXSUM");
 FORTRAN_NAME( COMPUTEMAPPEDFLUXSUM ,computemappedfluxsum )(
      CHFt_REAL(sum)
      ,CHFt_REAL(area)
      ,CHFt_CONST_FRA1(flux)
      ,CHFt_CONST_FRA1(J)
      ,CHFt_CONST_FIA1(mask)
      ,CHFt_BOX(stateBox)
      ,CHFt_CONST_REAL(areaScale) );
}
#endif  // GUARDCOMPUTEMAPPEDFLUXSUM 

}

#endif
