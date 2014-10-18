#ifndef _LEVELGEOMETRYF_F_H_
#define _LEVELGEOMETRYF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDFILLMAPPINGDISP 
#define GUARDFILLMAPPINGDISP 
// Prototype for Fortran procedure FILLMAPPINGDISP ...
//
void FORTRAN_NAME( FILLMAPPINGDISP ,fillmappingdisp )(
      CHFp_FRA(disp)
      ,CHFp_CONST_FRA(pos)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INTVECT(boxType) );

#define FORT_FILLMAPPINGDISP FORTRAN_NAME( inlineFILLMAPPINGDISP, inlineFILLMAPPINGDISP)
#define FORTNT_FILLMAPPINGDISP FORTRAN_NAME( FILLMAPPINGDISP, fillmappingdisp)

inline void FORTRAN_NAME(inlineFILLMAPPINGDISP, inlineFILLMAPPINGDISP)(
      CHFp_FRA(disp)
      ,CHFp_CONST_FRA(pos)
      ,CHFp_CONST_REALVECT(dXi)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INTVECT(boxType) )
{
 CH_TIMELEAF("FORT_FILLMAPPINGDISP");
 FORTRAN_NAME( FILLMAPPINGDISP ,fillmappingdisp )(
      CHFt_FRA(disp)
      ,CHFt_CONST_FRA(pos)
      ,CHFt_CONST_REALVECT(dXi)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INTVECT(boxType) );
}
#endif  // GUARDFILLMAPPINGDISP 

#ifndef GUARDCONTRACTMATRIXVECTORCC 
#define GUARDCONTRACTMATRIXVECTORCC 
// Prototype for Fortran procedure CONTRACTMATRIXVECTORCC ...
//
void FORTRAN_NAME( CONTRACTMATRIXVECTORCC ,contractmatrixvectorcc )(
      CHFp_FRA(vector)
      ,CHFp_CONST_FRA(matrix)
      ,CHFp_BOX(region) );

#define FORT_CONTRACTMATRIXVECTORCC FORTRAN_NAME( inlineCONTRACTMATRIXVECTORCC, inlineCONTRACTMATRIXVECTORCC)
#define FORTNT_CONTRACTMATRIXVECTORCC FORTRAN_NAME( CONTRACTMATRIXVECTORCC, contractmatrixvectorcc)

inline void FORTRAN_NAME(inlineCONTRACTMATRIXVECTORCC, inlineCONTRACTMATRIXVECTORCC)(
      CHFp_FRA(vector)
      ,CHFp_CONST_FRA(matrix)
      ,CHFp_BOX(region) )
{
 CH_TIMELEAF("FORT_CONTRACTMATRIXVECTORCC");
 FORTRAN_NAME( CONTRACTMATRIXVECTORCC ,contractmatrixvectorcc )(
      CHFt_FRA(vector)
      ,CHFt_CONST_FRA(matrix)
      ,CHFt_BOX(region) );
}
#endif  // GUARDCONTRACTMATRIXVECTORCC 

#ifndef GUARDCONTRACTWITHMETRIC 
#define GUARDCONTRACTWITHMETRIC 
// Prototype for Fortran procedure CONTRACTWITHMETRIC ...
//
void FORTRAN_NAME( CONTRACTWITHMETRIC ,contractwithmetric )(
      CHFp_FRA(dest)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(region)
      ,CHFp_CONST_FRA(metric) );

#define FORT_CONTRACTWITHMETRIC FORTRAN_NAME( inlineCONTRACTWITHMETRIC, inlineCONTRACTWITHMETRIC)
#define FORTNT_CONTRACTWITHMETRIC FORTRAN_NAME( CONTRACTWITHMETRIC, contractwithmetric)

inline void FORTRAN_NAME(inlineCONTRACTWITHMETRIC, inlineCONTRACTWITHMETRIC)(
      CHFp_FRA(dest)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(region)
      ,CHFp_CONST_FRA(metric) )
{
 CH_TIMELEAF("FORT_CONTRACTWITHMETRIC");
 FORTRAN_NAME( CONTRACTWITHMETRIC ,contractwithmetric )(
      CHFt_FRA(dest)
      ,CHFt_CONST_FRA(src)
      ,CHFt_BOX(region)
      ,CHFt_CONST_FRA(metric) );
}
#endif  // GUARDCONTRACTWITHMETRIC 

#ifndef GUARDCONTRACTVECTORS 
#define GUARDCONTRACTVECTORS 
// Prototype for Fortran procedure CONTRACTVECTORS ...
//
void FORTRAN_NAME( CONTRACTVECTORS ,contractvectors )(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA(vec1)
      ,CHFp_CONST_FRA(vec2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_FRA(gdn) );

#define FORT_CONTRACTVECTORS FORTRAN_NAME( inlineCONTRACTVECTORS, inlineCONTRACTVECTORS)
#define FORTNT_CONTRACTVECTORS FORTRAN_NAME( CONTRACTVECTORS, contractvectors)

inline void FORTRAN_NAME(inlineCONTRACTVECTORS, inlineCONTRACTVECTORS)(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA(vec1)
      ,CHFp_CONST_FRA(vec2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_FRA(gdn) )
{
 CH_TIMELEAF("FORT_CONTRACTVECTORS");
 FORTRAN_NAME( CONTRACTVECTORS ,contractvectors )(
      CHFt_FRA1(res)
      ,CHFt_CONST_FRA(vec1)
      ,CHFt_CONST_FRA(vec2)
      ,CHFt_BOX(region)
      ,CHFt_CONST_FRA(gdn) );
}
#endif  // GUARDCONTRACTVECTORS 

#ifndef GUARDSYMPROD 
#define GUARDSYMPROD 
// Prototype for Fortran procedure SYMPROD ...
//
void FORTRAN_NAME( SYMPROD ,symprod )(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA(vec1)
      ,CHFp_CONST_FRA(vec2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(adir)
      ,CHFp_CONST_INT(bdir) );

#define FORT_SYMPROD FORTRAN_NAME( inlineSYMPROD, inlineSYMPROD)
#define FORTNT_SYMPROD FORTRAN_NAME( SYMPROD, symprod)

inline void FORTRAN_NAME(inlineSYMPROD, inlineSYMPROD)(
      CHFp_FRA1(res)
      ,CHFp_CONST_FRA(vec1)
      ,CHFp_CONST_FRA(vec2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_INT(adir)
      ,CHFp_CONST_INT(bdir) )
{
 CH_TIMELEAF("FORT_SYMPROD");
 FORTRAN_NAME( SYMPROD ,symprod )(
      CHFt_FRA1(res)
      ,CHFt_CONST_FRA(vec1)
      ,CHFt_CONST_FRA(vec2)
      ,CHFt_BOX(region)
      ,CHFt_CONST_INT(adir)
      ,CHFt_CONST_INT(bdir) );
}
#endif  // GUARDSYMPROD 

}

#endif
