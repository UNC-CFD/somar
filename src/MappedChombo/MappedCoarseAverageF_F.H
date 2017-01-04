#ifndef _MAPPEDCOARSEAVERAGEF_F_H_
#define _MAPPEDCOARSEAVERAGEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDUNMAPPEDAVERAGE 
#define GUARDUNMAPPEDAVERAGE 
// Prototype for Fortran procedure UNMAPPEDAVERAGE ...
//
void FORTRAN_NAME( UNMAPPEDAVERAGE ,unmappedaverage )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INTVECT(refRatio)
      ,CHFp_BOX(bref) );

#define FORT_UNMAPPEDAVERAGE FORTRAN_NAME( inlineUNMAPPEDAVERAGE, inlineUNMAPPEDAVERAGE)
#define FORTNT_UNMAPPEDAVERAGE FORTRAN_NAME( UNMAPPEDAVERAGE, unmappedaverage)

inline void FORTRAN_NAME(inlineUNMAPPEDAVERAGE, inlineUNMAPPEDAVERAGE)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INTVECT(refRatio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_UNMAPPEDAVERAGE");
 FORTRAN_NAME( UNMAPPEDAVERAGE ,unmappedaverage )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INTVECT(refRatio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDUNMAPPEDAVERAGE 

#ifndef GUARDUNMAPPEDAVERAGEHARMONIC 
#define GUARDUNMAPPEDAVERAGEHARMONIC 
// Prototype for Fortran procedure UNMAPPEDAVERAGEHARMONIC ...
//
void FORTRAN_NAME( UNMAPPEDAVERAGEHARMONIC ,unmappedaverageharmonic )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INTVECT(refRatio)
      ,CHFp_BOX(bref) );

#define FORT_UNMAPPEDAVERAGEHARMONIC FORTRAN_NAME( inlineUNMAPPEDAVERAGEHARMONIC, inlineUNMAPPEDAVERAGEHARMONIC)
#define FORTNT_UNMAPPEDAVERAGEHARMONIC FORTRAN_NAME( UNMAPPEDAVERAGEHARMONIC, unmappedaverageharmonic)

inline void FORTRAN_NAME(inlineUNMAPPEDAVERAGEHARMONIC, inlineUNMAPPEDAVERAGEHARMONIC)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INTVECT(refRatio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_UNMAPPEDAVERAGEHARMONIC");
 FORTRAN_NAME( UNMAPPEDAVERAGEHARMONIC ,unmappedaverageharmonic )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INTVECT(refRatio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDUNMAPPEDAVERAGEHARMONIC 

#ifndef GUARDMAPPEDAVERAGE 
#define GUARDMAPPEDAVERAGE 
// Prototype for Fortran procedure MAPPEDAVERAGE ...
//
void FORTRAN_NAME( MAPPEDAVERAGE ,mappedaverage )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_CONST_FRA1(fineCCJ)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INTVECT(refRatio)
      ,CHFp_BOX(bref) );

#define FORT_MAPPEDAVERAGE FORTRAN_NAME( inlineMAPPEDAVERAGE, inlineMAPPEDAVERAGE)
#define FORTNT_MAPPEDAVERAGE FORTRAN_NAME( MAPPEDAVERAGE, mappedaverage)

inline void FORTRAN_NAME(inlineMAPPEDAVERAGE, inlineMAPPEDAVERAGE)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_CONST_FRA1(fineCCJ)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INTVECT(refRatio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_MAPPEDAVERAGE");
 FORTRAN_NAME( MAPPEDAVERAGE ,mappedaverage )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_CONST_FRA1(fineCCJ)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INTVECT(refRatio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDMAPPEDAVERAGE 

#ifndef GUARDMAPPEDAVERAGE2 
#define GUARDMAPPEDAVERAGE2 
// Prototype for Fortran procedure MAPPEDAVERAGE2 ...
//
void FORTRAN_NAME( MAPPEDAVERAGE2 ,mappedaverage2 )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_CONST_FRA1(fineCCJinv)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INTVECT(refRatio)
      ,CHFp_BOX(bref) );

#define FORT_MAPPEDAVERAGE2 FORTRAN_NAME( inlineMAPPEDAVERAGE2, inlineMAPPEDAVERAGE2)
#define FORTNT_MAPPEDAVERAGE2 FORTRAN_NAME( MAPPEDAVERAGE2, mappedaverage2)

inline void FORTRAN_NAME(inlineMAPPEDAVERAGE2, inlineMAPPEDAVERAGE2)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_CONST_FRA1(fineCCJinv)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INTVECT(refRatio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_MAPPEDAVERAGE2");
 FORTRAN_NAME( MAPPEDAVERAGE2 ,mappedaverage2 )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_CONST_FRA1(fineCCJinv)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INTVECT(refRatio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDMAPPEDAVERAGE2 

#ifndef GUARDUNMAPPEDAVERAGEFACE 
#define GUARDUNMAPPEDAVERAGEFACE 
// Prototype for Fortran procedure UNMAPPEDAVERAGEFACE ...
//
void FORTRAN_NAME( UNMAPPEDAVERAGEFACE ,unmappedaverageface )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INTVECT(refRatio)
      ,CHFp_BOX(refBox) );

#define FORT_UNMAPPEDAVERAGEFACE FORTRAN_NAME( inlineUNMAPPEDAVERAGEFACE, inlineUNMAPPEDAVERAGEFACE)
#define FORTNT_UNMAPPEDAVERAGEFACE FORTRAN_NAME( UNMAPPEDAVERAGEFACE, unmappedaverageface)

inline void FORTRAN_NAME(inlineUNMAPPEDAVERAGEFACE, inlineUNMAPPEDAVERAGEFACE)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INTVECT(refRatio)
      ,CHFp_BOX(refBox) )
{
 CH_TIMELEAF("FORT_UNMAPPEDAVERAGEFACE");
 FORTRAN_NAME( UNMAPPEDAVERAGEFACE ,unmappedaverageface )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(crseBox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INTVECT(refRatio)
      ,CHFt_BOX(refBox) );
}
#endif  // GUARDUNMAPPEDAVERAGEFACE 

#ifndef GUARDUNMAPPEDAVERAGEFACEHARMONIC 
#define GUARDUNMAPPEDAVERAGEFACEHARMONIC 
// Prototype for Fortran procedure UNMAPPEDAVERAGEFACEHARMONIC ...
//
void FORTRAN_NAME( UNMAPPEDAVERAGEFACEHARMONIC ,unmappedaveragefaceharmonic )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INTVECT(refRatio)
      ,CHFp_BOX(refBox) );

#define FORT_UNMAPPEDAVERAGEFACEHARMONIC FORTRAN_NAME( inlineUNMAPPEDAVERAGEFACEHARMONIC, inlineUNMAPPEDAVERAGEFACEHARMONIC)
#define FORTNT_UNMAPPEDAVERAGEFACEHARMONIC FORTRAN_NAME( UNMAPPEDAVERAGEFACEHARMONIC, unmappedaveragefaceharmonic)

inline void FORTRAN_NAME(inlineUNMAPPEDAVERAGEFACEHARMONIC, inlineUNMAPPEDAVERAGEFACEHARMONIC)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INTVECT(refRatio)
      ,CHFp_BOX(refBox) )
{
 CH_TIMELEAF("FORT_UNMAPPEDAVERAGEFACEHARMONIC");
 FORTRAN_NAME( UNMAPPEDAVERAGEFACEHARMONIC ,unmappedaveragefaceharmonic )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(crseBox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INTVECT(refRatio)
      ,CHFt_BOX(refBox) );
}
#endif  // GUARDUNMAPPEDAVERAGEFACEHARMONIC 

#ifndef GUARDMAPPEDAVERAGEFACE 
#define GUARDMAPPEDAVERAGEFACE 
// Prototype for Fortran procedure MAPPEDAVERAGEFACE ...
//
void FORTRAN_NAME( MAPPEDAVERAGEFACE ,mappedaverageface )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_CONST_FRA1(fineFCJ)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INTVECT(refRatio)
      ,CHFp_BOX(bref) );

#define FORT_MAPPEDAVERAGEFACE FORTRAN_NAME( inlineMAPPEDAVERAGEFACE, inlineMAPPEDAVERAGEFACE)
#define FORTNT_MAPPEDAVERAGEFACE FORTRAN_NAME( MAPPEDAVERAGEFACE, mappedaverageface)

inline void FORTRAN_NAME(inlineMAPPEDAVERAGEFACE, inlineMAPPEDAVERAGEFACE)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_CONST_FRA1(fineFCJ)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INTVECT(refRatio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_MAPPEDAVERAGEFACE");
 FORTRAN_NAME( MAPPEDAVERAGEFACE ,mappedaverageface )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_CONST_FRA1(fineFCJ)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INTVECT(refRatio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDMAPPEDAVERAGEFACE 

#ifndef GUARDINTERPTOCOARSENODES 
#define GUARDINTERPTOCOARSENODES 
// Prototype for Fortran procedure InterpToCoarseNodes ...
//
void FORTRAN_NAME( INTERPTOCOARSENODES ,interptocoarsenodes )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INTVECT(refRatio) );

#define FORT_INTERPTOCOARSENODES FORTRAN_NAME( inlineINTERPTOCOARSENODES, inlineINTERPTOCOARSENODES)
#define FORTNT_INTERPTOCOARSENODES FORTRAN_NAME( INTERPTOCOARSENODES, interptocoarsenodes)

inline void FORTRAN_NAME(inlineINTERPTOCOARSENODES, inlineINTERPTOCOARSENODES)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INTVECT(refRatio) )
{
 CH_TIMELEAF("FORT_INTERPTOCOARSENODES");
 FORTRAN_NAME( INTERPTOCOARSENODES ,interptocoarsenodes )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(crseBox)
      ,CHFt_CONST_INTVECT(refRatio) );
}
#endif  // GUARDINTERPTOCOARSENODES 

}

#endif
