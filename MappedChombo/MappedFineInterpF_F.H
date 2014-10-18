#ifndef _MAPPEDFINEINTERPF_F_H_
#define _MAPPEDFINEINTERPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDUNMAPPEDINTERPCONSTANT 
#define GUARDUNMAPPEDINTERPCONSTANT 
// Prototype for Fortran procedure unmappedinterpconstant ...
//
void FORTRAN_NAME( UNMAPPEDINTERPCONSTANT ,unmappedinterpconstant )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(coarseBox)
      ,CHFp_CONST_INTVECT(ref_ratio)
      ,CHFp_BOX(bref) );

#define FORT_UNMAPPEDINTERPCONSTANT FORTRAN_NAME( inlineUNMAPPEDINTERPCONSTANT, inlineUNMAPPEDINTERPCONSTANT)
#define FORTNT_UNMAPPEDINTERPCONSTANT FORTRAN_NAME( UNMAPPEDINTERPCONSTANT, unmappedinterpconstant)

inline void FORTRAN_NAME(inlineUNMAPPEDINTERPCONSTANT, inlineUNMAPPEDINTERPCONSTANT)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(coarseBox)
      ,CHFp_CONST_INTVECT(ref_ratio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_UNMAPPEDINTERPCONSTANT");
 FORTRAN_NAME( UNMAPPEDINTERPCONSTANT ,unmappedinterpconstant )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_BOX(coarseBox)
      ,CHFt_CONST_INTVECT(ref_ratio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDUNMAPPEDINTERPCONSTANT 

#ifndef GUARDUNMAPPEDINTERPLINEAR 
#define GUARDUNMAPPEDINTERPLINEAR 
// Prototype for Fortran procedure unmappedinterplinear ...
//
void FORTRAN_NAME( UNMAPPEDINTERPLINEAR ,unmappedinterplinear )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(slope)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INTVECT(ref_ratio)
      ,CHFp_BOX(bref) );

#define FORT_UNMAPPEDINTERPLINEAR FORTRAN_NAME( inlineUNMAPPEDINTERPLINEAR, inlineUNMAPPEDINTERPLINEAR)
#define FORTNT_UNMAPPEDINTERPLINEAR FORTRAN_NAME( UNMAPPEDINTERPLINEAR, unmappedinterplinear)

inline void FORTRAN_NAME(inlineUNMAPPEDINTERPLINEAR, inlineUNMAPPEDINTERPLINEAR)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(slope)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INTVECT(ref_ratio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_UNMAPPEDINTERPLINEAR");
 FORTRAN_NAME( UNMAPPEDINTERPLINEAR ,unmappedinterplinear )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(slope)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INTVECT(ref_ratio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDUNMAPPEDINTERPLINEAR 

#ifndef GUARDUNMAPPEDINTERPFACECONSTANT 
#define GUARDUNMAPPEDINTERPFACECONSTANT 
// Prototype for Fortran procedure UnmappedInterpFaceConstant ...
//
void FORTRAN_NAME( UNMAPPEDINTERPFACECONSTANT ,unmappedinterpfaceconstant )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INTVECT(ref_ratio)
      ,CHFp_BOX(bref)
      ,CHFp_CONST_INT(dir) );

#define FORT_UNMAPPEDINTERPFACECONSTANT FORTRAN_NAME( inlineUNMAPPEDINTERPFACECONSTANT, inlineUNMAPPEDINTERPFACECONSTANT)
#define FORTNT_UNMAPPEDINTERPFACECONSTANT FORTRAN_NAME( UNMAPPEDINTERPFACECONSTANT, unmappedinterpfaceconstant)

inline void FORTRAN_NAME(inlineUNMAPPEDINTERPFACECONSTANT, inlineUNMAPPEDINTERPFACECONSTANT)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(coarse)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INTVECT(ref_ratio)
      ,CHFp_BOX(bref)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_UNMAPPEDINTERPFACECONSTANT");
 FORTRAN_NAME( UNMAPPEDINTERPFACECONSTANT ,unmappedinterpfaceconstant )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(coarse)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INTVECT(ref_ratio)
      ,CHFt_BOX(bref)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDUNMAPPEDINTERPFACECONSTANT 

#ifndef GUARDUNMAPPEDINTERPLINEARFACE 
#define GUARDUNMAPPEDINTERPLINEARFACE 
// Prototype for Fortran procedure UnmappedInterpLinearFace ...
//
void FORTRAN_NAME( UNMAPPEDINTERPLINEARFACE ,unmappedinterplinearface )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(slope)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INTVECT(ref_ratio)
      ,CHFp_BOX(brefface) );

#define FORT_UNMAPPEDINTERPLINEARFACE FORTRAN_NAME( inlineUNMAPPEDINTERPLINEARFACE, inlineUNMAPPEDINTERPLINEARFACE)
#define FORTNT_UNMAPPEDINTERPLINEARFACE FORTRAN_NAME( UNMAPPEDINTERPLINEARFACE, unmappedinterplinearface)

inline void FORTRAN_NAME(inlineUNMAPPEDINTERPLINEARFACE, inlineUNMAPPEDINTERPLINEARFACE)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(slope)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INTVECT(ref_ratio)
      ,CHFp_BOX(brefface) )
{
 CH_TIMELEAF("FORT_UNMAPPEDINTERPLINEARFACE");
 FORTRAN_NAME( UNMAPPEDINTERPLINEARFACE ,unmappedinterplinearface )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(slope)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INTVECT(ref_ratio)
      ,CHFt_BOX(brefface) );
}
#endif  // GUARDUNMAPPEDINTERPLINEARFACE 

#ifndef GUARDUNMAPPEDINTERPLINEARINTERIORFACE 
#define GUARDUNMAPPEDINTERPLINEARINTERIORFACE 
// Prototype for Fortran procedure UnmappedInterpLinearInteriorFace ...
//
void FORTRAN_NAME( UNMAPPEDINTERPLINEARINTERIORFACE ,unmappedinterplinearinteriorface )(
      CHFp_FRA(fine)
      ,CHFp_BOX(bcoarse)
      ,CHFp_CONST_INTVECT(ref_ratio)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(interiorrefbox) );

#define FORT_UNMAPPEDINTERPLINEARINTERIORFACE FORTRAN_NAME( inlineUNMAPPEDINTERPLINEARINTERIORFACE, inlineUNMAPPEDINTERPLINEARINTERIORFACE)
#define FORTNT_UNMAPPEDINTERPLINEARINTERIORFACE FORTRAN_NAME( UNMAPPEDINTERPLINEARINTERIORFACE, unmappedinterplinearinteriorface)

inline void FORTRAN_NAME(inlineUNMAPPEDINTERPLINEARINTERIORFACE, inlineUNMAPPEDINTERPLINEARINTERIORFACE)(
      CHFp_FRA(fine)
      ,CHFp_BOX(bcoarse)
      ,CHFp_CONST_INTVECT(ref_ratio)
      ,CHFp_CONST_INT(facedir)
      ,CHFp_BOX(interiorrefbox) )
{
 CH_TIMELEAF("FORT_UNMAPPEDINTERPLINEARINTERIORFACE");
 FORTRAN_NAME( UNMAPPEDINTERPLINEARINTERIORFACE ,unmappedinterplinearinteriorface )(
      CHFt_FRA(fine)
      ,CHFt_BOX(bcoarse)
      ,CHFt_CONST_INTVECT(ref_ratio)
      ,CHFt_CONST_INT(facedir)
      ,CHFt_BOX(interiorrefbox) );
}
#endif  // GUARDUNMAPPEDINTERPLINEARINTERIORFACE 

#ifndef GUARDUNMAPPEDINTERPLIMITFACE 
#define GUARDUNMAPPEDINTERPLIMITFACE 
// Prototype for Fortran procedure UnmappedInterpLimitFace ...
//
void FORTRAN_NAME( UNMAPPEDINTERPLIMITFACE ,unmappedinterplimitface )(
      CHFp_FRA(islope)
      ,CHFp_FRA(jslope)
      ,CHFp_FRA(kslope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_BOX(bn)
      ,CHFp_BOX(validBox)
      ,CHFp_CONST_INT(normaldir) );

#define FORT_UNMAPPEDINTERPLIMITFACE FORTRAN_NAME( inlineUNMAPPEDINTERPLIMITFACE, inlineUNMAPPEDINTERPLIMITFACE)
#define FORTNT_UNMAPPEDINTERPLIMITFACE FORTRAN_NAME( UNMAPPEDINTERPLIMITFACE, unmappedinterplimitface)

inline void FORTRAN_NAME(inlineUNMAPPEDINTERPLIMITFACE, inlineUNMAPPEDINTERPLIMITFACE)(
      CHFp_FRA(islope)
      ,CHFp_FRA(jslope)
      ,CHFp_FRA(kslope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_BOX(bn)
      ,CHFp_BOX(validBox)
      ,CHFp_CONST_INT(normaldir) )
{
 CH_TIMELEAF("FORT_UNMAPPEDINTERPLIMITFACE");
 FORTRAN_NAME( UNMAPPEDINTERPLIMITFACE ,unmappedinterplimitface )(
      CHFt_FRA(islope)
      ,CHFt_FRA(jslope)
      ,CHFt_FRA(kslope)
      ,CHFt_CONST_FRA(state)
      ,CHFt_BOX(b)
      ,CHFt_BOX(bn)
      ,CHFt_BOX(validBox)
      ,CHFt_CONST_INT(normaldir) );
}
#endif  // GUARDUNMAPPEDINTERPLIMITFACE 

}

#endif
