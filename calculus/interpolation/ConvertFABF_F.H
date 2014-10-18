#ifndef _CONVERTFABF_F_H_
#define _CONVERTFABF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCONVERTFAB 
#define GUARDCONVERTFAB 
// Prototype for Fortran procedure CONVERTFAB ...
//
void FORTRAN_NAME( CONVERTFAB ,convertfab )(
      CHFp_FRA1(dest)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType)
      ,CHFp_CONST_FRA1(src)
      ,CHFp_CONST_INTVECT(srcBoxType) );

#define FORT_CONVERTFAB FORTRAN_NAME( inlineCONVERTFAB, inlineCONVERTFAB)
#define FORTNT_CONVERTFAB FORTRAN_NAME( CONVERTFAB, convertfab)

inline void FORTRAN_NAME(inlineCONVERTFAB, inlineCONVERTFAB)(
      CHFp_FRA1(dest)
      ,CHFp_BOX(destBox)
      ,CHFp_CONST_INTVECT(destBoxType)
      ,CHFp_CONST_FRA1(src)
      ,CHFp_CONST_INTVECT(srcBoxType) )
{
 CH_TIMELEAF("FORT_CONVERTFAB");
 FORTRAN_NAME( CONVERTFAB ,convertfab )(
      CHFt_FRA1(dest)
      ,CHFt_BOX(destBox)
      ,CHFt_CONST_INTVECT(destBoxType)
      ,CHFt_CONST_FRA1(src)
      ,CHFt_CONST_INTVECT(srcBoxType) );
}
#endif  // GUARDCONVERTFAB 

#ifndef GUARDCELLTOEDGE4TH 
#define GUARDCELLTOEDGE4TH 
// Prototype for Fortran procedure CellToEdge4th ...
//
void FORTRAN_NAME( CELLTOEDGE4TH ,celltoedge4th )(
      CHFp_FRA(Wface)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_CELLTOEDGE4TH FORTRAN_NAME( inlineCELLTOEDGE4TH, inlineCELLTOEDGE4TH)
#define FORTNT_CELLTOEDGE4TH FORTRAN_NAME( CELLTOEDGE4TH, celltoedge4th)

inline void FORTRAN_NAME(inlineCELLTOEDGE4TH, inlineCELLTOEDGE4TH)(
      CHFp_FRA(Wface)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_CELLTOEDGE4TH");
 FORTRAN_NAME( CELLTOEDGE4TH ,celltoedge4th )(
      CHFt_FRA(Wface)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDCELLTOEDGE4TH 

}

#endif
