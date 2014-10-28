#ifndef _EXTRAPOLATIONUTILSF_F_H_
#define _EXTRAPOLATIONUTILSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDEXTRAPOLATEFACENOEV 
#define GUARDEXTRAPOLATEFACENOEV 
// Prototype for Fortran procedure EXTRAPOLATEFACENOEV ...
//
void FORTRAN_NAME( EXTRAPOLATEFACENOEV ,extrapolatefacenoev )(
      CHFp_FRA(dest)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(toRegion)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(sidesign)
      ,CHFp_CONST_INT(order) );

#define FORT_EXTRAPOLATEFACENOEV FORTRAN_NAME( inlineEXTRAPOLATEFACENOEV, inlineEXTRAPOLATEFACENOEV)
#define FORTNT_EXTRAPOLATEFACENOEV FORTRAN_NAME( EXTRAPOLATEFACENOEV, extrapolatefacenoev)

inline void FORTRAN_NAME(inlineEXTRAPOLATEFACENOEV, inlineEXTRAPOLATEFACENOEV)(
      CHFp_FRA(dest)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(toRegion)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(sidesign)
      ,CHFp_CONST_INT(order) )
{
 CH_TIMELEAF("FORT_EXTRAPOLATEFACENOEV");
 FORTRAN_NAME( EXTRAPOLATEFACENOEV ,extrapolatefacenoev )(
      CHFt_FRA(dest)
      ,CHFt_CONST_FRA(src)
      ,CHFt_BOX(toRegion)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(sidesign)
      ,CHFt_CONST_INT(order) );
}
#endif  // GUARDEXTRAPOLATEFACENOEV 

#ifndef GUARDEXTRAPOLATEDIAGEDGEF 
#define GUARDEXTRAPOLATEDIAGEDGEF 
// Prototype for Fortran procedure ExtrapolateDiagEdgeF ...
//
void FORTRAN_NAME( EXTRAPOLATEDIAGEDGEF ,extrapolatediagedgef )(
      CHFp_FRA(dest)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(toRegion)
      ,CHFp_CONST_INT(dir1)
      ,CHFp_CONST_INT(sidesign1)
      ,CHFp_CONST_INT(dir2)
      ,CHFp_CONST_INT(sidesign2)
      ,CHFp_CONST_INT(order) );

#define FORT_EXTRAPOLATEDIAGEDGEF FORTRAN_NAME( inlineEXTRAPOLATEDIAGEDGEF, inlineEXTRAPOLATEDIAGEDGEF)
#define FORTNT_EXTRAPOLATEDIAGEDGEF FORTRAN_NAME( EXTRAPOLATEDIAGEDGEF, extrapolatediagedgef)

inline void FORTRAN_NAME(inlineEXTRAPOLATEDIAGEDGEF, inlineEXTRAPOLATEDIAGEDGEF)(
      CHFp_FRA(dest)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(toRegion)
      ,CHFp_CONST_INT(dir1)
      ,CHFp_CONST_INT(sidesign1)
      ,CHFp_CONST_INT(dir2)
      ,CHFp_CONST_INT(sidesign2)
      ,CHFp_CONST_INT(order) )
{
 CH_TIMELEAF("FORT_EXTRAPOLATEDIAGEDGEF");
 FORTRAN_NAME( EXTRAPOLATEDIAGEDGEF ,extrapolatediagedgef )(
      CHFt_FRA(dest)
      ,CHFt_CONST_FRA(src)
      ,CHFt_BOX(toRegion)
      ,CHFt_CONST_INT(dir1)
      ,CHFt_CONST_INT(sidesign1)
      ,CHFt_CONST_INT(dir2)
      ,CHFt_CONST_INT(sidesign2)
      ,CHFt_CONST_INT(order) );
}
#endif  // GUARDEXTRAPOLATEDIAGEDGEF 

#ifndef GUARDEXTRAPOLATEFACEFROMCC 
#define GUARDEXTRAPOLATEFACEFROMCC 
// Prototype for Fortran procedure EXTRAPOLATEFACEFROMCC ...
//
void FORTRAN_NAME( EXTRAPOLATEFACEFROMCC ,extrapolatefacefromcc )(
      CHFp_FRA(dest)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(toRegion)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(sidesign)
      ,CHFp_CONST_INT(order) );

#define FORT_EXTRAPOLATEFACEFROMCC FORTRAN_NAME( inlineEXTRAPOLATEFACEFROMCC, inlineEXTRAPOLATEFACEFROMCC)
#define FORTNT_EXTRAPOLATEFACEFROMCC FORTRAN_NAME( EXTRAPOLATEFACEFROMCC, extrapolatefacefromcc)

inline void FORTRAN_NAME(inlineEXTRAPOLATEFACEFROMCC, inlineEXTRAPOLATEFACEFROMCC)(
      CHFp_FRA(dest)
      ,CHFp_CONST_FRA(src)
      ,CHFp_BOX(toRegion)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(sidesign)
      ,CHFp_CONST_INT(order) )
{
 CH_TIMELEAF("FORT_EXTRAPOLATEFACEFROMCC");
 FORTRAN_NAME( EXTRAPOLATEFACEFROMCC ,extrapolatefacefromcc )(
      CHFt_FRA(dest)
      ,CHFt_CONST_FRA(src)
      ,CHFt_BOX(toRegion)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(sidesign)
      ,CHFt_CONST_INT(order) );
}
#endif  // GUARDEXTRAPOLATEFACEFROMCC 

}

#endif
