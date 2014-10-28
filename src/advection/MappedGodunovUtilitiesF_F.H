#ifndef _MAPPEDGODUNOVUTILITIESF_F_H_
#define _MAPPEDGODUNOVUTILITIESF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSECONDSLOPEDIFFSF 
#define GUARDSECONDSLOPEDIFFSF 
// Prototype for Fortran procedure SecondSlopeDiffsF ...
//
void FORTRAN_NAME( SECONDSLOPEDIFFSF ,secondslopediffsf )(
      CHFp_FRA(deltaWC)
      ,CHFp_FRA(deltaWL)
      ,CHFp_FRA(deltaWR)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_SECONDSLOPEDIFFSF FORTRAN_NAME( inlineSECONDSLOPEDIFFSF, inlineSECONDSLOPEDIFFSF)
#define FORTNT_SECONDSLOPEDIFFSF FORTRAN_NAME( SECONDSLOPEDIFFSF, secondslopediffsf)

inline void FORTRAN_NAME(inlineSECONDSLOPEDIFFSF, inlineSECONDSLOPEDIFFSF)(
      CHFp_FRA(deltaWC)
      ,CHFp_FRA(deltaWL)
      ,CHFp_FRA(deltaWR)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_SECONDSLOPEDIFFSF");
 FORTRAN_NAME( SECONDSLOPEDIFFSF ,secondslopediffsf )(
      CHFt_FRA(deltaWC)
      ,CHFt_FRA(deltaWL)
      ,CHFt_FRA(deltaWR)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDSECONDSLOPEDIFFSF 

#ifndef GUARDFOURTHSLOPEDIFFSF 
#define GUARDFOURTHSLOPEDIFFSF 
// Prototype for Fortran procedure FourthSlopeDiffsF ...
//
void FORTRAN_NAME( FOURTHSLOPEDIFFSF ,fourthslopediffsf )(
      CHFp_FRA(delta4WC)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(delta2W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_FOURTHSLOPEDIFFSF FORTRAN_NAME( inlineFOURTHSLOPEDIFFSF, inlineFOURTHSLOPEDIFFSF)
#define FORTNT_FOURTHSLOPEDIFFSF FORTRAN_NAME( FOURTHSLOPEDIFFSF, fourthslopediffsf)

inline void FORTRAN_NAME(inlineFOURTHSLOPEDIFFSF, inlineFOURTHSLOPEDIFFSF)(
      CHFp_FRA(delta4WC)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(delta2W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_FOURTHSLOPEDIFFSF");
 FORTRAN_NAME( FOURTHSLOPEDIFFSF ,fourthslopediffsf )(
      CHFt_FRA(delta4WC)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(delta2W)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDFOURTHSLOPEDIFFSF 

#ifndef GUARDVANLEERLIMITERF 
#define GUARDVANLEERLIMITERF 
// Prototype for Fortran procedure VanLeerLimiterF ...
//
void FORTRAN_NAME( VANLEERLIMITERF ,vanleerlimiterf )(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(dWleft)
      ,CHFp_CONST_FRA(dWright)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_BOX(box) );

#define FORT_VANLEERLIMITERF FORTRAN_NAME( inlineVANLEERLIMITERF, inlineVANLEERLIMITERF)
#define FORTNT_VANLEERLIMITERF FORTRAN_NAME( VANLEERLIMITERF, vanleerlimiterf)

inline void FORTRAN_NAME(inlineVANLEERLIMITERF, inlineVANLEERLIMITERF)(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(dWleft)
      ,CHFp_CONST_FRA(dWright)
      ,CHFp_CONST_INT(numslopes)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_VANLEERLIMITERF");
 FORTRAN_NAME( VANLEERLIMITERF ,vanleerlimiterf )(
      CHFt_FRA(dW)
      ,CHFt_CONST_FRA(dWleft)
      ,CHFt_CONST_FRA(dWright)
      ,CHFt_CONST_INT(numslopes)
      ,CHFt_BOX(box) );
}
#endif  // GUARDVANLEERLIMITERF 

#ifndef GUARDEXTPRESERVINGVANLEERLIMITERF 
#define GUARDEXTPRESERVINGVANLEERLIMITERF 
// Prototype for Fortran procedure ExtPreservingVanLeerLimiterF ...
//
void FORTRAN_NAME( EXTPRESERVINGVANLEERLIMITERF ,extpreservingvanleerlimiterf )(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(dWleft)
      ,CHFp_CONST_FRA(dWright)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox) );

#define FORT_EXTPRESERVINGVANLEERLIMITERF FORTRAN_NAME( inlineEXTPRESERVINGVANLEERLIMITERF, inlineEXTPRESERVINGVANLEERLIMITERF)
#define FORTNT_EXTPRESERVINGVANLEERLIMITERF FORTRAN_NAME( EXTPRESERVINGVANLEERLIMITERF, extpreservingvanleerlimiterf)

inline void FORTRAN_NAME(inlineEXTPRESERVINGVANLEERLIMITERF, inlineEXTPRESERVINGVANLEERLIMITERF)(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(dWleft)
      ,CHFp_CONST_FRA(dWright)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_EXTPRESERVINGVANLEERLIMITERF");
 FORTRAN_NAME( EXTPRESERVINGVANLEERLIMITERF ,extpreservingvanleerlimiterf )(
      CHFt_FRA(dW)
      ,CHFt_CONST_FRA(dWleft)
      ,CHFt_CONST_FRA(dWright)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDEXTPRESERVINGVANLEERLIMITERF 

#ifndef GUARDPLMNORMALPREDF 
#define GUARDPLMNORMALPREDF 
// Prototype for Fortran procedure PLMNormalPredF ...
//
void FORTRAN_NAME( PLMNORMALPREDF ,plmnormalpredf )(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_FRA(dW)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(nSlope)
      ,CHFp_BOX(box) );

#define FORT_PLMNORMALPREDF FORTRAN_NAME( inlinePLMNORMALPREDF, inlinePLMNORMALPREDF)
#define FORTNT_PLMNORMALPREDF FORTRAN_NAME( PLMNORMALPREDF, plmnormalpredf)

inline void FORTRAN_NAME(inlinePLMNORMALPREDF, inlinePLMNORMALPREDF)(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_FRA(dW)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(nSlope)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PLMNORMALPREDF");
 FORTRAN_NAME( PLMNORMALPREDF ,plmnormalpredf )(
      CHFt_FRA(dWMinus)
      ,CHFt_FRA(dWPlus)
      ,CHFt_CONST_FRA(dW)
      ,CHFt_CONST_FRA(lambda)
      ,CHFt_CONST_REAL(dtbydx)
      ,CHFt_CONST_INT(nSlope)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPLMNORMALPREDF 

#ifndef GUARDPPMLIMITERF 
#define GUARDPPMLIMITERF 
// Prototype for Fortran procedure PPMLimiterF ...
//
void FORTRAN_NAME( PPMLIMITERF ,ppmlimiterf )(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_BOX(box) );

#define FORT_PPMLIMITERF FORTRAN_NAME( inlinePPMLIMITERF, inlinePPMLIMITERF)
#define FORTNT_PPMLIMITERF FORTRAN_NAME( PPMLIMITERF, ppmlimiterf)

inline void FORTRAN_NAME(inlinePPMLIMITERF, inlinePPMLIMITERF)(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PPMLIMITERF");
 FORTRAN_NAME( PPMLIMITERF ,ppmlimiterf )(
      CHFt_FRA(dWMinus)
      ,CHFt_FRA(dWPlus)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPPMLIMITERF 

#ifndef GUARDPPMFACEVALUESF 
#define GUARDPPMFACEVALUESF 
// Prototype for Fortran procedure PPMFaceValuesF ...
//
void FORTRAN_NAME( PPMFACEVALUESF ,ppmfacevaluesf )(
      CHFp_FRA(WFace)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(dW)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_PPMFACEVALUESF FORTRAN_NAME( inlinePPMFACEVALUESF, inlinePPMFACEVALUESF)
#define FORTNT_PPMFACEVALUESF FORTRAN_NAME( PPMFACEVALUESF, ppmfacevaluesf)

inline void FORTRAN_NAME(inlinePPMFACEVALUESF, inlinePPMFACEVALUESF)(
      CHFp_FRA(WFace)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(dW)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_PPMFACEVALUESF");
 FORTRAN_NAME( PPMFACEVALUESF ,ppmfacevaluesf )(
      CHFt_FRA(WFace)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(dW)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDPPMFACEVALUESF 

#ifndef GUARDPPMNORMALPREDF 
#define GUARDPPMNORMALPREDF 
// Prototype for Fortran procedure PPMNormalPredF ...
//
void FORTRAN_NAME( PPMNORMALPREDF ,ppmnormalpredf )(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(nSlope)
      ,CHFp_BOX(box) );

#define FORT_PPMNORMALPREDF FORTRAN_NAME( inlinePPMNORMALPREDF, inlinePPMNORMALPREDF)
#define FORTNT_PPMNORMALPREDF FORTRAN_NAME( PPMNORMALPREDF, ppmnormalpredf)

inline void FORTRAN_NAME(inlinePPMNORMALPREDF, inlinePPMNORMALPREDF)(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_FRA(lambda)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_CONST_INT(nSlope)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PPMNORMALPREDF");
 FORTRAN_NAME( PPMNORMALPREDF ,ppmnormalpredf )(
      CHFt_FRA(dWMinus)
      ,CHFt_FRA(dWPlus)
      ,CHFt_CONST_FRA(lambda)
      ,CHFt_CONST_REAL(dtbydx)
      ,CHFt_CONST_INT(nSlope)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPPMNORMALPREDF 

#ifndef GUARDGETSECONDDIFF 
#define GUARDGETSECONDDIFF 
// Prototype for Fortran procedure GETSECONDDIFF ...
//
void FORTRAN_NAME( GETSECONDDIFF ,getseconddiff )(
      CHFp_FRA(d2W)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_GETSECONDDIFF FORTRAN_NAME( inlineGETSECONDDIFF, inlineGETSECONDDIFF)
#define FORTNT_GETSECONDDIFF FORTRAN_NAME( GETSECONDDIFF, getseconddiff)

inline void FORTRAN_NAME(inlineGETSECONDDIFF, inlineGETSECONDDIFF)(
      CHFp_FRA(d2W)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_GETSECONDDIFF");
 FORTRAN_NAME( GETSECONDDIFF ,getseconddiff )(
      CHFt_FRA(d2W)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDGETSECONDDIFF 

#ifndef GUARDFOURTHINTERPFACES 
#define GUARDFOURTHINTERPFACES 
// Prototype for Fortran procedure FourthInterpFaces ...
//
void FORTRAN_NAME( FOURTHINTERPFACES ,fourthinterpfaces )(
      CHFp_FRA(Wface)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_BOX(nextLoBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_BOX(nextHiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_FOURTHINTERPFACES FORTRAN_NAME( inlineFOURTHINTERPFACES, inlineFOURTHINTERPFACES)
#define FORTNT_FOURTHINTERPFACES FORTRAN_NAME( FOURTHINTERPFACES, fourthinterpfaces)

inline void FORTRAN_NAME(inlineFOURTHINTERPFACES, inlineFOURTHINTERPFACES)(
      CHFp_FRA(Wface)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_BOX(nextLoBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_BOX(nextHiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_FOURTHINTERPFACES");
 FORTRAN_NAME( FOURTHINTERPFACES ,fourthinterpfaces )(
      CHFt_FRA(Wface)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_BOX(nextLoBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_BOX(nextHiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDFOURTHINTERPFACES 

#ifndef GUARDCHECKCUBICLIMITERF 
#define GUARDCHECKCUBICLIMITERF 
// Prototype for Fortran procedure CheckCubicLimiterF ...
//
void FORTRAN_NAME( CHECKCUBICLIMITERF ,checkcubiclimiterf )(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(d2W)
      ,CHFp_CONST_FRA(dW2fcf)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_BOX(nextLoBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_BOX(nextHiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(innerCenterBox)
      ,CHFp_CONST_REAL(limitC)
      ,CHFp_CONST_REAL(C3)
      ,CHFp_CONST_REAL(eps) );

#define FORT_CHECKCUBICLIMITERF FORTRAN_NAME( inlineCHECKCUBICLIMITERF, inlineCHECKCUBICLIMITERF)
#define FORTNT_CHECKCUBICLIMITERF FORTRAN_NAME( CHECKCUBICLIMITERF, checkcubiclimiterf)

inline void FORTRAN_NAME(inlineCHECKCUBICLIMITERF, inlineCHECKCUBICLIMITERF)(
      CHFp_FRA(dWMinus)
      ,CHFp_FRA(dWPlus)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(d2W)
      ,CHFp_CONST_FRA(dW2fcf)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_BOX(nextLoBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_BOX(nextHiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(innerCenterBox)
      ,CHFp_CONST_REAL(limitC)
      ,CHFp_CONST_REAL(C3)
      ,CHFp_CONST_REAL(eps) )
{
 CH_TIMELEAF("FORT_CHECKCUBICLIMITERF");
 FORTRAN_NAME( CHECKCUBICLIMITERF ,checkcubiclimiterf )(
      CHFt_FRA(dWMinus)
      ,CHFt_FRA(dWPlus)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(d2W)
      ,CHFt_CONST_FRA(dW2fcf)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_BOX(nextLoBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_BOX(nextHiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(innerCenterBox)
      ,CHFt_CONST_REAL(limitC)
      ,CHFt_CONST_REAL(C3)
      ,CHFt_CONST_REAL(eps) );
}
#endif  // GUARDCHECKCUBICLIMITERF 

}

#endif
