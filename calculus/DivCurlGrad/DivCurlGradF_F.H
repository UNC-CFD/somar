#ifndef _DIVCURLGRADF_F_H_
#define _DIVCURLGRADF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMAPPEDMACGRAD 
#define GUARDMAPPEDMACGRAD 
// Prototype for Fortran procedure MAPPEDMACGRAD ...
//
void FORTRAN_NAME( MAPPEDMACGRAD ,mappedmacgrad )(
      CHFp_FRA1(edgeGrad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA(Jga)
      ,CHFp_BOX(edgeGrid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_INT(gradDir)
      ,CHFp_INT(edgeDir) );

#define FORT_MAPPEDMACGRAD FORTRAN_NAME( inlineMAPPEDMACGRAD, inlineMAPPEDMACGRAD)
#define FORTNT_MAPPEDMACGRAD FORTRAN_NAME( MAPPEDMACGRAD, mappedmacgrad)

inline void FORTRAN_NAME(inlineMAPPEDMACGRAD, inlineMAPPEDMACGRAD)(
      CHFp_FRA1(edgeGrad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA(Jga)
      ,CHFp_BOX(edgeGrid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_INT(gradDir)
      ,CHFp_INT(edgeDir) )
{
 CH_TIMELEAF("FORT_MAPPEDMACGRAD");
 FORTRAN_NAME( MAPPEDMACGRAD ,mappedmacgrad )(
      CHFt_FRA1(edgeGrad)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(extrap)
      ,CHFt_CONST_FRA(Jga)
      ,CHFt_BOX(edgeGrid)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_INT(gradDir)
      ,CHFt_INT(edgeDir) );
}
#endif  // GUARDMAPPEDMACGRAD 

#ifndef GUARDMAPPEDMACGRADORTHO 
#define GUARDMAPPEDMACGRADORTHO 
// Prototype for Fortran procedure MAPPEDMACGRADORTHO ...
//
void FORTRAN_NAME( MAPPEDMACGRADORTHO ,mappedmacgradortho )(
      CHFp_FRA1(edgeGrad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA(Jga)
      ,CHFp_BOX(edgeGrid)
      ,CHFp_CONST_REAL(dxDir)
      ,CHFp_INT(gradDir)
      ,CHFp_INT(edgeDir) );

#define FORT_MAPPEDMACGRADORTHO FORTRAN_NAME( inlineMAPPEDMACGRADORTHO, inlineMAPPEDMACGRADORTHO)
#define FORTNT_MAPPEDMACGRADORTHO FORTRAN_NAME( MAPPEDMACGRADORTHO, mappedmacgradortho)

inline void FORTRAN_NAME(inlineMAPPEDMACGRADORTHO, inlineMAPPEDMACGRADORTHO)(
      CHFp_FRA1(edgeGrad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA(Jga)
      ,CHFp_BOX(edgeGrid)
      ,CHFp_CONST_REAL(dxDir)
      ,CHFp_INT(gradDir)
      ,CHFp_INT(edgeDir) )
{
 CH_TIMELEAF("FORT_MAPPEDMACGRADORTHO");
 FORTRAN_NAME( MAPPEDMACGRADORTHO ,mappedmacgradortho )(
      CHFt_FRA1(edgeGrad)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(extrap)
      ,CHFt_CONST_FRA(Jga)
      ,CHFt_BOX(edgeGrid)
      ,CHFt_CONST_REAL(dxDir)
      ,CHFt_INT(gradDir)
      ,CHFt_INT(edgeDir) );
}
#endif  // GUARDMAPPEDMACGRADORTHO 

#ifndef GUARDMAPPEDMACGRADSCALE 
#define GUARDMAPPEDMACGRADSCALE 
// Prototype for Fortran procedure MappedMACGradScale ...
//
void FORTRAN_NAME( MAPPEDMACGRADSCALE ,mappedmacgradscale )(
      CHFp_FRA1(edgeGrad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA(Jga)
      ,CHFp_BOX(edgeGrid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(pdScale)
      ,CHFp_CONST_INT(gradDir)
      ,CHFp_CONST_INT(edgeDir) );

#define FORT_MAPPEDMACGRADSCALE FORTRAN_NAME( inlineMAPPEDMACGRADSCALE, inlineMAPPEDMACGRADSCALE)
#define FORTNT_MAPPEDMACGRADSCALE FORTRAN_NAME( MAPPEDMACGRADSCALE, mappedmacgradscale)

inline void FORTRAN_NAME(inlineMAPPEDMACGRADSCALE, inlineMAPPEDMACGRADSCALE)(
      CHFp_FRA1(edgeGrad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA(Jga)
      ,CHFp_BOX(edgeGrid)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(pdScale)
      ,CHFp_CONST_INT(gradDir)
      ,CHFp_CONST_INT(edgeDir) )
{
 CH_TIMELEAF("FORT_MAPPEDMACGRADSCALE");
 FORTRAN_NAME( MAPPEDMACGRADSCALE ,mappedmacgradscale )(
      CHFt_FRA1(edgeGrad)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(extrap)
      ,CHFt_CONST_FRA(Jga)
      ,CHFt_BOX(edgeGrid)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(pdScale)
      ,CHFt_CONST_INT(gradDir)
      ,CHFt_CONST_INT(edgeDir) );
}
#endif  // GUARDMAPPEDMACGRADSCALE 

#ifndef GUARDMAPPEDMACGRADORTHOSCALE 
#define GUARDMAPPEDMACGRADORTHOSCALE 
// Prototype for Fortran procedure MappedMACGradOrthoScale ...
//
void FORTRAN_NAME( MAPPEDMACGRADORTHOSCALE ,mappedmacgradorthoscale )(
      CHFp_FRA1(edgeGrad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA(Jga)
      ,CHFp_BOX(edgeGrid)
      ,CHFp_CONST_REAL(dxDir)
      ,CHFp_CONST_REAL(pdScaleDir)
      ,CHFp_CONST_INT(gradDir)
      ,CHFp_CONST_INT(edgeDir) );

#define FORT_MAPPEDMACGRADORTHOSCALE FORTRAN_NAME( inlineMAPPEDMACGRADORTHOSCALE, inlineMAPPEDMACGRADORTHOSCALE)
#define FORTNT_MAPPEDMACGRADORTHOSCALE FORTRAN_NAME( MAPPEDMACGRADORTHOSCALE, mappedmacgradorthoscale)

inline void FORTRAN_NAME(inlineMAPPEDMACGRADORTHOSCALE, inlineMAPPEDMACGRADORTHOSCALE)(
      CHFp_FRA1(edgeGrad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA(Jga)
      ,CHFp_BOX(edgeGrid)
      ,CHFp_CONST_REAL(dxDir)
      ,CHFp_CONST_REAL(pdScaleDir)
      ,CHFp_CONST_INT(gradDir)
      ,CHFp_CONST_INT(edgeDir) )
{
 CH_TIMELEAF("FORT_MAPPEDMACGRADORTHOSCALE");
 FORTRAN_NAME( MAPPEDMACGRADORTHOSCALE ,mappedmacgradorthoscale )(
      CHFt_FRA1(edgeGrad)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(extrap)
      ,CHFt_CONST_FRA(Jga)
      ,CHFt_BOX(edgeGrid)
      ,CHFt_CONST_REAL(dxDir)
      ,CHFt_CONST_REAL(pdScaleDir)
      ,CHFt_CONST_INT(gradDir)
      ,CHFt_CONST_INT(edgeDir) );
}
#endif  // GUARDMAPPEDMACGRADORTHOSCALE 

#ifndef GUARDMAPPEDCCGRADSCALE 
#define GUARDMAPPEDCCGRADSCALE 
// Prototype for Fortran procedure MappedCCGradScale ...
//
void FORTRAN_NAME( MAPPEDCCGRADSCALE ,mappedccgradscale )(
      CHFp_FRA1(grad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA(Jga)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(pdScale)
      ,CHFp_CONST_INT(gradDir) );

#define FORT_MAPPEDCCGRADSCALE FORTRAN_NAME( inlineMAPPEDCCGRADSCALE, inlineMAPPEDCCGRADSCALE)
#define FORTNT_MAPPEDCCGRADSCALE FORTRAN_NAME( MAPPEDCCGRADSCALE, mappedccgradscale)

inline void FORTRAN_NAME(inlineMAPPEDCCGRADSCALE, inlineMAPPEDCCGRADSCALE)(
      CHFp_FRA1(grad)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA(Jga)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(pdScale)
      ,CHFp_CONST_INT(gradDir) )
{
 CH_TIMELEAF("FORT_MAPPEDCCGRADSCALE");
 FORTRAN_NAME( MAPPEDCCGRADSCALE ,mappedccgradscale )(
      CHFt_FRA1(grad)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_FRA1(extrap)
      ,CHFt_CONST_FRA(Jga)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(pdScale)
      ,CHFt_CONST_INT(gradDir) );
}
#endif  // GUARDMAPPEDCCGRADSCALE 

#ifndef GUARDCRSEONESIDEGRAD 
#define GUARDCRSEONESIDEGRAD 
// Prototype for Fortran procedure CRSEONESIDEGRAD ...
//
void FORTRAN_NAME( CRSEONESIDEGRAD ,crseonesidegrad )(
      CHFp_FRA1(edgeGrad)
      ,CHFp_FIA1(mask)
      ,CHFp_BOX(loEdgeBox)
      ,CHFp_BOX(hiEdgeBox)
      ,CHFp_INT(dir)
      ,CHFp_INT(doLo)
      ,CHFp_INT(doHi) );

#define FORT_CRSEONESIDEGRAD FORTRAN_NAME( inlineCRSEONESIDEGRAD, inlineCRSEONESIDEGRAD)
#define FORTNT_CRSEONESIDEGRAD FORTRAN_NAME( CRSEONESIDEGRAD, crseonesidegrad)

inline void FORTRAN_NAME(inlineCRSEONESIDEGRAD, inlineCRSEONESIDEGRAD)(
      CHFp_FRA1(edgeGrad)
      ,CHFp_FIA1(mask)
      ,CHFp_BOX(loEdgeBox)
      ,CHFp_BOX(hiEdgeBox)
      ,CHFp_INT(dir)
      ,CHFp_INT(doLo)
      ,CHFp_INT(doHi) )
{
 CH_TIMELEAF("FORT_CRSEONESIDEGRAD");
 FORTRAN_NAME( CRSEONESIDEGRAD ,crseonesidegrad )(
      CHFt_FRA1(edgeGrad)
      ,CHFt_FIA1(mask)
      ,CHFt_BOX(loEdgeBox)
      ,CHFt_BOX(hiEdgeBox)
      ,CHFt_INT(dir)
      ,CHFt_INT(doLo)
      ,CHFt_INT(doHi) );
}
#endif  // GUARDCRSEONESIDEGRAD 

#ifndef GUARDSIMPLECRSEONESIDEGRAD 
#define GUARDSIMPLECRSEONESIDEGRAD 
// Prototype for Fortran procedure SIMPLECRSEONESIDEGRAD ...
//
void FORTRAN_NAME( SIMPLECRSEONESIDEGRAD ,simplecrseonesidegrad )(
      CHFp_FRA1(grad)
      ,CHFp_FIA1(mask)
      ,CHFp_BOX(loEdgeBox)
      ,CHFp_BOX(hiEdgeBox)
      ,CHFp_INT(dir)
      ,CHFp_INT(doLo)
      ,CHFp_INT(doHi) );

#define FORT_SIMPLECRSEONESIDEGRAD FORTRAN_NAME( inlineSIMPLECRSEONESIDEGRAD, inlineSIMPLECRSEONESIDEGRAD)
#define FORTNT_SIMPLECRSEONESIDEGRAD FORTRAN_NAME( SIMPLECRSEONESIDEGRAD, simplecrseonesidegrad)

inline void FORTRAN_NAME(inlineSIMPLECRSEONESIDEGRAD, inlineSIMPLECRSEONESIDEGRAD)(
      CHFp_FRA1(grad)
      ,CHFp_FIA1(mask)
      ,CHFp_BOX(loEdgeBox)
      ,CHFp_BOX(hiEdgeBox)
      ,CHFp_INT(dir)
      ,CHFp_INT(doLo)
      ,CHFp_INT(doHi) )
{
 CH_TIMELEAF("FORT_SIMPLECRSEONESIDEGRAD");
 FORTRAN_NAME( SIMPLECRSEONESIDEGRAD ,simplecrseonesidegrad )(
      CHFt_FRA1(grad)
      ,CHFt_FIA1(mask)
      ,CHFt_BOX(loEdgeBox)
      ,CHFt_BOX(hiEdgeBox)
      ,CHFt_INT(dir)
      ,CHFt_INT(doLo)
      ,CHFt_INT(doHi) );
}
#endif  // GUARDSIMPLECRSEONESIDEGRAD 

#ifndef GUARDCCGRADFC2D 
#define GUARDCCGRADFC2D 
// Prototype for Fortran procedure CCGRADFC2D ...
//
void FORTRAN_NAME( CCGRADFC2D ,ccgradfc2d )(
      CHFp_FRA1(sol)
      ,CHFp_CONST_FRA(CC)
      ,CHFp_CONST_FRA1(FC0)
      ,CHFp_CONST_FRA1(FC1)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(dir) );

#define FORT_CCGRADFC2D FORTRAN_NAME( inlineCCGRADFC2D, inlineCCGRADFC2D)
#define FORTNT_CCGRADFC2D FORTRAN_NAME( CCGRADFC2D, ccgradfc2d)

inline void FORTRAN_NAME(inlineCCGRADFC2D, inlineCCGRADFC2D)(
      CHFp_FRA1(sol)
      ,CHFp_CONST_FRA(CC)
      ,CHFp_CONST_FRA1(FC0)
      ,CHFp_CONST_FRA1(FC1)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_CCGRADFC2D");
 FORTRAN_NAME( CCGRADFC2D ,ccgradfc2d )(
      CHFt_FRA1(sol)
      ,CHFt_CONST_FRA(CC)
      ,CHFt_CONST_FRA1(FC0)
      ,CHFt_CONST_FRA1(FC1)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDCCGRADFC2D 

#ifndef GUARDCCGRADFC3D 
#define GUARDCCGRADFC3D 
// Prototype for Fortran procedure CCGRADFC3D ...
//
void FORTRAN_NAME( CCGRADFC3D ,ccgradfc3d )(
      CHFp_FRA1(sol)
      ,CHFp_CONST_FRA(CC)
      ,CHFp_CONST_FRA1(FC0)
      ,CHFp_CONST_FRA1(FC1)
      ,CHFp_CONST_FRA1(FC2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(dir) );

#define FORT_CCGRADFC3D FORTRAN_NAME( inlineCCGRADFC3D, inlineCCGRADFC3D)
#define FORTNT_CCGRADFC3D FORTRAN_NAME( CCGRADFC3D, ccgradfc3d)

inline void FORTRAN_NAME(inlineCCGRADFC3D, inlineCCGRADFC3D)(
      CHFp_FRA1(sol)
      ,CHFp_CONST_FRA(CC)
      ,CHFp_CONST_FRA1(FC0)
      ,CHFp_CONST_FRA1(FC1)
      ,CHFp_CONST_FRA1(FC2)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_INT(dir) )
{
 CH_TIMELEAF("FORT_CCGRADFC3D");
 FORTRAN_NAME( CCGRADFC3D ,ccgradfc3d )(
      CHFt_FRA1(sol)
      ,CHFt_CONST_FRA(CC)
      ,CHFt_CONST_FRA1(FC0)
      ,CHFt_CONST_FRA1(FC1)
      ,CHFt_CONST_FRA1(FC2)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_INT(dir) );
}
#endif  // GUARDCCGRADFC3D 

#ifndef GUARDDONEFORMCC 
#define GUARDDONEFORMCC 
// Prototype for Fortran procedure DONEFORMCC ...
//
void FORTRAN_NAME( DONEFORMCC ,doneformcc )(
      CHFp_FRA(alpha)
      ,CHFp_CONST_FRA(omega)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx) );

#define FORT_DONEFORMCC FORTRAN_NAME( inlineDONEFORMCC, inlineDONEFORMCC)
#define FORTNT_DONEFORMCC FORTRAN_NAME( DONEFORMCC, doneformcc)

inline void FORTRAN_NAME(inlineDONEFORMCC, inlineDONEFORMCC)(
      CHFp_FRA(alpha)
      ,CHFp_CONST_FRA(omega)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx) )
{
 CH_TIMELEAF("FORT_DONEFORMCC");
 FORTRAN_NAME( DONEFORMCC ,doneformcc )(
      CHFt_FRA(alpha)
      ,CHFt_CONST_FRA(omega)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx) );
}
#endif  // GUARDDONEFORMCC 

#ifndef GUARDMAPPEDFLUXDIVERGENCE1D 
#define GUARDMAPPEDFLUXDIVERGENCE1D 
// Prototype for Fortran procedure MAPPEDFLUXDIVERGENCE1D ...
//
void FORTRAN_NAME( MAPPEDFLUXDIVERGENCE1D ,mappedfluxdivergence1d )(
      CHFp_FRA(div)
      ,CHFp_CONST_FRA(flux0)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) );

#define FORT_MAPPEDFLUXDIVERGENCE1D FORTRAN_NAME( inlineMAPPEDFLUXDIVERGENCE1D, inlineMAPPEDFLUXDIVERGENCE1D)
#define FORTNT_MAPPEDFLUXDIVERGENCE1D FORTRAN_NAME( MAPPEDFLUXDIVERGENCE1D, mappedfluxdivergence1d)

inline void FORTRAN_NAME(inlineMAPPEDFLUXDIVERGENCE1D, inlineMAPPEDFLUXDIVERGENCE1D)(
      CHFp_FRA(div)
      ,CHFp_CONST_FRA(flux0)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_MAPPEDFLUXDIVERGENCE1D");
 FORTRAN_NAME( MAPPEDFLUXDIVERGENCE1D ,mappedfluxdivergence1d )(
      CHFt_FRA(div)
      ,CHFt_CONST_FRA(flux0)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDMAPPEDFLUXDIVERGENCE1D 

#ifndef GUARDMAPPEDFLUXDIVERGENCE2D 
#define GUARDMAPPEDFLUXDIVERGENCE2D 
// Prototype for Fortran procedure MAPPEDFLUXDIVERGENCE2D ...
//
void FORTRAN_NAME( MAPPEDFLUXDIVERGENCE2D ,mappedfluxdivergence2d )(
      CHFp_FRA(div)
      ,CHFp_CONST_FRA(flux0)
      ,CHFp_CONST_FRA(flux1)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx) );

#define FORT_MAPPEDFLUXDIVERGENCE2D FORTRAN_NAME( inlineMAPPEDFLUXDIVERGENCE2D, inlineMAPPEDFLUXDIVERGENCE2D)
#define FORTNT_MAPPEDFLUXDIVERGENCE2D FORTRAN_NAME( MAPPEDFLUXDIVERGENCE2D, mappedfluxdivergence2d)

inline void FORTRAN_NAME(inlineMAPPEDFLUXDIVERGENCE2D, inlineMAPPEDFLUXDIVERGENCE2D)(
      CHFp_FRA(div)
      ,CHFp_CONST_FRA(flux0)
      ,CHFp_CONST_FRA(flux1)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx) )
{
 CH_TIMELEAF("FORT_MAPPEDFLUXDIVERGENCE2D");
 FORTRAN_NAME( MAPPEDFLUXDIVERGENCE2D ,mappedfluxdivergence2d )(
      CHFt_FRA(div)
      ,CHFt_CONST_FRA(flux0)
      ,CHFt_CONST_FRA(flux1)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx) );
}
#endif  // GUARDMAPPEDFLUXDIVERGENCE2D 

#ifndef GUARDMAPPEDFLUXDIVERGENCE3D 
#define GUARDMAPPEDFLUXDIVERGENCE3D 
// Prototype for Fortran procedure MAPPEDFLUXDIVERGENCE3D ...
//
void FORTRAN_NAME( MAPPEDFLUXDIVERGENCE3D ,mappedfluxdivergence3d )(
      CHFp_FRA(div)
      ,CHFp_CONST_FRA(flux0)
      ,CHFp_CONST_FRA(flux1)
      ,CHFp_CONST_FRA(flux2)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx) );

#define FORT_MAPPEDFLUXDIVERGENCE3D FORTRAN_NAME( inlineMAPPEDFLUXDIVERGENCE3D, inlineMAPPEDFLUXDIVERGENCE3D)
#define FORTNT_MAPPEDFLUXDIVERGENCE3D FORTRAN_NAME( MAPPEDFLUXDIVERGENCE3D, mappedfluxdivergence3d)

inline void FORTRAN_NAME(inlineMAPPEDFLUXDIVERGENCE3D, inlineMAPPEDFLUXDIVERGENCE3D)(
      CHFp_FRA(div)
      ,CHFp_CONST_FRA(flux0)
      ,CHFp_CONST_FRA(flux1)
      ,CHFp_CONST_FRA(flux2)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx) )
{
 CH_TIMELEAF("FORT_MAPPEDFLUXDIVERGENCE3D");
 FORTRAN_NAME( MAPPEDFLUXDIVERGENCE3D ,mappedfluxdivergence3d )(
      CHFt_FRA(div)
      ,CHFt_CONST_FRA(flux0)
      ,CHFt_CONST_FRA(flux1)
      ,CHFt_CONST_FRA(flux2)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx) );
}
#endif  // GUARDMAPPEDFLUXDIVERGENCE3D 

#ifndef GUARDMAPPEDMACDIVSCALE2D 
#define GUARDMAPPEDMACDIVSCALE2D 
// Prototype for Fortran procedure MappedMACDivScale2D ...
//
void FORTRAN_NAME( MAPPEDMACDIVSCALE2D ,mappedmacdivscale2d )(
      CHFp_FRA(div)
      ,CHFp_CONST_FRA(flux0)
      ,CHFp_CONST_FRA(flux1)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(pdScale) );

#define FORT_MAPPEDMACDIVSCALE2D FORTRAN_NAME( inlineMAPPEDMACDIVSCALE2D, inlineMAPPEDMACDIVSCALE2D)
#define FORTNT_MAPPEDMACDIVSCALE2D FORTRAN_NAME( MAPPEDMACDIVSCALE2D, mappedmacdivscale2d)

inline void FORTRAN_NAME(inlineMAPPEDMACDIVSCALE2D, inlineMAPPEDMACDIVSCALE2D)(
      CHFp_FRA(div)
      ,CHFp_CONST_FRA(flux0)
      ,CHFp_CONST_FRA(flux1)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(pdScale) )
{
 CH_TIMELEAF("FORT_MAPPEDMACDIVSCALE2D");
 FORTRAN_NAME( MAPPEDMACDIVSCALE2D ,mappedmacdivscale2d )(
      CHFt_FRA(div)
      ,CHFt_CONST_FRA(flux0)
      ,CHFt_CONST_FRA(flux1)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(pdScale) );
}
#endif  // GUARDMAPPEDMACDIVSCALE2D 

#ifndef GUARDMAPPEDMACDIVSCALE3D 
#define GUARDMAPPEDMACDIVSCALE3D 
// Prototype for Fortran procedure MappedMACDivScale3D ...
//
void FORTRAN_NAME( MAPPEDMACDIVSCALE3D ,mappedmacdivscale3d )(
      CHFp_FRA(div)
      ,CHFp_CONST_FRA(flux0)
      ,CHFp_CONST_FRA(flux1)
      ,CHFp_CONST_FRA(flux2)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(pdScale) );

#define FORT_MAPPEDMACDIVSCALE3D FORTRAN_NAME( inlineMAPPEDMACDIVSCALE3D, inlineMAPPEDMACDIVSCALE3D)
#define FORTNT_MAPPEDMACDIVSCALE3D FORTRAN_NAME( MAPPEDMACDIVSCALE3D, mappedmacdivscale3d)

inline void FORTRAN_NAME(inlineMAPPEDMACDIVSCALE3D, inlineMAPPEDMACDIVSCALE3D)(
      CHFp_FRA(div)
      ,CHFp_CONST_FRA(flux0)
      ,CHFp_CONST_FRA(flux1)
      ,CHFp_CONST_FRA(flux2)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(pdScale) )
{
 CH_TIMELEAF("FORT_MAPPEDMACDIVSCALE3D");
 FORTRAN_NAME( MAPPEDMACDIVSCALE3D ,mappedmacdivscale3d )(
      CHFt_FRA(div)
      ,CHFt_CONST_FRA(flux0)
      ,CHFt_CONST_FRA(flux1)
      ,CHFt_CONST_FRA(flux2)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(pdScale) );
}
#endif  // GUARDMAPPEDMACDIVSCALE3D 

#ifndef GUARDMAPPEDCCDIVSCALE 
#define GUARDMAPPEDCCDIVSCALE 
// Prototype for Fortran procedure MappedCCDivScale ...
//
void FORTRAN_NAME( MAPPEDCCDIVSCALE ,mappedccdivscale )(
      CHFp_FRA1(div)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(pdScale) );

#define FORT_MAPPEDCCDIVSCALE FORTRAN_NAME( inlineMAPPEDCCDIVSCALE, inlineMAPPEDCCDIVSCALE)
#define FORTNT_MAPPEDCCDIVSCALE FORTRAN_NAME( MAPPEDCCDIVSCALE, mappedccdivscale)

inline void FORTRAN_NAME(inlineMAPPEDCCDIVSCALE, inlineMAPPEDCCDIVSCALE)(
      CHFp_FRA1(div)
      ,CHFp_CONST_FRA(flux)
      ,CHFp_CONST_FRA1(Jinv)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(pdScale) )
{
 CH_TIMELEAF("FORT_MAPPEDCCDIVSCALE");
 FORTRAN_NAME( MAPPEDCCDIVSCALE ,mappedccdivscale )(
      CHFt_FRA1(div)
      ,CHFt_CONST_FRA(flux)
      ,CHFt_CONST_FRA1(Jinv)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(pdScale) );
}
#endif  // GUARDMAPPEDCCDIVSCALE 

}

#endif
