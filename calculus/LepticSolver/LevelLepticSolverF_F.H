#ifndef _LEVELLEPTICSOLVERF_F_H_
#define _LEVELLEPTICSOLVERF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDLEPTICVERTHORIZGRAD 
#define GUARDLEPTICVERTHORIZGRAD 
// Prototype for Fortran procedure LepticVertHorizGrad ...
//
void FORTRAN_NAME( LEPTICVERTHORIZGRAD ,lepticverthorizgrad )(
      CHFp_FRA1(bcVals)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA(Jgz)
      ,CHFp_BOX(destFCBox)
      ,CHFp_CONST_INT(isign)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(scale) );

#define FORT_LEPTICVERTHORIZGRAD FORTRAN_NAME( inlineLEPTICVERTHORIZGRAD, inlineLEPTICVERTHORIZGRAD)
#define FORTNT_LEPTICVERTHORIZGRAD FORTRAN_NAME( LEPTICVERTHORIZGRAD, lepticverthorizgrad)

inline void FORTRAN_NAME(inlineLEPTICVERTHORIZGRAD, inlineLEPTICVERTHORIZGRAD)(
      CHFp_FRA1(bcVals)
      ,CHFp_CONST_FRA1(extrap)
      ,CHFp_CONST_FRA(Jgz)
      ,CHFp_BOX(destFCBox)
      ,CHFp_CONST_INT(isign)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REAL(scale) )
{
 CH_TIMELEAF("FORT_LEPTICVERTHORIZGRAD");
 FORTRAN_NAME( LEPTICVERTHORIZGRAD ,lepticverthorizgrad )(
      CHFt_FRA1(bcVals)
      ,CHFt_CONST_FRA1(extrap)
      ,CHFt_CONST_FRA(Jgz)
      ,CHFt_BOX(destFCBox)
      ,CHFt_CONST_INT(isign)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REAL(scale) );
}
#endif  // GUARDLEPTICVERTHORIZGRAD 

#ifndef GUARDLEPTICLAPACKVERTICALSOLVER 
#define GUARDLEPTICLAPACKVERTICALSOLVER 
// Prototype for Fortran procedure LepticLapackVerticalSolver ...
//
void FORTRAN_NAME( LEPTICLAPACKVERTICALSOLVER ,lepticlapackverticalsolver )(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(Jgzz)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dz)
      ,CHFp_CONST_REAL(dzCrse)
      ,CHFp_FRA1(DU)
      ,CHFp_FRA1(D)
      ,CHFp_FRA1(DL)
      ,CHFp_FRA1(B)
      ,CHFp_CONST_INT(loBCType)
      ,CHFp_CONST_INT(hiBCType) );

#define FORT_LEPTICLAPACKVERTICALSOLVER FORTRAN_NAME( inlineLEPTICLAPACKVERTICALSOLVER, inlineLEPTICLAPACKVERTICALSOLVER)
#define FORTNT_LEPTICLAPACKVERTICALSOLVER FORTRAN_NAME( LEPTICLAPACKVERTICALSOLVER, lepticlapackverticalsolver)

inline void FORTRAN_NAME(inlineLEPTICLAPACKVERTICALSOLVER, inlineLEPTICLAPACKVERTICALSOLVER)(
      CHFp_FRA1(phi)
      ,CHFp_CONST_FRA1(rhs)
      ,CHFp_CONST_FRA1(Jgzz)
      ,CHFp_BOX(region)
      ,CHFp_CONST_REAL(dz)
      ,CHFp_CONST_REAL(dzCrse)
      ,CHFp_FRA1(DU)
      ,CHFp_FRA1(D)
      ,CHFp_FRA1(DL)
      ,CHFp_FRA1(B)
      ,CHFp_CONST_INT(loBCType)
      ,CHFp_CONST_INT(hiBCType) )
{
 CH_TIMELEAF("FORT_LEPTICLAPACKVERTICALSOLVER");
 FORTRAN_NAME( LEPTICLAPACKVERTICALSOLVER ,lepticlapackverticalsolver )(
      CHFt_FRA1(phi)
      ,CHFt_CONST_FRA1(rhs)
      ,CHFt_CONST_FRA1(Jgzz)
      ,CHFt_BOX(region)
      ,CHFt_CONST_REAL(dz)
      ,CHFt_CONST_REAL(dzCrse)
      ,CHFt_FRA1(DU)
      ,CHFt_FRA1(D)
      ,CHFt_FRA1(DL)
      ,CHFt_FRA1(B)
      ,CHFt_CONST_INT(loBCType)
      ,CHFt_CONST_INT(hiBCType) );
}
#endif  // GUARDLEPTICLAPACKVERTICALSOLVER 

#ifndef GUARDLEPTICACCUMDIV 
#define GUARDLEPTICACCUMDIV 
// Prototype for Fortran procedure LepticAccumDiv ...
//
void FORTRAN_NAME( LEPTICACCUMDIV ,lepticaccumdiv )(
      CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA1(flux)
      ,CHFp_BOX(CCregion)
      ,CHFp_CONST_REAL(dxScale)
      ,CHFp_CONST_INT(FCdir) );

#define FORT_LEPTICACCUMDIV FORTRAN_NAME( inlineLEPTICACCUMDIV, inlineLEPTICACCUMDIV)
#define FORTNT_LEPTICACCUMDIV FORTRAN_NAME( LEPTICACCUMDIV, lepticaccumdiv)

inline void FORTRAN_NAME(inlineLEPTICACCUMDIV, inlineLEPTICACCUMDIV)(
      CHFp_FRA1(rhs)
      ,CHFp_CONST_FRA1(flux)
      ,CHFp_BOX(CCregion)
      ,CHFp_CONST_REAL(dxScale)
      ,CHFp_CONST_INT(FCdir) )
{
 CH_TIMELEAF("FORT_LEPTICACCUMDIV");
 FORTRAN_NAME( LEPTICACCUMDIV ,lepticaccumdiv )(
      CHFt_FRA1(rhs)
      ,CHFt_CONST_FRA1(flux)
      ,CHFt_BOX(CCregion)
      ,CHFt_CONST_REAL(dxScale)
      ,CHFt_CONST_INT(FCdir) );
}
#endif  // GUARDLEPTICACCUMDIV 

}

#endif
