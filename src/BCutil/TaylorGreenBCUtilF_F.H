#ifndef _TAYLORGREENBCUTILF_F_H_
#define _TAYLORGREENBCUTILF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDFILL_CARTVELSOLN 
#define GUARDFILL_CARTVELSOLN 
// Prototype for Fortran procedure FILL_CARTVELSOLN ...
//
void FORTRAN_NAME( FILL_CARTVELSOLN ,fill_cartvelsoln )(
      CHFp_FRA1(vel)
      ,CHFp_CONST_INT(solComp)
      ,CHFp_CONST_FRA1(xOfXi)
      ,CHFp_CONST_FRA1(yOfXi)
      ,CHFp_CONST_REAL(kx)
      ,CHFp_CONST_REAL(ky)
      ,CHFp_CONST_REAL(FofT)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(destBox) );

#define FORT_FILL_CARTVELSOLN FORTRAN_NAME( inlineFILL_CARTVELSOLN, inlineFILL_CARTVELSOLN)
#define FORTNT_FILL_CARTVELSOLN FORTRAN_NAME( FILL_CARTVELSOLN, fill_cartvelsoln)

inline void FORTRAN_NAME(inlineFILL_CARTVELSOLN, inlineFILL_CARTVELSOLN)(
      CHFp_FRA1(vel)
      ,CHFp_CONST_INT(solComp)
      ,CHFp_CONST_FRA1(xOfXi)
      ,CHFp_CONST_FRA1(yOfXi)
      ,CHFp_CONST_REAL(kx)
      ,CHFp_CONST_REAL(ky)
      ,CHFp_CONST_REAL(FofT)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(destBox) )
{
 CH_TIMELEAF("FORT_FILL_CARTVELSOLN");
 FORTRAN_NAME( FILL_CARTVELSOLN ,fill_cartvelsoln )(
      CHFt_FRA1(vel)
      ,CHFt_CONST_INT(solComp)
      ,CHFt_CONST_FRA1(xOfXi)
      ,CHFt_CONST_FRA1(yOfXi)
      ,CHFt_CONST_REAL(kx)
      ,CHFt_CONST_REAL(ky)
      ,CHFt_CONST_REAL(FofT)
      ,CHFt_CONST_REAL(time)
      ,CHFt_BOX(destBox) );
}
#endif  // GUARDFILL_CARTVELSOLN 

#ifndef GUARDFILL_PRESSURESOLN 
#define GUARDFILL_PRESSURESOLN 
// Prototype for Fortran procedure FILL_PRESSURESOLN ...
//
void FORTRAN_NAME( FILL_PRESSURESOLN ,fill_pressuresoln )(
      CHFp_FRA1(pressure)
      ,CHFp_CONST_FRA1(xOfXi)
      ,CHFp_CONST_FRA1(yOfXi)
      ,CHFp_CONST_REAL(kx)
      ,CHFp_CONST_REAL(ky)
      ,CHFp_CONST_REAL(FofT)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(destBox) );

#define FORT_FILL_PRESSURESOLN FORTRAN_NAME( inlineFILL_PRESSURESOLN, inlineFILL_PRESSURESOLN)
#define FORTNT_FILL_PRESSURESOLN FORTRAN_NAME( FILL_PRESSURESOLN, fill_pressuresoln)

inline void FORTRAN_NAME(inlineFILL_PRESSURESOLN, inlineFILL_PRESSURESOLN)(
      CHFp_FRA1(pressure)
      ,CHFp_CONST_FRA1(xOfXi)
      ,CHFp_CONST_FRA1(yOfXi)
      ,CHFp_CONST_REAL(kx)
      ,CHFp_CONST_REAL(ky)
      ,CHFp_CONST_REAL(FofT)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(destBox) )
{
 CH_TIMELEAF("FORT_FILL_PRESSURESOLN");
 FORTRAN_NAME( FILL_PRESSURESOLN ,fill_pressuresoln )(
      CHFt_FRA1(pressure)
      ,CHFt_CONST_FRA1(xOfXi)
      ,CHFt_CONST_FRA1(yOfXi)
      ,CHFt_CONST_REAL(kx)
      ,CHFt_CONST_REAL(ky)
      ,CHFt_CONST_REAL(FofT)
      ,CHFt_CONST_REAL(time)
      ,CHFt_BOX(destBox) );
}
#endif  // GUARDFILL_PRESSURESOLN 

#ifndef GUARDTAYLORGREENCCBC 
#define GUARDTAYLORGREENCCBC 
// Prototype for Fortran procedure TAYLORGREENCCBC ...
//
void FORTRAN_NAME( TAYLORGREENCCBC ,taylorgreenccbc )(
      CHFp_FRA(state)
      ,CHFp_BOX(ghostBox)
      ,CHFp_CONST_INTVECT(boxType)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(k)
      ,CHFp_CONST_INT(velDir)
      ,CHFp_CONST_REAL(FofT) );

#define FORT_TAYLORGREENCCBC FORTRAN_NAME( inlineTAYLORGREENCCBC, inlineTAYLORGREENCCBC)
#define FORTNT_TAYLORGREENCCBC FORTRAN_NAME( TAYLORGREENCCBC, taylorgreenccbc)

inline void FORTRAN_NAME(inlineTAYLORGREENCCBC, inlineTAYLORGREENCCBC)(
      CHFp_FRA(state)
      ,CHFp_BOX(ghostBox)
      ,CHFp_CONST_INTVECT(boxType)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REALVECT(dx)
      ,CHFp_CONST_REALVECT(k)
      ,CHFp_CONST_INT(velDir)
      ,CHFp_CONST_REAL(FofT) )
{
 CH_TIMELEAF("FORT_TAYLORGREENCCBC");
 FORTRAN_NAME( TAYLORGREENCCBC ,taylorgreenccbc )(
      CHFt_FRA(state)
      ,CHFt_BOX(ghostBox)
      ,CHFt_CONST_INTVECT(boxType)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_REALVECT(dx)
      ,CHFt_CONST_REALVECT(k)
      ,CHFt_CONST_INT(velDir)
      ,CHFt_CONST_REAL(FofT) );
}
#endif  // GUARDTAYLORGREENCCBC 

#ifndef GUARDTAYLORGREENFCBC 
#define GUARDTAYLORGREENFCBC 
// Prototype for Fortran procedure TAYLORGREENFCBC ...
//
void FORTRAN_NAME( TAYLORGREENFCBC ,taylorgreenfcbc )(
      CHFp_FRA(state)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REALVECT(k)
      ,CHFp_CONST_INT(velDir)
      ,CHFp_CONST_REAL(FofT) );

#define FORT_TAYLORGREENFCBC FORTRAN_NAME( inlineTAYLORGREENFCBC, inlineTAYLORGREENFCBC)
#define FORTNT_TAYLORGREENFCBC FORTRAN_NAME( TAYLORGREENFCBC, taylorgreenfcbc)

inline void FORTRAN_NAME(inlineTAYLORGREENFCBC, inlineTAYLORGREENFCBC)(
      CHFp_FRA(state)
      ,CHFp_BOX(faceBox)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(side)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REALVECT(k)
      ,CHFp_CONST_INT(velDir)
      ,CHFp_CONST_REAL(FofT) )
{
 CH_TIMELEAF("FORT_TAYLORGREENFCBC");
 FORTRAN_NAME( TAYLORGREENFCBC ,taylorgreenfcbc )(
      CHFt_FRA(state)
      ,CHFt_BOX(faceBox)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(side)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REALVECT(k)
      ,CHFt_CONST_INT(velDir)
      ,CHFt_CONST_REAL(FofT) );
}
#endif  // GUARDTAYLORGREENFCBC 

}

#endif
