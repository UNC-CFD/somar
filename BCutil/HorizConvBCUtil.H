/*******************************************************************************
 *    FILE: HorizConvBCUtil.H
 *    DESC: Horizontal convection problem a la Scotti & White.
 *    DATE: Wed 09 Apr 2014 10:37:56 AM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/
#ifndef __HorizConvBCUtil_H__INCLUDED__
#define __HorizConvBCUtil_H__INCLUDED__

#include "PhysBCUtil.H"


class HorizConvBCUtil: public PhysBCUtil
{
public:

    // Default constructor
    HorizConvBCUtil ();

    // Default destructor
    virtual ~HorizConvBCUtil ();

    // Factory
    virtual PhysBCUtil* newPhysBCUtil () const;

    // Fills a FAB with the initial scalars
    virtual void setScalarIC (FArrayBox&           a_scalarFAB,
		              const int            a_scalarComp,
			      const LevelGeometry& a_levGeo,
			      const DataIndex&     a_di) const;

    // Velocity BCs ------------------------------------------------------------

    // Used in single-compoent velocity TGA solves
    virtual BCMethodHolder viscousSolveFuncBC (int a_dir) const;


    // Scalar BCs --------------------------------------------------------------

    // used in scalar TGA solves
    virtual BCMethodHolder diffusiveSolveFuncBC () const;

    // Sets ghosts needed to calculate the diffusive term nu.L[scalar]
    virtual BCMethodHolder diffusiveSourceFuncBC () const;


protected:
    // The basic BC functions --------------------------------------------------

    /// Sets physical BCs on velocities
    virtual BCMethodHolder basicVelFuncBC (int a_veldir, bool a_isViscous) const;
};


#endif //!__HorizConvBCUtil_H__INCLUDED__

