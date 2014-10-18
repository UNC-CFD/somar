/*******************************************************************************
 *    FILE: InternalWaveBCUtil.H
 *    DESC: Simple internal wave demo.
 *    DATE: Mon 13 May 2013 01:10:23 AM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __InternalWaveBCUtil_H__INCLUDED__
#define __InternalWaveBCUtil_H__INCLUDED__

#include "PhysBCUtil.H"


class InternalWaveBCUtil: public PhysBCUtil
{
public:

    // Default constructor
    InternalWaveBCUtil ();

    // Default destructor
    virtual ~InternalWaveBCUtil ();

    // Factory
    virtual PhysBCUtil* newPhysBCUtil () const;

    // This is in case the BC's have an effect on the timestep.
    // Pass in currently computed dt, along with the cfl and dx. If the effect
    // of the BCs requires a decreased timestep, then the newly reduced timestep
    // is returned.
    virtual void computeBoundaryDt (Real&                a_dt,
                                    const Real           a_cfl,
                                    const LevelGeometry& a_levGeo) const;

    // Fills a FAB with the initial scalars
    virtual void setScalarIC (FArrayBox&           a_scalarFAB,
                              const int            a_scalarComp,
                              const LevelGeometry& a_levGeo,
                              const DataIndex&     a_di) const;

    // Fills a FAB with the background scalar
    virtual void setBackgroundScalar (FArrayBox&           a_scalarFAB,
                                      const int            a_scalarComp,
                                      const LevelGeometry& a_levGeo,
                                      const DataIndex&     a_di,
                                      const Real           a_time) const;
protected:
    static void bscalBCValues (Real*           a_pos,
                               int*            a_dir,
                               Side::LoHiSide* a_side,
                               Real*           a_value,
                               Real            a_derivScale,
                               Real            a_time);
};



#endif //!__InternalWaveBCUtil_H__INCLUDED__
