/*******************************************************************************
 *    FILE: VortexStreetBCUtil.H
 *    DESC: IBCs for the Von Karman vortex demonstration.
 *    DATE: Sun 03 Nov 2013 02:57:38 PM EST
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __VortexStreetBCUtil_H__INCLUDED__
#define __VortexStreetBCUtil_H__INCLUDED__

#include "PhysBCUtil.H"


class VortexStreetBCUtil: public PhysBCUtil
{
public:

    // Default constructor
    VortexStreetBCUtil ();

    // Default destructor
    virtual ~VortexStreetBCUtil ();

    // Factory
    virtual PhysBCUtil* newPhysBCUtil () const;

    // Fills a FAB with the initial velocity.
    // Locations are in mapped space, but components are Cartesian.
    virtual void setVelIC (FArrayBox&           a_velFAB,
                           const int            a_velComp,
                           const LevelGeometry& a_levGeo,
                           const DataIndex&     a_di) const;

    // Sets the Cartesian-based target velocity for the sponge layer.
    // By default, this function throws an error.
    virtual void fillVelSpongeLayerTarget (FArrayBox&           a_target,
                                           const int            a_velComp,
                                           const int            a_spongeDir,
                                           const Side::LoHiSide a_spongeSide,
                                           const LevelGeometry& a_levGeo,
                                           const DataIndex&     a_di,
                                           const Real           a_time);

    // Sets the target values for the scalar sponge layer. If we are using a
    // background scalar, then this function set the perturbation to zero.
    // Otherwise, an error is thrown and this function will need to be overridden.
    virtual void fillScalarSpongeLayerTarget (FArrayBox&           a_target,
                                              const int            a_scalarComp,
                                              const int            a_spongeDir,
                                              const Side::LoHiSide a_spongeSide,
                                              const LevelGeometry& a_levGeo,
                                              const DataIndex&     a_di,
                                              const Real           a_time);

protected:
    // Sets physical BCs on velocities
    virtual BCMethodHolder basicVelFuncBC (int a_veldir, bool a_isViscous) const;

    static bool s_paramsRead;
    static RealVect s_inflowVel;
};


#endif //!__VortexStreetBCUtil_H__INCLUDED__
