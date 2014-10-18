#ifndef __ADVECTIONTESTBCUTIL_HPP__INCLUDED__
#define __ADVECTIONTESTBCUTIL_HPP__INCLUDED__

#include "PhysBCUtil.H"


class AdvectionTestBCUtil: public PhysBCUtil
{
public:
    // Default constrictor
    AdvectionTestBCUtil ();

    // Default destructor
    virtual ~AdvectionTestBCUtil ();

    // This object is its own factory
    virtual PhysBCUtil* newPhysBCUtil () const;


    // ICs / background fields -------------------------------------------------

    // Fills a FAB with the initial velocity.
    // Locations are in mapped space, but components are Cartesian.
    virtual void setVelIC (FArrayBox&           a_velFAB,
                           const int            a_velComp,
                           const LevelGeometry& a_levGeo,
                           const DataIndex&     a_di) const;

    // Fills a FAB with the initial scalars
    virtual void setScalarIC (FArrayBox&           a_scalarFAB,
                              const int            a_scalarComp,
                              const LevelGeometry& a_levGeo,
                              const DataIndex&     a_di) const;
};


#endif //!__ADVECTIONTESTBCUTIL_HPP__INCLUDED__
