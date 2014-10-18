#ifndef __LOCKEXCHANGEBCUTIL_HPP__INCLUDED__
#define __LOCKEXCHANGEBCUTIL_HPP__INCLUDED__

#include "PhysBCUtil.H"


class LockExchangeBCUtil: public PhysBCUtil
{
public:
    // Default constrictor
    LockExchangeBCUtil ();

    // Default destructor
    virtual ~LockExchangeBCUtil ();

    // This object is its own factory
    virtual PhysBCUtil* newPhysBCUtil () const;

    // ICs / background fields -------------------------------------------------

    // Fills a FAB with the initial scalars
    virtual void setScalarIC (FArrayBox&           a_scalarFAB,
                              const int            a_scalarComp,
                              const LevelGeometry& a_levGeo,
                              const DataIndex&     a_di) const;

protected:
    // The basic BC functions --------------------------------------------------

    /// Sets physical BCs on a generic passive scalar
    virtual BCMethodHolder basicScalarFuncBC () const;

    // Sets physical BCs on velocities.
    virtual BCMethodHolder basicVelFuncBC (int a_veldir, bool a_isViscous) const;
};



#endif //!__LOCKEXCHANGEBCUTIL_HPP__INCLUDED__
