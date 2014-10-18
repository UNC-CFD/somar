#ifndef _AMRNAVIERSTOKESFACTORY_H_
#define _AMRNAVIERSTOKESFACTORY_H_

#include "REAL.H"
#include "MappedAMRLevelFactory.H"
class PhysBCUtil;


class AMRNavierStokesFactory : public MappedAMRLevelFactory
{
public:

    // Constructor
    AMRNavierStokesFactory ();

    // Destructor
    virtual ~AMRNavierStokesFactory ();

    // Virtual override of level factory
    virtual MappedAMRLevel* new_amrlevel () const;

    //
    void setPhysBC (const PhysBCUtil& a_physBC);

protected:

    //
    Real m_cfl;

    //
    PhysBCUtil* m_physBCPtr;
};


#endif
