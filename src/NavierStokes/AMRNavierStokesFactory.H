/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2014 University of North Carolina at Chapel Hill
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
 *  USA
 *
 *  For up-to-date contact information, please visit the repository homepage,
 *  https://github.com/somarhub.
 ******************************************************************************/
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
