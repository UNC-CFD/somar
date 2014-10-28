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
#ifndef __MappedLevelBackwardEuler_H__INCLUDED__
#define __MappedLevelBackwardEuler_H__INCLUDED__

#include "MappedBaseLevelHeatSolver.H"


// -----------------------------------------------------------------------------
// The public interface for solving heat problems via backward Euler.
// -----------------------------------------------------------------------------
class MappedLevelBackwardEuler: public MappedBaseLevelHeatSolver
{
    // These used to be template parameters before the great rewrite of 'fourteen.
    typedef LevelData<FArrayBox>    LevelDataType;
    typedef FluxBox                 FluxDataType;
    typedef MappedLevelFluxRegister FluxRegisterType;

public:
    // Full constructor
    MappedLevelBackwardEuler (const Vector<DisjointBoxLayout>&                                 a_grids,
                              const Vector<IntVect>&                                           a_refRat,
                              const ProblemDomain&                                             a_level0Domain,
                              RefCountedPtr<MappedAMRLevelOpFactory<LevelData<FArrayBox> > >&  a_opFact,
                              const RefCountedPtr<MappedAMRMultiGrid<LevelData<FArrayBox> > >& a_solver);

    // Destructor
    virtual ~MappedLevelBackwardEuler ();

    // Integrates the helmholtz equation represented by this object, placing
    // the new solution in a_phiNew.
    virtual void updateSoln (LevelDataType&           a_phiNew,
                             LevelDataType&           a_phiOld,
                             LevelDataType&           a_src,
                             LevelData<FluxDataType>& a_flux,
                             FluxRegisterType*        a_fineFluxRegPtr,
                             FluxRegisterType*        a_crseFluxRegPtr,
                             const LevelDataType*     a_crsePhiOldPtr,
                             const LevelDataType*     a_crsePhiNewPtr,
                             Real                     a_oldTime,
                             Real                     a_crseOldTime,
                             Real                     a_crseNewTime,
                             Real                     a_dt,
                             int                      a_level,
                             bool                     a_zeroPhi = true,
                             bool                     a_rhsAlreadyKappaWeighted = false,
                             int                      a_fluxStartComponent = 0);

private:
    // Disallowed operators.
    MappedLevelBackwardEuler& operator= (const MappedLevelBackwardEuler&);
    MappedLevelBackwardEuler (const MappedLevelBackwardEuler& a_opin);
    MappedLevelBackwardEuler ();
};



#endif //!__MappedLevelBackwardEuler_H__INCLUDED__
