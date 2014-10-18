/*******************************************************************************
 *    FILE: MappedLevelCrankNicolson.H
 *    DESC: The public interface for solving heat problems using the
 *          Crank-Nicolson algorithm.
 *    DATE: Tue 27 Aug 2013 09:28:17 AM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __MappedLevelCrankNicolson_H__INCLUDED__
#define __MappedLevelCrankNicolson_H__INCLUDED__

#include "MappedBaseLevelHeatSolver.H"


class MappedLevelCrankNicolson : public MappedBaseLevelHeatSolver
{
    // These used to be template parameters before the great rewrite of 'fourteen.
    typedef LevelData<FArrayBox>    LevelDataType;
    typedef FluxBox                 FluxDataType;
    typedef MappedLevelFluxRegister FluxRegisterType;

public:
    // Full constructor
    MappedLevelCrankNicolson (const Vector<DisjointBoxLayout>&                                 a_grids,
                              const Vector<IntVect>&                                           a_refRat,
                              const ProblemDomain&                                             a_level0Domain,
                              RefCountedPtr<MappedAMRLevelOpFactory<LevelData<FArrayBox> > >&  a_opFact,
                              const RefCountedPtr<MappedAMRMultiGrid<LevelData<FArrayBox> > >& a_solver);

    // Destructor
    virtual ~MappedLevelCrankNicolson ();

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
    MappedLevelCrankNicolson& operator= (const MappedLevelCrankNicolson&);
    MappedLevelCrankNicolson (const MappedLevelCrankNicolson& a_opin);
    MappedLevelCrankNicolson ();
};


#endif //!__MappedLevelCrankNicolson_H__INCLUDED__
