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
#ifndef __MappedLevelTGA_H__INCLUDED__
#define __MappedLevelTGA_H__INCLUDED__

#include "MappedBaseLevelHeatSolver.H"


// -----------------------------------------------------------------------------
// Implements TGA algorithm to update solution to heat equation
// The MappedLevelTGA class implements the Runge-Kutta-based approach
// to solving the heat equation due to Twizell, Gumel, and Arigu.
// -----------------------------------------------------------------------------
class MappedLevelTGA: public MappedBaseLevelHeatSolver
{
    // These used to be template parameters before the great rewrite of 'fourteen.
    typedef LevelData<FArrayBox>    LevelDataType;
    typedef FluxBox                 FluxDataType;
    typedef MappedLevelFluxRegister FluxRegisterType;

public:

    // Full constructor
    MappedLevelTGA (const Vector<DisjointBoxLayout>&                                 a_grids,
                    const Vector<IntVect>&                                           a_refRat,
                    const ProblemDomain&                                             a_level0Domain,
                    RefCountedPtr<MappedAMRLevelOpFactory<LevelData<FArrayBox> > >&  a_opFact,
                    const RefCountedPtr<MappedAMRMultiGrid<LevelData<FArrayBox> > >& a_solver);

    /// destructor
    virtual ~MappedLevelTGA ();

    // Computes the time-centered diffusion term L(phi). This can be used to
    // find contributions to the solution from diffusion. The diffusion term
    // is computed by computing a finite difference approximation for d phi/dt
    // using the updated and original values of phi and the time step. Most of
    // the arguments given here are passed along to updateSoln and retain their
    // significance therein.
    virtual void computeDiffusion (LevelDataType&           a_diffusiveTerm,
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

protected:
    // Update the solution assuming that the operator's coefficients are
    // independent of time. Same arguments as updateSoln.
    void updateSolnWithTimeIndependentOp (LevelDataType&           a_phiNew,
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
                                          int                      a_fluxStartComponent = 0);

    // Update the solution assuming that the operator's coefficients change
    // with time. Same arguments as updateSoln.
    void updateSolnWithTimeDependentOp (LevelDataType&           a_phiNew,
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
                                        int                      a_fluxStartComponent = 0);

    // The times within the integration step at which the operators are
    // evaluated.
    Real m_mu1, m_mu2, m_mu3, m_mu4, m_r1;

private:
    // Disallowed operators.
    MappedLevelTGA& operator= (const MappedLevelTGA&);
    MappedLevelTGA (const MappedLevelTGA& a_opin);
    MappedLevelTGA ();
};


#endif //!__MappedLevelTGA_H__INCLUDED__
