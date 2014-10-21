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
#ifndef __MappedBaseLevelHeatSolver_H__INCLUDED__
#define __MappedBaseLevelHeatSolver_H__INCLUDED__

#include "MappedAMRTGA.H"
#include "MappedLevelFluxRegister.H"


// -----------------------------------------------------------------------------
// The abstract base class for the Backwards Euler, Crank-Nicolson,
// and TGA heat solvers. This class is a copy of Chombo's
// BaseLevelHeatSolver, except this version allows anisotropic
// refinement.
// -----------------------------------------------------------------------------
class MappedBaseLevelHeatSolver
{
    // These used to be template parameters before the great rewrite of 'fourteen.
    typedef LevelData<FArrayBox>    LevelDataType;
    typedef FluxBox                 FluxDataType;
    typedef MappedLevelFluxRegister FluxRegisterType;

public:
    // Constructor
    MappedBaseLevelHeatSolver (const Vector<DisjointBoxLayout>&                         a_grids,
                               const Vector<IntVect>&                                   a_refRat,
                               const ProblemDomain&                                     a_level0Domain,
                               RefCountedPtr<MappedAMRLevelOpFactory<LevelDataType> >&  a_opFact,
                               const RefCountedPtr<MappedAMRMultiGrid<LevelDataType> >& a_solver);

    // Destructor
    virtual ~MappedBaseLevelHeatSolver ();

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
                             int                      a_fluxStartComponent = 0) = 0;

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

    // Calls set time and calls operator with given alpha and beta
    virtual void applyOperator (LevelDataType&          a_ans,
                                const LevelDataType&    a_phi,
                                const LevelDataType*    a_phiC,
                                int                     a_level,
                                Real                    a_alpha,
                                Real                    a_beta,
                                bool                    a_applyBC);

    // Applies the Helmholtz operator to the solution a_phi at the given
    // grid level. This will set a_ans to (I + a_mu * a_dt * L(a_phi).
    virtual void applyHelm (LevelDataType&          a_ans,
                            const LevelDataType&    a_phi,
                            const LevelDataType*    a_phiC,
                            int                     a_level,
                            Real                    a_mu,
                            Real                    a_dt,
                            bool                    a_homogeneousBC);

    // Adds flux contributions from the Helmholtz operator at the current
    // grid level.
    virtual void incrementFlux (LevelData<FluxDataType>& a_diffusiveFlux,
                                LevelDataType&           a_phi,
                                int                      a_level,
                                Real                     a_mu,
                                Real                     a_dt,
                                Real                     a_sign,
                                bool                     a_setToZero);

    // Solves the Helmholtz equation (I - a_mu * a_dt * L(a_phi) = a_rhs
    // for a_phi. Here it is assumed that a solution a_phiC exists
    // on a coarser grid level.
    virtual void solveHelm (LevelDataType& a_phi,
                            LevelDataType& a_phiC,
                            LevelDataType& a_rhs,
                            int            a_level,
                            Real           a_mu,
                            Real           a_dt,
                            bool           a_zeroPhi = true);

    // Sets the alpha and beta parameters in each Helmholtz operator to the
    // given values.
    virtual void resetSolverAlphaAndBeta (const Real& a_alpha,
                                          const Real& a_beta);

    // Creates a new Helmholtz operator for use by the TGA integrator.
    virtual MappedLevelTGAHelmOp<LevelDataType, FluxDataType>* newOp (
        const ProblemDomain& a_indexSpace,
        RefCountedPtr<MappedAMRLevelOpFactory<LevelDataType> >& a_opFact);

    // Returns the number of grid levels on which this integrator operates.
    virtual int size () const;

protected:
    // Interpolates a given quantity linearly in time using its beginning- and
    // end-of-step values and placing the result in \a a_data.
    virtual void timeInterp (LevelDataType&       a_data,
                             const LevelDataType& a_oldData,
                             const LevelDataType& a_newData,
                             Real                 a_time,
                             Real                 a_oldTime,
                             Real                 a_newTime,
                             int                  a_level);

    // Fills ghosts cells on source terms.
    virtual void setSourceGhostCells (LevelData<FArrayBox>&    a_src,
                                      const DisjointBoxLayout& a_grids,
                                      int a_lev);

    // The disjoint box layouts at every AMR grid level.
    Vector<DisjointBoxLayout> m_grids;

    // The refinement ratios between AMR grid levels.
    Vector<IntVect>           m_refRat;

    // The coarsest domain on which the Helmholtz equation is integrated.
    ProblemDomain             m_level0Domain;

    // An array of the solver's Helmholtz operators at each grid level,
    // casted to MappedLevelTGAHelmOp instances. These are owned by the solver,
    // so we shouldn't delete these.
    Vector<MappedLevelTGAHelmOp<LevelDataType, FluxDataType>*>  m_ops;

    // The multigrid solver used to solve the Helmholtz equation.
    RefCountedPtr<MappedAMRMultiGrid<LevelDataType> >           m_solver;
};



#endif //!__MappedBaseLevelHeatSolver_H__INCLUDED__
