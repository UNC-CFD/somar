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
#include "MappedBaseLevelHeatSolver.H"
#include "ExtrapolationUtilsF_F.H"


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
MappedBaseLevelHeatSolver::MappedBaseLevelHeatSolver (const Vector<DisjointBoxLayout>&                         a_grids,
                                                      const Vector<IntVect>&                                   a_refRat,
                                                      const ProblemDomain&                                     a_level0Domain,
                                                      RefCountedPtr<MappedAMRLevelOpFactory<LevelDataType> >&  a_opFact,
                                                      const RefCountedPtr<MappedAMRMultiGrid<LevelDataType> >& a_solver)
: m_grids(a_grids),
  m_refRat(a_refRat),
  m_level0Domain(a_level0Domain),
  m_ops(),
  m_solver(a_solver)
{
    m_ops.resize(a_grids.size());
    Vector< MappedAMRLevelOp<LevelDataType> * >& amrops =  m_solver->getAMROperators();

    for (int ilev = 0; ilev < m_ops.size(); ilev++) {
        m_ops[ilev] = dynamic_cast<MappedLevelTGAHelmOp<LevelDataType, FluxDataType>* >(amrops[ilev]);
        if (m_ops[ilev] == NULL) {
            MayDay::Error("dynamic cast failed---is that operator really a TGAHelmOp?");
        }
    }
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
MappedBaseLevelHeatSolver::~MappedBaseLevelHeatSolver ()
{;}


// -----------------------------------------------------------------------------
// Computes the time-centered diffusion term L(phi). This can be used to
// find contributions to the solution from diffusion. The diffusion term
// is computed by computing a finite difference approximation for d phi/dt
// using the updated and original values of phi and the time step. Most of
// the arguments given here are passed along to updateSoln and retain their
// significance therein.
// -----------------------------------------------------------------------------
void MappedBaseLevelHeatSolver::computeDiffusion (LevelDataType&           a_diffusiveTerm,
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
                                                  bool                     a_zeroPhi,
                                                  bool                     a_rhsAlreadyKappaWeighted,
                                                  int                      a_fluxStartComponent)
{
    // The operator has no time-dependent parameters. Life is easier.

    // first compute updated solution
    LevelDataType phiNew;

    m_ops[a_level]->create(phiNew, a_phiOld);
    m_ops[a_level]->setToZero(phiNew);
    if (!a_zeroPhi) {
        m_ops[a_level]->assign(phiNew, a_phiOld);
    }

    updateSoln(phiNew, a_phiOld, a_src, a_flux,
               a_fineFluxRegPtr, a_crseFluxRegPtr,
               a_crsePhiOldPtr, a_crsePhiNewPtr,
               a_oldTime, a_crseOldTime,
               a_crseNewTime, a_dt, a_level, a_zeroPhi,
               a_rhsAlreadyKappaWeighted, a_fluxStartComponent);

    // now subtract everything off to leave us with diffusive term
    m_ops[a_level]->incr(phiNew, a_phiOld, -1.0);
    m_ops[a_level]->scale(phiNew, 1.0 / a_dt);

    //now multiply by a if there is an a
    m_ops[a_level]->diagonalScale(phiNew, false);

    // and finally, subtract off a_src
    m_ops[a_level]->incr(phiNew, a_src, -1.0);

    // what's left should be the time-centered diffusive part of the update
    m_ops[a_level]->assign(a_diffusiveTerm, phiNew);
}


// -----------------------------------------------------------------------------
// Calls set time and calls operator with given alpha and beta
// -----------------------------------------------------------------------------
void MappedBaseLevelHeatSolver::applyOperator (LevelDataType&          a_ans,
                                               const LevelDataType&    a_phi,
                                               const LevelDataType*    a_phiC,
                                               int                     a_level,
                                               Real                    a_alpha,
                                               Real                    a_beta,
                                               bool                    a_applyBC)
{
    m_ops[a_level]->setAlphaAndBeta(a_alpha, a_beta);

    LevelDataType zero;
    m_ops[a_level]->create(zero, a_ans);
    m_ops[a_level]->setToZero(zero);

    // set a_ans = helm(a_phi)
    //           = (I + factor*laplacian)(a_phi)
    if (a_applyBC) {
        if ( (a_phiC == NULL) || (a_level == 0)) {
            m_ops[a_level]->applyOp(a_ans, a_phi, false);
        } else {
            m_ops[a_level]->AMROperatorNF(a_ans, a_phi, *a_phiC, false);
        }
    } else {
        m_ops[a_level]->applyOpNoBoundary(a_ans, a_phi, false);
    }
}


// -----------------------------------------------------------------------------
// Applies the Helmholtz operator to the solution a_phi at the given
// grid level. This will set a_ans to (I + a_mu * a_dt * L(a_phi).
// -----------------------------------------------------------------------------
void MappedBaseLevelHeatSolver::applyHelm (LevelDataType&          a_ans,
                                           const LevelDataType&    a_phi,
                                           const LevelDataType*    a_phiC,
                                           int                     a_level,
                                           Real                    a_mu,
                                           Real                    a_dt,
                                           bool                    a_homogeneousBC)
{
    // Set the operator's alpha and beta coefficients.
    Real factor  = a_mu * a_dt;
    m_ops[a_level]->setAlphaAndBeta(1.0, factor);

    LevelDataType zero;
    m_ops[a_level]->create(zero, a_ans);
    m_ops[a_level]->setToZero(zero);

    // set a_ans = helm(a_phi)
    //           = (I + factor*laplacian)(a_phi)
    if ( (a_phiC == NULL) || (a_level == 0)) {
        m_ops[a_level]->applyOp(a_ans, a_phi, a_homogeneousBC);
    } else {
        m_ops[a_level]->AMROperatorNF(a_ans, a_phi, *a_phiC, a_homogeneousBC);
    }
}


// -----------------------------------------------------------------------------
// Adds flux contributions from the Helmholtz operator at the current
// grid level.
// -----------------------------------------------------------------------------
void MappedBaseLevelHeatSolver::incrementFlux (LevelData<FluxDataType>& a_diffusiveFlux,
                                               LevelDataType&           a_phi,
                                               int                      a_level,
                                               Real                     a_mu,
                                               Real                     a_dt,
                                               Real                     a_sign,
                                               bool                     a_setToZero)
{
    Real factor  = a_sign * a_dt * a_mu;
    m_ops[a_level]->setAlphaAndBeta(1.0, factor);

    // increment flux
    m_ops[a_level]->fillGrad(a_phi);
    for (DataIterator  dit = a_phi.dataIterator(); dit.ok(); ++dit) {
        FluxDataType& thisFlux = a_diffusiveFlux[dit];
        FluxDataType tempFlux;
        tempFlux.define(thisFlux);

        tempFlux.setVal(0.0);
        if (a_setToZero) {
            thisFlux.setVal(0.0);
        }

        m_ops[a_level]->getFlux(tempFlux, a_phi, m_grids[a_level][dit], dit(), 1.0);
        thisFlux += tempFlux;
    }
}


// -----------------------------------------------------------------------------
// Solves the Helmholtz equation (I - a_mu * a_dt * L(a_phi) = a_rhs
// for a_phi. Here it is assumed that a solution a_phiC exists
// on a coarser grid level.
// -----------------------------------------------------------------------------
void MappedBaseLevelHeatSolver::solveHelm (LevelDataType& a_phi,
                                           LevelDataType& a_phiC,
                                           LevelDataType& a_rhs,
                                           int            a_level,
                                           Real           a_mu,
                                           Real           a_dt,
                                           bool           a_zeroPhi)
{
    if (a_zeroPhi) {
        m_ops[a_level]->setToZero(a_phi);
    }
    Vector<LevelDataType* > phi(m_grids.size(), NULL);
    Vector<LevelDataType* > rhs(m_grids.size(), NULL);
    phi[a_level] = &a_phi;
    rhs[a_level] = &a_rhs;
    if (a_level > 0) {
        phi[a_level - 1] = &a_phiC;
    }

    Real factor  = -a_dt * a_mu;
    resetSolverAlphaAndBeta(1.0, factor);

    m_solver->solve(phi, rhs, a_level, a_level, a_zeroPhi);
    int solverExitStatus = m_solver->m_exitStatus;
    if (solverExitStatus == 2 || solverExitStatus == 4 || solverExitStatus == 6) {
        // These status codes correspond to the cases in which
        // norm is neither small enough nor reduced enough.
        // Either we've reached the maximum number of iterations,
        // or we've hung.
        pout() << "BaseLevelTGA:: WARNING: solver exitStatus == "
               << solverExitStatus << std::endl;
    }
}


// -----------------------------------------------------------------------------
// Sets the alpha and beta parameters in each Helmholtz operator to the
// given values.
// -----------------------------------------------------------------------------
void MappedBaseLevelHeatSolver::resetSolverAlphaAndBeta (const Real& a_alpha,
                                                         const Real& a_beta)
{
    Vector<MappedMGLevelOp<LevelDataType>* > ops = m_solver->getAllOperators();
    for (int iop = 0; iop < ops.size(); iop++) {
        MappedLevelTGAHelmOp<LevelDataType, FluxDataType>* helmop =
            dynamic_cast<MappedLevelTGAHelmOp<LevelDataType, FluxDataType>*>(ops[iop]);
        helmop->setAlphaAndBeta(a_alpha, a_beta);
    }
    for (int iop = 0; iop < m_ops.size(); iop++) {
        m_ops[iop]->setAlphaAndBeta(a_alpha, a_beta);
    }
}


// -----------------------------------------------------------------------------
// Creates a new Helmholtz operator for use by the TGA integrator.
// -----------------------------------------------------------------------------
MappedLevelTGAHelmOp<MappedBaseLevelHeatSolver::LevelDataType, MappedBaseLevelHeatSolver::FluxDataType>*
MappedBaseLevelHeatSolver::newOp (const ProblemDomain&                                    a_indexSpace,
                                  RefCountedPtr<MappedAMRLevelOpFactory<LevelDataType> >& a_opFact)
{
    MappedLevelTGAHelmOp<LevelDataType, FluxDataType>* retval =
        dynamic_cast<MappedLevelTGAHelmOp<LevelDataType, FluxDataType>*>(a_opFact->AMRnewOp(a_indexSpace));

    return retval;
}


// -----------------------------------------------------------------------------
// Returns the number of grid levels on which this integrator operates.
// -----------------------------------------------------------------------------
int MappedBaseLevelHeatSolver::size () const
{
    return m_grids.size();
}


// -----------------------------------------------------------------------------
// Interpolates a given quantity linearly in time using its beginning- and
// end-of-step values and placing the result in \a a_data.
// -----------------------------------------------------------------------------
void MappedBaseLevelHeatSolver::timeInterp (LevelDataType&       a_data,
                                            const LevelDataType& a_oldData,
                                            const LevelDataType& a_newData,
                                            Real                 a_time,
                                            Real                 a_oldTime,
                                            Real                 a_newTime,
                                            int                  a_level)
{
    Real eps = 1.0e-10;
    CH_assert(a_newTime >= a_oldTime);
    Real diff = (a_newTime - a_oldTime);
    this->m_ops[a_level]->setToZero(a_data);
    if (diff < eps) {
        //no real time advance and don't want to divide by zero
        this->m_ops[a_level]->incr(a_data, a_oldData, 1.0);
    } else {
        Real factor = (a_time - a_oldTime) / (a_newTime - a_oldTime);
        this->m_ops[a_level]->incr(a_data, a_oldData, 1.0 - factor);
        this->m_ops[a_level]->incr(a_data, a_newData, factor);
    }
}


// -----------------------------------------------------------------------------
// Fills ghosts cells on source terms.
// -----------------------------------------------------------------------------
void MappedBaseLevelHeatSolver::setSourceGhostCells (LevelData<FArrayBox>&    a_src,
                                                     const DisjointBoxLayout& a_grids,
                                                     int a_lev)
{
    DataIterator dit = a_grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& srcFAB = a_src[dit];
        const Box& srcBox = srcFAB.box();
        const Box& valid  = a_grids[dit];

        for (int idir = 0; idir < SpaceDim; ++idir) {
            for (SideIterator sit; sit.ok(); ++sit) {
                const int isign = sign(sit());

                // Create the dest ghost box -- include corners too.
                Box destBox = adjCellBox(valid, idir, sit(), 1);
                for (int jdir = 0; jdir < SpaceDim; jdir++) {
                    if (jdir != idir) {
                        destBox.grow(jdir, 1);
                    }
                }

                // If this fails, srcFAB may not have ghost cells.
                CH_assert(srcBox.contains(destBox));

                // Use high order wherever possible.
                int extrapOrder = 3;
                if (valid.size(idir) < 4) --extrapOrder;

                // Do the extrapolation
                FORT_EXTRAPOLATEFACENOEV (
                    CHF_FRA(srcFAB),
                    CHF_FRA(srcFAB),
                    CHF_BOX(destBox),
                    CHF_CONST_INT(idir),
                    CHF_CONST_INT(isign),
                    CHF_CONST_INT(extrapOrder));

            } // end loop over sides (sit)
        } // end loop over directions (dir)
    } // end loop over grids (dit)
}
