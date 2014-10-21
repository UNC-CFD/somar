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
#include "AMRLepticSolver.H"


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
AMRLepticSolver::AMRLepticSolver ()
: m_eps(1E-6),
  m_hang(1E-15),
  m_normThresh(1E-30),
  m_imin(5),
  m_iterMax(20),
  m_verbosity(3),
  m_numMG(1),
  m_convergenceMetric(0.)
{
    this->clear();
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
AMRLepticSolver::~AMRLepticSolver ()
{
    this->clear();
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void AMRLepticSolver::setSolverParameters (const int  a_numMG,
                                           const int  a_iterMax,
                                           const Real a_eps,
                                           const Real a_hang,
                                           const Real a_normThresh)
{
    m_eps        =    a_eps;
    m_hang       =    a_hang;
    m_normThresh =    a_normThresh;
    m_iterMax    =    a_iterMax;
    m_numMG = a_numMG;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void AMRLepticSolver::define (const ProblemDomain&                            a_coarseDomain,
                              MappedAMRLevelOpFactory<LevelData<FArrayBox> >& a_factory,
                              const int                                       a_maxAMRLevels,
                              const int                                       a_verbosity)
{
    CH_TIME("AMRLepticSolver::define");

    m_verbosity = a_verbosity;

    this->clear();
    m_op.resize(a_maxAMRLevels, NULL);

    m_correction.resize(a_maxAMRLevels, NULL);
    m_residual.  resize(a_maxAMRLevels, NULL);
    m_resC.      resize(a_maxAMRLevels, NULL);
    m_resCopier. resize(a_maxAMRLevels);
    m_reverseCopier.resize(a_maxAMRLevels);

    ProblemDomain current = a_coarseDomain;
    for (int i = 0; i < a_maxAMRLevels; i++) {
        m_correction[i] = new LevelData<FArrayBox>();
        m_residual[i]   = new LevelData<FArrayBox>();
        m_resC[i]       = new LevelData<FArrayBox>();

        if (m_verbosity >= 5) {
            pout() << "calling a_factory.AMRnewOp(current) with lev = " << i
                   << " and index space = " << current.domainBox() << endl;
        }

        m_op[i] = a_factory.AMRnewOp(current);
        CH_assert(m_op[i] != NULL);

        // Only do this if it will be used (avoiding a reference to invalid
        // and/or unavailable refinement ratios)
        if (i < a_maxAMRLevels - 1) {
            ProblemDomain notSoCurrent = current;
            refine(current, notSoCurrent, a_factory.getFineRefRatio(current));
        }
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void AMRLepticSolver::clear ()
{
    CH_assert(m_correction.size() == m_op.size());
    CH_assert(m_residual.size() == m_op.size());
    CH_assert(m_resC.size() == m_op.size());

    // Free memory
    for (int i = 0; i < m_op.size(); i++) {
        m_op[i]->clear(*m_correction[i]);
        m_op[i]->clear(*m_residual[i]);
        m_op[i]->clear(*m_resC[i]);

        delete m_correction[i];
        delete m_residual[i];
        delete m_resC[i];
        delete m_op[i];

        m_correction[i] = NULL;
        m_residual[i] = NULL;
        m_resC[i] = NULL;
        m_op[i] = NULL;
    }
    m_correction.clear();
    m_residual.clear();
    m_resC.clear();
    m_op.clear();

    for (int l = 0; l < m_amrLepticSolver.size(); ++l) {
        delete m_amrLepticSolver[l];
        m_amrLepticSolver[l] = NULL;
    }
    m_amrLepticSolver.clear();
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void AMRLepticSolver::init (const Vector<LevelData<FArrayBox>*>& a_phi,
                            const Vector<LevelData<FArrayBox>*>& a_rhs,
                            const int                            l_max,
                            const int                            l_base)
{
    for (int i = l_base; i <= l_max; ++i) {
        MappedAMRLevelOp<LevelData<FArrayBox> >& op = *(m_op[i]);

        op.create(*m_correction[i], *a_phi[i]);
        op.create(*m_residual[i],   *a_rhs[i]);

        if (i != l_base) {
            const IntVect r = op.refToCoarser();
            op.createCoarsened(*m_resC[i], *a_rhs[i], r);
            op.buildCopier(m_resCopier[i], *a_rhs[i-1], *m_resC[i]);
            m_reverseCopier[i] = m_resCopier[i];
            m_reverseCopier[i].reverse();
        }
    }

    // Allocate and define the leptic solvers
    for (int l = 0; l < m_amrLepticSolver.size(); ++l) {
        delete m_amrLepticSolver[l];
        m_amrLepticSolver[l] = NULL;
    }
    m_amrLepticSolver = Vector<LevelLepticSolver*>(l_max+1, NULL);
    for (int l = l_base; l <= l_max; ++l) {
        m_amrLepticSolver[l] = new LevelLepticSolver;
        m_amrLepticSolver[l]->define(m_op[l], true);    // homogeneous
        m_amrLepticSolver[l]->setCrsePhiPtr(NULL);      // Just to be explicit
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void AMRLepticSolver::solve (Vector<LevelData<FArrayBox>*>&       a_phi,
                             const Vector<LevelData<FArrayBox>*>& a_rhs,
                             const int                            l_max,
                             const int                            l_base,
                             const bool                           a_zeroPhi,
                             const bool                           a_forceHomogeneous)
{
    CH_TIME("AMRLepticSolver::solve");
    init(a_phi, a_rhs, l_max, l_base);
    solveNoInit(a_phi, a_rhs, l_max, l_base, a_zeroPhi, a_forceHomogeneous);
    // AMRMultiGrid calls revert() here.
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void AMRLepticSolver::solveNoInit (Vector<LevelData<FArrayBox>*>&       a_phi,
                                   const Vector<LevelData<FArrayBox>*>& a_rhs,
                                   const int                            l_max,
                                   const int                            l_base,
                                   const bool                           a_zeroPhi,
                                   const bool                           a_forceHomogeneous)
{
    Vector<LevelData<FArrayBox>*> uberResidual(a_rhs.size());
    int lowlim = l_base;
    if (l_base > 0)  // we need an ubercorrection one level lower than l_base
        --lowlim;
    for (int ilev = lowlim; ilev <= l_max; ilev++) {
        uberResidual[ilev] = new LevelData<FArrayBox>();
        if (ilev >= l_base) {
            m_op[ilev]->create(*uberResidual[ilev], *a_rhs[ilev]);
        }
    }
    solveNoInitResid(a_phi, uberResidual, a_rhs, l_max, l_base, a_zeroPhi, a_forceHomogeneous);
    for (int i = lowlim; i <= l_max; i++) {
        m_op[i]->clear(*uberResidual[i]);
        delete uberResidual[i];
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void AMRLepticSolver::solveNoInitResid (Vector<LevelData<FArrayBox>*>&       a_phi,
                                        Vector<LevelData<FArrayBox>*>&       uberResidual,
                                        const Vector<LevelData<FArrayBox>*>& a_rhs,
                                        const int                            l_max,
                                        const int                            l_base,
                                        const bool                           a_zeroPhi,
                                        const bool                           a_forceHomogeneous)
{
    CH_TIMERS("AMRLepticSolver::solveNoInit");
    CH_TIMER("AMRLepticSolver::AMRVcycle", vtimer);

    CH_assert(l_base <= l_max);
    CH_assert(a_rhs.size() == a_phi.size());

    //these correspond to the residual and correction
    //that live in AMRSolver
    Vector<LevelData<FArrayBox>*> uberCorrection(a_rhs.size());

    // This stores the best solution
    Vector<LevelData<FArrayBox>*> bestPhi(a_phi.size());

    int lowlim = l_base;
    if (l_base > 0)  // we need an ubercorrection one level lower than l_base
        lowlim--;
    bool outputIntermediates = false;

    for (int ilev = lowlim; ilev <= l_max; ilev++) {
        uberCorrection[ilev] = new LevelData<FArrayBox>();
        m_op[ilev]->create(*uberCorrection[ilev], *a_phi[ilev]);
        if (ilev >= l_base) {
            m_op[ilev]->create(*uberResidual[ilev], *a_rhs[ilev]);
            m_op[ilev]->setToZero(*(uberResidual[ilev]));
        }
        m_op[ilev]->setToZero(*(uberCorrection[ilev]));

        // Allocate bestPhi
        bestPhi[ilev] = new LevelData<FArrayBox>();
        m_op[ilev]->create(*bestPhi[ilev], *a_phi[ilev]);
    }
    //  m_op[0]->dumpStuff(uberResidual, string("initialRes.hdf5"));
    if (a_zeroPhi)
        for (int ilev = l_base; ilev <= l_max; ++ilev) {
            m_op[ilev]->setToZero(*(a_phi[ilev]));
        }
    //compute initial residual and initialize internal residual to it

    // Set bestPhi equal to initial guess
    for (int ilev = lowlim; ilev <= l_max; ilev++) {
        m_op[ilev]->assign(*bestPhi[ilev], *a_phi[ilev]);
    }

    Real initial_rnorm = 0;
    {
        CH_TIME("Initial AMR Residual");
        initial_rnorm = computeAMRResidual(uberResidual, a_phi, a_rhs, l_max, l_base, a_forceHomogeneous, true);
    }

    if (m_convergenceMetric != 0.) {
        initial_rnorm = m_convergenceMetric;
    }

    Real rnorm = initial_rnorm;
    Real norm_last = 2 * initial_rnorm;
    Real best_rnorm = rnorm;
    bool useBestPhi = false;
    bool somethingConverged = false;  // This remains false if no iters converge to a sol'n better than what was given to us.

    int iter = 0;
    if (m_verbosity == 2) {
        pout() << "    AMRLepticSolver:: " << std::flush;
    } else if (m_verbosity > 2) {
        pout() << "    AMRLepticSolver:: iteration = "
               << std::fixed << iter << ", residual norm = "
               << std::scientific << rnorm << std::endl;
    }

    bool goNorm = rnorm > m_normThresh;                        //iterate if norm is not small enough
    bool goRedu = rnorm > m_eps * initial_rnorm;               //iterate if initial norm is not reduced enough
    bool goIter = iter < m_iterMax;                            //iterate if iter < max iteration count
    bool goHang = iter < m_imin || rnorm < (1 - m_hang) * norm_last; //iterate if we didn't hang

    while (goIter && goRedu && goHang && goNorm) {
        if (outputIntermediates) {
            char strresname[100];
            sprintf(strresname, "amrleptic.res.iter.%03d.hdf5", iter);
            string nameres(strresname);
            outputAMR(uberResidual, nameres, l_max, l_base);
        }

        norm_last = rnorm;

        //this generates a correction from the current residual
        CH_START(vtimer);
        AMRVCycle(uberCorrection, uberResidual, l_max, l_max, l_base);
        CH_STOP(vtimer);

        // Do post VCycle stuff
        rnorm = postVCycleOps(uberResidual, uberCorrection, a_phi, a_rhs, l_max, l_base, a_forceHomogeneous);
        iter++;

        // Check if residual has dropped.
        if (rnorm <= best_rnorm) {
            // OK. Save this as the best solution and move on.
            best_rnorm = rnorm;
            for (int ilev = l_base; ilev <= l_max; ++ilev) {
                m_op[ilev]->assign(*bestPhi[ilev], *a_phi[ilev]);
            }
            useBestPhi = false;
            somethingConverged = true;
        } else {
            // Bad. VCycle did not produce a better solution.
            useBestPhi = true;
        }

        // Let the user know what happened
        if (m_verbosity >= 3) { ////
            if (useBestPhi) {
                pout() << "   [D] ";   // Mark iterations that have diverged.
            } else {
                pout() << "       ";
            }
            pout() << "AMRLepticSolver:: iteration = "
                   << std::fixed << iter << ", residual norm = "
                   << std::scientific << rnorm;
            if (rnorm > 0.0) {
                pout() << ", rate = " << std::scientific << norm_last / rnorm;
            }
            pout() << std::endl;
        }

        goNorm = rnorm > m_normThresh;                        //keep iterating if norm is not small enough
        goRedu = rnorm > m_eps * initial_rnorm;               //keep iterating if initial norm is not reduced enough
        goIter = iter < m_iterMax;                            //keep iterating if iter < max iteration count
        goHang = iter < m_imin || rnorm < (1 - m_hang) * norm_last; //keep iterating if we didn't hang
    }

    // Use the best solution available
    if (useBestPhi) {
        rnorm = best_rnorm;
        for (int ilev = l_base; ilev <= l_max; ++ilev) {
            m_op[ilev]->assign(*a_phi[ilev], *bestPhi[ilev]);
        }
    }

    // Did we blow up?
    if ((rnorm > 10.*initial_rnorm) && (rnorm > 10.*m_eps)) {
        pout() << "solver seems to have blown up" << endl;
        MayDay::Error("kaboom");

        // TODO: Maybe we shouldn't just give up.
    }

    if (!somethingConverged && rnorm >= initial_rnorm && rnorm >= m_eps) {
        pout() << "AMRLepticSolver solver blew up." << endl;
        std::cout << "AMRLepticSolver solver blew up." << endl;
        // MayDay::Error("AMRLepticSolver solver blew up");
    }

    // The solver has finished. Figure out the final state of the solution.
    m_exitStatus = int(!goRedu) + int(!goIter) * 2 + int(!goHang) * 4 + int(!goNorm) * 8;

    if (m_verbosity == 2) {
        if (initial_rnorm != 0.0) {
            pout() << std::fixed << iter << " iters\t rel res = "
                   << std::scientific << rnorm/initial_rnorm << std::endl;
        } else {
            pout() << std::fixed << iter << " iters\t rel res = "
                   << std::scientific << rnorm << " / " << initial_rnorm << std::endl;
        }
    } else if (m_verbosity > 2) {
        pout() << "    AMRLepticSolver:: iteration = " << std::fixed << iter
               << ", residual norm = " << std::scientific << rnorm << std::endl;
    }

    if (m_verbosity > 1) {
        if (!goIter && goRedu && goNorm) { // goRedu=T, goIter=F, goHang=?, goNorm=T
            // m_exitStatus == 0 + 2 + 0|4 + 0 = 2|6
            pout() << "    AMRLepticSolver:: WARNING: Exit because max iteration count exceeded" << std::endl;
        }
        if (!goHang && goRedu && goNorm) { // goRedu=T, goIter=?, goHang=F, goNorm=T
            // m_exitStatus == 0 + 0|2 + 4 + 0 = 4|6
            pout() << "    AMRLepticSolver:: WARNING: Exit because of solver hang" << std::endl;
        }
        if (m_verbosity > 4) {
            pout() << "    AMRLepticSolver:: exitStatus = " << m_exitStatus << std::endl;
        }
    }

    // Clean up after ourselves
    for (int i = lowlim; i <= l_max; i++) {
        m_op[i]->clear(*uberCorrection[i]);
        delete uberCorrection[i];
        delete bestPhi[i];
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void AMRLepticSolver::AMRVCycle (Vector<LevelData<FArrayBox>*>&       a_uberCorrection,
                                 const Vector<LevelData<FArrayBox>*>& a_uberResidual,
                                 const int                            ilev,
                                 const int                            l_max,
                                 const int                            l_base)
{
    if (ilev == l_max) {
        for (int level = l_base; level <= l_max; level++) {
            m_op[level]->assignLocal(*m_residual[level], *a_uberResidual[level]);
            m_op[level]->setToZero(*m_correction[level]);
        }
    }

    if (l_max == l_base) {
        CH_assert(ilev == l_base);
        if (m_verbosity >= 5) pout() << "levelLepticSolver " << ilev << "\n";
        m_amrLepticSolver[l_base]->solve(*(a_uberCorrection[ilev]), *(a_uberResidual[ilev]));
    } else if (ilev == l_base) {
        if (m_verbosity >= 5) pout() << "levelLepticSolver " << ilev << "\n";
        m_amrLepticSolver[l_base]->solve(*(a_uberCorrection[ilev]), *(a_uberResidual[ilev]));

        m_op[ilev]->incr(*(a_uberCorrection[ilev]), *(m_correction[ilev]), 1.0);
    } else {
        //============= Downsweep ========================

        if (m_verbosity >= 5) pout() << "Calling LevelLepticSolver on level " << ilev << "...\n";
        m_amrLepticSolver[ilev]->solve(*(m_correction[ilev]), *(m_residual[ilev]));
        m_op[ilev]->incr(*(a_uberCorrection[ilev]), *(m_correction[ilev]), 1.0);

        // Set next coarser level correction to zero
        m_op[ilev - 1]->setToZero(*(m_correction[ilev - 1]));

        // Recompute residual on next coarser level
        //  for the valid region NOT covered by this level.
        computeAMRResidualLevel(m_residual,
                                a_uberCorrection,
                                a_uberResidual,
                                l_max, l_base, ilev - 1,
                                true);

        // Compute the restriction of the residual to the coarser level resC.
        if (m_verbosity >= 5) pout() << "restrict residual " << ilev << " to " << ilev - 1 << "\n";
        m_op[ilev]->AMRRestrictS(*(m_resC[ilev]),
                                 *(m_residual[ilev]),
                                 *(m_correction[ilev]),
                                 *(m_correction[ilev - 1]),
                                 *(a_uberCorrection[ilev]));

        // Overwrite residual on the valid region of the next coarser level
        //  with coarsened residual from this level
        m_op[ilev - 1]->assignCopier(*m_residual[ilev - 1], *(m_resC[ilev]), m_resCopier[ilev]);

        //============finish Compute residual for the next coarser level======

        for (int img = 0; img < m_numMG; img++) {
            AMRVCycle(a_uberCorrection, a_uberResidual, ilev - 1, l_max, l_base);
        }

        //================= Upsweep ======================
        //increment the correction with coarser version
        if (m_verbosity >= 5) pout() << "Prolong correction " << ilev - 1 << " to " << ilev << "\n";
        m_op[ilev]->AMRProlongS(*(m_correction[ilev]), *(m_correction[ilev - 1]),
                                *m_resC[ilev], m_reverseCopier[ilev]);

        //recompute residual
        m_op[ilev]->AMRUpdateResidual(*(m_residual[ilev]), *(m_correction[ilev]), *(m_correction[ilev - 1]));

        //compute correction to the correction
        if (m_verbosity >= 5) pout() << "Calling relax on level " << ilev << "...\n";
        LevelData<FArrayBox>& dCorr = *(a_uberCorrection[ilev]); // user uberCorrection as holder for correction to correction
        m_op[ilev]->setToZero(dCorr);
        m_amrLepticSolver[ilev]->solve(dCorr, *(m_residual[ilev]));

        //correct the correction with the correction to the correction
        if (m_verbosity >= 5) {
            LevelData<FArrayBox> newRes;
            m_op[ilev]->create(newRes, *(m_correction[ilev]));
            m_op[ilev]->residual(newRes, *(m_correction[ilev]), *(m_residual[ilev]), true);

            LevelData<FArrayBox> fineRes;
            IntVect ref = IntVect::Unit;
            Real rnorm = m_op[ilev]->AMRNorm(newRes, fineRes, ref, 0);

            pout() << "Pre residual max norm = " << rnorm << endl;
        }
        m_op[ilev]->incr(*(m_correction[ilev]), dCorr, 1.0);
        if (m_verbosity >= 5) {
            LevelData<FArrayBox> newRes;
            m_op[ilev]->create(newRes, *(m_correction[ilev]));
            m_op[ilev]->residual(newRes, *(m_correction[ilev]), *(m_residual[ilev]), true);

            LevelData<FArrayBox> fineRes;
            IntVect ref = IntVect::Unit;
            Real rnorm = m_op[ilev]->AMRNorm(newRes, fineRes, ref, 0);

            pout() << "Post residual max norm = " << rnorm << endl;
        }

        m_op[ilev]->assignLocal(*(a_uberCorrection[ilev]), *(m_correction[ilev]));
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
Real AMRLepticSolver::postVCycleOps (Vector<LevelData<FArrayBox>*>&       a_uberResidual,
                                     Vector<LevelData<FArrayBox>*>&       a_uberCorrection,
                                     Vector<LevelData<FArrayBox>*>&       a_phi,
                                     const Vector<LevelData<FArrayBox>*>& a_rhs,
                                     const int                            l_max,
                                     const int                            l_base,
                                     const bool                           a_forceHomogeneous)
{
    // Increment phi by correction and reset correction to zero
    for (int ilev = l_base; ilev <= l_max; ilev++) {
        m_op[ilev]->incr(*(a_phi[ilev]), *(a_uberCorrection[ilev]), 1.0);
        m_op[ilev]->setToZero(*(a_uberCorrection[ilev]));
    }


    // For solvers with accuracy higher than 2nd order
    //  consistency between levels has to be explicitly enforced.
    // Qinghai Zhang
    if (m_op[0]->orderOfAccuracy() > 2) {
        for (int ilev = l_max; ilev > l_base; ilev--) {
            m_op[ilev]->enforceCFConsistency(*a_phi[ilev - 1], *a_phi[ilev]);
        }
    }
    //------end enforcing consistency.

    // recompute residual
    return computeAMRResidual(a_uberResidual, a_phi, a_rhs, l_max, l_base, a_forceHomogeneous, true);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void AMRLepticSolver::computeAMRResidualLevel (Vector<LevelData<FArrayBox>*>&       a_resid,
                                               const Vector<LevelData<FArrayBox>*>& a_phi,
                                               const Vector<LevelData<FArrayBox>*>& a_rhs,
                                               const int                            l_max,
                                               const int                            l_base,
                                               const int                            ilev,
                                               const bool                           a_homogeneousBC)
{
    CH_TIME("AMRLepticSolver::computeAMRResidualLevel");

    //m_op[ilev]->setToZero(*(a_resid[l_max]));
    if (l_max != l_base) {
        if (ilev == l_max) {
            m_op[l_max]->AMRResidualNF(*(a_resid[l_max]), *(a_phi[l_max]),
                                       *(a_phi[l_max - 1]), *(a_rhs[l_max]),
                                       a_homogeneousBC);
        } else if (ilev == l_base) {
            if (l_base == 0) {
                m_op[l_base]->AMRResidualNC(*(a_resid[l_base]), *(a_phi[l_base + 1]),
                                            *(a_phi[l_base]),  *(a_rhs[l_base]),
                                            a_homogeneousBC, m_op[l_base + 1]);
            } else {
                m_op[l_base]->AMRResidual(*a_resid[l_base], *a_phi[l_base + 1], *a_phi[l_base],
                                          *a_phi[l_base - 1], *a_rhs[l_base],
                                          a_homogeneousBC, m_op[l_base + 1]);
            }
        } else {
            m_op[ilev]->AMRResidual(*a_resid[ilev], *a_phi[ilev + 1], *a_phi[ilev],
                                    *a_phi[ilev - 1], *a_rhs[ilev],
                                    a_homogeneousBC, m_op[ilev + 1]);
        }
    } else {
        CH_assert(ilev == l_base);
        if (l_base == 0) {
            m_op[l_max]->residual(*a_resid[l_max], *a_phi[l_max], *a_rhs[l_max], a_homogeneousBC);
        } else {
            m_op[l_max]->AMRResidualNF(*(a_resid[l_max]), *(a_phi[l_max]),
                                       *(a_phi[l_max - 1]), *(a_rhs[l_max]),
                                       a_homogeneousBC);
        }
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
Real AMRLepticSolver::computeAMRResidual (Vector<LevelData<FArrayBox>*>&       a_resid,
                                          const Vector<LevelData<FArrayBox>*>& a_phi,
                                          const Vector<LevelData<FArrayBox>*>& a_rhs,
                                          const int                            l_max,
                                          const int                            l_base,
                                          const bool                           a_homogeneousBC,
                                          const bool                           a_computeNorm)
{
    CH_TIME("AMRLepticSolver::computeAMRResidual");

    Real rnorm = 0;
    Real localNorm = 0;
    for (int ilev = l_base; ilev <= l_max; ilev++) {
        //always used at top level where bcs are inhomogeneous
        computeAMRResidualLevel(a_resid,
                                a_phi,
                                a_rhs,
                                l_max, l_base, ilev, a_homogeneousBC);
        if (a_computeNorm) {
            if (ilev == l_max) {
                localNorm = m_op[ilev]->localMaxNorm(*a_resid[ilev]);
            } else {
                m_op[ilev]->zeroCovered(*a_resid[ilev], *m_resC[ilev + 1], m_resCopier[ilev + 1]);
                localNorm = m_op[ilev]->localMaxNorm(*a_resid[ilev]);
            }
            rnorm = Max(localNorm, rnorm);
        }
    }
#ifdef CH_MPI
    if (a_computeNorm) {
        CH_TIME("MPI_Allreduce");
        Real recv;
        int result = MPI_Allreduce(&rnorm, &recv, 1, MPI_CH_REAL,
                                   MPI_MAX, Chombo_MPI::comm);
        if (result != MPI_SUCCESS) {
            //bark!!!
            MayDay::Error("sorry, but I had a communcation error on norm");
        }
        rnorm = recv;
    }
#endif

    return rnorm; // if a_computeNorm is false, then this just returns zero.
}


// -----------------------------------------------------------------------------
// Write a_data to HDF5.
// -----------------------------------------------------------------------------
void AMRLepticSolver::outputAMR (const Vector<LevelData<FArrayBox>*>& a_data,
                                 const std::string                    a_name,
                                 const int                            a_lmax,
                                 const int                            a_lbase)
{
    Vector<LevelData<FArrayBox>*> outputData;
    for (int ilev = a_lbase; ilev <= a_lmax; ilev++) {
        outputData.push_back(a_data[ilev]);
    }
    m_op[a_lbase]->outputAMR(outputData, a_name);
}
