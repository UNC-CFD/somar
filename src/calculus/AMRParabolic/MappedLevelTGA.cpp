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
#include "MappedLevelTGA.H"


// -----------------------------------------------------------------------------
// Full constructor
// -----------------------------------------------------------------------------
MappedLevelTGA::MappedLevelTGA (const Vector<DisjointBoxLayout>&                                 a_grids,
                                const Vector<IntVect>&                                           a_refRat,
                                const ProblemDomain&                                             a_level0Domain,
                                RefCountedPtr<MappedAMRLevelOpFactory<LevelData<FArrayBox> > >&  a_opFact,
                                const RefCountedPtr<MappedAMRMultiGrid<LevelData<FArrayBox> > >& a_solver)
: MappedBaseLevelHeatSolver(a_grids, a_refRat, a_level0Domain, a_opFact, a_solver),
  m_mu1(0.0),
  m_mu2(0.0),
  m_mu3(0.0),
  m_mu4(0.0),
  m_r1(0.0)
{
    Real tgaEpsilon = 1.e-12;
#ifdef CH_USE_FLOAT
    tgaEpsilon = sqrt(tgaEpsilon);
#endif
    Real a = 2.0 - sqrt(2.0) - tgaEpsilon;
    m_mu1 = (a - sqrt( a * a - 4.0 * a + 2.0)) / 2.0 ;
    m_mu2 = (a + sqrt( a * a - 4.0 * a + 2.0)) / 2.0 ;
    m_mu3 = (1.0 - a);
    m_mu4 = 0.5 - a;

    Real discr = sqrt(a * a - 4.0 * a + 2.0);
    m_r1 = (2.0 * a - 1.0) / (a + discr);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
MappedLevelTGA::~MappedLevelTGA ()
{;}


// -----------------------------------------------------------------------------
// Computes the time-centered diffusion term L(phi). This can be used to
// find contributions to the solution from diffusion. The diffusion term
// is computed by computing a finite difference approximation for d phi/dt
// using the updated and original values of phi and the time step. Most of
// the arguments given here are passed along to updateSoln and retain their
// significance therein.
// -----------------------------------------------------------------------------
void MappedLevelTGA::computeDiffusion (LevelDataType&           a_diffusiveTerm,
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
    // in this function, we call the updateSoln function, compute
    // the update, then subtract off the extra pieces to return the
    // diffusive part of the update

    if (!this->m_ops[a_level]->isTimeDependent()) {
        // The operator has no time-dependent parameters. Life is easier.

        // first compute updated solution
        LevelDataType phiNew;

        this->m_ops[a_level]->create(phiNew, a_phiOld);
        this->m_ops[a_level]->setToZero(phiNew);
        if (!a_zeroPhi) {
            this->m_ops[a_level]->assign(phiNew, a_phiOld);
        }

        updateSoln(phiNew, a_phiOld, a_src, a_flux,
                   a_fineFluxRegPtr, a_crseFluxRegPtr,
                   a_crsePhiOldPtr, a_crsePhiNewPtr,
                   a_oldTime, a_crseOldTime,
                   a_crseNewTime, a_dt, a_level, a_zeroPhi,
                   a_rhsAlreadyKappaWeighted, a_fluxStartComponent);

        // now subtract everything off to leave us with diffusive term
        this->m_ops[a_level]->incr(phiNew, a_phiOld, -1.0);
        this->m_ops[a_level]->scale(phiNew, 1.0 / a_dt);

        //now multiply by a if there is an a
        this->m_ops[a_level]->diagonalScale(phiNew, false);

        // and finally, subtract off a_src
        this->m_ops[a_level]->incr(phiNew, a_src, -1.0);

        // what's left should be the time-centered diffusive part of the update
        this->m_ops[a_level]->assign(a_diffusiveTerm, phiNew);
    } else {
        // The operator has time-dependent coefficients. We must be more careful!
        MayDay::Error("Is this ever called???");
        MayDay::Error("This hasn't been updated to deal with a_fluxStartComponent");

        //                                                      n+1    n
        //                                                n   (a    - a )
        // There's an extra source term (phi (da/dt) = phi  * -----------
        //                                                       dt
        // from the time-changing density that we need to subtract from the RHS.
        LevelDataType phidadt, aOld, aNew, rhs;
        this->m_ops[a_level]->create(phidadt, a_phiOld);
        this->m_ops[a_level]->create(aOld, a_phiOld);
        this->m_ops[a_level]->create(aNew, a_phiOld);
        this->m_ops[a_level]->create(rhs, a_phiOld);
        for (DataIterator dit = a_phiOld.disjointBoxLayout().dataIterator(); dit.ok(); ++dit) {
            // Set the old and new a coefficients.
            this->m_ops[a_level]->setTime(a_oldTime, 0.0, a_dt);
            aOld[dit()].setVal(1.);
            this->m_ops[a_level]->diagonalScale(aOld);
            this->m_ops[a_level]->setTime(a_oldTime, 1.0, a_dt);
            aNew[dit()].setVal(1.);
            this->m_ops[a_level]->diagonalScale(aNew);

            // Compute the above term.
            phidadt[dit()].axby(aNew[dit()], aOld[dit()], 1.0 / a_dt, -1.0 / a_dt);
            phidadt[dit()] *= a_phiOld[dit()];

            // Make a new right hand side out of the difference between the
            // source term and phidadt.
            rhs[dit()].axby(a_src[dit()], phidadt[dit()], 1.0, -1.0);
        }

        // Compute the updated solution.
        LevelDataType phiNew;
        this->m_ops[a_level]->create(phiNew, a_phiOld);
        this->m_ops[a_level]->setToZero(phiNew);
        if (!a_zeroPhi) {
            this->m_ops[a_level]->assign(phiNew, a_phiOld);
        }
        updateSoln(phiNew, a_phiOld, rhs, a_flux,
                   a_fineFluxRegPtr, a_crseFluxRegPtr,
                   a_crsePhiOldPtr, a_crsePhiNewPtr,
                   a_oldTime, a_crseOldTime,
                   a_crseNewTime, a_dt, a_level, a_zeroPhi);

        //         n+1         n
        // ([a phi]   - [a phi] )
        // ---------------------- - a_src  -> a_diffusiveTerm.
        //          dt
        for (DataIterator dit = a_phiOld.disjointBoxLayout().dataIterator(); dit.ok(); ++dit) {
            aNew[dit()] *= phiNew[dit()];
            aOld[dit()] *= a_phiOld[dit()];
            a_diffusiveTerm[dit()].axby(aNew[dit()], aOld[dit()], 1.0 / a_dt, -1.0 / a_dt);
            a_diffusiveTerm[dit()] -= a_src[dit()];
        }
    }
}


// -----------------------------------------------------------------------------
// Integrates the helmholtz equation represented by this object, placing
// the new solution in a_phiNew.
// -----------------------------------------------------------------------------
void MappedLevelTGA::updateSoln (LevelDataType&           a_phiNew,
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
    // If our operators are time-independent, do the "easy" thing.
    if (!this->m_ops[a_level]->isTimeDependent()) {
        this->updateSolnWithTimeIndependentOp(a_phiNew, a_phiOld, a_src, a_flux,
                                              a_fineFluxRegPtr, a_crseFluxRegPtr,
                                              a_crsePhiOldPtr, a_crsePhiNewPtr,
                                              a_oldTime, a_crseOldTime, a_crseNewTime,
                                              a_dt, a_level, a_zeroPhi,
                                              a_fluxStartComponent);
    } else {
        MayDay::Error("Is this ever called???");
        // We have to assume that the operator and its coefficients are
        // time-dependent.
        this->updateSolnWithTimeDependentOp(a_phiNew, a_phiOld, a_src, a_flux,
                                            a_fineFluxRegPtr, a_crseFluxRegPtr,
                                            a_crsePhiOldPtr, a_crsePhiNewPtr,
                                            a_oldTime, a_crseOldTime, a_crseNewTime,
                                            a_dt, a_level, a_zeroPhi,
                                            a_fluxStartComponent);
    }
}


// -----------------------------------------------------------------------------
// Update the solution assuming that the operator's coefficients are
// independent of time. Same arguments as updateSoln.
// -----------------------------------------------------------------------------
void MappedLevelTGA::updateSolnWithTimeIndependentOp (LevelDataType&           a_phiNew,
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
                                                      int                      a_fluxStartComponent)
{
    CH_assert(!this->m_ops[a_level]->isTimeDependent());
    int ncomp = a_phiNew.nComp();
    Interval intervBase(0, ncomp - 1);
    Interval intervFlux(a_fluxStartComponent, a_fluxStartComponent + ncomp - 1);

    CH_assert(a_level >= 0);
    CH_assert(a_level <  this->m_grids.size());
    CH_assert((a_level == 0) || (a_crsePhiOldPtr != NULL));
    CH_assert((a_level == 0) || (a_crsePhiNewPtr != NULL));
    CH_assert(a_crseNewTime >= a_crseOldTime);
    CH_assert(a_dt >= 0.);

    LevelDataType rhst, srct, phis;
    this->m_ops[a_level]->create(rhst, a_src);
    this->m_ops[a_level]->create(srct, a_phiNew);
    this->m_ops[a_level]->create(phis, a_phiNew);

    this->m_ops[a_level]->setToZero(srct);
    this->m_ops[a_level]->setToZero(rhst);
    //this makes srct = a_src*dt
    this->m_ops[a_level]->incr(srct, a_src, a_dt);

    //save input guess if we are not using zero
    if (!a_zeroPhi) {
        this->m_ops[a_level]->setToZero(phis);
        this->m_ops[a_level]->incr(phis, a_phiNew, 1.0);
    }

    // Divide the source S by the identity coefficient a. srct = rhs*dt/a
    this->m_ops[a_level]->divideByIdentityCoef(srct);

    LevelDataType coarseData;

    if ((a_crsePhiOldPtr != NULL) && (a_level > 0)) {
        this->m_ops[a_level - 1]->create(coarseData, *a_crsePhiOldPtr);
        setSourceGhostCells(srct, this->m_grids[a_level], a_level);
    }

    //from here on k is kappa and L is kappa L
    //this makes rhs hold       (k*a I + mu4 k L) (S/a) dt
    this->applyHelm(rhst, srct, NULL, a_level, m_mu4, a_dt, true); //true means the application is homogeneous
    this->incrementFlux(a_flux, srct, a_level, m_mu4, a_dt, -1.0, true);

    // first need to compute coarse-level BC at this level's old time
    if (a_level > 0) {
        this->timeInterp(coarseData, *a_crsePhiOldPtr, *a_crsePhiNewPtr,
                         a_oldTime, a_crseOldTime, a_crseNewTime, a_level - 1);
    }

    //this makes a_phiNew hold (k*a I + mu3 k L) phi^n
    //'true' apply CF and domain BC
    this->applyHelm(a_phiNew,   a_phiOld, &coarseData, a_level, m_mu3, a_dt, false);
    this->incrementFlux(a_flux, a_phiOld,              a_level, m_mu3, a_dt, -1., false);

    //this makes rhs hold (k*a I + mu3 L) phi^n + dt(k*a I +mu4 L) S/a
    this->m_ops[a_level]->incr(rhst, a_phiNew, 1.0);

    // now construct coarse-level BC for intermediate solve
    // intermediate solution will be at time = oldTime + (1-r1)dt
    if (a_level > 0) {
        Real intermediateTime = a_oldTime + (1.0 - m_r1) * a_dt;

        this->timeInterp(coarseData, *a_crsePhiOldPtr, *a_crsePhiNewPtr,
                         intermediateTime, a_crseOldTime, a_crseNewTime, a_level - 1);
    }

    //if user thought her original guess was good, use it.
    if (!a_zeroPhi) {
        this->m_ops[a_level]->setToZero(a_phiNew);
        this->m_ops[a_level]->incr(a_phiNew, phis,  1.0);
    }

    //by now rhs =  ((k*a I + mu3 L) phi^n + dt(k*a I +mu4 L) S/a )
    //this makes phinew = (k*a I - mu2 k L)^-1 (rhs)
    this->solveHelm(a_phiNew, coarseData, rhst, a_level, m_mu2, a_dt, a_zeroPhi);
    this->incrementFlux(a_flux, a_phiNew,       a_level, m_mu2, a_dt, -1.0, false);

    //this puts the answer into rhst so we can do the final solve
    this->m_ops[a_level]->assign(rhst, a_phiNew);


    // now construct CF-BC for final solve
    if (a_level > 0) {
        Real newTime = a_oldTime + a_dt;
        this->timeInterp(coarseData, *a_crsePhiOldPtr, *a_crsePhiNewPtr,
                         newTime, a_crseOldTime, a_crseNewTime, a_level - 1);
    }

    //by now rhs =  ((k*a I + mu3 L) phi^n + dt(k*a I +mu4 L) S/a )
    //this makes rhs hold k*a[(k*a I - mu2 L)^-1 (rhs)]
    this->m_ops[a_level]->diagonalScale(rhst, true);

    //if user thought her original guess was good, use it again.
    if (!a_zeroPhi) {
        this->m_ops[a_level]->setToZero(a_phiNew);
        this->m_ops[a_level]->incr(a_phiNew, phis,  1.0);
    }

    //this makes phinew = (k*a I - mu1 L)^-1 [ka ((k*a I - mu2 L)^-1 (orig rhs))]
    this->solveHelm(a_phiNew, coarseData, rhst, a_level, m_mu1, a_dt, a_zeroPhi);
    this->incrementFlux(a_flux, a_phiNew, a_level, m_mu1, a_dt, -1.0, false);

    // now increment flux registers -- note that because of the way
    // we defined the fluxes, the dt multiplier is already in the
    // flux
    if ((a_fineFluxRegPtr != NULL) && (a_level < this->m_grids.size() - 1)) {
        const RealVect& dx = this->m_ops[a_level]->getDx();

        for (DataIterator dit = this->m_grids[a_level].dataIterator(); dit.ok(); ++dit) {
            FluxDataType& thisFlux = a_flux[dit];
            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real fluxMult = -a_dt / dx[dir];
                a_fineFluxRegPtr->incrementCoarse(thisFlux[dir],
                                                  fluxMult, dit(),
                                                  intervBase, // source
                                                  intervFlux, // dest
                                                  dir);
            }
        }
    } // end if there is a finer-level

    if ((a_crseFluxRegPtr != NULL) && (a_level > 0)) {
        const RealVect& dx = this->m_ops[a_level-1]->getDx();

        for (DataIterator dit = this->m_grids[a_level].dataIterator(); dit.ok(); ++dit) {
            FluxDataType& thisFlux = a_flux[dit];
            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real fluxMult = -a_dt / dx[dir];
                a_crseFluxRegPtr->incrementFine(thisFlux[dir], fluxMult, dit(),
                                                intervBase, // source
                                                intervFlux, // dest
                                                dir);
            }
        }
    } // end if there is a coarser level
}


// -----------------------------------------------------------------------------
// Update the solution assuming that the operator's coefficients change
// with time. Same arguments as updateSoln.
// -----------------------------------------------------------------------------
void MappedLevelTGA::updateSolnWithTimeDependentOp (LevelDataType&           a_phiNew,
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
                                                    int                      a_fluxStartComponent)
{
    MayDay::Error("Is this ever called???");

    int ncomp = a_phiNew.nComp();
    Interval intervBase(0, ncomp - 1);
    Interval intervFlux(a_fluxStartComponent, a_fluxStartComponent + ncomp - 1);

    CH_assert(a_level >= 0);
    CH_assert(a_level <  this->m_grids.size());
    CH_assert((a_level == 0) || (a_crsePhiOldPtr != NULL));
    CH_assert((a_level == 0) || (a_crsePhiNewPtr != NULL));
    CH_assert(a_crseNewTime >= a_crseOldTime);
    CH_assert(a_dt >= 0.);

    LevelDataType rhst, srct, phis;
    this->m_ops[a_level]->create(rhst, a_src);
    this->m_ops[a_level]->create(srct, a_phiNew);
    this->m_ops[a_level]->create(phis, a_phiNew);

    this->m_ops[a_level]->setToZero(srct);
    this->m_ops[a_level]->setToZero(rhst);
    this->m_ops[a_level]->incr(srct, a_src, 1.0);

    //save input guess if we are not using zero
    if (!a_zeroPhi) {
        this->m_ops[a_level]->setToZero(phis);
        this->m_ops[a_level]->incr(phis, a_phiNew, 1.0);
    }

    // In the comments, we use superscripts to denote different time centerings
    // in the operator L, the identity coefficient a, and the solution phi.
    // A '0' superscript means the quantity at time n.
    // A '1/2' superscript means the quantity at time n + 1/2.
    // A '*' superscript means the quantity at the "intermediate" time used
    //       by the TGA scheme.
    // A '1' superscript means the quantity at time n + 1.
    //
    // In embedded boundary problems we use k for kappa and L means kappa L.

    //-----------------------------------------------------
    //                                         1/2
    //          dt       1/2          1/2     S
    // Compute ----  x (a    I + mu4 L   ) x ---    -> rhst
    //           1/2                           1/2
    //          a                             a
    //-----------------------------------------------------

    // First, feed the half-step time to the operator to set its coefficients.
    this->m_ops[a_level]->setTime(a_oldTime, 0.5, a_dt);

    //                                                  1/2
    // Divide the source S by the identity coefficient a   .
    this->m_ops[a_level]->divideByIdentityCoef(srct);

    // Set the ghost cells in the source at the coarser level.
    LevelDataType coarseData;
    if ((a_crsePhiOldPtr != NULL) && (a_level > 0)) {
        this->m_ops[a_level - 1]->create(coarseData, *a_crsePhiOldPtr);
        setSourceGhostCells(srct, this->m_grids[a_level], a_level);
    }

    //     1/2                    1/2
    // (k a    I + mu4 k L) (S / a   ) -> rhst
    // (The false parameter tells the operator to use the source's ghost data
    //  instead of applying boundary conditions).
    //homogeneous BC for L(rhs)
    this->applyHelm(rhst, srct, NULL, a_level, m_mu4, a_dt, true);
    this->incrementFlux(a_flux, srct, a_level, a_dt * m_mu4, a_dt, -1.0, true);

    //--------------------------------------------------------
    //          1        0           0       0
    // Compute --- x (k a I + mu4 k L ) x phi  and add it to rhst.
    //           0
    //          a
    //--------------------------------------------------------

    // First we compute the coarse-level boundary condition data at time n.
    if (a_level > 0) {
        this->timeInterp(coarseData, *a_crsePhiOldPtr, *a_crsePhiNewPtr,
                         a_oldTime, a_crseOldTime, a_crseNewTime, a_level - 1);
    }

    // Feed the old time to the operator.
    this->m_ops[a_level]->setTime(a_oldTime, 0.0, a_dt);

    //             0         0     0
    // Compute (k a I + mu3 L ) phi and store it in a_phiNew for now.
    // (The true parameter means apply the coarse-fine and domain boundary
    //  conditions).
    //inhommogeneous bc for solution
    this->applyHelm(a_phiNew, a_phiOld, &coarseData, a_level, m_mu3, a_dt, false);
    this->incrementFlux(a_flux, a_phiOld, a_level, m_mu3, a_dt, -1., false);

    //                     n
    // Divide a_phiNew by a .
    this->m_ops[a_level]->divideByIdentityCoef(a_phiNew);

    // Add a_phiNew to rhst to obtain the right hand side of our
    // TGA update equation.
    this->m_ops[a_level]->incr(rhst, a_phiNew, 1.0);

    //-----------------------------------------------------------------
    //           *           *     n    *
    // Solve (k a I + mu4 k L ) phi  = a rhst and place the answer in rhst,
    //-----------------------------------------------------------------

    // Construct the coarse-level boundary condition data for our TGA
    //                                               *   n
    // intermediate solve. The intermediate time is t = t + (1-r1)dt.
    if (a_level > 0) {
        Real intermediateTime = a_oldTime + (1.0 - m_r1) * a_dt;

        this->timeInterp(coarseData, *a_crsePhiOldPtr, *a_crsePhiNewPtr,
                         intermediateTime, a_crseOldTime, a_crseNewTime, a_level - 1);
    }

    // If the user thought her original guess was good, use it.
    if (!a_zeroPhi) {
        this->m_ops[a_level]->setToZero(a_phiNew);
        this->m_ops[a_level]->incr(a_phiNew, phis,  1.0);
    }

    // Feed the intermediate time to the operator to get the time-centered
    // coefficients.
    this->m_ops[a_level]->setTime(a_oldTime, 1.0 - m_r1, a_dt);

    //                   *
    // Multiply rhst by a .
    this->m_ops[a_level]->diagonalScale(rhst, false);

    //                   *         * -1
    // This computes (k a I - mu2 L )   (rhst) and stores it in a_phiNew.
    this->solveHelm(a_phiNew, coarseData, rhst, a_level, m_mu2, a_dt, a_zeroPhi);
    this->incrementFlux(a_flux, a_phiNew, a_level, m_mu2, a_dt, -1.0, false);

    // This puts the answer into rhst.
    this->m_ops[a_level]->assign(rhst, a_phiNew);

    //--------------------------------------------------
    //           1           1     1                      1
    // Solve (k a I + k mu4 L ) phi  = rhst and scale by a .
    //--------------------------------------------------

    // Compute the coarse-fine boundary condition data for the final solve.
    if (a_level > 0) {
        Real newTime = a_oldTime + a_dt;
        this->timeInterp(coarseData, *a_crsePhiOldPtr, *a_crsePhiNewPtr,
                         newTime, a_crseOldTime, a_crseNewTime, a_level - 1);
    }

    // Feed the new time to the operator.
    this->m_ops[a_level]->setTime(a_oldTime, 1.0, a_dt);

    //                  1
    // Scale rhst by k a .
    this->m_ops[a_level]->diagonalScale(rhst);

    // If the user thought her original guess was good, use it again.
    if (!a_zeroPhi) {
        this->m_ops[a_level]->setToZero(a_phiNew);
        this->m_ops[a_level]->incr(a_phiNew, phis, 1.0);
    }

    //             1           1  -1    1     *           * -1  *
    // Compute [k a I - mu1 k L )]   k a [(k a I - mu2 k L )   a (orig rhs)].
    this->solveHelm(a_phiNew, coarseData, rhst, a_level, m_mu1, a_dt, a_zeroPhi);
    this->incrementFlux(a_flux, a_phiNew, a_level, m_mu1, a_dt, -1.0, false);

    // Now increment the flux registers -- note that because of the way
    // we defined the fluxes, the dt multiplier is already in the
    // flux.
    if ((a_fineFluxRegPtr != NULL) && (a_level < this->m_grids.size() - 1)) {
        const RealVect& dx = this->m_ops[a_level]->getDx();

        for (DataIterator dit = this->m_grids[a_level].dataIterator(); dit.ok(); ++dit) {
            FluxDataType& thisFlux = a_flux[dit];
            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real fluxMult = -a_dt / dx[dir];
                a_fineFluxRegPtr->incrementCoarse(thisFlux[dir],
                                                  fluxMult, dit(),
                                                  intervBase, // source
                                                  intervFlux, // dest
                                                  dir);
            }
        }
    } // end if there is a finer-level

    if ((a_crseFluxRegPtr != NULL) && (a_level > 0)) {
        const RealVect& dx = this->m_ops[a_level-1]->getDx();

        for (DataIterator dit = this->m_grids[a_level].dataIterator(); dit.ok(); ++dit) {
            FluxDataType& thisFlux = a_flux[dit];
            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real fluxMult = -a_dt / dx[dir];
                a_crseFluxRegPtr->incrementFine(thisFlux[dir], fluxMult, dit(),
                                                intervBase, // source
                                                intervFlux, // dest
                                                dir);
            }
        }
    } // end if there is a coarser level
}
