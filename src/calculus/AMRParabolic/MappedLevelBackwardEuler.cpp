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
#include "MappedLevelBackwardEuler.H"
#include "ExtrapolationUtilsF_F.H"
#include "FluxBox.H"


// -----------------------------------------------------------------------------
// Full constructor
// -----------------------------------------------------------------------------
MappedLevelBackwardEuler::MappedLevelBackwardEuler (const Vector<DisjointBoxLayout>&                                 a_grids,
                                                    const Vector<IntVect>&                                           a_refRat,
                                                    const ProblemDomain&                                             a_level0Domain,
                                                    RefCountedPtr<MappedAMRLevelOpFactory<LevelData<FArrayBox> > >&  a_opFact,
                                                    const RefCountedPtr<MappedAMRMultiGrid<LevelData<FArrayBox> > >& a_solver)
: MappedBaseLevelHeatSolver(a_grids, a_refRat, a_level0Domain, a_opFact, a_solver)
{;}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
MappedLevelBackwardEuler::~MappedLevelBackwardEuler ()
{;}


// -----------------------------------------------------------------------------
// Integrates the helmholtz equation represented by this object, placing
// the new solution in a_phiNew.
// -----------------------------------------------------------------------------
void MappedLevelBackwardEuler::updateSoln (LevelDataType&           a_phiNew,
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

    LevelDataType rhst, phit;
    this->m_ops[a_level]->create(rhst, a_src);
    this->m_ops[a_level]->create(phit, a_phiNew);

    this->m_ops[a_level]->setToZero(phit);
    this->m_ops[a_level]->setToZero(rhst);
    if (a_zeroPhi) {
        this->m_ops[a_level]->setToZero(a_phiNew);
    }

    //set phit to a_phiOld, rhst to a_src * a_dt
    this->m_ops[a_level]->incr(phit, a_phiOld, 1.0);
    this->m_ops[a_level]->incr(rhst, a_src   , a_dt);

    //multiply phi old by kappa*acoef
    this->m_ops[a_level]->diagonalScale(phit, true);

    //multiply rhs by kappa (but NOT by a)
    if (!a_rhsAlreadyKappaWeighted)
        this->m_ops[a_level]->kappaScale(rhst);

    //add both together to make rhs for euler solve = kappa a phiold  + kappa rhs dt
    this->m_ops[a_level]->incr(rhst, phit, 1.0);

    //solve for phi new = (kappa a I - kappa dt L) phinew = kappa a phiold  + kappa rhs
    //this makes phinew = (k*a I - mu2 L)^-1 (rhs)
    LevelDataType coarseData;
    if ((a_crsePhiOldPtr != NULL) && (a_level > 0)) {
        this->m_ops[a_level - 1]->create(coarseData, *a_crsePhiOldPtr);
        this->m_ops[a_level - 1]->setToZero(coarseData);

        this->timeInterp(coarseData, *a_crsePhiOldPtr, *a_crsePhiNewPtr,
                         a_oldTime, a_crseOldTime, a_crseNewTime, a_level - 1);
    }

    this->solveHelm(a_phiNew, coarseData, rhst, a_level, 1.0, a_dt, a_zeroPhi);
    this->incrementFlux(a_flux, a_phiNew,       a_level, 1.0, a_dt, -1.0, true);

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
