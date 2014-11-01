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
#include "MappedAMRPoissonOp.H"
#include "MappedAMRPoissonOpF_F.H"
#include "MappedAMRPoissonOpOrthoF_F.H"
#include "DivCurlGradF_F.H"
#include "MappedCoarseAverage.H"
#include "MappedCoarseAverageF_F.H"
#include "HomogeneousCFInterp.H"
#include "computeMappedSum.H"
#include "AnisotropicRefinementTools.H"
#include "LevelGeometry.H"
#include "ExtrapolationUtils.H"
#include "Constants.H"
#include "MiscUtils.H"
#include "Debug.H"
#include "ProblemContext.H"
#include "Printing.H"
#include "CornerCopier.H"
#include "AMRIO.H"
#include "InterpF_F.H" // Only needed for 2 InterpHomo functions.


#define HEADER() \
//    if(procID() == 0) {cout << m_levid << " [" << __LINE__ << "]: " << __PRETTY_FUNCTION__ << endl;}

static const bool considerCellSizes = false;

int MappedAMRPoissonOp::s_maxCoarse = 4;    // The smallest allowable grid size


// -----------------------------------------------------------------------------
// Creates a simple LinearOp for a LinearSolver.
// If you have a coarser level, this will only handle homogeneous CFBCs.
// Take care of your CFBCs, then solve the residual equation.
// If a_JinvPtr.isNull(), this op will just use Jinv=1.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::linearOpDefine (BCMethodHolder&                             a_bc,
                                         const RealVect&                             a_dx,
                                         const RealVect&                             a_dxCrse,
                                         const Real                                  a_alpha,
                                         const Real                                  a_beta,
                                         const bool                                  a_isDiagonal,
                                         const bool                                  a_horizontalOp,
                                         const RefCountedPtr<LevelData<FluxBox> >&   a_JgupPtr,
                                         const RefCountedPtr<LevelData<FArrayBox> >& a_JinvPtr,
                                         const int                                   a_precondMode)
{
    CH_assert(a_JgupPtr != NULL);
    CH_assert(a_JgupPtr->isDefined());

    const DisjointBoxLayout& grids = a_JgupPtr->getBoxes();
    DataIterator dit = grids.dataIterator();
    const ProblemDomain& domain = grids.physDomain();
    Copier exCopier(grids, grids, domain, IntVect::Unit, true);
    CFRegion cfregion(grids, domain);

    this->define(grids, a_dx, domain, a_bc, exCopier, cfregion);

    m_preCondSmoothIters = 0;
    m_FCJgup = a_JgupPtr;
    if (a_JinvPtr.isNull()) {
        m_CCJinv = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>);
        m_CCJinv->define(grids, 1);
        m_levelOps.setVal(*m_CCJinv, 1.0);
    } else {
        m_CCJinv = (RefCountedPtr<LevelData<FArrayBox> >&)a_JinvPtr;
    }

    m_isDiagonal = a_isDiagonal;
    m_horizontalOp = a_horizontalOp;
    m_time = quietNAN;

    m_precondMode = a_precondMode;
    m_precondRelaxPtr = NULL;
    m_relaxPtr = NULL;
    m_restrictPtr = NULL;
    m_prolongPtr = NULL;

    m_dxCrse = a_dxCrse;

    m_aCoef = a_alpha;
    m_bCoef = a_beta;
    m_alpha = a_alpha;
    m_beta = a_beta;

    // Cache the diagonal
    m_lapDiag = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>);
    m_lapDiag->define(grids, 1);

    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& lapDiagFAB = (*m_lapDiag)[dit];
        const FluxBox& JgupFB = (*m_FCJgup)[dit];
        const FArrayBox& JinvFAB = (*m_CCJinv)[dit];

#       if CH_SPACEDIM == 2
            if (!m_horizontalOp) {
                FORT_FILLMAPPEDLAPDIAG2D(CHF_FRA1(lapDiagFAB,0),
                                         CHF_CONST_FRA(JgupFB[0]),
                                         CHF_CONST_FRA(JgupFB[1]),
                                         CHF_CONST_FRA1(JinvFAB,0),
                                         CHF_BOX(lapDiagFAB.box()),
                                         CHF_CONST_REALVECT(m_dx));
            } else {
                FORT_FILLMAPPEDLAPDIAG1D(CHF_FRA1(lapDiagFAB,0),
                                         CHF_CONST_FRA1(JgupFB[0],0),
                                         CHF_CONST_FRA1(JinvFAB,0),
                                         CHF_BOX(lapDiagFAB.box()),
                                         CHF_CONST_REALVECT(m_dx));
            }
#       elif CH_SPACEDIM == 3
            if (!m_horizontalOp) {
                FORT_FILLMAPPEDLAPDIAG3D(CHF_FRA1(lapDiagFAB,0),
                                         CHF_CONST_FRA(JgupFB[0]),
                                         CHF_CONST_FRA(JgupFB[1]),
                                         CHF_CONST_FRA(JgupFB[2]),
                                         CHF_CONST_FRA1(JinvFAB,0),
                                         CHF_BOX(lapDiagFAB.box()),
                                         CHF_CONST_REALVECT(m_dx));
            } else {
                FORT_FILLMAPPEDLAPDIAG2D(CHF_FRA1(lapDiagFAB,0),
                                         CHF_CONST_FRA(JgupFB[0]),
                                         CHF_CONST_FRA(JgupFB[1]),
                                         CHF_CONST_FRA1(JinvFAB,0),
                                         CHF_BOX(lapDiagFAB.box()),
                                         CHF_CONST_REALVECT(m_dx));
            }
#       else
#           error Bad CH_SPACEDIM
#       endif
    }
}


// -----------------------------------------------------------------------------
// Full define function for an operator with both finer and coarser levels.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::define(const DisjointBoxLayout& a_grids,
                                const DisjointBoxLayout& a_gridsFiner,
                                const DisjointBoxLayout& a_gridsCoarser,
                                const RealVect&          a_dxLevel,
                                const IntVect&           a_refRatio,
                                const IntVect&           a_refRatioFiner,
                                const ProblemDomain&     a_domain,
                                BCMethodHolder&          a_bc,
                                const Copier&            a_exchange,
                                const CFRegion&          a_cfregion)

{
    CH_TIME("MappedAMRPoissonOp::define1");

    this->define(a_grids, a_gridsCoarser, a_dxLevel, a_refRatio, a_domain, a_bc,
                 a_exchange, a_cfregion);
    m_refToFiner = a_refRatioFiner;

    if (a_gridsFiner.isClosed()) {
        ProblemDomain fineDomain = refine(m_domain, m_refToFiner);
        m_levfluxreg.define(a_gridsFiner,
                            a_grids,
                            fineDomain,
                            m_refToFiner,
                            1 /* ncomp*/);
    }
}


// -----------------------------------------------------------------------------
// Full define function for an operator with finer levels, but no coarser.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::define(const DisjointBoxLayout& a_grids,
                                const DisjointBoxLayout& a_gridsFiner,
                                const RealVect&          a_dxLevel,
                                const IntVect&           a_refRatio, // dummy arg
                                const IntVect&           a_refRatioFiner,
                                const ProblemDomain&     a_domain,
                                BCMethodHolder&          a_bc,
                                const Copier&            a_exchange,
                                const CFRegion&          a_cfregion)
{
    CH_TIME("MappedAMRPoissonOp::define2");

    CH_assert(a_refRatio == IntVect::Unit);

    // Calls the MG version of define
    this->define(a_grids, a_dxLevel, a_domain, a_bc, a_exchange, a_cfregion);
    m_refToFiner = a_refRatioFiner;

    ProblemDomain fineDomain = refine(m_domain, m_refToFiner);
    m_levfluxreg.define(a_gridsFiner, a_grids, fineDomain, m_refToFiner,
                        1); // ncomp
}


// -----------------------------------------------------------------------------
// Full define function for an operator with coarser levels, but no finer.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::define(const DisjointBoxLayout& a_grids,
                                const DisjointBoxLayout& a_coarse,
                                const RealVect&          a_dxLevel,
                                const IntVect&           a_refRatio,
                                const ProblemDomain&     a_domain,
                                BCMethodHolder&          a_bc,
                                const Copier&            a_exchange,
                                const CFRegion&          a_cfregion)
{
    CH_TIME("MappedAMRPoissonOp::define3");

    this->define(a_grids, a_dxLevel, a_domain, a_bc,
                 a_exchange, a_cfregion);

    m_refToCoarser = a_refRatio;
    m_dxCrse = RealVect(a_refRatio) * a_dxLevel;
    m_refToFiner = IntVect::Unit;

    m_interpWithCoarser.define(a_grids, &a_coarse, a_dxLevel,
                               m_refToCoarser, 1, m_domain);
}


// -----------------------------------------------------------------------------
// Full define function for an operator with no coarser or finer levels.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::define(const DisjointBoxLayout& a_grids,
                                const RealVect&          a_dx,
                                const ProblemDomain&     a_domain,
                                BCMethodHolder&          a_bc,
                                const Copier&            a_exchange,
                                const CFRegion&          a_cfregion)
{
    CH_TIME("MappedAMRPoissonOp::define4");

    m_bc     = a_bc;
    m_domain = a_domain;
    m_dx     = a_dx;
    m_dxCrse = RealVect(D_DECL(quietNAN, quietNAN, quietNAN));

    // Redefined in AMRLevelOp<LevelData<FArrayBox> >::define virtual function.
    m_refToCoarser = IntVect::Unit;
    m_refToFiner   = IntVect::Unit;

    // These get set again after define is called
    m_alpha = 0.0;
    m_beta  = 1.0;
    m_aCoef  = quietNAN;
    m_bCoef  = quietNAN;

    m_exchangeCopier = a_exchange;
    m_cfregion = a_cfregion;

    m_validDomain = a_domain.domainBox();
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (a_domain.isPeriodic(dir)) {
            m_validDomain.grow(dir, 1);
        }
    }
}


// -----------------------------------------------------------------------------
// Full define function for an operator with no coarser or finer levels.
// This version is not used by the factory.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::define(const DisjointBoxLayout& a_grids,
                                const RealVect&          a_dx,
                                const ProblemDomain&     a_domain,
                                BCMethodHolder&          a_bc)
{
    CH_TIME("MappedAMRPoissonOp::define5");
    MayDay::Error("MappedAMRPoissonOp::define5 function may be outdated");

    Copier exchangeCopier;
    exchangeCopier.define(a_grids, a_grids, IntVect::Unit, true);
    exchangeCopier.trimEdges(a_grids, IntVect::Unit);

    CFRegion cfregion(a_grids, a_domain);

    this->define(a_grids, a_dx, a_domain, a_bc, exchangeCopier, cfregion);
}


// -----------------------------------------------------------------------------
// Full define function that mimics the old PoissonOp.
// Deprecated, throws an error.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::define(const DisjointBoxLayout& a_grids,
                                const DisjointBoxLayout* a_baseBAPtr,
                                const RealVect&          a_dx,
                                const IntVect&           a_refRatio,
                                const ProblemDomain&     a_domain,
                                BCMethodHolder&          a_bc)
{
    TODO();
    UNDEFINED_FUNCTION();
    return;
}



// -----------------------------------------------------------------------------
// Defines this op to simply calculate Laplacians and residuals.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::define(BCMethodHolder&                             a_bc,
                                const LevelGeometry&                        a_levGeo,
                                const Real                                  a_alpha,
                                const Real                                  a_beta,
                                const RefCountedPtr<LevelData<FArrayBox> >& a_lapDiagsPtr)
{
    CH_TIME("MappedAMRPoissonOp::define7");

    // Get this level's data from levGeo
    const ProblemDomain& domain = a_levGeo.getDomain();
    const DisjointBoxLayout& grids = a_levGeo.getBoxes();
    const RealVect& dx = a_levGeo.getDx();

    // Define default coarser data
    const DisjointBoxLayout* crseGridsPtr = NULL;
    RealVect crseDx(D_DECL(-1.0, -1.0, -1.0));
    IntVect crseRefRatio = IntVect::Unit;

    // Get coarser level's data from levGeo
    const LevelGeometry* crseLevGeoPtr = a_levGeo.getCoarserPtr();
    if (crseLevGeoPtr != NULL) {
        crseGridsPtr = &(crseLevGeoPtr->getBoxes());
        crseRefRatio = crseLevGeoPtr->getFineRefRatio();
        crseDx = crseLevGeoPtr->getDx();

        CH_assert(crseRefRatio == a_levGeo.getCrseRefRatio());
        CH_assert(crseDx == dx * RealVect(crseRefRatio));
    }

    // Define default finer data
    const DisjointBoxLayout* fineGridsPtr = NULL;
    RealVect fineDx(D_DECL(-1.0, -1.0, -1.0));
    IntVect fineRefRatio = IntVect::Unit;

    // Get finer level's data from levGeo
    const LevelGeometry* fineLevGeoPtr = a_levGeo.getFinerPtr();
    if (fineLevGeoPtr != NULL) {
        fineGridsPtr = &(fineLevGeoPtr->getBoxes());
        fineRefRatio = fineLevGeoPtr->getCrseRefRatio();
        fineDx = fineLevGeoPtr->getDx();

        CH_assert(fineRefRatio == a_levGeo.getFineRefRatio());
        CH_assert(fineDx * RealVect(fineRefRatio) == dx);
    }

    // Define Copier and CFRegion
    Copier exchangeCopier;
    exchangeCopier.define(grids, grids, domain, IntVect::Unit, true);
    if (LevelGeometry::isDiagonal()) {
        exchangeCopier.trimEdges(grids, IntVect::Unit);
    }

    CFRegion cfregion;
    cfregion.define(grids, domain);

    // Call the appropriate define function
    if (crseGridsPtr != NULL) {
        if (fineGridsPtr != NULL) {
            // Intermediate level
            this->define(grids, *fineGridsPtr, *crseGridsPtr, dx, crseRefRatio, fineRefRatio,
                         domain, a_bc, exchangeCopier, cfregion);
        } else {
            // Top level
            this->define(grids, *crseGridsPtr, dx, crseRefRatio,
                         domain, a_bc, exchangeCopier, cfregion);
        }
    } else {
        if (fineGridsPtr != NULL) {
            if (fineGridsPtr->size() > 0) {
                // Bottom level with a non-empty finer level
                const IntVect dummyRef = IntVect::Unit;
                this->define(grids, *fineGridsPtr, dx, dummyRef, fineRefRatio,
                             domain, a_bc, exchangeCopier, cfregion);
            } else {
                // Bottom level with an empty finer level (ie, only level)
                // This is needed to define the op when the finer levGeos
                // have not yet been regridded.
                this->define(grids, dx,
                             domain, a_bc, exchangeCopier, cfregion);
            }
        } else {
            // Only level
            this->define(grids, dx,
                         domain, a_bc, exchangeCopier, cfregion);
        }
    }

    // Define members that AMRLevelOp does not know about.
    m_alpha  = a_alpha;
    m_beta   = a_beta;
    m_aCoef  = a_alpha;
    m_bCoef  = a_beta;
    m_dxCrse = crseDx;

    // Just set to some default value.
    // The preconditioner should not be used by this op anyway.
    m_preCondSmoothIters = 2;

    // Grab pointers to the metric
    m_FCJgup  = a_levGeo.getFCJgupPtr();
    m_CCJinv  = a_levGeo.getCCJinvPtr();

    // Calculate the operator diagonals
    if (a_lapDiagsPtr.isNull()) {
        m_lapDiag = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>);
        m_lapDiag->define(grids, 1, IntVect::Zero);

        DataIterator dit = m_lapDiag->dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& lapDiagFAB = (*m_lapDiag)[dit];
            const FluxBox& JgupFB = (*m_FCJgup)[dit];
            const FArrayBox& JinvFAB = (*m_CCJinv)[dit];

#           if CH_SPACEDIM == 2
                FORT_FILLMAPPEDLAPDIAG2D(CHF_FRA1(lapDiagFAB,0),
                                         CHF_CONST_FRA(JgupFB[0]),
                                         CHF_CONST_FRA(JgupFB[1]),
                                         CHF_CONST_FRA1(JinvFAB,0),
                                         CHF_BOX(lapDiagFAB.box()),
                                         CHF_CONST_REALVECT(dx));
#           elif CH_SPACEDIM == 3
                FORT_FILLMAPPEDLAPDIAG3D(CHF_FRA1(lapDiagFAB,0),
                                         CHF_CONST_FRA(JgupFB[0]),
                                         CHF_CONST_FRA(JgupFB[1]),
                                         CHF_CONST_FRA(JgupFB[2]),
                                         CHF_CONST_FRA1(JinvFAB,0),
                                         CHF_BOX(lapDiagFAB.box()),
                                         CHF_CONST_REALVECT(dx));
#           else
#               error Bad CH_SPACEDIM
#           endif
        }
    } else {
        // We were given a prepared set of diagonals.
        m_lapDiag = a_lapDiagsPtr;
    }

    m_validDomain = domain.domainBox();
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (domain.isPeriodic(dir)) {
            m_validDomain.grow(dir, 1);
        }
    }

    // TEMPORARY!!!
    m_levid = std::string("(-, 0)");

    CH_assert(m_dx == dx);
    CH_assert(m_dxCrse == crseDx);
    CH_assert(m_refToCoarser == crseRefRatio);
    CH_assert(m_refToFiner == fineRefRatio);
}


// -----------------------------------------------------------------------------
// Installs a new metric
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::setJgup (RefCountedPtr<LevelData<FluxBox> >& a_JgupPtr)
{
    const DisjointBoxLayout& grids = a_JgupPtr->getBoxes();
    CH_assert(grids.compatible(m_FCJgup->getBoxes()));
    CH_assert(grids.compatible(m_CCJinv->getBoxes()));
    CH_assert(grids.compatible(m_lapDiag->getBoxes()));

    // Collect the user's homegrown Jgup ptr.
    m_FCJgup = a_JgupPtr;

    // Calculate the operator's new diagonals
    m_lapDiag = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>);
    m_lapDiag->define(grids, 1, IntVect::Zero);

    DataIterator dit = m_lapDiag->dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& lapDiagFAB = (*m_lapDiag)[dit];
        const FluxBox& JgupFB = (*m_FCJgup)[dit];
        const FArrayBox& JinvFAB = (*m_CCJinv)[dit];

#       if CH_SPACEDIM == 2
            FORT_FILLMAPPEDLAPDIAG2D(CHF_FRA1(lapDiagFAB,0),
                                     CHF_CONST_FRA(JgupFB[0]),
                                     CHF_CONST_FRA(JgupFB[1]),
                                     CHF_CONST_FRA1(JinvFAB,0),
                                     CHF_BOX(lapDiagFAB.box()),
                                     CHF_CONST_REALVECT(m_dx));
#       elif CH_SPACEDIM == 3
            FORT_FILLMAPPEDLAPDIAG3D(CHF_FRA1(lapDiagFAB,0),
                                     CHF_CONST_FRA(JgupFB[0]),
                                     CHF_CONST_FRA(JgupFB[1]),
                                     CHF_CONST_FRA(JgupFB[2]),
                                     CHF_CONST_FRA1(JinvFAB,0),
                                     CHF_BOX(lapDiagFAB.box()),
                                     CHF_CONST_REALVECT(m_dx));
#       else
#           error Bad CH_SPACEDIM
#       endif
    }
}


// -----------------------------------------------------------------------------
// Installs a new metric
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::setJgup (const FillJgupInterface* a_customFillJgupPtr)
{
    const DisjointBoxLayout& grids = m_FCJgup->getBoxes();
    DataIterator dit = grids.dataIterator();

    // It's important that we create a new pointer. If we just overwrite the
    // data at the old address, we'll disturb everyone using that pointer!
    RefCountedPtr<LevelData<FluxBox> > newFCJgup(new LevelData<FluxBox>);
    newFCJgup->define(grids, SpaceDim);

    // Fill new holder.
    for (dit.reset(); dit.ok(); ++dit) {
        for (int mu = 0; mu < SpaceDim; ++mu) {
            FArrayBox& JgupFAB = (*newFCJgup)[dit][mu];

            for (int nu = 0; nu < SpaceDim; ++nu) {
                a_customFillJgupPtr->fill_Jgup(JgupFAB, // FAB
                                               nu,      // FAB comp
                                               mu,      // first index
                                               nu,      // second index
                                               m_dx,    // dx
                                               1.0);    // scale
            }
        }
    }

    // Update diagonals and such.
    this->setJgup(newFCJgup);
}


// -----------------------------------------------------------------------------
// Sets the time used to evaluate BCs.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::setTime(Real a_time) {
    m_time = a_time;
}


// -----------------------------------------------------------------------------
// Returns the time used to evaluate BCs.
// -----------------------------------------------------------------------------
Real MappedAMRPoissonOp::getTime() const {
    return m_time;
}


// -----------------------------------------------------------------------------
// Access function. The aCoef and bCoef stuff is for the TGA solvers.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::setAlphaAndBeta (const Real& a_alpha,
                                          const Real& a_beta)
{
    m_alpha = a_alpha * m_aCoef;
    m_beta  = a_beta  * m_bCoef;

    CH_assert(m_relaxPtr != NULL);
    m_relaxPtr->setAlphaAndBeta(m_alpha, m_beta);
    if (m_precondRelaxPtr != m_relaxPtr) {
        m_precondRelaxPtr->setAlphaAndBeta(m_alpha, m_beta);
    }

    DataIterator dit = m_lapDiag->dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& lapDiagFAB = (*m_lapDiag)[dit];
        const FluxBox& JgupFB = (*m_FCJgup)[dit];
        const FArrayBox& JinvFAB = (*m_CCJinv)[dit];

#       if CH_SPACEDIM == 2
            FORT_FILLMAPPEDLAPDIAG2D(CHF_FRA1(lapDiagFAB,0),
                                     CHF_CONST_FRA(JgupFB[0]),
                                     CHF_CONST_FRA(JgupFB[1]),
                                     CHF_CONST_FRA1(JinvFAB,0),
                                     CHF_BOX(lapDiagFAB.box()),
                                     CHF_CONST_REALVECT(m_dx));
#       elif CH_SPACEDIM == 3
            FORT_FILLMAPPEDLAPDIAG3D(CHF_FRA1(lapDiagFAB,0),
                                     CHF_CONST_FRA(JgupFB[0]),
                                     CHF_CONST_FRA(JgupFB[1]),
                                     CHF_CONST_FRA(JgupFB[2]),
                                     CHF_CONST_FRA1(JinvFAB,0),
                                     CHF_BOX(lapDiagFAB.box()),
                                     CHF_CONST_REALVECT(m_dx));
#       else
#           error Bad CH_SPACEDIM
#       endif
    }
}


// LinearOp functions ----------------------------------------------------------

// -----------------------------------------------------------------------------
// Calculates a_lhs = a_rhs - L[a_phi] and interpolates ghosts from a
// coarser grid of zeros.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::residual(LevelData<FArrayBox>&       a_lhs,
                                  const LevelData<FArrayBox>& a_phi,
                                  const LevelData<FArrayBox>& a_rhs,
                                  bool                        a_homogeneous)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::residual");

    LevelData<FArrayBox>& phiRef = (LevelData<FArrayBox>&)a_phi;
    this->interpCFGhosts(phiRef, NULL, true);

    this->residualI(a_lhs, a_phi, a_rhs, a_homogeneous);
}


// -----------------------------------------------------------------------------
// Calculate a_lhs = a_rhs - L[a_phi] and ignore the coarse-fine boundary
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::residualI(LevelData<FArrayBox>&       a_lhs,
                                   const LevelData<FArrayBox>& a_phi,
                                   const LevelData<FArrayBox>& a_rhs,
                                   bool                        a_homogeneous)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::residualI");
    CH_assert(&a_lhs != &a_rhs); // Does not work in place!

    const DisjointBoxLayout& grids = a_lhs.disjointBoxLayout();
    DataIterator dit = a_lhs.dataIterator();

    // Sanity check
    DBLCHECK(m_FCJgup->getBoxes(), grids);
    // CH_assert(m_FCJgup->getBoxes().compatible(grids));

    // Calculate L[phi]
    this->applyOpI(a_lhs, a_phi, a_homogeneous);

    // Loop through grids and calculate rhs-L[phi].
    for (dit.begin(); dit.ok(); ++dit) {
        FArrayBox& lhsFAB = a_lhs[dit];
        const FArrayBox& phiFAB = a_phi[dit];
        const FArrayBox& rhsFAB = a_rhs[dit];
        const Box& region = grids[dit];

        // Subtract from rhs
        FORT_SUBTRACTOP(CHF_FRA(lhsFAB),
                        CHF_CONST_FRA(rhsFAB),
                        CHF_CONST_FRA(lhsFAB),
                        CHF_BOX(region));
    }
}


// -----------------------------------------------------------------------------
// this preconditioner first initializes phihat to (DiagOfOp)phihat = rhshat
// then smooths with a couple of passes of levelGSRB
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::preCond (LevelData<FArrayBox>&       a_phi,
                                  const LevelData<FArrayBox>& a_rhs)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::preCond");

    CH_assert(a_phi.nComp() == a_rhs.nComp());
    const DisjointBoxLayout& grids = a_rhs.getBoxes();
    DataIterator dit = a_phi.dataIterator();

    if (m_preCondSmoothIters == 0) {
        // Do nothing
        a_rhs.copyTo(a_phi);

    } else {
        switch (m_precondMode) {
        case ProblemContext::PrecondMode::DiagRelax:
        case ProblemContext::PrecondMode::DiagLineRelax:
            // Divide by operator's diagonals
            for (dit.begin(); dit.ok(); ++dit)
            {
                FArrayBox& phiFAB = a_phi[dit];
                const FArrayBox& rhsFAB = a_rhs[dit];
                const FArrayBox& lapDiagFAB = (*m_lapDiag)[dit];
                const Box& region = grids[dit];

                FORT_DIAGPRECOND(CHF_FRA(phiFAB),
                                 CHF_CONST_FRA(rhsFAB),
                                 CHF_CONST_FRA1(lapDiagFAB,0),
                                 CHF_BOX(region),
                                 CHF_CONST_REAL(m_alpha),
                                 CHF_CONST_REAL(m_beta));
            }

            // Smooth the result with lineGSRB
            CH_assert(m_precondRelaxPtr != NULL);
            for (int iter = 0; iter < m_preCondSmoothIters; ++iter) {
                m_precondRelaxPtr->relax(a_phi, a_rhs);
            }
            break;

        case ProblemContext::PrecondMode::None:
            // Do nothing
            a_rhs.copyTo(a_phi);
            break;

        default:
            MayDay::Error("Bad m_precondMode");
        };
    }
}


// -----------------------------------------------------------------------------
// Fill in all ghosts then apply the operator
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::applyOp(LevelData<FArrayBox>&       a_lhs,
                                 const LevelData<FArrayBox>& a_phi,
                                 bool                        a_homogeneous)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::applyOp");

#ifndef NDEBUG
    if (!(a_lhs.getBoxes().physDomain() == m_FCJgup->getBoxes().physDomain()) ||
        !(a_phi.getBoxes().physDomain() == m_FCJgup->getBoxes().physDomain())) {

        pout() << "this = (MappedAMRPoissonOp*)" << (void*)this
               << "\nthis->m_levid = " << m_levid
               << "\nthis->m_FCJgup->getBoxes().physDomain() = " << m_FCJgup->getBoxes().physDomain()
               << "\na_lhs.getBoxes().physDomain() = " << a_lhs.getBoxes().physDomain()
               << "\na_phi.getBoxes().physDomain() = " << a_lhs.getBoxes().physDomain()
               << endl;

        MayDay::Error("Operator defined on a different level than the fields given");
    }
#endif
    LevelData<FArrayBox>& phiRef = (LevelData<FArrayBox>&)a_phi;
    this->interpCFGhosts(phiRef, NULL, true);

    this->applyOpI(a_lhs, a_phi, a_homogeneous);
}


// -----------------------------------------------------------------------------
// Apply the operator on this level only - ignore the bdry to the coarser level.
// The ghosts at the CF boundary need to be filled before calling this function.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::applyOpI (LevelData<FArrayBox>&       a_lhs,
                                   const LevelData<FArrayBox>& a_phi,
                                   bool                        a_homogeneous)
{
    HEADER();
    CH_TIMERS("MappedAMRPoissonOp::applyOpI");
    CH_TIMER("MappedAMRPoissonOp::applyOpI::exchange", exchangeTmr);
    CH_TIMER("MappedAMRPoissonOp::applyOpI::extrap", extrapTmr);
    CH_TIMER("MappedAMRPoissonOp::applyOpI::grad", gradTmr);
    CH_TIMER("MappedAMRPoissonOp::applyOpI::div", divTmr);

    const DisjointBoxLayout& grids = a_lhs.disjointBoxLayout();
    DataIterator dit = a_phi.dataIterator();
    const int ncomps = a_phi.nComp();

    // Sanity check
    DBLCHECK(m_FCJgup->getBoxes(), grids);
    // CH_assert(m_FCJgup->getBoxes().compatible(grids));

    // This is OK if we use it to only change ghost cells
    LevelData<FArrayBox>& phiRef = const_cast< LevelData<FArrayBox>& >(a_phi);

    // Begin exchange of valid data
    CH_START(exchangeTmr);
    this->exchangeComplete(phiRef);
    CH_STOP(exchangeTmr);

    // Loop through grids and perform calculations that do not need ghosts.
    for (dit.reset(); dit.ok(); ++dit) {
        const FArrayBox& phiFAB    = a_phi[dit];
        FArrayBox&       phiRefFAB = phiRef[dit];
        FArrayBox&       lhsFAB    = a_lhs[dit];
        const FluxBox&   JgupFB    = (*m_FCJgup)[dit];
        const FArrayBox& JinvFAB   = (*m_CCJinv)[dit];
        const Box&       valid     = grids[dit];
        const Box&       phiBox    = phiFAB.box();

        // Create temp holders
        FluxBox fluxFB(valid, ncomps);
#ifndef NDEBUG
        fluxFB.setVal(quietNAN);
#endif

        // Extrapolate ghosts for non-diagonal derivatives
        FArrayBox extrapFAB(phiFAB.box(), ncomps);
        CH_START(extrapTmr);
        this->fillExtrap(extrapFAB, phiFAB, 2);
        CH_STOP(extrapTmr);

        // Set physical BCs
        m_bc.setGhosts(phiRefFAB, &extrapFAB, valid, m_domain, m_dx, dit(), &JgupFB, a_homogeneous, m_time);

        // Compute fluxes
        CH_START(gradTmr);
        if (!m_horizontalOp) {
            D_TERM(this->getFluxComplete(fluxFB[0], phiFAB, extrapFAB, surroundingNodes(valid,0), dit(), 0);,
                   this->getFluxComplete(fluxFB[1], phiFAB, extrapFAB, surroundingNodes(valid,1), dit(), 1);,
                   this->getFluxComplete(fluxFB[2], phiFAB, extrapFAB, surroundingNodes(valid,2), dit(), 2);)
        } else {
            D_TERM(,
                   this->getFluxComplete(fluxFB[0], phiFAB, extrapFAB, surroundingNodes(valid,0), dit(), 0);,
                   this->getFluxComplete(fluxFB[1], phiFAB, extrapFAB, surroundingNodes(valid,1), dit(), 1);)
        }
        CH_STOP(gradTmr);

        // Set boundary fluxes
        m_bc.setFluxes(fluxFB, &extrapFAB, valid, m_domain, m_dx, dit(), &JgupFB, a_homogeneous, m_time);

        // Compute divergence
        CH_START(divTmr);
#if CH_SPACEDIM == 2
        if (!m_horizontalOp) {
            FORT_MAPPEDFLUXDIVERGENCE2D(
                CHF_FRA(lhsFAB),
                CHF_CONST_FRA(fluxFB[0]),
                CHF_CONST_FRA(fluxFB[1]),
                CHF_CONST_FRA1(JinvFAB,0),
                CHF_BOX(valid),
                CHF_CONST_REALVECT(m_dx));
        } else {
            FORT_MAPPEDFLUXDIVERGENCE1D(
                CHF_FRA(lhsFAB),
                CHF_CONST_FRA(fluxFB[0]),
                CHF_CONST_FRA1(JinvFAB,0),
                CHF_BOX(valid),
                CHF_CONST_REAL(m_dx[0]));
        }
#elif CH_SPACEDIM == 3
        if (!m_horizontalOp) {
            FORT_MAPPEDFLUXDIVERGENCE3D(
                CHF_FRA(lhsFAB),
                CHF_CONST_FRA(fluxFB[0]),
                CHF_CONST_FRA(fluxFB[1]),
                CHF_CONST_FRA(fluxFB[2]),
                CHF_CONST_FRA1(JinvFAB,0),
                CHF_BOX(valid),
                CHF_CONST_REALVECT(m_dx));
        } else {
            FORT_MAPPEDFLUXDIVERGENCE2D(
                CHF_FRA(lhsFAB),
                CHF_CONST_FRA(fluxFB[0]),
                CHF_CONST_FRA(fluxFB[1]),
                CHF_CONST_FRA1(JinvFAB,0),
                CHF_BOX(valid),
                CHF_CONST_REALVECT(m_dx));
        }
#else
#error Bad CH_SPACEDIM
#endif
        CH_STOP(divTmr);

        // Right now, lhsFAB = beta*L[phiFAB].
        // We want lhsFAB = (alpha + beta*L)[phiFAB].
        // This utility function does just that for us.
        if (m_alpha != 0.0) {
            const Real betaDummy = 1.0;
            FORT_AXBYIP(
                CHF_FRA(lhsFAB),
                CHF_CONST_FRA(phiFAB),
                CHF_CONST_REAL(m_alpha),
                CHF_CONST_REAL(betaDummy),
                CHF_BOX(valid));
        }
    }
}


// -----------------------------------------------------------------------------
// Same as applyOp, but assumes ALL of phi's ghosts are filled.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::applyOpNoBoundary(LevelData<FArrayBox>&       a_lhs,
                                           const LevelData<FArrayBox>& a_phi,
                                           const bool                  a_homogeneous)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::applyOpNoBoundary");

    MayDay::Error("MappedAMRPoissonOp::applyOpNoBoundary should be updated");

    const DisjointBoxLayout& grids = a_lhs.disjointBoxLayout();
    DataIterator dit = a_phi.dataIterator();
    const int ncomps = a_phi.nComp();
    const int extrapOrder = 2;

    // Set up dimensional info
    IntVect ghostVect(D_DECL(1,1,1));
    int dims = SpaceDim;
    if (m_horizontalOp) {
        ghostVect[CH_SPACEDIM-1] = 0;
        --dims;
    }

    // Sanity check
    DBLCHECK(m_FCJgup->getBoxes(), grids);
    // CH_assert(m_FCJgup->getBoxes().compatible(grids));

    // Allocate scratch space.
    LevelData<FArrayBox> extrap(grids, ncomps, ghostVect);
    LevelData<FluxBox> flux(grids, ncomps, IntVect::Zero);

    // Loop through grids and perform calculations that do not need ghosts.
    for (dit.reset(); dit.ok(); ++dit) {
        FluxBox&         fluxFB    = flux[dit];
        const FArrayBox& phiFAB    = a_phi[dit];
        FArrayBox&       extrapFAB = extrap[dit];
        const FluxBox&   JgupFB    = (*m_FCJgup)[dit];
        Box              valid     = grids[dit];

        // Extrapolate ghosts for non-diagonal derivatives
        if (!m_isDiagonal) {
            extrapFAB.copy(phiFAB);

            for (int fdir = 0; fdir < dims; ++fdir) {
                SideIterator fsit;
                for (fsit.reset(); fsit.ok(); ++fsit) {
                    ExtrapolateFaceAndCopy(extrapFAB, extrapFAB, valid,
                                           fdir, fsit(), extrapOrder);
                }
                valid.grow(fdir, 1);
            }
            valid = grids[dit];
        }

        // Compute fluxes
#ifndef NDEBUG
        fluxFB.setVal(quietNAN);
#endif

        if (!m_horizontalOp) {
            D_TERM(this->getFluxDuringExchange(fluxFB[0], phiFAB, extrapFAB, surroundingNodes(valid,0), dit(), 0);,
                   this->getFluxDuringExchange(fluxFB[1], phiFAB, extrapFAB, surroundingNodes(valid,1), dit(), 1);,
                   this->getFluxDuringExchange(fluxFB[2], phiFAB, extrapFAB, surroundingNodes(valid,2), dit(), 2);)
        } else {
            D_TERM(,
                   this->getFluxDuringExchange(fluxFB[0], phiFAB, extrapFAB, surroundingNodes(valid,0), dit(), 0);,
                   this->getFluxDuringExchange(fluxFB[1], phiFAB, extrapFAB, surroundingNodes(valid,1), dit(), 1);)
        }
    }

    // Loop through grids and perform remaining calculations.
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox&       lhsFAB    = a_lhs[dit];
        FluxBox&         fluxFB    = flux[dit];
        const FArrayBox& phiFAB    = a_phi[dit];
        const FArrayBox& extrapFAB = extrap[dit];
        const FluxBox&   JgupFB    = (*m_FCJgup)[dit];
        const FArrayBox& JinvFAB   = (*m_CCJinv)[dit];
        const Box&       valid     = grids[dit];

        // Compute remaining fluxes
        if (!m_horizontalOp) {
            D_TERM(this->getFluxAfterExchange(fluxFB[0], phiFAB, extrapFAB, surroundingNodes(valid,0), dit(), 0);,
                   this->getFluxAfterExchange(fluxFB[1], phiFAB, extrapFAB, surroundingNodes(valid,1), dit(), 1);,
                   this->getFluxAfterExchange(fluxFB[2], phiFAB, extrapFAB, surroundingNodes(valid,2), dit(), 2);)
        } else {
            D_TERM(,
                   this->getFluxAfterExchange(fluxFB[0], phiFAB, extrapFAB, surroundingNodes(valid,0), dit(), 0);,
                   this->getFluxAfterExchange(fluxFB[1], phiFAB, extrapFAB, surroundingNodes(valid,1), dit(), 1);)
        }

        // Set boundary fluxes
        m_bc.setFluxes(fluxFB, &extrapFAB, valid, m_domain, m_dx, dit(), &JgupFB, a_homogeneous, m_time);

        // Compute divergence
#if CH_SPACEDIM == 2
        if (!m_horizontalOp) {
            FORT_MAPPEDFLUXDIVERGENCE2D(
                CHF_FRA(lhsFAB),
                CHF_CONST_FRA(fluxFB[0]),
                CHF_CONST_FRA(fluxFB[1]),
                CHF_CONST_FRA1(JinvFAB,0),
                CHF_BOX(valid),
                CHF_CONST_REALVECT(m_dx));
        } else {
            FORT_MAPPEDFLUXDIVERGENCE1D(
                CHF_FRA(lhsFAB),
                CHF_CONST_FRA(fluxFB[0]),
                CHF_CONST_FRA1(JinvFAB,0),
                CHF_BOX(valid),
                CHF_CONST_REAL(m_dx[0]));
        }
#elif CH_SPACEDIM == 3
        if (!m_horizontalOp) {
            FORT_MAPPEDFLUXDIVERGENCE3D(
                CHF_FRA(lhsFAB),
                CHF_CONST_FRA(fluxFB[0]),
                CHF_CONST_FRA(fluxFB[1]),
                CHF_CONST_FRA(fluxFB[2]),
                CHF_CONST_FRA1(JinvFAB,0),
                CHF_BOX(valid),
                CHF_CONST_REALVECT(m_dx));
        } else {
            FORT_MAPPEDFLUXDIVERGENCE2D(
                CHF_FRA(lhsFAB),
                CHF_CONST_FRA(fluxFB[0]),
                CHF_CONST_FRA(fluxFB[1]),
                CHF_CONST_FRA1(JinvFAB,0),
                CHF_BOX(valid),
                CHF_CONST_REALVECT(m_dx));
        }
#else
#error Bad CH_SPACEDIM
#endif

        // Right now, lhsFAB = beta*L[phiFAB].
        // We want lhsFAB = (alpha + beta*L)[phiFAB].
        // This utility function does just that for us.
        if (m_alpha != 0.0) {
            const Real betaDummy = 1.0;
            FORT_AXBYIP(CHF_FRA(lhsFAB),
                        CHF_CONST_FRA(phiFAB),
                        CHF_CONST_REAL(m_alpha),
                        CHF_CONST_REAL(betaDummy),
                        CHF_BOX(valid));
        }
    }
}


// -----------------------------------------------------------------------------
// Define a_lhs to be a coarsened version of a_rhs.
// This function does not fill a_lhs with data.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::createCoarsened(LevelData<FArrayBox>&       a_lhs,
                                         const LevelData<FArrayBox>& a_rhs,
                                         const IntVect&              a_refRat)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::createCoarsened");

    DisjointBoxLayout crseGrids;
    const DisjointBoxLayout& fineGrids = a_rhs.disjointBoxLayout();
    CH_assert(coarsenable(fineGrids, a_refRat));

    if (a_refRat == IntVect(D_DECL(2,2,2))) {
        if (m_coarsenedMGrids.size() == 0)
            coarsen(m_coarsenedMGrids, fineGrids, 2);
        crseGrids = m_coarsenedMGrids;
    } else {
        coarsen(crseGrids, fineGrids, a_refRat);
    }

    const int ncomp = a_rhs.nComp();
    const IntVect ghostVect = a_rhs.ghostVect();
    CH_assert( !(m_horizontalOp && ghostVect[CH_SPACEDIM-1] > 0) );
    a_lhs.define(crseGrids, ncomp, ghostVect);
}


// -----------------------------------------------------------------------------
// Just calls this->norm
// -----------------------------------------------------------------------------
Real MappedAMRPoissonOp::localMaxNorm(const LevelData<FArrayBox>& a_x)
{
    CH_TIME("MappedAMRPoissonOp::localMaxNorm");

    Real localMax = 0;
    int nComp = a_x.nComp();
    DataIterator dit = a_x.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        localMax = Max(localMax, a_x[dit].norm(a_x.box(dit()), 0, 0, nComp));
    }
    return localMax;
}


// MGLevelOp functions ---------------------------------------------------------

// MappedAMRPoissonOp::relax moved to be with its friends (levelJacobi, etc).

// -----------------------------------------------------------------------------
// MG function to define (but not fill) the coarser level data holder.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::createCoarser(LevelData<FArrayBox>&       a_coarse,
                                       const LevelData<FArrayBox>& a_fine,
                                       bool                        a_ghosted)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::createCoarser");

    CH_assert(m_domain == a_fine.getBoxes().physDomain());
    CH_assert(coarsenable(a_fine.getBoxes(), m_mgCrseRefRatio));

    if (m_coarsenedMGrids.size() == 0) {
        coarsen(m_coarsenedMGrids, a_fine.disjointBoxLayout(), m_mgCrseRefRatio);
    }

    const IntVect& ghost = a_fine.ghostVect();
    CH_assert( !(m_horizontalOp && ghost[CH_SPACEDIM-1] > 0) );
    a_coarse.define(m_coarsenedMGrids, a_fine.nComp(), ghost);
}


// -----------------------------------------------------------------------------
// Computes restriction of residual to a coarser multgrid level
// res[2h] = I[h->2h] (rhs[h] - L(phi[h]))
// This op is on the finer grid.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                          LevelData<FArrayBox>&       a_phiFine,
                                          const LevelData<FArrayBox>& a_rhsFine)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::restrictResidual");

    // Sanity checks
    CH_assert(a_resCoarse.nComp() == a_phiFine.nComp());
    CH_assert(a_resCoarse.nComp() == a_rhsFine.nComp());
    CH_assert(m_domain == a_phiFine.getBoxes().physDomain());
    CH_assert(m_domain == a_rhsFine.getBoxes().physDomain());
    DBLCHECK(m_FCJgup->getBoxes(), a_rhsFine.getBoxes());

    // Compute the fine-level residual.
    const DisjointBoxLayout& fineGrids = a_rhsFine.getBoxes();
    const int ncomps = a_phiFine.nComp();
    LevelData<FArrayBox> resFine(fineGrids, ncomps, IntVect::Zero);
    this->residual(resFine, a_phiFine, a_rhsFine, true);

    // Restrict the residual to the coarser grids.
    CH_assert(m_restrictPtr != NULL);
    m_restrictPtr->restrict(a_resCoarse, resFine);
}


// AMRLevelOp overrides --------------------------------------------------------

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::AMRResidual(LevelData<FArrayBox>&                    a_residual,
                                     const LevelData<FArrayBox>&              a_phiFine,
                                     const LevelData<FArrayBox>&              a_phi,
                                     const LevelData<FArrayBox>&              a_phiCoarse,
                                     const LevelData<FArrayBox>&              a_rhs,
                                     const bool                               a_homogeneousPhysBC,
                                     MappedAMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
    CH_TIME("MappedAMRPoissonOp::AMRResidual");

    this->AMROperator(a_residual, a_phiFine, a_phi, a_phiCoarse,
                      a_homogeneousPhysBC, a_finerOp);

    this->axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
}


// -----------------------------------------------------------------------------
// Calculate the residual assuming there are no coarser AMR levels.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::AMRResidualNC(LevelData<FArrayBox>&                    a_residual,
                                       const LevelData<FArrayBox>&              a_phiFine,
                                       const LevelData<FArrayBox>&              a_phi,
                                       const LevelData<FArrayBox>&              a_rhs,
                                       const bool                               a_homogeneousPhysBC,
                                       MappedAMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::AMRResidualNC");

    this->AMROperatorNC(a_residual, a_phiFine, a_phi,
                        a_homogeneousPhysBC, a_finerOp);

    this->axby(a_residual, a_residual, a_rhs, -1.0, 1.0);
}


// -----------------------------------------------------------------------------
// Calculate the residual assuming there are no finer AMR levels.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::AMRResidualNF(LevelData<FArrayBox>&       a_residual,
                                       const LevelData<FArrayBox>& a_phi,
                                       const LevelData<FArrayBox>& a_phiCoarse,
                                       const LevelData<FArrayBox>& a_rhs,
                                       const bool                  a_homogeneousPhysBC)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::AMRResidualNF");

    // Fill ghosts on CF interface using quadratic interpolation
    if (a_phiCoarse.isDefined()) {
        LevelData<FArrayBox>& phiRef = const_cast<LevelData<FArrayBox>&>(a_phi);
        this->interpCFGhosts(phiRef, &a_phiCoarse, false);
    }

    // Apply the physical BCs and calculate the residual on this level.
    this->residualI(a_residual, a_phi, a_rhs, a_homogeneousPhysBC);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::AMROperator(LevelData<FArrayBox>&                    a_LofPhi,
                                     const LevelData<FArrayBox>&              a_phiFine,
                                     const LevelData<FArrayBox>&              a_phi,
                                     const LevelData<FArrayBox>&              a_phiCoarse,
                                     const bool                               a_homogeneousPhysBC,
                                     MappedAMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::AMROperator");

    CH_assert(a_phi.isDefined());

    // Fill ghosts on CF interface using quadratic interpolation
    if (a_phiCoarse.isDefined()) {
        LevelData<FArrayBox>& phiRef = const_cast<LevelData<FArrayBox>&>(a_phi);
        this->interpCFGhosts(phiRef, &a_phiCoarse, false);
    }

    // Set BCs and apply the operator on this level.
    this->applyOpI(a_LofPhi, a_phi, a_homogeneousPhysBC);

    // Reflux
    if (a_phiFine.isDefined()) {
        CH_assert(a_finerOp != NULL);
        this->reflux(a_phiFine, a_phi, a_LofPhi, a_finerOp);
    }
}


// -----------------------------------------------------------------------------
// Apply the AMR op, including CF matching, and assume no coarser AMR levels.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::AMROperatorNC(LevelData<FArrayBox>&                    a_LofPhi,
                                       const LevelData<FArrayBox>&              a_phiFine,
                                       const LevelData<FArrayBox>&              a_phi,
                                       const bool                               a_homogeneousPhysBC,
                                       MappedAMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::AMROperatorNC");

    CH_assert(a_phi.isDefined());

    // Set BCs and apply the operator on this level.
    this->applyOpI(a_LofPhi, a_phi, a_homogeneousPhysBC);

    // Reflux
    if (a_phiFine.isDefined()) {
        CH_assert(a_finerOp != NULL);
        this->reflux(a_phiFine, a_phi, a_LofPhi, a_finerOp);
    }
}


// -----------------------------------------------------------------------------
// Apply the AMR operator, including CF matching. Assume no finer AMR level.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::AMROperatorNF(LevelData<FArrayBox>&       a_LofPhi,
                                       const LevelData<FArrayBox>& a_phi,
                                       const LevelData<FArrayBox>& a_phiCoarse,
                                       const bool                  a_homogeneousPhysBC)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::AMROperatorNF");

    CH_assert(a_phi.isDefined());

    if (a_phiCoarse.isDefined()) {
        LevelData<FArrayBox>& phiRef = (LevelData<FArrayBox>&)a_phi;
        this->interpCFGhosts(phiRef, &a_phiCoarse, false);
    }

    // apply physical boundary conditions in applyOpI
    this->applyOpI(a_LofPhi, a_phi, a_homogeneousPhysBC);
}


// -----------------------------------------------------------------------------
// Calculate the residual and average the result to the next coarser AMR level.
//   a_resCoarse = I[h->2h]( a_residual - L(a_correction, a_coarseCorrection))
// It is assumed that a_resCoarse has already been filled in with the coarse
// version of AMRResidualNF and that this operation is free to overwrite in the
// overlap regions. This function operates on the finer level.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::AMRRestrict(LevelData<FArrayBox>&       a_resCoarse,
                                     const LevelData<FArrayBox>& a_residual,
                                     const LevelData<FArrayBox>& a_correction,
                                     const LevelData<FArrayBox>& a_coarseCorrection)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::AMRRestrict");

    // Create scratch space
    LevelData<FArrayBox> r;
    this->create(r, a_residual);

    // Do the restriction
    this->AMRRestrictS(a_resCoarse, a_residual, a_correction, a_coarseCorrection, r);
}


// -----------------------------------------------------------------------------
// Calculate the residual and average the result to the next coarser AMR level.
//   a_resCoarse = I[h->2h]( a_residual - L(a_correction, a_coarseCorrection))
// It is assumed that a_resCoarse has already been filled in with the coarse
// version of AMRResidualNF and that this operation is free to overwrite in the
// overlap regions. This function operates on the finer level. The 'S' means
// that scratch space is passed in so that extra allocation is unnecessary.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::AMRRestrictS(LevelData<FArrayBox>&       a_resCoarse,
                                      const LevelData<FArrayBox>& a_residual,
                                      const LevelData<FArrayBox>& a_correction,
                                      const LevelData<FArrayBox>& a_coarseCorrection,
                                      LevelData<FArrayBox>&       a_scratch)
{
    HEADER();
    CH_TIME("MappedAMRPoissonOp::AMRRestrictS");

    // Calculate the residual as if there is no finer AMR level.
    this->AMRResidualNF(a_scratch, a_correction,
                        a_coarseCorrection, a_residual,
                        true);

    // Restrict the residual to the coarser grids.
    CH_assert(m_restrictPtr != NULL);
    m_restrictPtr->restrict(a_resCoarse, a_scratch);
}


// -----------------------------------------------------------------------------
// a_correction += I[2h->h](a_coarseCorrection)
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::AMRProlong(LevelData<FArrayBox>&       a_correction,
                                    const LevelData<FArrayBox>& a_coarseCorrection)
{
    CH_TIME("MappedAMRPoissonOp::AMRProlong");

    // Unlike the MG vesion, we cannot be sure that the coarse and fine grids
    // are compatible. So, we create a version of the coarse correction that
    // is compatible with the fine grids, then use the MG prolongation function.

    DisjointBoxLayout coarsenedGrids;
    coarsen(coarsenedGrids, a_correction.disjointBoxLayout(), m_refToCoarser);

    LevelData<FArrayBox> eCoar(coarsenedGrids, a_correction.nComp(), a_coarseCorrection.ghostVect());
    a_coarseCorrection.copyTo(eCoar.interval(), eCoar, eCoar.interval());

    this->prolongIncrement(a_correction, eCoar);
}


// -----------------------------------------------------------------------------
// Performs a correction of fine level data from coarse data.
// The 'S' means that scratch space (a_temp) is provided so extra allocation is
// unnecessary. a_temp must be defined with coarse level grids prior to call.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::AMRProlongS(LevelData<FArrayBox>&       a_correction,
                                     const LevelData<FArrayBox>& a_coarseCorrection,
                                     LevelData<FArrayBox>&       a_temp,
                                     const Copier&               a_copier)
{
    CH_TIME("MappedAMRPoissonOp::AMRProlong");

    // Unlike the MG version, we cannot be sure that the coarse and fine grids
    // are compatible. So, we create a version of the coarse correction that
    // is compatible with the fine grids, then use the MG prolongation function.

    DisjointBoxLayout coarsenedGrids;
    coarsen(coarsenedGrids, a_correction.disjointBoxLayout(), m_refToCoarser);

    LevelData<FArrayBox> eCoar(coarsenedGrids, a_correction.nComp(), a_coarseCorrection.ghostVect());
    a_coarseCorrection.copyTo(eCoar.interval(), eCoar, eCoar.interval());

    this->prolongIncrement(a_correction, eCoar);
}


// -----------------------------------------------------------------------------
// This function just calls AMRResidualNF
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::AMRUpdateResidual(LevelData<FArrayBox>&       a_residual,
                                           const LevelData<FArrayBox>& a_correction,
                                           const LevelData<FArrayBox>& a_coarseCorrection)
{
    LevelData<FArrayBox> oldRes;
    this->create(oldRes, a_residual);
    this->assign(oldRes, a_residual);
    this->AMRResidualNF(a_residual,
                        a_correction,
                        a_coarseCorrection,
                        oldRes,
                        true);
}


// -----------------------------------------------------------------------------
// Compute norm over all cells on coarse not covered by finer
// -----------------------------------------------------------------------------
Real MappedAMRPoissonOp::AMRNorm(const LevelData<FArrayBox>& a_coarResid,
                                 const LevelData<FArrayBox>& a_fineResid,
                                 const IntVect&              a_refRat,
                                 const int                   a_ord)

{
    CH_TIME("MappedAMRPoissonOp::AMRNorm");

    //create temp and zero out under finer grids
    LevelData<FArrayBox> coarTemp;
    m_levelOps.create(coarTemp, a_coarResid);
    m_levelOps.assign(coarTemp, a_coarResid);

    if (a_fineResid.isDefined()) {
        const DisjointBoxLayout& coarGrids = a_coarResid.disjointBoxLayout();
        const DisjointBoxLayout& fineGrids = a_fineResid.disjointBoxLayout();

        int ncomp = coarTemp.nComp();

        for (DataIterator dit = coarGrids.dataIterator(); dit.ok(); ++dit) {
            FArrayBox& coarTempFAB = coarTemp[dit];
            LayoutIterator litFine = fineGrids.layoutIterator();

            for (litFine.reset(); litFine.ok(); ++litFine) {
                Box overlayBox = coarTempFAB.box();
                Box coarsenedGrid = coarsen(fineGrids[litFine], a_refRat);

                overlayBox &= coarsenedGrid;

                if (!overlayBox.isEmpty()) {
                    coarTempFAB.setVal(0.0, overlayBox, 0, ncomp);
                }
            }
        }
    }

    // return norm of temp
    return norm(coarTemp, a_ord);
}


// -----------------------------------------------------------------------------
// Perform refluxing between this and the next finer level.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::reflux(const LevelData<FArrayBox>&              a_phiFine,
                                const LevelData<FArrayBox>&              a_phi,
                                LevelData<FArrayBox>&                    a_residual,
                                MappedAMRLevelOp<LevelData<FArrayBox> >* a_finerOp)
{
    HEADER();
    CH_TIMERS("MappedAMRPoissonOp::reflux");
    CH_TIMER("MappedAMRPoissonOp::reflux::incrementCoarse", t2);
    CH_TIMER("MappedAMRPoissonOp::reflux::incrementFine", t3);

    // Collect data and perform sanity checks
    MappedAMRPoissonOp* fineOp = dynamic_cast<MappedAMRPoissonOp*>(a_finerOp);
    CH_assert(fineOp != NULL);

    const DisjointBoxLayout& crseGrids = a_phi.getBoxes();
    CH_assert(this->m_FCJgup->getBoxes() == crseGrids);

    const DisjointBoxLayout& fineGrids = a_phiFine.getBoxes();
    CH_assert(fineOp->m_FCJgup->getBoxes() == fineGrids);

    const Interval interv(0, a_phi.nComp()-1);
    const int ncomps = a_phi.nComp();
    CH_assert(a_phiFine.nComp() == ncomps);

    // Set up dimensional info
    int dims = SpaceDim;
    if (m_horizontalOp) --dims;

    // 1. Initialize the level flux register to zero...
    m_levfluxreg.setToZero();


    // 2. Increment the coarse side...
    CH_START(t2);
    for (DataIterator ditc(crseGrids); ditc.ok(); ++ditc) {
        const FArrayBox& crsePhi = a_phi[ditc];
        const Box& crseRegion = crseGrids[ditc];

        if (m_levfluxreg.hasCF(ditc())) {
            for (int idir = 0; idir < dims; ++idir) {
                // Calculate the fluxes
                FArrayBox crseFlux;
                this->getFlux(crseFlux, crsePhi, ditc(), idir);

                // Increment the register
                const Real scale = 1.0 / m_dx[idir];
                m_levfluxreg.incrementCoarse(crseFlux, scale, ditc(),
                                             interv, interv, idir);
            }
        }
    }
    CH_STOP(t2);

    // 3. Increment the fine side...

    // Fill the fine level ghosts at the CF interface by
    // quadratic interpolation from the coarser level
    {
        LevelData<FArrayBox>& phiFineRef = const_cast<LevelData<FArrayBox>&>(a_phiFine);
        fineOp->interpCFGhosts(phiFineRef, &a_phi, false);
    }

    CH_START(t3);
    for (DataIterator ditf(fineGrids); ditf.ok(); ++ditf) {
        const FArrayBox& finePhi = a_phiFine[ditf];
        const Box& fineRegion = fineGrids[ditf];

        for (int idir = 0; idir < dims; ++idir) {
            SideIterator sit;
            for (sit.begin(); sit.ok(); sit.next()) {
                if (m_levfluxreg.hasCF(ditf(), sit())) {
                    const Side::LoHiSide hiorlo = sit();
                    const Box fluxBox = bdryBox(fineRegion, idir, hiorlo, 1);

                    // Calculate the fluxes
                    FArrayBox fineFlux(fluxBox, ncomps);
                    fineOp->getFlux(fineFlux, finePhi, fluxBox,
                                    ditf(), idir,
                                    1); // refToFiner (1 since we are using the finer op)

                    // Increment the register
                    const Real scale = 1.0 / m_dx[idir];
                    m_levfluxreg.incrementFine(fineFlux, scale, ditf(),
                                               interv, interv, idir, hiorlo);
                }
            }
        }
    }
    CH_STOP(t3);

    // 4. Do the refluxing...
    CH_assert(m_CCJinv->getBoxes() == a_residual.getBoxes());
    m_levfluxreg.reflux(a_residual, *m_CCJinv);
}


// -----------------------------------------------------------------------------
// Debugging aid for solvers.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::write(const LevelData<FArrayBox>* a_data,
                               const char*                 a_filename)
{
    CH_TIME("MappedAMRPoissonOp::write");

#   ifdef CH_USE_HDF5
        writeLevelname(a_data, a_filename);
#   else
        MayDay::Warning("MappedAMRPoissonOp::write unimplemented since CH_USE_HDF5 undefined");
#   endif
}


// -----------------------------------------------------------------------------
// Use your relaxation operator to remove the high frequency wave numbers from
// the correction so that it may be averaged to a coarser refinement.
// A point relaxtion scheme, for example takes the form
// a_correction -= lambda*(L(a_correction) - a_residual).
// This function is part of the MGLevelOp interface.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::relax(LevelData<FArrayBox>&       a_e,
                               const LevelData<FArrayBox>& a_residual,
                               int                         a_iterations)
{
    CH_assert(m_relaxPtr != NULL);
    CH_assert(m_relaxPtr->isDefined());

    for (int i = 0; i < a_iterations; ++i) {
        m_relaxPtr->relax(a_e, a_residual);
    }
}


// -----------------------------------------------------------------------------
// Computes fluxes of a_data at faces between valid cells of a_data.
// Only does one direction. a_flux gets redefined.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::getFlux(FArrayBox&       a_flux,
                                 const FArrayBox& a_data,
                                 const DataIndex& a_di,
                                 int              a_dir,
                                 int              a_ref) const
{
    CH_TIME("MappedAMRPoissonOp::getFlux1");

#ifndef NDEBUG
    // Set up dimensional info
    int dims = SpaceDim;
    if (m_horizontalOp) --dims;

    // Sanity checks
    CH_assert(a_dir >= 0);
    CH_assert(a_dir < dims);
    CH_assert(!a_data.box().isEmpty());
#endif

    // Calculate destination flux box
    const DisjointBoxLayout& grids = m_FCJgup->getBoxes();
    const Box valid = grow(a_data.box(), 1) & grids[a_di];
    const Box edgebox = surroundingNodes(valid, a_dir);

    // We want to compute the fluxes at faces between valid cells of a_data.
    // If the assert fails, the data box was too small (one cell wide, in fact).
    CH_assert(!edgebox.isEmpty());

    // Redefine the output data holder.
    a_flux.resize(edgebox, a_data.nComp());

    // Caculate the fluxes over the entire box
    this->getFlux(a_flux, a_data, edgebox, a_di, a_dir, a_ref);
}


// -----------------------------------------------------------------------------
// Computes fluxes of a_data at faces between valid cells of a_data.
// Only does one direction, the edgebox is passed in, and the flux array
// must already be defined.
//
// This can be used as a gradient funciton since it calculates
//   J * Grad[a_data]^{a_dir} over the FC a_edgebox.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::getFlux (FArrayBox&       a_flux,
                                  const FArrayBox& a_data,
                                  const Box&       a_edgebox,
                                  const DataIndex& a_di,
                                  int              a_dir,
                                  int              a_ref) const
{
    CH_TIME("MappedAMRPoissonOp::getFlux2");

    // Set up dimensional info
    int dims = SpaceDim;
    if (m_horizontalOp) --dims;

    // Sanity checks
    CH_assert(a_dir >= 0);
    CH_assert(a_dir < dims);
    CH_assert(!a_data.box().isEmpty());
    CH_assert(!a_edgebox.isEmpty());
    CH_assert( a_flux.contains(a_edgebox) );

    // Create references for convenience
    const DisjointBoxLayout& grids = m_FCJgup->getBoxes();
    const FArrayBox& JgupFAB = (*m_FCJgup)[a_di][a_dir];
    const int ncomps = a_data.nComp();

#ifndef NDEBUG
    {
        // Initialize fluxes
        a_flux.setVal(quietNAN);

        // If this fails, we won't have enough Jgup data
        const Box validFCBox = surroundingNodes(grids[a_di], a_dir);
        CH_assert( validFCBox.contains(a_edgebox) );
    }
#endif

    // Extrapolate ghosts for non-diagonal derivatives
    // This is equivalent to using one-sided derivatives
    FArrayBox extrap(a_data.box(), ncomps);
    this->fillExtrap(extrap, a_data, 2);

    // Compute all fluxes
    this->getFluxComplete(a_flux, a_data, extrap, a_edgebox, a_di, a_dir, a_ref);
}


// -----------------------------------------------------------------------------
// Computes fluxes of a_data at faces between valid cells of a_data.
// Only does one direction, the edgebox is passed in, and the flux array
// must already be defined.
//
// This can be used as a gradient funciton since it calculates
//   J * Grad[a_data]^{a_dir} over the FC a_edgebox.
//
// One sided differencing in the off-diagonal derivatives is done via
// extrapolation of a_data's ghosts (stored in extrap) and then using a
// standard gradient stencil. extrap must have it's ghosts extrapolated and
// a_data must have it's BCs set prior to call.
//
// This version only updates cells that do not require ghosts.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::getFluxDuringExchange (FArrayBox&       a_flux,
                                                const FArrayBox& a_data,
                                                const FArrayBox& a_extrap,
                                                const Box&       a_edgebox,
                                                const DataIndex& a_di,
                                                int              a_dir,
                                                int              a_ref) const
{
    CH_TIME("MappedAMRPoissonOp::getFluxDuringExchange");
    TODO();
    UNDEFINED_FUNCTION();

#ifndef NDEBUG
    // Set up dimensional info
    int dims = SpaceDim;
    if (m_horizontalOp) --dims;

    // Sanity checks
    CH_assert(a_dir >= 0);
    CH_assert(a_dir < dims);
    CH_assert(!a_data.box().isEmpty());
    CH_assert(!a_edgebox.isEmpty());
    CH_assert( a_flux.contains(a_edgebox) );
    CH_assert(m_FCJgup != NULL);
#endif

    // Create references for convenience
    const DisjointBoxLayout& grids = m_FCJgup->getBoxes();
    const Box& validCCBox = grids[a_di];
    const FArrayBox& JgupFAB = (*m_FCJgup)[a_di][a_dir];

    // Compute fluxes over interior
    Box FCBox = grow(validCCBox, -BASISV(a_dir));               // TODO: Swap order to prevent empty boxes.
    FCBox.surroundingNodes(a_dir);
    FCBox &= a_edgebox;

    if (!FCBox.isEmpty()) {
        if (m_isDiagonal) {
            // Don't need to do anything special for ortho + horizontalOp.
            const Real scale = m_beta * a_ref / m_dx[a_dir];
            if (!m_horizontalOp) {
                FORT_MAPPEDGETFLUXORTHO(
                    CHF_FRA(a_flux),
                    CHF_CONST_FRA(a_data),
                    CHF_CONST_FRA1(JgupFAB,a_dir),
                    CHF_BOX(FCBox),
                    CHF_CONST_REAL(scale),
                    CHF_CONST_INT(a_dir));
            } else {
                MayDay::Error("FORT_MAPPEDGETFLATFLUXORTHO is not yet written");
            }
        } else {
            const Real scale = m_beta * a_ref;
            if (!m_horizontalOp) {
                FORT_MAPPEDGETFLUX(
                    CHF_FRA(a_flux),
                    CHF_CONST_FRA(a_data),
                    CHF_CONST_FRA(a_extrap),
                    CHF_CONST_FRA(JgupFAB),
                    CHF_BOX(FCBox),
                    CHF_CONST_REAL(scale),
                    CHF_CONST_REALVECT(m_dx),
                    CHF_CONST_INT(a_dir));
            } else {
#if CH_SPACEDIM == 3
                FORT_MAPPEDGETFLATFLUX3D(
                    CHF_FRA(a_flux),
                    CHF_CONST_FRA(a_data),
                    CHF_CONST_FRA(a_extrap),
                    CHF_CONST_FRA(JgupFAB),
                    CHF_BOX(FCBox),
                    CHF_CONST_REAL(scale),
                    CHF_CONST_REALVECT(m_dx),
                    CHF_CONST_INT(a_dir));
#else
                MayDay::Error("getFluxDuringExchange - is this needed?");
                CH_assert(a_dir == 0);
                FORT_MAPPEDGETFLATFLUX2D(
                    CHF_FRA(a_flux),
                    CHF_CONST_FRA(a_data),
                    CHF_CONST_FRA1(JgupFAB,0),
                    CHF_BOX(FCBox),
                    CHF_CONST_REAL(scale),
                    CHF_CONST_REAL(m_dx[0]));
#endif
            }
        } // end if diagonal or not
    } // end if FCBox is not empty
}


// -----------------------------------------------------------------------------
// Computes fluxes of a_data at faces between valid cells of a_data.
// Only does one direction, the edgebox is passed in, and the flux array
// must already be defined.
//
// This can be used as a gradient funciton since it calculates
//   J * Grad[a_data]^{a_dir} over the FC a_edgebox.
//
// One sided differencing in the off-diagonal derivatives is done via
// extrapolation of a_data's ghosts (stored in extrap) and then using a
// standard gradient stencil. extrap must have it's ghosts extrapolated and
// a_data must have it's BCs set prior to call.
//
// This version only updates cells that need ghosts -- finish exchange first!
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::getFluxAfterExchange(FArrayBox&       a_flux,
                                              const FArrayBox& a_data,
                                              const FArrayBox& a_extrap,
                                              const Box&       a_edgebox,
                                              const DataIndex& a_di,
                                              int              a_dir,
                                              int              a_ref) const
{
    CH_TIME("MappedAMRPoissonOp::getFluxAfterExchange");
    TODO();
    UNDEFINED_FUNCTION();

#ifndef NDEBUG
    // Set up dimensional info
    int dims = SpaceDim;
    if (m_horizontalOp) --dims;

    // Sanity checks
    CH_assert(a_dir >= 0);
    CH_assert(a_dir < dims);
    CH_assert(!a_data.box().isEmpty());
    CH_assert(!a_edgebox.isEmpty());
    CH_assert( a_flux.contains(a_edgebox) );
#endif

    // Create references for convenience
    const DisjointBoxLayout& grids = m_FCJgup->getBoxes();
    const Box& validCCBox = grids[a_di];
    const FArrayBox& JgupFAB = (*m_FCJgup)[a_di][a_dir];

    // Compute fluxes over faces
    SideIterator fsit;
    for (fsit.reset(); fsit.ok(); ++fsit) {
        Box FCBox = bdryBox(validCCBox, a_dir, fsit(), 1);
        FCBox &= a_edgebox;

        if (!FCBox.isEmpty()) {
            if (m_isDiagonal) {
                // Don't need to do anything special for ortho + horizontalOp.
                const Real scale = m_beta * a_ref / m_dx[a_dir];
                FORT_MAPPEDGETFLUXORTHO(CHF_FRA(a_flux),
                                        CHF_CONST_FRA(a_data),
                                        CHF_CONST_FRA1(JgupFAB,a_dir),
                                        CHF_BOX(FCBox),
                                        CHF_CONST_REAL(scale),
                                        CHF_CONST_INT(a_dir));
            } else {
                const Real scale = m_beta * a_ref;
                if (!m_horizontalOp) {
                    FORT_MAPPEDGETFLUX(CHF_FRA(a_flux),
                                       CHF_CONST_FRA(a_data),
                                       CHF_CONST_FRA(a_extrap),
                                       CHF_CONST_FRA(JgupFAB),
                                       CHF_BOX(FCBox),
                                       CHF_CONST_REAL(scale),
                                       CHF_CONST_REALVECT(m_dx),
                                       CHF_CONST_INT(a_dir));
                } else {
#if CH_SPACEDIM == 3
                    FORT_MAPPEDGETFLATFLUX3D(CHF_FRA(a_flux),
                                             CHF_CONST_FRA(a_data),
                                             CHF_CONST_FRA(a_extrap),
                                             CHF_CONST_FRA(JgupFAB),
                                             CHF_BOX(FCBox),
                                             CHF_CONST_REAL(scale),
                                             CHF_CONST_REALVECT(m_dx),
                                             CHF_CONST_INT(a_dir));
#else
                    MayDay::Error("FORT_MAPPEDGETFLATFLUX2D is not yet written");
#endif
                }
            }
        }
    }
}


// -----------------------------------------------------------------------------
// Computes fluxes of a_data at faces between valid cells of a_data.
// Only does one direction, the edgebox is passed in, and the flux array
// must already be defined.
//
// This can be used as a gradient funciton since it calculates
//   J * Grad[a_data]^{a_dir} over the FC a_edgebox.
//
// One sided differencing in the off-diagonal derivatives is done via
// extrapolation of a_data's ghosts (stored in extrap) and then using a
// standard gradient stencil. extrap must have it's ghosts extrapolated and
// a_data must have it's BCs set prior to call.
//
// This version only updates all faces.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::getFluxComplete (FArrayBox&       a_flux,
                                          const FArrayBox& a_data,
                                          const FArrayBox& a_extrap,
                                          const Box&       a_edgebox,
                                          const DataIndex& a_di,
                                          int              a_dir,
                                          int              a_ref) const
{
    CH_TIME("MappedAMRPoissonOp::getFluxComplete");

#ifndef NDEBUG
    // Set up dimensional info
    int dims = SpaceDim;
    if (m_horizontalOp) --dims;

    // Sanity checks
    CH_assert(a_dir >= 0);
    CH_assert(a_dir < dims);
    CH_assert(!a_data.box().isEmpty());
    CH_assert(!a_edgebox.isEmpty());
    CH_assert( a_flux.contains(a_edgebox) );
    CH_assert(m_FCJgup != NULL);
#endif

    // Create references for convenience
    const FArrayBox& JgupFAB = (*m_FCJgup)[a_di][a_dir];

    if (!a_edgebox.isEmpty()) {
        if (m_isDiagonal) {
            // Don't need to do anything special for ortho + horizontalOp.
            const Real scale = m_beta * a_ref / m_dx[a_dir];
            FORT_MAPPEDGETFLUXORTHO(
                CHF_FRA(a_flux),
                CHF_CONST_FRA(a_data),
                CHF_CONST_FRA1(JgupFAB,a_dir),
                CHF_BOX(a_edgebox),
                CHF_CONST_REAL(scale),
                CHF_CONST_INT(a_dir));
        } else {
            const Real scale = m_beta * a_ref;
            if (!m_horizontalOp) {
                FORT_MAPPEDGETFLUX(
                    CHF_FRA(a_flux),
                    CHF_CONST_FRA(a_data),
                    CHF_CONST_FRA(a_extrap),
                    CHF_CONST_FRA(JgupFAB),
                    CHF_BOX(a_edgebox),
                    CHF_CONST_REAL(scale),
                    CHF_CONST_REALVECT(m_dx),
                    CHF_CONST_INT(a_dir));
            } else {
#if CH_SPACEDIM == 3
                FORT_MAPPEDGETFLATFLUX3D(
                    CHF_FRA(a_flux),
                    CHF_CONST_FRA(a_data),
                    CHF_CONST_FRA(a_extrap),
                    CHF_CONST_FRA(JgupFAB),
                    CHF_BOX(a_edgebox),
                    CHF_CONST_REAL(scale),
                    CHF_CONST_REALVECT(m_dx),
                    CHF_CONST_INT(a_dir));
#else
                CH_assert(a_dir == 0);
                FORT_MAPPEDGETFLATFLUX2D(
                    CHF_FRA(a_flux),
                    CHF_CONST_FRA(a_data),
                    CHF_CONST_FRA1(JgupFAB,0),
                    CHF_BOX(a_edgebox),
                    CHF_CONST_REAL(scale),
                    CHF_CONST_REAL(m_dx[0]));
#endif
            }
        } // end if diagonal or not
    } // end if FCBox is not empty
}


// -----------------------------------------------------------------------------
// For TGA stuff
// Since fillGrad is a noop, this function must perform the flux calculation.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::getFlux (FluxBox&                    a_flux,
                                  const LevelData<FArrayBox>& a_data,
                                  const Box&                  a_grid,
                                  const DataIndex&            a_dit,
                                  Real                        a_scale)
{
    // Set up dimensional info
    int dims = SpaceDim;
    if (m_horizontalOp) --dims;

    const FArrayBox& dataFAB = a_data[a_dit];
    a_flux.define(a_grid, dataFAB.nComp());

    for (int idir = 0; idir < dims; ++idir) {
        this->getFlux(a_flux[idir],
                      dataFAB,
                      a_flux[idir].box(),
                      a_dit,
                      idir,
                      1.0);
        a_flux[idir] *= a_scale;
    }
}


// -----------------------------------------------------------------------------
// Debugging function
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::outputAMR (const Vector<LevelData<FArrayBox>*>& a_rhs,
                                    const std::string a_name)
{
    Vector<string> vNames(a_rhs[0]->nComp());
    for (int comp = 0; comp < a_rhs[0]->nComp(); ++comp) {
        char compname[100];
        sprintf(compname, "comp_%03d", comp);
        string compstr(compname);
        vNames[comp] = compstr;
    }

    Vector<DisjointBoxLayout> vGrids(a_rhs.size());
    Vector<IntVect> refRatios(a_rhs.size());
    for (int lev = 0; lev < a_rhs.size(); ++lev) {
        vGrids[lev] = a_rhs[lev]->getBoxes();
        if (lev > 0) {
            Box fineBox = a_rhs[lev]->getBoxes().physDomain().domainBox();
            Box crseBox = a_rhs[lev-1]->getBoxes().physDomain().domainBox();
            refRatios[lev-1] = fineBox.size() / crseBox.size();
        }
    }

    Real dt = 0.0;
    Real dummyTime = 0.0;
    WriteAnisotropicAMRHierarchyHDF5(a_name,
                                     vGrids, a_rhs, vNames,
                                     vGrids[0].physDomain().domainBox(),
                                     m_dx, dt, dummyTime,
                                     refRatios, a_rhs.size());
}


// -----------------------------------------------------------------------------
// Calls MappedQuadCFInterp or homogeneousCFInterp. Also fills corner ghosts.
// This operator is the finer level.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::interpCFGhosts (LevelData<FArrayBox>&       a_phi,
                                         const LevelData<FArrayBox>* a_phiCoarsePtr,
                                         const bool                  a_isHomogeneous) const
{
    if (a_isHomogeneous) {
        const IntVect activeDirs = this->getActiveDirs();
        homogeneousCFInterp(a_phi, m_dx, m_dxCrse, m_cfregion, activeDirs);
        ExtrapolateCFEV(a_phi, m_cfregion, 2, activeDirs);

    } else {
        CH_assert(a_phiCoarsePtr != NULL); // If inhomogeneous, we need coarse level data.
        CH_assert(a_phi.getBoxes().physDomain() == m_domain); // This level must be the fine level.
        CH_assert(!m_horizontalOp); // Because we don't have the coarse phi available.

        const DisjointBoxLayout& fineGrids = a_phi.getBoxes();
        const DisjointBoxLayout& crseGrids = a_phiCoarsePtr->getBoxes();
        const ProblemDomain& crseDomain = crseGrids.physDomain();
        const int ncomps = a_phi.nComp();
        const IntVect refRatio = this->refToCoarser();
        const IntVect activeDirs = this->getActiveDirs();

        // This replaces valid-adjacent ghosts with unlimited, 2nd-order values.
        MappedQuadCFInterp& castQuadInterp = (MappedQuadCFInterp&)m_interpWithCoarser;
        castQuadInterp.coarseFineInterp(a_phi, *a_phiCoarsePtr);

        ExtrapolateCFEV(a_phi, m_cfregion, 2, activeDirs);
    }
}


// -----------------------------------------------------------------------------
// Exchanges all necessary ghosts.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::exchangeComplete (LevelData<FArrayBox>& a_phi) const
{
    // phiRef.exchange(m_exchangeCopier);

    const DisjointBoxLayout& grids = a_phi.getBoxes();
    const IntVect& ghostVect = a_phi.ghostVect();

    Copier exCopier(grids, grids, m_domain, ghostVect, true);
    a_phi.exchange(exCopier);

    if (!m_horizontalOp) {
        CornerCopier exCornerCopier(grids, grids, m_domain, ghostVect, true);
        a_phi.exchange(exCornerCopier);
    } else {
        CH_assert(ghostVect == IntVect::Unit - BASISV(SpaceDim-1));
    }
}


// -----------------------------------------------------------------------------
// Copies a_phi and extrapolates ghosts at physical boundaries.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::fillExtrap (FArrayBox&       a_extrap,
                                     const FArrayBox& a_phi,
                                     const int        a_order) const
{
#ifndef NDEBUG
    a_extrap.setVal(quietNAN);
#endif

    // If metric is diagonal, we shouldn't need a_extrap.
    if (m_isDiagonal) return;

    // Calculate the region that contains valid data.
    const Box& stateBox = a_phi.box();
    Box validPhi = stateBox & m_validDomain; // Used to be grids[dit]
    const int dims = m_horizontalOp? SpaceDim-1: SpaceDim;

    // Copy the valid values
    a_extrap.copy(a_phi, validPhi);

    // Extrapolate the invalid data.
    for (int fdir = 0; fdir < dims; ++fdir) {
        ExtrapolateFaceAndCopy(a_extrap, a_extrap, validPhi, fdir, Side::Lo, a_order);
        ExtrapolateFaceAndCopy(a_extrap, a_extrap, validPhi, fdir, Side::Hi, a_order);
        validPhi.grow(fdir, 1);
        validPhi &= stateBox;
    }
}
