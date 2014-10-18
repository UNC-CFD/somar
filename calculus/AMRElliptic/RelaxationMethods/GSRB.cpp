/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Copyright (C) 2014 Edward Santilli & Alberto Scotti
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
#include "GSRB.H"
#include "GSRBF_F.H"

// NOTE: Since I usually edit the relax() functions more often than the BaseGSRB
// functions, I placed the messy BaseGSRB functions at the end of this file.


// -----------------------------------------------------------------------------
// The complete GSRB method with no shortcuts for speed.
// -----------------------------------------------------------------------------
void LevelGSRB::relax (LevelData<FArrayBox>&       a_phi,
                       const LevelData<FArrayBox>& a_rhs)
{
    CH_TIME("LevelGSRB::relax");

    // Sanity checks
    CH_assert(isDefined());
    CH_assert(a_phi.isDefined());
    CH_assert(a_rhs.isDefined());
    CH_assert(a_phi.ghostVect() >= m_activeDirs);
    CH_assert(a_phi.nComp() == a_rhs.nComp());
    CH_assert(a_phi.getBoxes().compatible(a_rhs.getBoxes()));
    CH_assert(a_phi.getBoxes().compatible(m_FCJgup->getBoxes()));

    // Gather needed info
    const DisjointBoxLayout& grids = a_phi.disjointBoxLayout();
    DataIterator dit = a_phi.dataIterator();
    const ProblemDomain& domain = grids.physDomain();
    const Box domInterior = grow(domain.domainBox(), -m_activeDirs);

    // Loop over red and black passes
    for (int whichPass = 0; whichPass < 2; ++whichPass) {
        // Synchronous communication
        a_phi.exchange(m_exchangeCopier);

        // Fill boundary ghosts and extrapolate.
        this->fillGhostsAndExtrapolate(a_phi);

        // Interior relaxation
        for (dit.reset(); dit.ok(); ++dit) {
            Box validInterior = grids[dit];
            validInterior &= domInterior;

            this->fullStencilGSRB(a_phi[dit], a_rhs[dit], validInterior, dit(), whichPass);
        }

        // Boundary relaxation
        this->boundaryGSRB(a_phi, a_rhs, whichPass, false);

    } // end loop over red and black passes (whichPass)
}


// -----------------------------------------------------------------------------
// GSRB that only performs one exchange instead of two, asynchronously.
// -----------------------------------------------------------------------------
void LooseGSRB::relax (LevelData<FArrayBox>&       a_phi,
                       const LevelData<FArrayBox>& a_rhs)
{
    CH_TIME("LooseGSRB::relax");

    // Sanity checks
    CH_assert(isDefined());
    CH_assert(a_phi.isDefined());
    CH_assert(a_rhs.isDefined());
    CH_assert(a_phi.ghostVect() >= m_activeDirs);
    CH_assert(a_phi.nComp() == a_rhs.nComp());
    CH_assert(a_phi.getBoxes().compatible(a_rhs.getBoxes()));
    CH_assert(a_phi.getBoxes().compatible(m_FCJgup->getBoxes()));

    // Gather needed info
    const DisjointBoxLayout& grids = a_phi.disjointBoxLayout();
    DataIterator dit = a_phi.dataIterator();
    const ProblemDomain& domain = grids.physDomain();

    // Begin communication
    a_phi.exchangeBegin(m_exchangeCopier);

    // Fill boundary ghosts and extrapolate.
    this->fillGhostsAndExtrapolate(a_phi);

    // Loop through grids and perform calculations that do not need ghosts.
    for (dit.reset(); dit.ok(); ++dit) {
        const Box interior = grow(grids[dit], -m_activeDirs);
        this->fullStencilGSRB(a_phi[dit], a_rhs[dit], interior, dit(), 0);
        this->fullStencilGSRB(a_phi[dit], a_rhs[dit], interior, dit(), 1);
    }

    // End communication
    a_phi.exchangeEnd();

    // Boundary sweeps
    this->boundaryGSRB(a_phi, a_rhs, 0, true);
    this->boundaryGSRB(a_phi, a_rhs, 1, true);
}


// -----------------------------------------------------------------------------
// GSRB that uses a tridiagonal solver for vertical derivatives.
// WARINING: You better not be vertically decomposing your grids!
// -----------------------------------------------------------------------------
void LineGSRB::relax (LevelData<FArrayBox>&       a_phi,
                      const LevelData<FArrayBox>& a_rhs)
{
    CH_TIME("LineGSRB::relax");

    // Sanity checks
    CH_assert(isDefined());
    CH_assert(a_phi.isDefined());
    CH_assert(a_rhs.isDefined());
    CH_assert(a_phi.ghostVect() >= m_activeDirs);
    CH_assert(a_phi.nComp() == a_rhs.nComp());
    CH_assert(a_phi.getBoxes().compatible(a_rhs.getBoxes()));
    CH_assert(a_phi.getBoxes().compatible(m_FCJgup->getBoxes()));

    // It doesn't make sense to use this in a horizontal solve.
    CH_assert(m_activeDirs == IntVect::Unit);

    // Gather geometric info
    const DisjointBoxLayout& grids = a_phi.getBoxes();
    DataIterator dit = grids.dataIterator();
    const ProblemDomain& domain = grids.physDomain();
    const Box& domBox = domain.domainBox();
    const BCDescriptor& fluxDesc = m_bc.getFluxDescriptor();
    int stencil[CH_SPACEDIM][2];

    // Loop over red and black passes.
    for (int whichPass = 0; whichPass < 2; ++whichPass) {
        // Fill ghosts.
        a_phi.exchange(m_exchangeCopier);
        this->fillGhostsAndExtrapolate(a_phi); // TODO: Try not extrapolating on black

        // Loop over each grid in the domain.
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& phiFAB = a_phi[dit];
            const FArrayBox& rhsFAB = a_rhs[dit];

            // This is the box over which we need a solution.
            // phiFAB is defined over valid + 1 ghost layer.
            // rhsFAB is defined over valid with no ghosts.
            const Box& valid = grids[dit];

            // These contain all of the metric info you will ever need.
            // In 2D:
            //  JgupFB[0] is FC in the x-dir. It's 2 components are J*g^{0 0} and J*g^{0 1}.
            //  JgupFB[1] is FC in the y-dir. It's 2 components are J*g^{1 0} and J*g^{1 1}.
            //
            // In 3D:
            //  JgupFB[0] is FC in the x-dir. It's 3 components are J*g^{0 0}, J*g^{0 1}, and J*g^{0 2}.
            //  JgupFB[1] is FC in the y-dir. It's 3 components are J*g^{1 0}, J*g^{1 1}, and J*g^{1 2}.
            //  JgupFB[2] is FC in the z-dir. It's 3 components are J*g^{2 0}, J*g^{2 1}, and J*g^{2 2}.
            //
            // In any number of dimensions:
            //  JinvFAB is CC and has only 1 component equal to 1.0 / the Jacobian determinant.
            const FluxBox&   JgupFB     = (*m_FCJgup)[dit];
            const FArrayBox& JinvFAB    = (*m_CCJinv)[dit];

            // extrapFAB is exactly equal to phiFAB inside valid. However, it's ghost cells have been
            // extrapolated. If you use extrapFAB to compute a derivative at the boundary, the result
            // will be the same as if you used a biased stencil to avoid the ghost cell. You MUST use
            // extrapFAB to compute cross-derivative terms such as
            //   J^{-1} \partial_0 Jg^{0 1} \partial_1 \phi
            // and use phiFAB for the diagonal entries
            //   J^{-1} \partial_0 Jg^{0 0} \partial_0 \phi
            const FArrayBox& extrapFAB  = m_extrap[dit];

            // More sanity checks
            CH_assert(phiFAB    .box().contains(valid));
            CH_assert(rhsFAB    .box().contains(valid));
            CH_assert(JgupFB    .box().contains(valid));
            CH_assert(JinvFAB   .box().contains(valid));
            CH_assert(extrapFAB .box().contains(valid));

            // Gather BC types. Once in fortran, the stencil values will be one of these:
            //  BCType_None.......Use the ghost. It came from an exchange.
            //  BCType_Neum.......Set the zero Neumann BC explicitly. DO NOT use the ghost!
            //  BCType_Diri.......Use the ghost. It is interpolated from this level and a zero Dirichlet BC.
            //  BCType_Periodic...Use the ghost. It came from an neighboring grid on this level.
            //  BCType_CF.........Use the ghost. It is interpolated from this level and a zeroed-out, coarser level.
            //  BCType_Undefined..Throw an error. You should never see this value.
            //
            // If you get a BCType_None in the vertical, throw and error. This means the grids
            // have been decomposed improperly for a vertical tridiagonal solve.
            //
            // If you get BCType_Periodic in the vertical, throw an error for now. There's no
            // reason why you should have a periodic vertical BC for internal wave problems!
            //
            // Once we know the line solver is working, we can speed things up by
            // creating these stencil values in the lineGSRB constructor.
            D_TERM(
            stencil[0][Side::Lo] = fluxDesc.stencil(valid, domain, 0, Side::Lo);
            stencil[0][Side::Hi] = fluxDesc.stencil(valid, domain, 0, Side::Hi);,
            stencil[1][Side::Lo] = fluxDesc.stencil(valid, domain, 1, Side::Lo);
            stencil[1][Side::Hi] = fluxDesc.stencil(valid, domain, 1, Side::Hi);,
            stencil[2][Side::Lo] = fluxDesc.stencil(valid, domain, 2, Side::Lo);
            stencil[2][Side::Hi] = fluxDesc.stencil(valid, domain, 2, Side::Hi);)

            // Examples can be found in various places...
            // LevelLepticSolver.cpp: verticalLineSolver
            //   Sets up a vertical line solve and calls fortran to do the heavy lifting.
            // LevelLepticSolverF.ChF: LepticLapackVerticalSolver
            //   Sets up a vertical tridiagonal matrix equation and solves using Lapack.
            // TridiagUtils.cpp: PoissonNN1D
            //   Sets up a Neumann-Neumann tridiagonal probelem and solves using a
            //   modified version of the Thomas algorithm to deal with the zero eigenvalue.
            //   The modified solver is in TriDiagUtilsF.ChF: TriDiagPoissonNN1DFAB.
            // ThomasAlgorithm.f90
            //   Contains two fortran 90 functions that solve tridiagonal systems.

            // Create workspace for Fortran
            int Nz = valid.size(SpaceDim-1);
            // CH_assert(Nz == domain.size(SpaceDim-1)); // Will fail with AMR!!!

            Box workspace(IntVect::Unit, IntVect(D_DECL(Nz,1,1)));
            BaseFab<Real> D(workspace, 1);
            BaseFab<Real> B(workspace, 1);

            workspace.setSmall(0, 2);
            BaseFab<Real> DU(workspace, 1);

            workspace.setSmall(0, 1);
            workspace.setBig(0, Nz-1);
            BaseFab<Real> DL(workspace, 1);

            // We only want to preserve the horizontal index to compute
            // the red-black ordering.
            IntVect validShift = IntVect::Zero;
            validShift[SpaceDim-1] = valid.smallEnd(SpaceDim-1);

            if (SpaceDim != 2) MayDay::Error("Need to write 3D line solver");

            // Solve over valid region.
            for (int comp = 0; comp < a_phi.nComp(); ++comp) {
                FORT_LINEGSRBITER2D(
                    CHF_FRA1_SHIFT(phiFAB, comp, validShift),
                    CHF_CONST_FRA1_SHIFT(extrapFAB, comp, validShift),
                    CHF_CONST_FRA1_SHIFT(rhsFAB, comp, validShift),
                    CHF_CONST_FRA_SHIFT(JgupFB[0], validShift),
                    CHF_CONST_FRA_SHIFT(JgupFB[1], validShift),
                    CHF_CONST_FRA1_SHIFT(JinvFAB, 0, validShift),
                    CHF_BOX_SHIFT(valid, validShift),
                    CHF_CONST_REALVECT(m_dx),
                    CHF_CONST_REAL(m_crseDx[SpaceDim-1]),
                    CHF_CONST_REAL(m_alpha),
                    CHF_CONST_REAL(m_beta),
                    CHF_CONST_INT(whichPass),
                    CHF_FRA1(DU,0),
                    CHF_FRA1(D,0),
                    CHF_FRA1(DL,0),
                    CHF_FRA1(B,0),
                    CHF_CONST_INT(stencil[0][Side::Lo]),
                    CHF_CONST_INT(stencil[0][Side::Hi]),
                    CHF_CONST_INT(stencil[1][Side::Lo]),
                    CHF_CONST_INT(stencil[1][Side::Hi]));
            } // end loop over components (comp)
        } // end loop over grids (dit)
    } // end loop over red and black passes (whichPass)
}


// -----------------------------------------------------------------------------
// Applies one sweep of Red-Black Gauss-Seidel to a_phi on just one box.
// This version uses the full stencil and does not take special care at the
// box boundaries.
//
// NOTE: At the moment, m_activeDirs must be (1,1) or (1,0) in 2D and
// (1,1,1) or (1,1,0) in 3D.
// -----------------------------------------------------------------------------
void BaseGSRB::fullStencilGSRB (FArrayBox&       a_phi,
                                const FArrayBox& a_rhs,
                                const Box&       a_region,
                                const DataIndex& a_index,
                                const int        a_whichPass) const
{
    CH_TIME("BaseGSRB::fullStencilGSRB");

    // Sanity checks
    CH_assert(isDefined());
    CH_assert(a_phi.nComp() == m_extrap.nComp()); // If this trips, send comps into define.
    CH_assert(a_phi.nComp() == a_rhs.nComp());
    CH_assert(a_phi.box().contains(a_region));
    CH_assert(a_rhs.box().contains(a_region));
    CH_assert(a_whichPass == 0 || a_whichPass == 1);

    // Grab references to metric and diagonals.
    const FluxBox&   JgupFB     = (*m_FCJgup)[a_index];
    const FArrayBox& JinvFAB    = (*m_CCJinv)[a_index];
    const FArrayBox& lapDiagFAB = (*m_lapDiag)[a_index];
    const FArrayBox& extrapFAB  = m_extrap[a_index];

    // More sanity checks
    CH_assert(JgupFB.box().contains(a_region));
    CH_assert(JinvFAB.box().contains(a_region));
    CH_assert(lapDiagFAB.box().contains(a_region));
    CH_assert(extrapFAB.box().contains(a_region));

    // For now, m_activeDirs can't be just anything.
    const bool horizontalOp = (m_activeDirs[SpaceDim-1] == 0);
    D_TERM(,
           CH_assert(m_activeDirs[0] == 1);,
           CH_assert(m_activeDirs[1] == 1);)

    // Call the appropriate Fortran routine.
    if (CH_SPACEDIM == 2) {
        if (!horizontalOp) {
            if (m_isDiagonal) {
                FORT_GSRBITER2DORTHO(CHF_FRA(a_phi),
                                     CHF_CONST_FRA(a_rhs),
                                     CHF_CONST_FRA1(JgupFB[0],0),
                                     CHF_CONST_FRA1(JgupFB[1],1),
                                     CHF_CONST_FRA1(JinvFAB,0),
                                     CHF_CONST_FRA1(lapDiagFAB,0),
                                     CHF_BOX(a_region),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REAL(m_alpha),
                                     CHF_CONST_REAL(m_beta),
                                     CHF_CONST_INT(a_whichPass));
            } else {
                FORT_GSRBITER2D(CHF_FRA(a_phi),
                                CHF_CONST_FRA(extrapFAB),
                                CHF_CONST_FRA(a_rhs),
                                CHF_CONST_FRA(JgupFB[0]),
                                CHF_CONST_FRA(JgupFB[1]),
                                CHF_CONST_FRA1(JinvFAB,0),
                                CHF_CONST_FRA1(lapDiagFAB,0),
                                CHF_BOX(a_region),
                                CHF_CONST_REALVECT(m_dx),
                                CHF_CONST_REAL(m_alpha),
                                CHF_CONST_REAL(m_beta),
                                CHF_CONST_INT(a_whichPass));
            }
        } else {
            const int lineDir = 0;
            FORT_GSRBITER1D(CHF_FRA(a_phi),
                            CHF_CONST_FRA(a_rhs),
                            CHF_CONST_FRA(JgupFB[0]),
                            CHF_CONST_FRA1(JinvFAB,0),
                            CHF_CONST_FRA1(lapDiagFAB,0),
                            CHF_BOX(a_region),
                            CHF_CONST_REALVECT(m_dx),
                            CHF_CONST_REAL(m_alpha),
                            CHF_CONST_REAL(m_beta),
                            CHF_CONST_INT(lineDir),
                            CHF_CONST_INT(a_whichPass));
        }
    } else if (CH_SPACEDIM == 3) {
        if (!horizontalOp) {
            if (m_isDiagonal) {
                FORT_GSRBITER3DORTHO(CHF_FRA(a_phi),
                                     CHF_CONST_FRA(a_rhs),
                                     CHF_CONST_FRA1(JgupFB[0],0),
                                     CHF_CONST_FRA1(JgupFB[1],1),
                                     CHF_CONST_FRA1(JgupFB[2],2),
                                     CHF_CONST_FRA1(JinvFAB,0),
                                     CHF_CONST_FRA1(lapDiagFAB,0),
                                     CHF_BOX(a_region),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REAL(m_alpha),
                                     CHF_CONST_REAL(m_beta),
                                     CHF_CONST_INT(a_whichPass));
            } else {
                FORT_GSRBITER3D(CHF_FRA(a_phi),
                                CHF_CONST_FRA(extrapFAB),
                                CHF_CONST_FRA(a_rhs),
                                CHF_CONST_FRA(JgupFB[0]),
                                CHF_CONST_FRA(JgupFB[1]),
                                CHF_CONST_FRA(JgupFB[2]),
                                CHF_CONST_FRA1(JinvFAB,0),
                                CHF_CONST_FRA1(lapDiagFAB,0),
                                CHF_BOX(a_region),
                                CHF_CONST_REALVECT(m_dx),
                                CHF_CONST_REAL(m_alpha),
                                CHF_CONST_REAL(m_beta),
                                CHF_CONST_INT(a_whichPass));
            }
        } else {
            if (m_isDiagonal) {
                FORT_GSRBITER2DORTHO(CHF_FRA(a_phi),
                                     CHF_CONST_FRA(a_rhs),
                                     CHF_CONST_FRA1(JgupFB[0],0),
                                     CHF_CONST_FRA1(JgupFB[1],1),
                                     CHF_CONST_FRA1(JinvFAB,0),
                                     CHF_CONST_FRA1(lapDiagFAB,0),
                                     CHF_BOX(a_region),
                                     CHF_CONST_REALVECT(m_dx),
                                     CHF_CONST_REAL(m_alpha),
                                     CHF_CONST_REAL(m_beta),
                                     CHF_CONST_INT(a_whichPass));
            } else {
                FORT_GSRBITER2D(CHF_FRA(a_phi),
                                CHF_CONST_FRA(extrapFAB),
                                CHF_CONST_FRA(a_rhs),
                                CHF_CONST_FRA(JgupFB[0]),
                                CHF_CONST_FRA(JgupFB[1]),
                                CHF_CONST_FRA1(JinvFAB,0),
                                CHF_CONST_FRA1(lapDiagFAB,0),
                                CHF_BOX(a_region),
                                CHF_CONST_REALVECT(m_dx),
                                CHF_CONST_REAL(m_alpha),
                                CHF_CONST_REAL(m_beta),
                                CHF_CONST_INT(a_whichPass));
            }
        }
    } else {
        MayDay::Error("Bad SpaceDim");
    }
}


// -----------------------------------------------------------------------------
// Relaxes only on cells that abut boundaries and special care is taken when
// choosing the stencils. If the a_doAllBoundaries flag is set, we iterate over
// every boundary of every grid in the layout, not just at physical boundaries.
//
// NOTE: At the moment, m_activeDirs must be (1,1) or (1,0) in 2D and
// (1,1,1) or (1,1,0) in 3D.
// -----------------------------------------------------------------------------
void BaseGSRB::boundaryGSRB (LevelData<FArrayBox>&       a_phi,
                             const LevelData<FArrayBox>& a_rhs,
                             const int                   a_whichPass,
                             const bool                  a_doAllBoundaries) const
{
    CH_TIME("BaseGSRB::boundaryGSRB");

    // Sanity checks
    CH_assert(isDefined());
    CH_assert(a_phi.nComp() == m_extrap.nComp()); // If this trips, send comps into define.
    CH_assert(a_phi.nComp() == a_rhs.nComp());
    CH_assert(a_phi.getBoxes().compatible(m_extrap.getBoxes()));
    CH_assert(a_phi.getBoxes().compatible(m_FCJgup->getBoxes()));
    CH_assert(a_phi.getBoxes().physDomain() == m_extrap.getBoxes().physDomain());
    CH_assert(a_phi.getBoxes().physDomain() == m_FCJgup->getBoxes().physDomain());
    CH_assert(a_whichPass == 0 || a_whichPass == 1);

    // For now, m_activeDirs can't be just anything.
    const bool horizontalOp = (m_activeDirs[SpaceDim-1] == 0);
    D_TERM(,
           CH_assert(m_activeDirs[0] == 1);,
           CH_assert(m_activeDirs[1] == 1);)

    const Vector<BoundaryBoxData>& boundaryBoxDataRef
        = (a_doAllBoundaries? m_simpleBoundaryBoxData: m_boundaryBoxData);

    for (int idx = 0; idx < boundaryBoxDataRef.size(); ++idx) {
        const BoundaryBoxData data = boundaryBoxDataRef[idx];

        // Grab references to metric and diagonals.
        const FluxBox&   JgupFB     = (*m_FCJgup)[data.index];
        const FArrayBox& JinvFAB    = (*m_CCJinv)[data.index];
        const FArrayBox& lapDiagFAB = (*m_lapDiag)[data.index];

        // Grab references to fields
        FArrayBox&       phiFAB     = a_phi[data.index];
        const FArrayBox& extrapFAB  = m_extrap[data.index];
        const FArrayBox& rhsFAB     = a_rhs[data.index];

#if CH_SPACEDIM == 2
        if (!horizontalOp) {
            if (m_isDiagonal) {
                FORT_GSRBBOUNDARYITER2DORTHO(CHF_FRA(phiFAB),
                                             CHF_CONST_FRA(rhsFAB),
                                             CHF_CONST_FRA1(JgupFB[0],0),
                                             CHF_CONST_FRA1(JgupFB[1],1),
                                             CHF_CONST_FRA1(JinvFAB,0),
                                             CHF_BOX(data.validBdry),
                                             CHF_CONST_REALVECT(m_dx),
                                             CHF_CONST_REAL(m_alpha),
                                             CHF_CONST_REAL(m_beta),
                                             CHF_CONST_INT(data.stencil[0][Side::Lo]),
                                             CHF_CONST_INT(data.stencil[0][Side::Hi]),
                                             CHF_CONST_INT(data.stencil[1][Side::Lo]),
                                             CHF_CONST_INT(data.stencil[1][Side::Hi]),
                                             CHF_CONST_INT(a_whichPass));
            } else {
                FORT_GSRBBOUNDARYITER2D(CHF_FRA(phiFAB),
                                        CHF_CONST_FRA(extrapFAB),
                                        CHF_CONST_FRA(rhsFAB),
                                        CHF_CONST_FRA(JgupFB[0]),
                                        CHF_CONST_FRA(JgupFB[1]),
                                        CHF_CONST_FRA1(JinvFAB,0),
                                        CHF_BOX(data.validBdry),
                                        CHF_CONST_REALVECT(m_dx),
                                        CHF_CONST_REAL(m_alpha),
                                        CHF_CONST_REAL(m_beta),
                                        CHF_CONST_INT(data.stencil[0][Side::Lo]),
                                        CHF_CONST_INT(data.stencil[0][Side::Hi]),
                                        CHF_CONST_INT(data.stencil[1][Side::Lo]),
                                        CHF_CONST_INT(data.stencil[1][Side::Hi]),
                                        CHF_CONST_INT(a_whichPass));
            }
        } else {
            const int edir = 0;
            FORT_GSRBBOUNDARYITER1D(CHF_FRA(phiFAB),
                                    CHF_CONST_FRA(rhsFAB),
                                    CHF_CONST_FRA(JgupFB[0]),
                                    CHF_CONST_FRA1(JinvFAB,0),
                                    CHF_BOX(data.validBdry),
                                    CHF_CONST_REALVECT(m_dx),
                                    CHF_CONST_REAL(m_alpha),
                                    CHF_CONST_REAL(m_beta),
                                    CHF_CONST_INT(data.stencil[edir][Side::Lo]),
                                    CHF_CONST_INT(data.stencil[edir][Side::Hi]),
                                    CHF_CONST_INT(edir),
                                    CHF_CONST_INT(a_whichPass));
        }
#elif CH_SPACEDIM == 3
        if (!horizontalOp) {
            if (m_isDiagonal) {
                FORT_GSRBBOUNDARYITER3DORTHO(CHF_FRA(phiFAB),
                                             CHF_CONST_FRA(rhsFAB),
                                             CHF_CONST_FRA1(JgupFB[0],0),
                                             CHF_CONST_FRA1(JgupFB[1],1),
                                             CHF_CONST_FRA1(JgupFB[2],2),
                                             CHF_CONST_FRA1(JinvFAB,0),
                                             CHF_BOX(data.validBdry),
                                             CHF_CONST_REALVECT(m_dx),
                                             CHF_CONST_REAL(m_alpha),
                                             CHF_CONST_REAL(m_beta),
                                             CHF_CONST_INT(data.stencil[0][Side::Lo]),
                                             CHF_CONST_INT(data.stencil[0][Side::Hi]),
                                             CHF_CONST_INT(data.stencil[1][Side::Lo]),
                                             CHF_CONST_INT(data.stencil[1][Side::Hi]),
                                             CHF_CONST_INT(data.stencil[2][Side::Lo]),
                                             CHF_CONST_INT(data.stencil[2][Side::Hi]),
                                             CHF_CONST_INT(a_whichPass));
            } else {
                FORT_GSRBBOUNDARYITER3D(CHF_FRA(phiFAB),
                                        CHF_CONST_FRA(extrapFAB),
                                        CHF_CONST_FRA(rhsFAB),
                                        CHF_CONST_FRA(JgupFB[0]),
                                        CHF_CONST_FRA(JgupFB[1]),
                                        CHF_CONST_FRA(JgupFB[2]),
                                        CHF_CONST_FRA1(JinvFAB,0),
                                        CHF_BOX(data.validBdry),
                                        CHF_CONST_REALVECT(m_dx),
                                        CHF_CONST_REAL(m_alpha),
                                        CHF_CONST_REAL(m_beta),
                                        CHF_CONST_INT(data.stencil[0][Side::Lo]),
                                        CHF_CONST_INT(data.stencil[0][Side::Hi]),
                                        CHF_CONST_INT(data.stencil[1][Side::Lo]),
                                        CHF_CONST_INT(data.stencil[1][Side::Hi]),
                                        CHF_CONST_INT(data.stencil[2][Side::Lo]),
                                        CHF_CONST_INT(data.stencil[2][Side::Hi]),
                                        CHF_CONST_INT(a_whichPass));
            }
        } else {
            if (m_isDiagonal) {
                FORT_GSRBBOUNDARYITER2DORTHO(CHF_FRA(phiFAB),
                                             CHF_CONST_FRA(rhsFAB),
                                             CHF_CONST_FRA1(JgupFB[0],0),
                                             CHF_CONST_FRA1(JgupFB[1],1),
                                             CHF_CONST_FRA1(JinvFAB,0),
                                             CHF_BOX(data.validBdry),
                                             CHF_CONST_REALVECT(m_dx),
                                             CHF_CONST_REAL(m_alpha),
                                             CHF_CONST_REAL(m_beta),
                                             CHF_CONST_INT(data.stencil[0][Side::Lo]),
                                             CHF_CONST_INT(data.stencil[0][Side::Hi]),
                                             CHF_CONST_INT(data.stencil[1][Side::Lo]),
                                             CHF_CONST_INT(data.stencil[1][Side::Hi]),
                                             CHF_CONST_INT(a_whichPass));
            } else {
                FORT_GSRBBOUNDARYITER2D(CHF_FRA(phiFAB),
                                        CHF_CONST_FRA(extrapFAB),
                                        CHF_CONST_FRA(rhsFAB),
                                        CHF_CONST_FRA(JgupFB[0]),
                                        CHF_CONST_FRA(JgupFB[1]),
                                        CHF_CONST_FRA1(JinvFAB,0),
                                        CHF_BOX(data.validBdry),
                                        CHF_CONST_REALVECT(m_dx),
                                        CHF_CONST_REAL(m_alpha),
                                        CHF_CONST_REAL(m_beta),
                                        CHF_CONST_INT(data.stencil[0][Side::Lo]),
                                        CHF_CONST_INT(data.stencil[0][Side::Hi]),
                                        CHF_CONST_INT(data.stencil[1][Side::Lo]),
                                        CHF_CONST_INT(data.stencil[1][Side::Hi]),
                                        CHF_CONST_INT(a_whichPass));
            }
        }
#else
#   error Bad SpaceDim
#endif
    }
}
