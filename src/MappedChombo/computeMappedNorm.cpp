#include "computeMappedNorm.H"
#include "computeMappedNormF_F.H"
#include "computeNorm.H"
#include "LevelGeometry.H"
#include "EdgeToCell.H"


// FArrayBox versions...

// -----------------------------------------------------------------------------
// Private function.
// Performs most of the computation for the norm functions. This version is not
// public and does not perform MPI communication or the final root.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
static Real computeLocalMappedNormPow (const LevelData<FArrayBox>& a_phi,
                                       const DisjointBoxLayout*    a_finerGridsPtr,
                                       const LevelGeometry&        a_levGeo,
                                       const int                   a_p,
                                       const int                   a_comp)
{
    // Sanity checks
    CH_assert(0 <= a_p);
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < a_phi.nComp());
    CH_assert(a_levGeo.getBoxes().compatible(a_phi.getBoxes()));

    // Gather geometric structures
    const Real dxScale = a_levGeo.getDx().product();
    const LevelData<FArrayBox>& JRef = a_levGeo.getCCJ();
    const IntVect& fineRefRatio = a_levGeo.getFineRefRatio();
    const DisjointBoxLayout& grids = a_phi.getBoxes();
    DataIterator dit = grids.dataIterator();

    // Loop over this level's grids, adding to norm^p.
    Real localNorm = 0.0;
    for (dit.reset(); dit.ok(); ++dit) {
        // Create references for convenience / allocate calucation space.
        const FArrayBox& phiFAB = a_phi[dit];
        const FArrayBox& JFAB = JRef[dit];
        const Box& valid = grids[dit];
        FArrayBox tempFAB(valid, 1);

        // If needed, we can adopt other centerings later.
        CH_assert(phiFAB.box().type() == IntVect::Zero);

        // Phi should be equal to or larger than the valid region.
        CH_assert(phiFAB.box().contains(valid));

        // Copy the state data.
        const Interval srcInterval(a_comp, a_comp);
        const Interval destInterval(0, 0);
        tempFAB.copy(valid,
                     destInterval,
                     valid,
                     phiFAB,
                     srcInterval);

        // Zero out the cells covered by a finer grid.
        if (a_finerGridsPtr != NULL) {
            LayoutIterator litFine = a_finerGridsPtr->layoutIterator();
            for (litFine.reset(); litFine.ok(); ++litFine) {
                // Calculate covered region
                Box coveredBox(a_finerGridsPtr->get(litFine()));
                coveredBox.coarsen(fineRefRatio);
                coveredBox &= valid;

                // Set data to zero
                if (!coveredBox.isEmpty()) {
                    tempFAB.setVal(0.0, coveredBox, 0, 1);
                }
            }
        } // end if there is a finer level

        // Add to localNorm
        FORT_COMPUTEMAPPEDNORMPOW(CHF_REAL(localNorm),
                                  CHF_CONST_FRA1(tempFAB,0),
                                  CHF_CONST_FRA1(JFAB,0),
                                  CHF_BOX(valid),
                                  CHF_CONST_REAL(dxScale),
                                  CHF_CONST_INT(a_p));
    } // end loop over this level's grids

    return localNorm;
}


// -----------------------------------------------------------------------------
// Returns the volume-weighted norm of phi over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real computeMappedNorm (const Vector<LevelData<FArrayBox>*>& a_phi,
                        const LevelGeometry&                 a_levGeo,
                        const int                            a_p,
                        const int                            a_comp,
                        const int                            a_lBase)
{
    // Sanity check on a_lBase
    const int vectorSize = a_phi.size();
    CH_assert(0 <= a_lBase);
    CH_assert(a_lBase < vectorSize);

    // We need at least one defined level.
    CH_assert(a_phi[a_lBase] != NULL);

    // Find highest defined level.
    int lev;
    for (lev = a_lBase; lev < vectorSize; ++lev) {
        if (a_phi[lev] == NULL) break;
        if (!a_phi[lev]->isDefined()) break;
    }
    const int topLevel = lev - 1;

    // Get every level's levGeo.
    Vector<const LevelGeometry*> vLevGeos = a_levGeo.getAMRLevGeos();

    // Loop over levels and add to total norm^p.
    Real localNormPow = 0.0;
    for (lev = a_lBase; lev <= topLevel; ++lev) {
        const LevelData<FArrayBox>& levelPhi = *a_phi[lev];

        const DisjointBoxLayout* finerGridsPtr = NULL;
        if (lev < topLevel) {
            finerGridsPtr = &(a_phi[lev+1]->getBoxes());
        }

        const LevelGeometry& levGeoRef = *vLevGeos[lev];
        CH_assert(levGeoRef.getBoxes() == levelPhi.getBoxes()); // This may be overdoing it.

        // Add norm^p to running total.
        if (a_p == 0) {
            localNormPow = Max(localNormPow, computeLocalMappedNormPow(levelPhi,
                                                                       finerGridsPtr,
                                                                       levGeoRef,
                                                                       a_p,
                                                                       a_comp));
        } else {
            localNormPow += computeLocalMappedNormPow(levelPhi,
                                                      finerGridsPtr,
                                                      levGeoRef,
                                                      a_p,
                                                      a_comp);
        }
    }

    // Compute global norm (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalNorm = 0.0;
    int result;
    if (a_p == 0) {
        result = MPI_Allreduce(&localNormPow, &globalNorm, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
    } else {
        result = MPI_Allreduce(&localNormPow, &globalNorm, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
    }

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedNorm");
    }

#else
    Real globalNorm = localNormPow;
#endif

    // Take the root, if necessary, and return the result.
    if (a_p != 0) {
        const Real invp = 1.0 / Real(a_p);
        globalNorm = pow(globalNorm, invp);
    }

    return globalNorm;
}


// -----------------------------------------------------------------------------
// Returns the volume-weighted norm of phi over the valid region.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real computeMappedNorm (const LevelData<FArrayBox>& a_phi,
                        const DisjointBoxLayout*    a_finerGridsPtr,
                        const LevelGeometry&        a_levGeo,
                        const int                   a_p,
                        const int                   a_comp)
{
    // Compute the local norm^p.
    // Assertions are done in this function call.
    Real localNormPow = computeLocalMappedNormPow(a_phi,
                                                  a_finerGridsPtr,
                                                  a_levGeo,
                                                  a_p,
                                                  a_comp);

    // Compute global norm (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalNorm = 0.0;
    int result;
    if (a_p == 0) {
        result = MPI_Allreduce(&localNormPow, &globalNorm, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
    } else {
        result = MPI_Allreduce(&localNormPow, &globalNorm, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
    }

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedNorm");
    }

#else
    Real globalNorm = localNormPow;
#endif

    // Take the root, if necessary, and return the result.
    if (a_p != 0) {
        const Real invp = 1.0 / Real(a_p);
        globalNorm = pow(globalNorm, invp);
    }

    return globalNorm;
}


// Unmapped FArrayBox versions...

// -----------------------------------------------------------------------------
// Private function.
// Performs most of the computation for the norm functions. This version is not
// public and does not perform MPI communication or the final root.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
static Real computeLocalUnmappedNormPow (const LevelData<FArrayBox>& a_phi,
                                         const DisjointBoxLayout*    a_finerGridsPtr,
                                         const LevelGeometry&        a_levGeo,
                                         const int                   a_p,
                                         const int                   a_comp)
{
    // Sanity checks
    CH_assert(0 <= a_p);
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < a_phi.nComp());

    // Gather geometric structures
    const Real dxScale = a_levGeo.getDx().product();
    const IntVect& fineRefRatio = a_levGeo.getFineRefRatio();
    const DisjointBoxLayout& grids = a_phi.getBoxes();
    DataIterator dit = grids.dataIterator();

    // Loop over this level's grids, adding to norm^p.
    Real localNorm = 0.0;
    for (dit.reset(); dit.ok(); ++dit) {
        // Create references for convenience / allocate calucation space.
        const FArrayBox& phiFAB = a_phi[dit];
        const Box& valid = grids[dit];
        FArrayBox tempFAB(valid, 1);

        // If needed, we can adopt other centerings later.
        CH_assert(phiFAB.box().type() == IntVect::Zero);

        // Phi should be equal to or larger than the valid region.
        CH_assert(phiFAB.box().contains(valid));

        // Copy the state data.
        const Interval srcInterval(a_comp, a_comp);
        const Interval destInterval(0, 0);
        tempFAB.copy(valid,
                     destInterval,
                     valid,
                     phiFAB,
                     srcInterval);

        // Zero out the cells covered by a finer grid.
        if (a_finerGridsPtr != NULL) {
            LayoutIterator litFine = a_finerGridsPtr->layoutIterator();
            for (litFine.reset(); litFine.ok(); ++litFine) {
                // Calculate covered region
                Box coveredBox(a_finerGridsPtr->get(litFine()));
                coveredBox.coarsen(fineRefRatio);
                coveredBox &= valid;

                // Set data to zero
                if (!coveredBox.isEmpty()) {
                    tempFAB.setVal(0.0, coveredBox, 0, 1);
                }
            }
        } // end if there is a finer level

        // Add to localNorm
        FORT_COMPUTEUNMAPPEDNORMPOW(CHF_REAL(localNorm),
                                    CHF_CONST_FRA1(tempFAB,0),
                                    CHF_BOX(valid),
                                    CHF_CONST_REAL(dxScale),
                                    CHF_CONST_INT(a_p));
    } // end loop over this level's grids

    return localNorm;
}


// -----------------------------------------------------------------------------
// Returns the norm of phi over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real computeUnmappedNorm (const Vector<LevelData<FArrayBox>*>& a_phi,
                          const LevelGeometry&                 a_levGeo,
                          const int                            a_p,
                          const int                            a_comp,
                          const int                            a_lBase)
{
    // Sanity check on a_lBase
    const int vectorSize = a_phi.size();
    CH_assert(0 <= a_lBase);
    CH_assert(a_lBase < vectorSize);

    // We need at least one defined level.
    CH_assert(a_phi[a_lBase] != NULL);

    // Find highest defined level.
    int lev;
    for (lev = a_lBase; lev < vectorSize; ++lev) {
        if (a_phi[lev] == NULL) break;
        if (!a_phi[lev]->isDefined()) break;
    }
    const int topLevel = lev - 1;

    // Get every level's levGeo.
    Vector<const LevelGeometry*> vLevGeos = a_levGeo.getAMRLevGeos();

    // Loop over levels and add to total norm^p.
    Real localNormPow = 0.0;
    for (lev = a_lBase; lev <= topLevel; ++lev) {
        const LevelData<FArrayBox>& levelPhi = *a_phi[lev];
        const LevelGeometry& levGeoRef = *vLevGeos[lev];
        const DisjointBoxLayout* finerGridsPtr = NULL;

        if (lev < topLevel) {
            finerGridsPtr = &(a_phi[lev+1]->getBoxes());
        }

        // Add norm^p to running total.
        if (a_p == 0) {
            localNormPow = Max(localNormPow, computeLocalUnmappedNormPow(levelPhi,
                                                                         finerGridsPtr,
                                                                         levGeoRef,
                                                                         a_p,
                                                                         a_comp));
        } else {
            localNormPow += computeLocalUnmappedNormPow(levelPhi,
                                                        finerGridsPtr,
                                                        levGeoRef,
                                                        a_p,
                                                        a_comp);
        }
    }

    // Compute global norm (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalNorm = 0.0;
    int result;
    if (a_p == 0) {
        result = MPI_Allreduce(&localNormPow, &globalNorm, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
    } else {
        result = MPI_Allreduce(&localNormPow, &globalNorm, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
    }

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedNorm");
    }

#else
    Real globalNorm = localNormPow;
#endif

    // Take the root, if necessary, and return the result.
    if (a_p != 0) {
        const Real invp = 1.0 / Real(a_p);
        globalNorm = pow(globalNorm, invp);
    }

    return globalNorm;
}


// -----------------------------------------------------------------------------
// Returns the norm of phi over the valid region.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real computeUnmappedNorm (const LevelData<FArrayBox>& a_phi,
                          const DisjointBoxLayout*    a_finerGridsPtr,
                          const LevelGeometry&        a_levGeo,
                          const int                   a_p,
                          const int                   a_comp)
{
    // Compute the local norm^p.
    // Assertions are done in this function call.
    Real localNormPow = computeLocalUnmappedNormPow(a_phi,
                                                    a_finerGridsPtr,
                                                    a_levGeo,
                                                    a_p,
                                                    a_comp);

    // Compute global norm (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalNorm = 0.0;
    int result;
    if (a_p == 0) {
        result = MPI_Allreduce(&localNormPow, &globalNorm, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
    } else {
        result = MPI_Allreduce(&localNormPow, &globalNorm, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
    }

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedNorm");
    }

#else
    Real globalNorm = localNormPow;
#endif

    // Take the root, if necessary, and return the result.
    if (a_p != 0) {
        const Real invp = 1.0 / Real(a_p);
        globalNorm = pow(globalNorm, invp);
    }

    return globalNorm;
}


// FluxBox versions...

// -----------------------------------------------------------------------------
// Private function.
// Performs most of the computation for the norm functions. This version is not
// public and does not perform MPI communication or the final root.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
static RealVect computeLocalMappedNormPow (const LevelData<FluxBox>& a_phi,
                                           const DisjointBoxLayout*  a_finerGridsPtr,
                                           const LevelGeometry&      a_levGeo,
                                           const int                 a_p,
                                           const int                 a_comp)
{
    // Sanity checks
    CH_assert(0 <= a_p);
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < a_phi.nComp());
    CH_assert(a_levGeo.getBoxes().compatible(a_phi.getBoxes()));

    // Gather geometric structures
    const Real dxScale = a_levGeo.getDx().product();
    const LevelData<FArrayBox>& JRef = a_levGeo.getCCJ();
    const IntVect& fineRefRatio = a_levGeo.getFineRefRatio();
    const DisjointBoxLayout& grids = a_phi.getBoxes();
    DataIterator dit = grids.dataIterator();

    // Loop over this level's grids and add to norm
    RealVect localNormPow(D_DECL(0.,0.,0.));
    for (dit.reset(); dit.ok(); ++dit) {
        // Create references / allocate work space
        const FluxBox& phiFlux = a_phi[dit];
        const FArrayBox& JFAB = JRef[dit];
        const Box& valid = grids[dit];
        FArrayBox tempFAB(valid, 1);

        // Loop over directions.
        for (int dir = 0; dir < SpaceDim; ++dir) {
            // Send data to cell-centered holder.
            EdgeToCell(phiFlux, a_comp, tempFAB, 0, valid, dir);

            // Zero out the cells covered by a finer grid.
            if (a_finerGridsPtr != NULL) {
                LayoutIterator litFine = a_finerGridsPtr->layoutIterator();
                for (litFine.reset(); litFine.ok(); ++litFine) {
                    // Calculate covered region
                    Box coveredBox(a_finerGridsPtr->get(litFine()));
                    coveredBox.coarsen(fineRefRatio);
                    coveredBox &= valid;

                    // Set data to zero
                    if (!coveredBox.isEmpty()) {
                        tempFAB.setVal(0.0, coveredBox, 0, 1);
                    }
                }
            } // end if there is a finer level

            // Add to norm^p.
            Real& localNormPowDir = localNormPow[dir];
            FORT_COMPUTEMAPPEDNORMPOW(CHF_REAL(localNormPowDir),
                                      CHF_CONST_FRA1(tempFAB,0),
                                      CHF_CONST_FRA1(JFAB,0),
                                      CHF_BOX(valid),
                                      CHF_CONST_REAL(dxScale),
                                      CHF_CONST_INT(a_p));
        }
    } // end loop over this level's grids

    return localNormPow;
}


// -----------------------------------------------------------------------------
// Returns the volume-weighted norm of phi over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// LIMITATIONS: This function can only handle face-centered data.
// -----------------------------------------------------------------------------
RealVect computeMappedNorm (const Vector<LevelData<FluxBox>*>& a_phi,
                            const LevelGeometry&               a_levGeo,
                            const int                          a_p,
                            const int                          a_comp,
                            const int                          a_lBase)
{
    // Sanity check on a_lBase
    const int vectorSize = a_phi.size();
    CH_assert(0 <= a_lBase);
    CH_assert(a_lBase < vectorSize);

    // We need at least one defined level.
    CH_assert(a_phi[a_lBase] != NULL);

    // Find highest defined level.
    int lev;
    for (lev = a_lBase; lev < vectorSize; ++lev) {
        if (a_phi[lev] == NULL) break;
        if (!a_phi[lev]->isDefined()) break;
    }
    const int topLevel = lev - 1;

    // Get every level's levGeo.
    Vector<const LevelGeometry*> vLevGeos = a_levGeo.getAMRLevGeos();

    // Loop over levels and add to total norm^p.
    RealVect localNormPow(D_DECL(0.,0.,0.));
    for (lev = a_lBase; lev <= topLevel; ++lev) {
        const LevelData<FluxBox>& levelPhi = *a_phi[lev];

        const DisjointBoxLayout* finerGridsPtr = NULL;
        if (lev < topLevel) {
            finerGridsPtr = &(a_phi[lev+1]->getBoxes());
        }

        const LevelGeometry& levGeoRef = *vLevGeos[lev];
        CH_assert(levGeoRef.getBoxes() == levelPhi.getBoxes()); // This may be overdoing it.

        // Add norm^p to running total.
        if (a_p == 0) {
            localNormPow = Max(localNormPow, computeLocalMappedNormPow(levelPhi,
                                                                       finerGridsPtr,
                                                                       levGeoRef,
                                                                       a_p,
                                                                       a_comp));
        } else {
            localNormPow += computeLocalMappedNormPow(levelPhi,
                                                      finerGridsPtr,
                                                      levGeoRef,
                                                      a_p,
                                                      a_comp);
        }
    }

    // Compute global norm (this is where the MPI communication happens)
#ifdef CH_MPI
    RealVect globalNorm(D_DECL(0.,0.,0.));
    for (int dir = 0; dir < SpaceDim; ++dir) {
        int result;
        if (a_p == 0) {
            result = MPI_Allreduce(&(localNormPow[dir]), &(globalNorm[dir]), 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
        } else {
            result = MPI_Allreduce(&(localNormPow[dir]), &(globalNorm[dir]), 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
        }

        if (result != MPI_SUCCESS) {
            MayDay::Error("Sorry, but I had a communication error in computeMappedNorm");
        }
    }
#else
    RealVect& globalNorm = localNormPow;
#endif

    // Take the root, if necessary, and return the result.
    if (a_p != 0) {
        const Real invp = 1.0 / Real(a_p);
        D_TERM(globalNorm[0] = pow(globalNorm[0], invp);,
               globalNorm[1] = pow(globalNorm[1], invp);,
               globalNorm[2] = pow(globalNorm[2], invp);)
    }

    return globalNorm;
}


// -----------------------------------------------------------------------------
// Returns the volume-weighted norm of phi over the valid region.
// LIMITATIONS: This function can only handle face-centered data.
// -----------------------------------------------------------------------------
RealVect computeMappedNorm (const LevelData<FluxBox>& a_phi,
                            const DisjointBoxLayout*  a_finerGridsPtr,
                            const LevelGeometry&      a_levGeo,
                            const int                 a_p,
                            const int                 a_comp)
{
    // Compute the local norm^p.
    // Assertions are done in this function call.
    RealVect localNormPow = computeLocalMappedNormPow(a_phi,
                                                      a_finerGridsPtr,
                                                      a_levGeo,
                                                      a_p,
                                                      a_comp);

    // Compute global norm (this is where the MPI communication happens)
#ifdef CH_MPI
    RealVect globalNorm(D_DECL(0.,0.,0.));
    for (int dir = 0; dir < SpaceDim; ++dir) {
        int result;
        if (a_p == 0) {
            result = MPI_Allreduce(&(localNormPow[dir]), &(globalNorm[dir]), 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
        } else {
            result = MPI_Allreduce(&(localNormPow[dir]), &(globalNorm[dir]), 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
        }

        if (result != MPI_SUCCESS) {
            MayDay::Error("Sorry, but I had a communication error in computeMappedNorm");
        }
    }
#else
    RealVect& globalNorm = localNormPow;
#endif

    // Take the root, if necessary, and return the result.
    if (a_p != 0) {
        const Real invp = 1.0 / Real(a_p);
        D_TERM(globalNorm[0] = pow(globalNorm[0], invp);,
               globalNorm[1] = pow(globalNorm[1], invp);,
               globalNorm[2] = pow(globalNorm[2], invp);)
    }

    return globalNorm;
}


// Ghost cell versions (FArrayBox only)...

// -----------------------------------------------------------------------------
// Private function.
// Performs most of the computation for the norm functions. This version is not
// public and does not perform MPI communication or the final root.
// NOTE: This function can only handle cell-centered data.
// NOTE: This does not check for finer data coverings, nor should it since,
// technically, all ghosts are invalid data. All ghosts in the mode are used.
// -----------------------------------------------------------------------------
static Real computeLocalMappedNormPowGhost (const LevelData<FArrayBox>& a_phi,
                                            const LevelGeometry&        a_levGeo,
                                            const GhostNormMode         a_mode,
                                            const int                   a_p,
                                            const int                   a_comp)
{
    // Sanity checks
    CH_assert(0 <= a_p);
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < a_phi.nComp());
    CH_assert(a_levGeo.getBoxes().compatible(a_phi.getBoxes()));

    CH_assert(a_phi.ghostVect().sum() != 0);
    CH_assert(a_levGeo.getCCJ().ghostVect() >= a_phi.ghostVect());

    // Gather geometric structures
    const Real dxScale = a_levGeo.getDx().product();
    const LevelData<FArrayBox>& JRef = a_levGeo.getCCJ();
    const ProblemDomain& domain = a_levGeo.getDomain();
    const Box& domBox = domain.domainBox();
    const DisjointBoxLayout& grids = a_phi.getBoxes();
    DataIterator dit = grids.dataIterator();

    IntVect ghostVect = a_phi.ghostVect();
    D_TERM(if (ghostVect[0] > 1) ghostVect[0] = 1;,
           if (ghostVect[1] > 1) ghostVect[1] = 1;,
           if (ghostVect[2] > 1) ghostVect[2] = 1;)

    // Loop over this level's grids, adding to norm^p.
    Real localNorm = 0.0;
    for (dit.reset(); dit.ok(); ++dit) {
        // Create references for convenience / allocate calucation space.
        const FArrayBox& phiFAB = a_phi[dit];
        const FArrayBox& JFAB = JRef[dit];
        const Box& valid = grids[dit];

        // Phi should be cell-centered with at least one ghost layer.
        CH_assert(phiFAB.box().type() == IntVect::Zero);
        CH_assert(phiFAB.box().contains(grow(valid, ghostVect)));

        // Loop over the boundary directions and sides, skipping periodic dirs
        for (int fdir = 0; fdir < SpaceDim; ++fdir) {
            if (domain.isPeriodic(fdir)) continue;
            if (ghostVect[fdir] == 0) continue;

            SideIterator sit;
            for (sit.reset(); sit.ok(); ++sit) {

                // Find ghost region
                Box ghostBox = adjCellBox(valid, fdir, sit(), 1);

                // Do we want to include these ghosts?
                bool doCalc;
                switch (a_mode) {
                case PHYS_GHOSTS:
                    doCalc = (ghostBox & domBox).isEmpty();
                    break;
                case INTERP_GHOSTS:
                    doCalc = !(ghostBox & domBox).isEmpty();
                    break;
                case ALL_GHOSTS:
                    doCalc = true;
                    break;
                default:
                    MayDay::Error("computeMappedNormGhost: Bad mode");
                }

                if (!doCalc) continue;

                // This may be overkill. We may just want to adjust the
                // ghostBox to its intersection with phiFAB.
                CH_assert(phiFAB.box().contains(ghostBox));

                // Add to localNorm
                FORT_COMPUTEMAPPEDNORMPOW(CHF_REAL(localNorm),
                                          CHF_CONST_FRA1(phiFAB,0),
                                          CHF_CONST_FRA1(JFAB,0),
                                          CHF_BOX(ghostBox),
                                          CHF_CONST_REAL(dxScale),
                                          CHF_CONST_INT(a_p));
            } // end loop over boundary sides (sit)
        } // end loop over boundary directions (fdir)
    } // end loop over this level's grids

    return localNorm;
}


// -----------------------------------------------------------------------------
// Returns the volume-weighted norm of phi's ghosts over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// NOTE: This function can only handle cell-centered data.
// NOTE: This does not check for finer data coverings, nor should it since,
// technically, all ghosts are invalid data. All ghosts in the mode are used.
// -----------------------------------------------------------------------------
Real computeMappedNormGhost (const Vector<LevelData<FArrayBox>*>& a_phi,
                             const LevelGeometry&                 a_levGeo,
                             const GhostNormMode                  a_mode,
                             const int                            a_p,
                             const int                            a_comp,
                             const int                            a_lBase)
{
    // Sanity check on a_lBase
    const int vectorSize = a_phi.size();
    CH_assert(0 <= a_lBase);
    CH_assert(a_lBase < vectorSize);

    // We need at least one defined level.
    CH_assert(a_phi[a_lBase] != NULL);

    // Find highest defined level.
    int lev;
    for (lev = a_lBase; lev < vectorSize; ++lev) {
        if (a_phi[lev] == NULL) break;
        if (!a_phi[lev]->isDefined()) break;
    }
    const int topLevel = lev - 1;

    // Get every level's levGeo.
    Vector<const LevelGeometry*> vLevGeos = a_levGeo.getAMRLevGeos();

    // Loop over levels and add to total norm^p.
    Real localNormPow = 0.0;
    for (lev = a_lBase; lev <= topLevel; ++lev) {
        const LevelData<FArrayBox>& levelPhi = *a_phi[lev];
        const LevelGeometry& levGeoRef = *vLevGeos[lev];
        CH_assert(levGeoRef.getBoxes() == levelPhi.getBoxes()); // This may be overdoing it.

        // Add norm^p to running total.
        if (a_p == 0) {
            localNormPow = Max(localNormPow, computeLocalMappedNormPowGhost(levelPhi,
                                                                            levGeoRef,
                                                                            a_mode,
                                                                            a_p,
                                                                            a_comp));
        } else {
            localNormPow += computeLocalMappedNormPowGhost(levelPhi,
                                                           levGeoRef,
                                                           a_mode,
                                                           a_p,
                                                           a_comp);
        }
    }

    // Compute global norm (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalNorm = 0.0;
    int result;
    if (a_p == 0) {
        result = MPI_Allreduce(&localNormPow, &globalNorm, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
    } else {
        result = MPI_Allreduce(&localNormPow, &globalNorm, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
    }

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedNormGhost");
    }

#else
    Real globalNorm = localNormPow;
#endif

    // Take the root, if necessary, and return the result.
    if (a_p != 0) {
        const Real invp = 1.0 / Real(a_p);
        globalNorm = pow(globalNorm, invp);
    }

    return globalNorm;
}


// -----------------------------------------------------------------------------
// Returns the volume-weighted norm of phi over the first ghost layer.
// NOTE: This function can only handle cell-centered data.
// NOTE: This does not check for finer data coverings, nor should it since,
// technically, all ghosts are invalid data. All ghosts in the mode are used.
// -----------------------------------------------------------------------------
Real computeMappedNormGhost (const LevelData<FArrayBox>& a_phi,
                             const LevelGeometry&        a_levGeo,
                             const GhostNormMode         a_mode,
                             const int                   a_p,
                             const int                   a_comp)
{
    // Compute the local norm^p.
    // Assertions are done in this function call.
    Real localNormPow = computeLocalMappedNormPowGhost(a_phi,
                                                       a_levGeo,
                                                       a_mode,
                                                       a_p,
                                                       a_comp);

    // Compute global norm (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalNorm = 0.0;
    int result;
    if (a_p == 0) {
        result = MPI_Allreduce(&localNormPow, &globalNorm, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
    } else {
        result = MPI_Allreduce(&localNormPow, &globalNorm, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
    }

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedNormGhost");
    }

#else
    Real globalNorm = localNormPow;
#endif

    // Take the root, if necessary, and return the result.
    if (a_p != 0) {
        const Real invp = 1.0 / Real(a_p);
        globalNorm = pow(globalNorm, invp);
    }

    return globalNorm;
}


// -----------------------------------------------------------------------------
// Returns the max value of phi over all valid regions.
// -----------------------------------------------------------------------------
Real computeMax(const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<IntVect>&                a_nRefFine,
                const Interval                        a_comps,
                const int                             a_lBase)
{
    int numLevels = a_phi.size();

    // it is often the case that while a_phi has many possible
    // levels, only a subset of them are defined -- check that
    // just to be sure
    if (a_phi[numLevels - 1] == NULL) {
        int lev = numLevels - 1;

        while (a_phi[lev] == NULL) {
            lev--;
        }

        numLevels = lev + 1;
    }

    Real max;
    Real maxOnLevel;

    max = -HUGE_VAL;

    // loop over levels
    for (int lev = a_lBase; lev < numLevels; lev++) {
        //in case there are extra levels which are not defined
        if (a_phi[lev] != NULL) {
            CH_assert(a_phi[lev]->isDefined());

            LevelData<FArrayBox>& thisPhi = *(a_phi[lev]);
            const DisjointBoxLayout* finerGridsPtr = NULL;

            if (lev < numLevels - 1) {
                finerGridsPtr = &(a_phi[lev + 1]->getBoxes());

                maxOnLevel = computeMax(thisPhi, finerGridsPtr, a_nRefFine[lev],
                                        a_comps);
            } else {
                int bogusRefRatio = 100000;
                maxOnLevel = computeMax(thisPhi, finerGridsPtr, bogusRefRatio,
                                        a_comps);
            }

            if (maxOnLevel > max) {
                max = maxOnLevel;
            }
        }
    }

    // shouldn't need to do broadcast/gather thing
    return max;
}


// -----------------------------------------------------------------------------
// Returns the maximum value of phi over the valid region.
// -----------------------------------------------------------------------------
Real computeMax(const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGridsPtr,
                const IntVect&              a_nRefFine,
                const Interval              a_comps)
{
    Real levelMax = -99999999999.9;
    Real thisMax;

    const DisjointBoxLayout& levelGrids = a_phi.getBoxes();
    DataIterator dit = a_phi.dataIterator();

    LevelData<FArrayBox> temp(levelGrids, a_comps.size());
    Interval tempComps(0, a_comps.size() - 1);

    for (dit.reset(); dit.ok(); ++dit) {
        temp[dit()].copy(temp[dit()].box(), tempComps,
                         temp[dit()].box(), a_phi[dit()], a_comps);

        if (a_finerGridsPtr != NULL) {
            LayoutIterator litFine = a_finerGridsPtr->layoutIterator();

            // now loop over fine boxes and set covered regions to 0
            for (litFine.reset(); litFine.ok(); ++litFine) {
                Box coveredBox(a_finerGridsPtr->get(litFine()));
                coveredBox.coarsen(a_nRefFine);
                coveredBox &= temp[dit()].box();

                if (!coveredBox.isEmpty()) {
                    temp[dit()].setVal(levelMax, coveredBox, 0, tempComps.size());
                }
            }  // end loop over fine-grid boxes
        } // end if there is a finer level

        // while we're looping over the grids, get max as well
        // need to loop over comps
        for (int comp = tempComps.begin(); comp <= tempComps.end(); ++comp) {
            thisMax = temp[dit()].max(comp);
            if (thisMax > levelMax) {
                levelMax = thisMax;
            }
        }
    } // end loop over this level's grids

    // need to do MPI gather here
#ifdef CH_MPI
    Real recv;
    int result = MPI_Allreduce(&levelMax, &recv, 1, MPI_CH_REAL,
                               MPI_MAX, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMax");
    }

    levelMax = recv;
#endif

    return levelMax;
}


// -----------------------------------------------------------------------------
// Returns the min value of phil over all valid regions.
// -----------------------------------------------------------------------------
Real computeMin(const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<IntVect>&                a_nRefFine,
                const Interval                        a_comps,
                const int                             a_lBase)
{
    int numLevels = a_phi.size();

    // it is often the case that while a_phi has many possible
    // levels, only a subset of them are defined -- check that
    // just to be sure
    if (a_phi[numLevels - 1] == NULL) {
        int lev = numLevels - 1;

        while (a_phi[lev] == NULL) {
            lev--;
        }

        numLevels = lev + 1;
    }

    Real min;
    Real minOnLevel;

    min = HUGE_VAL;

    // loop over levels
    for (int lev = a_lBase; lev < numLevels; lev++) {
        //in case there are extra levels which are not defined
        if (a_phi[lev] != NULL) {
            CH_assert(a_phi[lev]->isDefined());

            LevelData<FArrayBox>& thisPhi = *(a_phi[lev]);
            const DisjointBoxLayout* finerGridsPtr = NULL;

            if (lev < numLevels - 1) {
                finerGridsPtr = &(a_phi[lev + 1]->getBoxes());

                minOnLevel = computeMin(thisPhi, finerGridsPtr, a_nRefFine[lev],
                                        a_comps);
            } else {
                int bogusRefRatio = 1000000;
                minOnLevel = computeMin(thisPhi, finerGridsPtr, bogusRefRatio,
                                        a_comps);
            }

            if (minOnLevel < min) {
                min = minOnLevel;
            }
        }
    }

    // shouldn't need to do broadcast/gather thing
    return min;
}


// -----------------------------------------------------------------------------
// Returns the minimum value of phi over the valid region.
// -----------------------------------------------------------------------------
Real computeMin(const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGridsPtr,
                const IntVect&              a_nRefFine,
                const Interval              a_comps)
{
    Real levelMin = 99999999999.9;
    Real thisMin;

    const DisjointBoxLayout& levelGrids = a_phi.getBoxes();
    DataIterator dit = a_phi.dataIterator();

    LevelData<FArrayBox> temp(levelGrids, a_comps.size());
    Interval tempComps(0, a_comps.size() - 1);

    for (dit.reset(); dit.ok(); ++dit) {
        temp[dit()].copy(temp[dit()].box(), tempComps,
                         temp[dit()].box(), a_phi[dit()], a_comps);

        if (a_finerGridsPtr != NULL) {
            LayoutIterator litFine = a_finerGridsPtr->layoutIterator();

            // now loop over fine boxes and set covered regions to 0
            for (litFine.reset(); litFine.ok(); ++litFine) {
                Box coveredBox(a_finerGridsPtr->get(litFine()));
                coveredBox.coarsen(a_nRefFine);
                coveredBox &= temp[dit()].box();

                if (!coveredBox.isEmpty()) {
                    temp[dit()].setVal(levelMin, coveredBox, 0, tempComps.size());
                }
            }  // end loop over fine-grid boxes
        } // end if there is a finer level

        // while we're looping over the grids, get min as well
        // need to loop over comps
        for (int comp = tempComps.begin(); comp <= tempComps.end(); ++comp) {
            thisMin = temp[dit()].min(comp);
            if (thisMin < levelMin) {
                levelMin = thisMin;
            }
        }
    } // end loop over this level's grids

#ifdef CH_MPI
    Real recv;
    int result = MPI_Allreduce(&levelMin, &recv, 1, MPI_CH_REAL,
                               MPI_MIN, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMin");
    }

    levelMin = recv;
#endif

    return levelMin;
}
