#include "computeMappedSum.H"
#include "computeMappedSumF_F.H"
#include "LevelGeometry.H"
#include "EllipticBCInterface.H"


// FArrayBox versions...

// -----------------------------------------------------------------------------
// Create a grid mask. 0 = covered by finer grid, 1 = valid data.
// WARNING: a_maskFAB.box() must equal grids[dit]. No ghosts!!!
// -----------------------------------------------------------------------------
void createGridMask (BaseFab<int>&            a_maskFAB,
                     const DisjointBoxLayout* a_finerGridsPtr,
                     const IntVect&           a_refRatio)
{
    a_maskFAB.setVal(1);

    if (a_finerGridsPtr != NULL) {
        LayoutIterator litFine = a_finerGridsPtr->layoutIterator();
        for (litFine.reset(); litFine.ok(); ++litFine) {
            // Calculate covered region
            Box coveredBox(a_finerGridsPtr->get(litFine()));
            coveredBox.coarsen(a_refRatio);
            coveredBox.convert(a_maskFAB.box().type());
            coveredBox &= a_maskFAB.box();

            // Set mask to zero on covered region
            if (!coveredBox.isEmpty()) {
                a_maskFAB.setVal(0, coveredBox, 0, 1);
            }
        }
    } // end if there is a finer level
}


// -----------------------------------------------------------------------------
// This is a version of computeMappedSum that works without a LevelGeometry
// object. This is useful in the MappedAMRPoissonOp, where only the metric
// is cached. This function performs data exchanges.
// NOTE: This asks for Jinv, not J!
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real computeMappedSum (Real&                       a_vol,
                       const LevelData<FArrayBox>& a_phi,
                       const DisjointBoxLayout*    a_finerGridsPtr,
                       const IntVect&              a_fineRefRatio,
                       const RealVect&             a_dx,
                       const LevelData<FArrayBox>& a_CCJinv,
                       const int                   a_comp)
{
    // Sanity checks
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < a_phi.nComp());
    CH_assert(a_CCJinv.getBoxes().compatible(a_phi.getBoxes()));

    // Gather geometric structures
    const Real dxScale = a_dx.product();
    const DisjointBoxLayout& grids = a_phi.getBoxes();
    DataIterator dit = grids.dataIterator();

    // Loop over this level's grids, adding to sum.
    Real localSum = 0.0;
    Real localVol = 0.0;
    for (dit.reset(); dit.ok(); ++dit) {
        // Create references for convenience / allocate calucation space.
        const FArrayBox& phiFAB = a_phi[dit];
        const FArrayBox& JinvFAB = a_CCJinv[dit];
        const Box& valid = grids[dit];

        // If needed, we can adopt other centerings later.
        CH_assert(phiFAB.box().type() == IntVect::Zero);

        // Phi should be equal to or larger than the valid region.
        CH_assert(phiFAB.box().contains(valid));

        // Create a grid mask. 0 = covered by finer grid, 1 = valid data.
        BaseFab<int> maskFAB(valid, 1);
        createGridMask(maskFAB, a_finerGridsPtr, a_fineRefRatio);

        // Add to localSum
        FORT_COMPUTEMAPPEDSUMINV(
            CHF_REAL(localSum),
            CHF_REAL(localVol),
            CHF_CONST_FRA1(phiFAB, a_comp),
            CHF_CONST_FRA1(JinvFAB,0),
            CHF_CONST_FIA1(maskFAB,0),
            CHF_BOX(valid),
            CHF_CONST_REAL(dxScale));
    } // end loop over this level's grids

    // Compute global sum (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalSum = 0.0;
    int result = MPI_Allreduce(&localSum, &globalSum, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedSum");
    }

    Real globalVol = 0.0;
    result = MPI_Allreduce(&localVol, &globalVol, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedSum");
    }

#else
    Real globalSum = localSum;
    Real globalVol = localVol;
#endif

    a_vol += globalVol;
    return globalSum;
}


// -----------------------------------------------------------------------------
// Private function.
// Performs most of the computation for the sum functions. This version is not
// public and does not perform MPI communication.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
static Real computeLocalMappedSum (Real&                       a_vol,
                                   const LevelData<FArrayBox>& a_phi,
                                   const DisjointBoxLayout*    a_finerGridsPtr,
                                   const LevelGeometry&        a_levGeo,
                                   const int                   a_comp)
{
    // Sanity checks
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < a_phi.nComp());
    CH_assert(a_levGeo.getBoxes().compatible(a_phi.getBoxes()));
    if (a_finerGridsPtr != NULL) {
        CH_assert(a_levGeo.getFinerPtr() != NULL);
        CH_assert(a_levGeo.getFinerPtr()->getBoxes().compatible(*a_finerGridsPtr));
    }

    // Gather geometric structures
    const Real dxScale = a_levGeo.getDx().product();
    const LevelData<FArrayBox>& JRef = a_levGeo.getCCJ();
    const IntVect& fineRefRatio = a_levGeo.getFineRefRatio();
    const DisjointBoxLayout& grids = a_phi.getBoxes();
    DataIterator dit = grids.dataIterator();

    // Loop over this level's grids, adding to sum.
    Real localSum = 0.0;
    for (dit.reset(); dit.ok(); ++dit) {
        // Create references for convenience / allocate calucation space.
        const FArrayBox& phiFAB = a_phi[dit];
        const FArrayBox& JFAB = JRef[dit];
        const Box& valid = grids[dit];

        // If needed, we can adopt other centerings later.
        CH_assert(phiFAB.box().type() == IntVect::Zero);

        // Phi should be equal to or larger than the valid region.
        CH_assert(phiFAB.box().contains(valid));

        // Create a grid mask. 0 = covered by finer grid, 1 = valid data.
        BaseFab<int> maskFAB(valid, 1);
        createGridMask(maskFAB, a_finerGridsPtr, fineRefRatio);

        // Add to localSum
        FORT_COMPUTEMAPPEDSUM(
            CHF_REAL(localSum),
            CHF_REAL(a_vol),
            CHF_CONST_FRA1(phiFAB,a_comp),
            CHF_CONST_FRA1(JFAB,0),
            CHF_CONST_FIA1(maskFAB,0),
            CHF_BOX(valid),
            CHF_CONST_REAL(dxScale));
    } // end loop over this level's grids

    return localSum;
}


// -----------------------------------------------------------------------------
// Returns the integral of phi over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real computeMappedSum (const Vector<LevelData<FArrayBox>*>& a_phi,
                       const LevelGeometry&                 a_levGeo,
                       const int                            a_comp,
                       const int                            a_lBase)
{
    Real vol;
    return computeMappedSum(vol, a_phi, a_levGeo, a_comp, a_lBase);
}


// -----------------------------------------------------------------------------
// Returns the integral of phi over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real computeMappedSum (Real&                                a_vol,
                       const Vector<LevelData<FArrayBox>*>& a_phi,
                       const LevelGeometry&                 a_levGeo,
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

    // Loop over levels and add to total sum.
    Real localSum = 0.0;
    Real localVol = 0.0;
    for (lev = a_lBase; lev <= topLevel; ++lev) {
        const LevelData<FArrayBox>& levelPhi = *a_phi[lev];

        const DisjointBoxLayout* finerGridsPtr = NULL;
        if (lev < topLevel) {
            finerGridsPtr = &(a_phi[lev+1]->getBoxes());
        }

        const LevelGeometry& levGeoRef = *vLevGeos[lev];
        CH_assert(levGeoRef.getBoxes() == levelPhi.getBoxes()); // This may be overdoing it.

        // Add sum to running total.
        localSum += computeLocalMappedSum(localVol,
                                          levelPhi,
                                          finerGridsPtr,
                                          levGeoRef,
                                          a_comp);
    }

    // Compute global sum (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalSum = 0.0;
    int result = MPI_Allreduce(&localSum, &globalSum, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedSum");
    }

    Real globalVol = 0.0;
    result = MPI_Allreduce(&localVol, &globalVol, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedSum");
    }

#else
    Real globalSum = localSum;
    Real globalVol = localVol;
#endif

    a_vol += globalVol;
    return globalSum;
}


// -----------------------------------------------------------------------------
// Returns the integral of phi over the valid region.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real computeMappedSum (const LevelData<FArrayBox>& a_phi,
                       const DisjointBoxLayout*    a_finerGridsPtr,
                       const LevelGeometry&        a_levGeo,
                       const int                   a_comp)
{
    Real vol;
    return computeMappedSum(vol, a_phi, a_finerGridsPtr, a_levGeo, a_comp);
}


// -----------------------------------------------------------------------------
// Returns the integral of phi over the valid region.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real computeMappedSum (Real&                       a_vol,
                       const LevelData<FArrayBox>& a_phi,
                       const DisjointBoxLayout*    a_finerGridsPtr,
                       const LevelGeometry&        a_levGeo,
                       const int                   a_comp)
{
    // Compute the sum.
    // Assertions are done in this function call.
    Real localVol = 0.0;
    Real localSum = computeLocalMappedSum(localVol,
                                          a_phi,
                                          a_finerGridsPtr,
                                          a_levGeo,
                                          a_comp);

    // Compute global sum (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalSum = 0.0;
    int result = MPI_Allreduce(&localSum, &globalSum, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedSum");
    }

    Real globalVol = 0.0;
    result = MPI_Allreduce(&localVol, &globalVol, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedSum");
    }

#else
    Real globalSum = localSum;
    Real globalVol = localVol;
#endif

    a_vol += globalVol;
    return globalSum;
}


// -----------------------------------------------------------------------------
// This is a version of computeUnmappedSum that works without a LevelGeometry
// object. This is useful in the MappedAMRPoissonOp, where only the metric
// is cached. This function performs data exchanges.
// NOTE: This asks for Jinv, not J!
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real computeUnmappedSum (Real&                       a_vol,
                         const LevelData<FArrayBox>& a_phi,
                         const DisjointBoxLayout*    a_finerGridsPtr,
                         const IntVect&              a_fineRefRatio,
                         const RealVect&             a_dx,
                         const LevelData<FArrayBox>& a_CCJinv,
                         const int                   a_comp)
{
    // Sanity checks
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < a_phi.nComp());
    CH_assert(a_CCJinv.getBoxes().compatible(a_phi.getBoxes()));

    // Gather geometric structures
    const Real dxScale = a_dx.product();
    const DisjointBoxLayout& grids = a_phi.getBoxes();
    DataIterator dit = grids.dataIterator();

    // Loop over this level's grids, adding to sum.
    Real localSum = 0.0;
    Real localVol = 0.0;
    for (dit.reset(); dit.ok(); ++dit) {
        // Create references for convenience / allocate calucation space.
        const FArrayBox& phiFAB = a_phi[dit];
        const FArrayBox& JinvFAB = a_CCJinv[dit];
        const Box& valid = grids[dit];

        // If needed, we can adopt other centerings later.
        CH_assert(phiFAB.box().type() == IntVect::Zero);

        // Phi should be equal to or larger than the valid region.
        CH_assert(phiFAB.box().contains(valid));

        // Create a grid mask. 0 = covered by finer grid, 1 = valid data.
        BaseFab<int> maskFAB(valid, 1);
        createGridMask(maskFAB, a_finerGridsPtr, a_fineRefRatio);

        // Add to localSum
        FORT_COMPUTEUNMAPPEDSUMINV(
            CHF_REAL(localSum),
            CHF_REAL(localVol),
            CHF_CONST_FRA1(phiFAB,a_comp),
            CHF_CONST_FRA1(JinvFAB,0),
            CHF_CONST_FIA1(maskFAB,0),
            CHF_BOX(valid),
            CHF_CONST_REAL(dxScale));
    } // end loop over this level's grids

    // Compute global sum (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalSum = 0.0;
    int result = MPI_Allreduce(&localSum, &globalSum, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedSum");
    }

    Real globalVol = 0.0;
    result = MPI_Allreduce(&localVol, &globalVol, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedSum");
    }

#else
    Real globalSum = localSum;
    Real globalVol = localVol;
#endif

    a_vol += globalVol;
    return globalSum;
}


// -----------------------------------------------------------------------------
// Private function.
// Performs most of the computation for the sum functions. This version is not
// public and does not perform MPI communication.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
static Real computeLocalUnmappedSum (Real&                       a_vol,
                                     const LevelData<FArrayBox>& a_phi,
                                     const DisjointBoxLayout*    a_finerGridsPtr,
                                     const LevelGeometry&        a_levGeo,
                                     const int                   a_comp)
{
    // Sanity checks
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < a_phi.nComp());
    CH_assert(a_levGeo.getBoxes().compatible(a_phi.getBoxes()));

    // Gather geometric structures
    const Real dxScale = a_levGeo.getDx().product();
    const LevelData<FArrayBox>& JRef = a_levGeo.getCCJ();
    const IntVect& fineRefRatio = a_levGeo.getFineRefRatio();
    const DisjointBoxLayout& grids = a_phi.getBoxes();
    DataIterator dit = grids.dataIterator();

    // Loop over this level's grids, adding to sum.
    Real localSum = 0.0;
    for (dit.reset(); dit.ok(); ++dit) {
        // Create references for convenience / allocate calucation space.
        const FArrayBox& phiFAB = a_phi[dit];
        const FArrayBox& JFAB = JRef[dit];
        const Box& valid = grids[dit];

        // If needed, we can adopt other centerings later.
        CH_assert(phiFAB.box().type() == IntVect::Zero);

        // Phi should be equal to or larger than the valid region.
        CH_assert(phiFAB.box().contains(valid));

        // Create a grid mask. 0 = covered by finer grid, 1 = valid data.
        BaseFab<int> maskFAB(valid, 1);
        createGridMask(maskFAB, a_finerGridsPtr, fineRefRatio);

        // Add to localSum
        FORT_COMPUTEUNMAPPEDSUM(
            CHF_REAL(localSum),
            CHF_REAL(a_vol),
            CHF_CONST_FRA1(phiFAB,a_comp),
            CHF_CONST_FRA1(JFAB,0),
            CHF_CONST_FIA1(maskFAB,0),
            CHF_BOX(valid),
            CHF_CONST_REAL(dxScale));
    } // end loop over this level's grids

    return localSum;
}


// -----------------------------------------------------------------------------
// Returns the integral of phi over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real computeUnmappedSum (const Vector<LevelData<FArrayBox>*>& a_phi,
                         const LevelGeometry&                 a_levGeo,
                         const int                            a_comp,
                         const int                            a_lBase)
{
    Real vol;
    return computeUnmappedSum(vol, a_phi, a_levGeo, a_comp, a_lBase);
}


// -----------------------------------------------------------------------------
// Returns the integral of phi over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real computeUnmappedSum (Real&                                a_vol,
                         const Vector<LevelData<FArrayBox>*>& a_phi,
                         const LevelGeometry&                 a_levGeo,
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

    // Loop over levels and add to total sum.
    Real localSum = 0.0;
    Real localVol = 0.0;
    for (lev = a_lBase; lev <= topLevel; ++lev) {
        const LevelData<FArrayBox>& levelPhi = *a_phi[lev];

        const DisjointBoxLayout* finerGridsPtr = NULL;
        if (lev < topLevel) {
            finerGridsPtr = &(a_phi[lev+1]->getBoxes());
        }

        const LevelGeometry& levGeoRef = *vLevGeos[lev];

        // Add sum to running total.
        localSum += computeLocalUnmappedSum(localVol,
                                            levelPhi,
                                            finerGridsPtr,
                                            levGeoRef,
                                            a_comp);
    }

    // Compute global sum (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalSum = 0.0;
    int result = MPI_Allreduce(&localSum, &globalSum, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedSum");
    }

    Real globalVol = 0.0;
    result = MPI_Allreduce(&localVol, &globalVol, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedSum");
    }

#else
    Real globalSum = localSum;
    Real globalVol = localVol;
#endif

    a_vol += globalVol;
    return globalSum;
}


// -----------------------------------------------------------------------------
// Returns the integral of phi over the valid region.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real computeUnmappedSum (const LevelData<FArrayBox>& a_phi,
                         const DisjointBoxLayout*    a_finerGridsPtr,
                         const LevelGeometry&        a_levGeo,
                         const int                   a_comp)
{
    Real vol;
    return computeUnmappedSum(vol, a_phi, a_finerGridsPtr, a_levGeo, a_comp);
}


// -----------------------------------------------------------------------------
// Returns the integral of phi over the valid region.
// LIMITATIONS: This function can only handle cell-centered data.
// -----------------------------------------------------------------------------
Real computeUnmappedSum (Real&                       a_vol,
                         const LevelData<FArrayBox>& a_phi,
                         const DisjointBoxLayout*    a_finerGridsPtr,
                         const LevelGeometry&        a_levGeo,
                         const int                   a_comp)
{
    // Compute the sum.
    // Assertions are done in this function call.
    Real localVol = 0.0;
    Real localSum = computeLocalUnmappedSum(localVol,
                                            a_phi,
                                            a_finerGridsPtr,
                                            a_levGeo,
                                            a_comp);

    // Compute global sum (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalSum = 0.0;
    int result = MPI_Allreduce(&localSum, &globalSum, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedSum");
    }

    Real globalVol = 0.0;
    result = MPI_Allreduce(&localVol, &globalVol, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedSum");
    }

#else
    Real globalSum = localSum;
    Real globalVol = localVol;
#endif

    a_vol += globalVol;
    return globalSum;
}


// FluxBox versions (boundary integrals only)...

// -----------------------------------------------------------------------------
// Private function.
// Performs most of the computation for the bdry sum functions. This version is
// not public and does not perform MPI communication.
// -----------------------------------------------------------------------------
static Real computeLocalMappedBdrySum (Real&                     a_area,
                                       const LevelData<FluxBox>& a_flux,
                                       const DisjointBoxLayout*  a_finerGridsPtr,
                                       const LevelGeometry&      a_levGeo,
                                       const int                 a_comp)
{
    // Sanity checks
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < a_flux.nComp());
    CH_assert(a_levGeo.getBoxes().compatible(a_flux.getBoxes()));

    // Gather geometric structures
    const RealVect dx = a_levGeo.getDx();
    const RealVect dArea = dx.product() / dx;
    const IntVect& fineRefRatio = a_levGeo.getFineRefRatio();
    const ProblemDomain& domain = a_levGeo.getDomain();
    const Box& domBox = domain.domainBox();
    const DisjointBoxLayout& grids = a_flux.getBoxes();
    DataIterator dit = grids.dataIterator();

    // Loop over this level's grids, adding to the integral.
    Real localSum = 0.0;
    for (dit.reset(); dit.ok(); ++dit) {
        const Box& valid = grids[dit];

        // Loop over boundary faces. Do not include periodic directions.
        for (int dir = 0; dir < SpaceDim; ++dir) {
            if (domain.isPeriodic(dir)) continue;

            const Box domFaceBox = surroundingNodes(domBox, dir);
            const Real dAreaDir = dArea[dir];

            SideIterator sit;
            for (sit.reset(); sit.ok(); ++sit) {
                // Is this face at the domain boundary?
                Box faceBox = bdryBox(valid, dir, sit(), 1);
                if ((faceBox & domFaceBox).isEmpty()) continue;

                // Create references / allocate work space
                const FArrayBox& fluxFAB = a_flux[dit][dir];
                // FArrayBox tempFAB(faceBox, 1);

                // Sanity checks
                CH_assert(domFaceBox.contains(faceBox));
                CH_assert(fluxFAB.box().contains(faceBox));

                BaseFab<int> maskFAB(faceBox, 1);
                createGridMask(maskFAB, a_finerGridsPtr, fineRefRatio);

                // Get a face-centered J.
                FArrayBox JFAB(faceBox, 1);
                a_levGeo.fill_J(JFAB);

                // Add to localSum
                FORT_COMPUTEMAPPEDFLUXSUM(
                    CHF_REAL(localSum),
                    CHF_REAL(a_area),
                    CHF_CONST_FRA1(fluxFAB,a_comp),
                    CHF_CONST_FRA1(JFAB,0),
                    CHF_CONST_FIA1(maskFAB,0),
                    CHF_BOX(faceBox),
                    CHF_CONST_REAL(dAreaDir));

            } // end loop over boundary sides (sit)
        } // end loop over boundary directions (dir)
    } // end loop over this level's grids

    return localSum;
}


// -----------------------------------------------------------------------------
// Returns the integral of a_flux around the boundary over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// -----------------------------------------------------------------------------
Real computeMappedBdrySum (Real&                              a_area,
                           const Vector<LevelData<FluxBox>*>& a_flux,
                           const LevelGeometry&               a_levGeo,
                           const int                          a_comp,
                           const int                          a_lBase)
{
    // Sanity check on a_lBase
    const int vectorSize = a_flux.size();
    CH_assert(0 <= a_lBase);
    CH_assert(a_lBase < vectorSize);

    // We need at least one defined level.
    CH_assert(a_flux[a_lBase] != NULL);

    // Find highest defined level.
    int lev;
    for (lev = a_lBase; lev < vectorSize; ++lev) {
        if (a_flux[lev] == NULL) break;
        if (!a_flux[lev]->isDefined()) break;
    }
    const int topLevel = lev - 1;

    // Get every level's levGeo.
    Vector<const LevelGeometry*> vLevGeos = a_levGeo.getAMRLevGeos();

    // Loop over levels and add to total sum.
    Real localSum = 0.0;
    Real localArea = 0.0;
    for (lev = a_lBase; lev <= topLevel; ++lev) {
        const LevelData<FluxBox>& levelFlux = *a_flux[lev];

        const DisjointBoxLayout* finerGridsPtr = NULL;
        if (lev < topLevel) {
            finerGridsPtr = &(a_flux[lev+1]->getBoxes());
        }

        const LevelGeometry& levGeoRef = *vLevGeos[lev];
        CH_assert(levGeoRef.getBoxes() == levelFlux.getBoxes()); // This may be overdoing it.

        // Add sum to running total.
        localSum += computeLocalMappedBdrySum(localArea,
                                              levelFlux,
                                              finerGridsPtr,
                                              levGeoRef,
                                              a_comp);
    }

    // Compute global sum (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalSum = 0.0;
    int result = MPI_Allreduce(&localSum, &globalSum, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedBdrySum");
    }

    Real globalArea = 0.0;
    result = MPI_Allreduce(&localArea, &globalArea, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedBdrySum");
    }

#else
    Real globalSum = localSum;
    Real globalArea = localArea;
#endif

    a_area += globalArea;
    return globalSum;
}


// -----------------------------------------------------------------------------
// Returns the integral of phi over the valid physical boundary.
// -----------------------------------------------------------------------------
Real computeMappedBdrySum (Real&                     a_area,
                           const LevelData<FluxBox>& a_flux,
                           const DisjointBoxLayout*  a_finerGridsPtr,
                           const LevelGeometry&      a_levGeo,
                           const int                 a_comp)
{
    // Compute the sum.
    // Assertions are done in this function call.
    Real localArea = 0.0;
    Real localSum = computeLocalMappedBdrySum(localArea,
                                              a_flux,
                                              a_finerGridsPtr,
                                              a_levGeo,
                                              a_comp);

    // Compute global sum (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalSum = 0.0;
    int result = MPI_Allreduce(&localSum, &globalSum, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedBdrySum");
    }

    Real globalArea = 0.0;
    result = MPI_Allreduce(&localArea, &globalArea, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedBdrySum");
    }

#else
    Real globalSum = localSum;
    Real globalArea = localArea;
#endif

    a_area += globalArea;
    return globalSum;
}


// -----------------------------------------------------------------------------
// Returns the integral of the results of a_bc.setFluxes() over the valid
// physical boundary.
// -----------------------------------------------------------------------------
Real computeMappedBdrySum (Real&                 a_area,
                           const BCMethodHolder& a_bc,
                           const Real            a_time,
                           const LevelGeometry&  a_levGeo,
                           const int             a_comp,
                           const int             a_lBase)
{
    // Get every level's levGeo.
    Vector<const LevelGeometry*> vLevGeos = a_levGeo.getAMRLevGeos();

    // Sanity check on a_lBase
    const int vectorSize = vLevGeos.size();
    CH_assert(0 <= a_lBase);
    CH_assert(a_lBase < vectorSize);

    // Find highest defined level.
    int lev;
    for (lev = a_lBase; lev < vectorSize; ++lev) {
        if (!vLevGeos[lev]->isDefined()) break;
        if (vLevGeos[lev]->getBoxes().size() == 0) break;
    }
    const int topLevel = lev - 1;

    // Loop over levels and add to total sum.
    Real localSum = 0.0;
    Real localArea = 0.0;
    for (lev = a_lBase; lev <= topLevel; ++lev) {
        // Get this level's geometry
        const LevelGeometry& levGeoRef = *vLevGeos[lev];
        const ProblemDomain& domain = levGeoRef.getDomain();
        const Box& domBox = domain.domainBox();
        const RealVect& dx = levGeoRef.getDx();
        const RealVect dArea = dx.product() / dx;
        const IntVect& fineRefRatio = levGeoRef.getFineRefRatio();
        const DisjointBoxLayout& grids = levGeoRef.getBoxes();
        DataIterator dit = grids.dataIterator();

        // Get the finer grids and zero out cells covered by finer data.
        const DisjointBoxLayout* finerGridsPtr = NULL;
        if (lev < topLevel) {
            finerGridsPtr = &(vLevGeos[lev+1]->getBoxes());
        }

        for (dit.reset(); dit.ok(); ++dit) {
            const Box& valid = grids[dit];

            for (int dir = 0; dir < SpaceDim; ++dir) {
                if (domain.isPeriodic(dir)) continue;

                const Box domFaceBox = surroundingNodes(domBox, dir);
                const Real dAreaDir = dArea[dir];

                SideIterator sit;
                for (sit.reset(); sit.ok(); ++sit) {
                    // Is this face at the domain boundary?
                    Box faceBox = bdryBox(valid, dir, sit(), 1);
                    if ((faceBox & domFaceBox).isEmpty()) continue;

                    // Create calculation space and fill boundary values.
                    FArrayBox fluxFAB(faceBox, 1);
                    fluxFAB.setVal(1.0e300);

                    const FluxBox* FCJgupPtr = &(levGeoRef.getFCJgup()[dit]);
                    a_bc.setFluxes(fluxFAB,
                                   NULL,
                                   valid,
                                   domain,
                                   dx,
                                   dit(),
                                   FCJgupPtr,
                                   dir,
                                   false,
                                   a_time,
                                   Interval(a_comp, a_comp));

                    BaseFab<int> maskFAB(faceBox, 1);
                    createGridMask(maskFAB, finerGridsPtr, fineRefRatio);

                    // Get a face-centered J.
                    FArrayBox JFAB(faceBox, 1);
                    levGeoRef.fill_J(JFAB);

                    // Add to localSum
                    FORT_COMPUTEMAPPEDFLUXSUM(
                        CHF_REAL(localSum),
                        CHF_REAL(localArea),
                        CHF_CONST_FRA1(fluxFAB,0),
                        CHF_CONST_FRA1(JFAB,0),
                        CHF_CONST_FIA1(maskFAB,0),
                        CHF_BOX(faceBox),
                        CHF_CONST_REAL(dAreaDir));

                } // end loop over boundary sides (sit)
            } // end loop over boundary directions (dir)
        } // end loop over this level's grids (dit)
    } // end loop over levels (lev)

    // Compute global sum (this is where the MPI communication happens)
#ifdef CH_MPI
    Real globalSum = 0.0;
    int result = MPI_Allreduce(&localSum, &globalSum, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedBdrySum");
    }

    Real globalArea = 0.0;
    result = MPI_Allreduce(&localArea, &globalArea, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedBdrySum");
    }

#else
    Real globalSum = localSum;
    Real globalArea = localArea;
#endif

    a_area += globalArea;
    return globalSum;
}

