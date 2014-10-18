#include "HomogeneousCFInterp.H"
#include "InterpF_F.H" // Only needed for 2 InterpHomo functions.


// -----------------------------------------------------------------------------
// Interpolate ghosts at CF interface using zeros on coarser grids.
// This version acts on a set of IntVects.
// -----------------------------------------------------------------------------
void interpOnIVSHomo (LevelData<FArrayBox>& a_phif,
                      const DataIndex&      a_index,
                      const int             a_dir,
                      const Side::LoHiSide  a_side,
                      const IntVectSet&     a_interpIVS,
                      const Real            a_fineDxDir,
                      const Real            a_crseDxDir)
{
    CH_TIME("interpOnIVSHomo");

    // Sanity checks
    CH_assert((a_dir >= 0) && (a_dir < SpaceDim));
    CH_assert(a_phif.ghostVect()[a_dir] >= 1);

    IVSIterator fine_ivsit(a_interpIVS);
    FArrayBox& a_phi = a_phif[a_index];
    const int isign = sign(a_side);

    if (a_phi.box().size(a_dir) == 3) {
        // Linear interpolation of fine ghosts assuming
        // all zeros on coarser level.

        // we are in a 1-wide box
        IntVect iv;
        Real pa;
        Real factor = 1.0 - 2.0 * a_fineDxDir / (a_fineDxDir + a_crseDxDir);
        for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit) {
            iv = fine_ivsit();
            iv[a_dir] -= isign;
            // Use linear interpolation
            for (int ivar = 0; ivar < a_phif.nComp(); ivar++) {
                pa = a_phi(iv, ivar);
                a_phi(fine_ivsit(), ivar) = factor * pa;
            }
        }
    } else {
        // Quadratic interpolation of fine ghosts assuming
        // all zeros on coarser level.

        // Symbolic reduced version of CF quadratic stencil
        Real pa, pb;
        Real c1 = 2.0*(a_crseDxDir-a_fineDxDir)/(a_crseDxDir+    a_fineDxDir); //first inside point
        Real c2 =    -(a_crseDxDir-a_fineDxDir)/(a_crseDxDir+3.0*a_fineDxDir); // next point inward
        IntVect ivf;
        for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit) {
            ivf = fine_ivsit();
            // Use quadratic interpolation
            for (int ivar = 0; ivar < a_phif.nComp(); ++ivar) {
                ivf[a_dir]-=2*isign;
                pa = a_phi(ivf, ivar);
                ivf[a_dir]+=isign;
                pb = a_phi(ivf, ivar);

                ivf[a_dir]+=isign;
                a_phi(fine_ivsit(), ivar) = c1*pb + c2*pa;
            } //end loop over components
        } //end loop over fine intvects
    }
}


// -----------------------------------------------------------------------------
// Interpolate ghosts at CF interface using zeros on coarser grids.
// Only do one face.
// -----------------------------------------------------------------------------
void homogeneousCFInterp (LevelData<FArrayBox>& a_phif,
                          const DataIndex&      a_index,
                          const int             a_dir,
                          const Side::LoHiSide  a_side,
                          const Real            a_fineDxDir,
                          const Real            a_crseDxDir,
                          const CFRegion&       a_cfRegion)
{
    CH_TIME("homogeneousCFInterp (single grid)");

    // Sanity checks
    CH_assert((a_dir >= 0) && (a_dir < SpaceDim));
    CH_assert(a_phif.ghostVect()[a_dir] >= 1);

    // Get the C/F region
    const CFIVS* cfivs_ptr = NULL;
    if (a_side == Side::Lo) {
        CFRegion& castCFRegion = (CFRegion&)a_cfRegion;
        cfivs_ptr = &(castCFRegion.loCFIVS(a_index, a_dir));
    } else {
        CFRegion& castCFRegion = (CFRegion&)a_cfRegion;
        cfivs_ptr = &(castCFRegion.hiCFIVS(a_index, a_dir));
    }

    if (cfivs_ptr->isPacked()) {
        // The C/F region is a box.
        // Iterate over the box and interpolate.
        int ihiorlo = sign(a_side);
        FArrayBox& phiFab = a_phif[a_index];
        const Box region = cfivs_ptr->packedBox();

        CH_assert(a_crseDxDir > 0.0 || region.isEmpty());

        if (phiFab.box().size(a_dir) == 3) {
            // Linear interpolation of fine ghosts assuming
            // all zeros on coarser level.
            FORTNT_INTERPHOMOLINEAR(
                CHF_FRA(phiFab),
                CHF_BOX(region),
                CHF_CONST_REAL(a_fineDxDir),
                CHF_CONST_REAL(a_crseDxDir),
                CHF_CONST_INT(a_dir),
                CHF_CONST_INT(ihiorlo));
        } else {
            // Quadratic interpolation of fine ghosts assuming
            // all zeros on coarser level.
            FORTNT_INTERPHOMO(
                CHF_FRA(phiFab),
                CHF_BOX(region),
                CHF_CONST_REAL(a_fineDxDir),
                CHF_CONST_REAL(a_crseDxDir),
                CHF_CONST_INT(a_dir),
                CHF_CONST_INT(ihiorlo));
        }
    } else {
        // The C/F region is sparse.
        // Iterate over the IVS and interpolate.
        const IntVectSet& interp_ivs = cfivs_ptr->getFineIVS();

        if (!interp_ivs.isEmpty()) {
            CH_assert(a_crseDxDir > 0.0);
            interpOnIVSHomo(a_phif,
                            a_index,
                            a_dir,
                            a_side,
                            interp_ivs,
                            a_fineDxDir,
                            a_crseDxDir);
        }
    }
}


// -----------------------------------------------------------------------------
// Interpolate ghosts at CF interface using zeros on coarser grids.
// -----------------------------------------------------------------------------
void homogeneousCFInterp (LevelData<FArrayBox>& a_phif,
                          const RealVect&       a_fineDx,
                          const RealVect&       a_crseDx,
                          const CFRegion&       a_cfRegion,
                          const IntVect&        a_applyDirs)
{
    CH_TIME("homogeneousCFInterp (full level)");

    // Loop over grids, directions, and sides and call the worker function.
    DataIterator dit = a_phif.dataIterator();
    for (dit.begin(); dit.ok(); ++dit) {
        if (a_phif[dit].box().isEmpty()) continue;

        for (int dir = 0; dir < SpaceDim; ++dir) {
            if (a_applyDirs[dir] == 0) continue;

            SideIterator sit;
            for (sit.begin(); sit.ok(); sit.next()) {
                homogeneousCFInterp(a_phif,
                                    dit(),
                                    dir,
                                    sit(),
                                    a_fineDx[dir],
                                    a_crseDx[dir],
                                    a_cfRegion);
            }
        }
    }
}

