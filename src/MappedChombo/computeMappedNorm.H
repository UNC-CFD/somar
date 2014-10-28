/*******************************************************************************
 *    FILE: computeMappedNorm.H
 *    DESC: Computes p-norms in mapped spaces with anisotropic refinement.
 *    DATE: Fri 07 Jun 2013 06:51:03 PM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __computeMappedNorm_H__INCLUDED__
#define __computeMappedNorm_H__INCLUDED__

#include "LevelData.H"
class LevelGeometry;
class RealVect;


// FArrayBox versions...

// Returns the volume-weighted norm of phi over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// LIMITATIONS: This function can only handle cell-centered data.
Real computeMappedNorm (const Vector<LevelData<FArrayBox>*>& a_phi,
                        const LevelGeometry&                 a_levGeo,
                        const int                            a_p = 2,
                        const int                            a_comp = 0,
                        const int                            a_lBase = 0);

// Returns the volume-weighted norm of phi over the valid region
// LIMITATIONS: This function can only handle cell-centered data.
Real computeMappedNorm (const LevelData<FArrayBox>& a_phi,
                        const DisjointBoxLayout*    a_finerGridsPtr,
                        const LevelGeometry&        a_levGeo,
                        const int                   a_p = 2,
                        const int                   a_comp = 0);


// Returns the norm of phi over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// LIMITATIONS: This function can only handle cell-centered data.
Real computeUnmappedNorm (const Vector<LevelData<FArrayBox>*>& a_phi,
                          const LevelGeometry&                 a_levGeo,
                          const int                            a_p = 2,
                          const int                            a_comp = 0,
                          const int                            a_lBase = 0);

// Returns the norm of phi over the valid region.
// LIMITATIONS: This function can only handle cell-centered data.
Real computeUnmappedNorm (const LevelData<FArrayBox>& a_phi,
                          const DisjointBoxLayout*    a_finerGridsPtr,
                          const LevelGeometry&        a_levGeo,
                          const int                   a_p = 2,
                          const int                   a_comp = 0);


// FluxBox versions...

// Returns the volume-weighted norm of phi over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// LIMITATIONS: This function can only handle face-centered data.
RealVect computeMappedNorm (const Vector<LevelData<FluxBox>*>& a_phi,
                            const LevelGeometry&               a_levGeo,
                            const int                          a_p = 2,
                            const int                          a_comp = 0,
                            const int                          a_lBase = 0);

// Returns the volume-weighted norm of phi over the valid region
// LIMITATIONS: This function can only handle face-centered data.
RealVect computeMappedNorm (const LevelData<FluxBox>& a_phi,
                            const DisjointBoxLayout*  a_finerGridsPtr,
                            const LevelGeometry&      a_levGeo,
                            const int                 a_p = 2,
                            const int                 a_comp = 0);


// Ghost cell versions (FArrayBox only)...

// This describes what ghosts to use in the norm calculations.
enum GhostNormMode {
    PHYS_GHOSTS = 0,
    INTERP_GHOSTS,
    ALL_GHOSTS
};

// Returns the volume-weighted norm of phi's ghosts over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// NOTE: This function can only handle cell-centered data.
// NOTE: This does not check for finer data coverings, nor should it since,
// technically, all ghosts are invalid data. All ghosts in the mode are used.
Real computeMappedNormGhost (const Vector<LevelData<FArrayBox>*>& a_phi,
                             const LevelGeometry&                 a_levGeo,
                             const GhostNormMode                  a_mode,
                             const int                            a_p = 2,
                             const int                            a_comp = 0,
                             const int                            a_lBase = 0);

// Returns the volume-weighted norm of phi over the first ghost layer.
// NOTE: This function can only handle cell-centered data.
// NOTE: This does not check for finer data coverings, nor should it since,
// technically, all ghosts are invalid data. All ghosts in the mode are used.
Real computeMappedNormGhost (const LevelData<FArrayBox>& a_phi,
                             const LevelGeometry&        a_levGeo,
                             const GhostNormMode         a_mode,
                             const int                   a_p = 2,
                             const int                   a_comp = 0);

// Returns the max value of phi over all valid regions.
Real computeMax(const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<IntVect>&                a_nRefFine,
                const Interval                        a_comps = Interval(0,0),
                const int                             a_lBase = 0);

// Returns the maximum value of phi over the valid region.
Real computeMax(const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGrids,
                const IntVect&              a_nRefFine,
                const Interval              a_comps = Interval(0,0));

// Returns the min value of phil over all valid regions.
Real computeMin(const Vector<LevelData<FArrayBox>* >& a_phi,
                const Vector<IntVect>&                a_nRefFine,
                const Interval                        a_comps = Interval(0,0),
                const int                             a_lBase = 0);

// Returns the minimum value of phi over the valid region.
Real computeMin(const LevelData<FArrayBox>& a_phi,
                const DisjointBoxLayout*    a_finerGrids,
                const IntVect&              a_nRefFine,
                const Interval              a_comps = Interval(0,0));


#endif //!__computeMappedNorm_H__INCLUDED__
