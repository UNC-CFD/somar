/*******************************************************************************
 *    FILE: computeMappedSum.H
 *    DESC: Computes integrals in mapped spaces with anisotropic refinement.
 *    DATE: Sat 08 Jun 2013 08:55:38 PM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __computeMappedSum_H__INCLUDED__
#define __computeMappedSum_H__INCLUDED__

#include "LevelData.H"
class RealVect;
class LevelGeometry;
class BCMethodHolder;


// FArrayBox versions...

// This is a version of computeMappedSum that works without a LevelGeometry
// object. This is useful in the MappedAMRPoissonOp, where only the metric
// is cached. This function performs data exchanges.
// NOTE: This asks for Jinv, not J!
// LIMITATIONS: This function can only handle cell-centered data.
Real computeMappedSum (Real&                       a_vol,
                       const LevelData<FArrayBox>& a_phi,
                       const DisjointBoxLayout*    a_finerGridsPtr,
                       const IntVect&              a_fineRefRatio,
                       const RealVect&             a_dx,
                       const LevelData<FArrayBox>& a_CCJinv,
                       const int                   a_comp);

// Returns the integral of a_phi over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// LIMITATIONS: This function can only handle cell-centered data.
Real computeMappedSum (const Vector<LevelData<FArrayBox>*>& a_phi,
                       const LevelGeometry&                 a_levGeo,
                       const int                            a_comp = 0,
                       const int                            a_lBase = 0);

// Returns the integral of a_phi over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// LIMITATIONS: This function can only handle cell-centered data.
Real computeMappedSum (Real&                                a_vol,
                       const Vector<LevelData<FArrayBox>*>& a_phi,
                       const LevelGeometry&                 a_levGeo,
                       const int                            a_comp = 0,
                       const int                            a_lBase = 0);

// Returns the integral of phi over the valid region.
// LIMITATIONS: This function can only handle cell-centered data.
Real computeMappedSum (const LevelData<FArrayBox>& a_phi,
                       const DisjointBoxLayout*    a_finerGridsPtr,
                       const LevelGeometry&        a_levGeo,
                       const int                   a_comp = 0);

// Returns the integral of phi over the valid region.
// LIMITATIONS: This function can only handle cell-centered data.
Real computeMappedSum (Real&                       a_vol,
                       const LevelData<FArrayBox>& a_phi,
                       const DisjointBoxLayout*    a_finerGridsPtr,
                       const LevelGeometry&        a_levGeo,
                       const int                   a_comp = 0);

// This is a version of computeUnmappedSum that works without a LevelGeometry
// object. This is useful in the MappedAMRPoissonOp, where only the metric
// is cached. This function performs data exchanges.
// NOTE: This asks for Jinv, not J!
// LIMITATIONS: This function can only handle cell-centered data.
Real computeUnmappedSum (Real&                       a_vol,
                         const LevelData<FArrayBox>& a_phi,
                         const DisjointBoxLayout*    a_finerGridsPtr,
                         const IntVect&              a_fineRefRatio,
                         const RealVect&             a_dx,
                         const LevelData<FArrayBox>& a_CCJinv,
                         const int                   a_comp);

// Returns the integral of a_phi over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// LIMITATIONS: This function can only handle cell-centered data.
Real computeUnmappedSum (const Vector<LevelData<FArrayBox>*>& a_phi,
                         const LevelGeometry&                 a_levGeo,
                         const int                            a_comp = 0,
                         const int                            a_lBase = 0);

// Returns the integral of a_phi over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
// LIMITATIONS: This function can only handle cell-centered data.
Real computeUnmappedSum (Real&                                a_vol,
                         const Vector<LevelData<FArrayBox>*>& a_phi,
                         const LevelGeometry&                 a_levGeo,
                         const int                            a_comp = 0,
                         const int                            a_lBase = 0);

// Returns the integral of phi over the valid region.
// LIMITATIONS: This function can only handle cell-centered data.
Real computeUnmappedSum (const LevelData<FArrayBox>& a_phi,
                         const DisjointBoxLayout*    a_finerGridsPtr,
                         const LevelGeometry&        a_levGeo,
                         const int                   a_comp = 0);

// Returns the integral of phi over the valid region.
// LIMITATIONS: This function can only handle cell-centered data.
Real computeUnmappedSum (Real&                       a_vol,
                         const LevelData<FArrayBox>& a_phi,
                         const DisjointBoxLayout*    a_finerGridsPtr,
                         const LevelGeometry&        a_levGeo,
                         const int                   a_comp = 0);


// FluxBox versions (boundary integrals only)...

// Returns the integral of a_flux around the boundary over an AMR hierarchy.
// a_levGeo can be any levGeo in the hierarchy.
Real computeMappedBdrySum (Real&                              a_area,
                           const Vector<LevelData<FluxBox>*>& a_flux,
                           const LevelGeometry&               a_levGeo,
                           const int                          a_comp = 0,
                           const int                          a_lBase = 0);

// Returns the integral of a_flux over the valid physical boundary.
Real computeMappedBdrySum (Real&                     a_area,
                           const LevelData<FluxBox>& a_flux,
                           const DisjointBoxLayout*  a_finerGridsPtr,
                           const LevelGeometry&      a_levGeo,
                           const int                 a_comp = 0);

// Returns the integral of the results of a_bc.setFluxes() over the valid
// physical boundary.
Real computeMappedBdrySum (Real&                 a_area,
                           const BCMethodHolder& a_bc,
                           const Real            a_time,
                           const LevelGeometry&  a_levGeo,
                           const int             a_comp = 0,
                           const int             a_lBase = 0);


#endif //!__computeMappedSum_H__INCLUDED__


