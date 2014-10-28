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
#ifndef AMRPressureSolver_H__INCLUDED__
#define AMRPressureSolver_H__INCLUDED__

#include "EllipticBCInterface.H"
#include "AMREllipticSolver.H"
#include "LinearSolver.H"
#include "LevelGeometry.H"


// -----------------------------------------------------------------------------
// Encapsulates any solver that will do the projecting.
//
// This class serves the following purposes:
// 1. Makes switching solvers used by the projectors easy.
// 2. Simplifies the complex AMR solver construction process.
// 3. Eliminates messy code duplication in the projectors.
//    (No more 4000 line file!)
// -----------------------------------------------------------------------------
class AMRPressureSolver
{
public:
    // Default constructor -- leaves object unusable.
    AMRPressureSolver ();

    // Destructor
    virtual ~AMRPressureSolver ();

    // Overrides the default AMRMG settings.
    // This can only be called _before_ define.
    virtual void setAMRMGParameters (const int  a_AMRMG_imin,
                                     const int  a_AMRMG_imax,
                                     const Real a_AMRMG_eps,
                                     const int  a_AMRMG_maxDepth,
                                     const int  a_AMRMG_num_precond_iters,
                                     const int  a_AMRMG_num_smooth_down,
                                     const int  a_AMRMG_num_smooth_up,
                                     const int  a_AMRMG_num_smooth_bottom,
                                     const int  a_AMRMG_precondMode,
                                     const int  a_AMRMG_relaxMode,
                                     const int  a_AMRMG_numMG,
                                     const Real a_AMRMG_hang,
                                     const Real a_AMRMG_norm_thresh,
                                     const int  a_AMRMG_verbosity);

    // Overrides the default bottom solver settings.
    // This can only be called _before_ define.
    virtual void setBottomParameters (const Real a_bottom_imax,
                                      const int  a_bottom_numRestarts,
                                      const Real a_bottom_eps,
                                      const Real a_bottom_reps,
                                      const Real a_bottom_hang,
                                      const Real a_bottom_small,
                                      const int  a_bottom_normType,
                                      const int  a_bottom_verbosity);

    // Allocates everything and leaves object useable.
    // This will not erase the solver parameters.
    virtual void levelDefine (BCMethodHolder           a_bc,
                              const LevelGeometry&     a_levGeo,
                              const int                a_numLevels,
                              const FillJgupInterface* a_customFillJgupPtr = NULL);

    // Allocates everything and leaves object useable.
    // This will not erase the solver parameters.
    virtual void define (BCMethodHolder           a_bc,
                         const LevelGeometry&     a_levGeo,
                         const Box&               a_lminDomBox,
                         const int                a_numLevels,
                         const FillJgupInterface* a_customFillJgupPtr = NULL);

    // Frees memory and leaves object unuseable.
    // This will not erase the solver parameters.
    virtual void undefine ();

    // Solves the Poisson problem.
    virtual void solve (Vector<LevelData<FArrayBox>*>&       a_phi,
                        const Vector<LevelData<FArrayBox>*>& a_rhs,
                        const int                            a_lmin,
                        const int                            a_lmax,
                        const bool                           a_zeroPhi = true,
                        const bool                           a_forceHomogeneous = false);

    // Solves the Poisson problem.
    virtual void levelSolve (LevelData<FArrayBox>&       a_phi,
                             const LevelData<FArrayBox>* a_crsePhiPtr,
                             const LevelData<FArrayBox>& a_rhs,
                             const bool                  a_zeroPhi = true,
                             const bool                  a_forceHomogeneous = false);

    // Is this object in a useable state?
    virtual inline bool isDefined () const;

    // Even if this is a level solver, we may need to apply CF-BCs.
    virtual inline int getNumLevels () const;

    // This collects levgeos and puts them into the correct vector index.
    static void gatherAMRLevGeos (Vector<const LevelGeometry*>& a_amrLevGeos,
                                  const LevelGeometry&          a_levGeo,
                                  const Box&                    a_lminDomBox,
                                  const int                     a_lmin,
                                  const int                     a_lmax);

protected:
    static bool s_useLevelLepticSolver;
    static bool s_useLevelMGSolver;

    static bool s_useAMRLepticSolver;
    static bool s_useAMRMGSolver;

    bool m_isDefined;
    bool m_isLevelSolver;
    int  m_numLevels; // Even if this is a level solver, we may need to apply CF-BCs.

    // AMRMG settings
    int  m_AMRMG_imin;
    int  m_AMRMG_imax;
    Real m_AMRMG_eps;
    int  m_AMRMG_maxDepth;
    int  m_AMRMG_num_precond_iters;
    int  m_AMRMG_num_smooth_down;
    int  m_AMRMG_num_smooth_up;
    int  m_AMRMG_num_smooth_bottom;
    int  m_AMRMG_precondMode;
    int  m_AMRMG_relaxMode;
    int  m_AMRMG_numMG;
    Real m_AMRMG_hang;
    Real m_AMRMG_norm_thresh;
    int  m_AMRMG_verbosity;

    // Bottom solver settings
    Real m_bottom_imax;
    int  m_bottom_numRestarts;
    Real m_bottom_reps;
    Real m_bottom_eps;
    Real m_bottom_hang;
    Real m_bottom_small;
    int  m_bottom_normType;
    int  m_bottom_verbosity;

    // The solvers
    AMREllipticSolver<LevelData<FArrayBox> >*  m_lepticSolverPtr;
    AMREllipticSolver<LevelData<FArrayBox> >*  m_amrmgSolverPtr;
    LinearSolver<LevelData<FArrayBox> >*       m_bottomSolverPtr;

private:
    // Copy and assignment not allowed
    AMRPressureSolver (const AMRPressureSolver&);
    void operator= (const AMRPressureSolver&);
};


// -----------------------------------------------------------------------------
// Is this object in a useable state?
// -----------------------------------------------------------------------------
inline bool AMRPressureSolver::isDefined () const
{
    return m_isDefined;
}


// -----------------------------------------------------------------------------
// Even if this is a level solver, we may need to apply CF-BCs.
// -----------------------------------------------------------------------------
inline int AMRPressureSolver::getNumLevels () const
{
    CH_assert(isDefined());
    return m_numLevels;
}



#endif //!AMRPressureSolver_H__INCLUDED__
