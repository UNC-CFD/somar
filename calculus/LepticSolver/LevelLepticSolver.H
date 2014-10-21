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
#ifndef __LevelLepticSolver_H__INCLUDED__
#define __LevelLepticSolver_H__INCLUDED__

#include "LevelData.H"
#include "FluxBox.H"
#include "LinearSolver.H"
#include "Tuple.H"
#include "CFRegion.H"
#include "CornerCopier.H"
#include "EllipticBCInterface.H"
#include "RealVect.H"
#include "MappedAMRPoissonOpFactory.H"
#include "BiCGStabSolver.H"
#include "BoundaryData.H"


// -----------------------------------------------------------------------------
// Solves L[phi] = Div[u*] on a single level using the leptic iterative method.
// See: http://dx.doi.org/10.1016/j.jcp.2011.06.022
//
// You are looking at the prototype of the leptic method. This solver has
// been working well for us, but the class itself is not 100% user-friendly.
// There are a lot of parameters to set -- many of which are hard coded. In the
// next version of this code, we hope to create a more versatile solver that is
// a bit more automatic. That is, if you just need to solve equations without
// getting under the hood, we should be able to provide a leptic solver with
// less knobs to turn. -ES
// -----------------------------------------------------------------------------
class LevelLepticSolver: public LinearSolver<LevelData<FArrayBox> >
{
    friend void test_AMRLepticSolver ();

public:
    // Constructor -- leaves object unusable.
    LevelLepticSolver ();

    // Destructor
    virtual ~LevelLepticSolver ();

    // Deallocates all memory and brings this solver back to an undefined state.
    virtual void undefine ();

    // Full define -- leave object in a useable state.
    // This is an override of the pure virtual LinearSolver function.
    //
    // The operator is needed to calculate residuals, to provide metric data,
    // to provide BC types, etc. This solver asks more of the op than the
    // LinearOp interface can provide. Therefore, the op MUST be able to be
    // cast into a LepticOperator. MappedAMRPoissonOp can do this. This function
    // asks for a LinearOp and not a LepticOperator so that it can be used as a
    // bottom solver.
    virtual void define (LinearOp<LevelData<FArrayBox> >* a_operator,
                         bool                             a_homogeneous = false);

    // Reset whether the solver is homogeneous.
    // This is an override of the pure virtual LinearSolver function.
    virtual void setHomogeneous (bool a_homogeneous);

    // Used for setting CFBCs. Set to NULL for homogeneous CFBCs.
    virtual void setCrsePhiPtr (const LevelData<FArrayBox>* a_crsePhiPtr);

    // Initializes all parameters with default values.
    virtual void setDefaultParameters ();

    // Sets the solver's overall parameters.
    // a_horizRhsTol is set in this function and not in one of the horiz param
    // functions because it is a statement that controls the overall performance
    // of the leptic solver, not the horizontal solver. If |horizontal rhs| is less
    // than a_horizRhsTol, then a horizontal solve will not be performed.
    virtual void setParameters (const int  a_maxOrder,
                                const Real a_eps,
                                const Real a_hang,
                                const int  a_normType,
                                const int  a_verbosity,
                                const Real a_horizRhsTol);

    // Sets the horizontal MG solver's parameters.
    virtual void setHorizMGParameters (const int  a_imin,
                                       const int  a_imax,
                                       const int  a_numSmoothDown,
                                       const int  a_numSmoothBottom,
                                       const int  a_numSmoothUp,
                                       const int  a_numSmoothPreCond,
                                       const int  a_relaxType,
                                       const int  a_maxDepth,
                                       const Real a_eps,
                                       const Real a_hang,
                                       const Real a_normThresh,
                                       const int  a_verbosity);

    // Sets the horizontal bottom solver's parameters.
    virtual void setHorizBottomParameters (const int  a_imax,
                                           const Real a_eps,
                                           const int  a_numRestarts,
                                           const Real a_hang,
                                           const int  a_verbosity);

    // Solve L[phi] = rhs.
    // This is an override of the pure virtual LinearSolver function.
    virtual void solve (LevelData<FArrayBox>&       a_phi,
                        const LevelData<FArrayBox>& a_rhs);

    // Set a convergence metric, along with solver tolerance, if desired.
    // This is an override of the non-pure virtual LinearSolver function.
    virtual void setConvergenceMetrics (Real a_metric,
                                        Real a_tolerance)
    {;}

    // Inform the user why we stopped solving.
    struct ExitStatus {
        enum {
            NONE     = -1,  // We haven't solved anything yet.
            CONVERGE =  0,  // Converged to tolerance
            ITER     =  1,  // Exited prematurely -- max iters reached.
            HANG     =  2,  // Exited prematurely -- convergence slowed or stopped.
            DIVERGE  =  3,  // Exited prematurely -- solver made some progress, then diverged.
            KABOOM   =  4   // Exited prematurely -- solver immediately diverged without making progress.
        };
    };
    inline virtual int exitStatusInt () const;
    virtual char* exitStatusStr () const;

protected:
    // Computing utilities...

    // Compute -\partial_m \bar{Jgup^{mj} \partial_j \phi}.
    // This does set BCs and exchanges phi. This does not subtract excess/H.
    // Careful! a_rhs will be defined over m_flatGrids, not m_horizGrids!
    virtual void computeHorizRHS (BoxLayoutData<FArrayBox>& a_rhs,
                                  LevelData<FArrayBox>&     a_phi,
                                  const BoundaryData<Real>& a_bdryData) const;

    // Computes scale*Jgup^{SpaceDim-1,m}*D[phi]/dx^[m], m=[0,SpaceDim-1).
    // Used to generate BCs for the vertical problems.
    // No ghosts need to be set. All ghosts will be extrapolated.
    // It only makes sense to use this function on Neum-Neum problems.
    // It will be up to you to add the excess if needed.
    virtual void levelVertHorizGradient (BoundaryData<Real>&   a_vertBdryData,
                                         LevelData<FArrayBox>& a_phi,
                                         const Real            a_scale) const;

    // Computes the excess function, excess = hiNeumBC - loNeumBC - Integral[rhs].
    virtual void computeVerticalExcess (BoxLayoutData<FArrayBox>&   a_excess,
                                        const LevelData<FArrayBox>& a_rhs,
                                        const BoundaryData<Real>&   a_bdryData) const;

    // Computes the vertical solutions.
    // vertRhs will be temporarily altered, but then restored.
    // No BCs will be set on vertPhi and this is a non-blocking function.
    virtual void verticalLineSolver (LevelData<FArrayBox>&            a_vertPhi,
                                     LevelData<FArrayBox>&            a_vertRhs,
                                     const LayoutData<Tuple<int,2> >& a_vertBCTypes,
                                     const BoundaryData<Real>&        a_bdryData) const;

    // Computes the horizontal solutions.
    virtual void horizontalSolver (LevelData<FArrayBox>&       a_phi,
                                   const LevelData<FArrayBox>& a_rhs);


    // Miscellaneous utilities...

    // Adds a vertical correction to the vertical solution.
    virtual void addVerticalCorrection (LevelData<FArrayBox>&       a_phiTotal,
                                        const LevelData<FArrayBox>& a_vertCor) const;

    // Adds a horizontal correction to the vertical solution.
    virtual void addHorizontalCorrection (LevelData<FArrayBox>&       a_vertPhi,
                                          const LevelData<FArrayBox>& a_horizCor) const;

    // Static utility
    // This fills ghost cells with the adjacent BCType.
    // If all BCTypes are Neum or periodic, a_doHorizSolve will be set to true.
    static void gatherVerticalBCTypes (LayoutData<Tuple<int,2> >& a_BCTypes,
                                       bool&                      a_doHorizSolve,
                                       const ProblemDomain&       a_domain,
                                       CFRegion&                  a_CFRegion,
                                       const BCMethodHolder&      a_bc);

    // Static utility
    // Sets all data referenced by the vector of data indices to zero.
    static void setToZero (BoxLayoutData<FArrayBox>& a_data,
                           const Vector<DataIndex>&  a_vdi);

    // Brings average[a_phi] to zero.
    static void setZeroAvg (LevelData<FArrayBox>& a_phi);


    // Member variables...
    bool         m_isDefined;
    bool         m_homogeneous;

    int          m_exitStatus;
    Vector<Real> m_resNorms;

    int          m_maxOrder;
    Real         m_eps;
    int          m_verbosity;
    int          m_normType;
    Real         m_hang;
    Real         m_horiz_rhsTol;

    // Horizontal MG solver parameters
    int          m_horiz_imin;
    int          m_horiz_imax;
    int          m_horiz_numSmoothDown;
    int          m_horiz_numSmoothBottom;
    int          m_horiz_numSmoothUp;
    int          m_horiz_numSmoothPreCond;
    int          m_horiz_relaxType;
    int          m_horiz_numMG;
    int          m_horiz_maxDepth;
    Real         m_horiz_eps;
    Real         m_horiz_hang;
    Real         m_horiz_normThresh;
    int          m_horiz_verbosity;

    int          m_horizBottom_imax;
    Real         m_horizBottom_eps;
    int          m_horizBottom_numRestarts;
    Real         m_horizBottom_hang;
    int          m_horizBottom_verbosity;

    // Get these from the LepticOperator
    LinearOp<LevelData<FArrayBox> >*     m_origOpPtr;   // DO NOT CALL DELETE ON THIS!
    RealVect                             m_dx;
    RealVect                             m_dxCrse;
    RefCountedPtr<LevelData<FluxBox> >   m_origJgupPtr;
    RefCountedPtr<LevelData<FArrayBox> > m_origJinvPtr;
    const LevelData<FArrayBox>*          m_crsePhiPtr;

    // Vertical grid stuff
    Copier                                                  m_origToVertCopier;
    Copier                                                  m_vertToOrigCopier;
    ProblemDomain                                           m_domain;
    DisjointBoxLayout                                       m_grids;
    Copier                                                  m_exCopier;
    CornerCopier                                            m_exCornerCopier;
    CFRegion                                                m_CFRegion;
    RefCountedPtr<MappedAMRLevelOp<LevelData<FArrayBox> > > m_opPtr;
    RefCountedPtr<LevelData<FluxBox> >                      m_JgupPtr;
    RefCountedPtr<LevelData<FArrayBox> >                    m_JinvPtr;

    // Flat grid stuff
    BoxLayout                                               m_flatGrids;
    Vector<DataIndex>                                       m_flatDI;
    Vector<DataIndex>                                       m_flatDIComplement;
    LayoutData<Tuple<int,2> >                               m_vertBCTypes;

    // Horizontal grid stuff
    bool                                        m_doHorizSolve;
    ProblemDomain                               m_horizDomain;
    DisjointBoxLayout                           m_horizGrids;
    Copier                                      m_flatToHorizCopier;
    Copier                                      m_horizToFlatCopier;
    MappedAMRPoissonOpFactory*                  m_horizOpFactoryPtr;
    MappedAMRMultiGrid<LevelData<FArrayBox> >*  m_horizSolverPtr;
    BiCGStabSolver<LevelData<FArrayBox> >*      m_horizBottomSolverPtr;
    bool                                        m_horizRemoveAvg;

    static IntVect s_vmask;  // (0,0,1)
    static IntVect s_hmask;  // (1,1,0)
};



// -----------------------------------------------------------------------------
// Inform the user why we stopped solving.
// This just returns the integer code.
// -----------------------------------------------------------------------------
int LevelLepticSolver::exitStatusInt () const
{
    return m_exitStatus;
}


#endif //!__LevelLepticSolver_H__INCLUDED__
