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
#include "AMRPressureSolver.H"
#include "ProblemContext.H"
#include "MappedAMRPoissonOpFactory.H"
#include "AMRLepticSolver.H"
#include "AMRMultiGrid.H"
#include "BiCGStabSolver.H"
#include "RelaxSolver.H"
#include "SetValLevel.H"


// NOTE: The level solver switches are no longer static.

// To use just one solver, set the bools accordingly.
// To use both solvers in tandem, set both bools to true and fiddle with the
// AMRPressureSolver::solve function. I'll make this more automated another day.
bool AMRPressureSolver::s_useAMRLepticSolver = false;
bool AMRPressureSolver::s_useAMRMGSolver = true;


// -----------------------------------------------------------------------------
// Default constructor -- leaves object unusable.
// -----------------------------------------------------------------------------
AMRPressureSolver::AMRPressureSolver ()
: m_useLevelLepticSolver(false),
  m_useLevelMGSolver(false),
  m_isDefined(false),
  m_numLevels(-1),
  m_lepticSolverPtr(NULL),
  m_amrmgSolverPtr(NULL),
  m_bottomSolverPtr(NULL)
{
    // Set the default parameters
    const ProblemContext* ctx = ProblemContext::getInstance();

    setAMRMGParameters(ctx->AMRMG_imin,
                       ctx->AMRMG_imax,
                       ctx->AMRMG_eps,
                       ctx->AMRMG_maxDepth,
                       ctx->AMRMG_num_smooth_precond,
                       ctx->AMRMG_num_smooth_down,
                       ctx->AMRMG_num_smooth_up,
                       ctx->AMRMG_num_smooth_bottom,
                       ctx->AMRMG_precondMode,
                       ctx->AMRMG_relaxMode,
                       ctx->AMRMG_numMG,
                       ctx->AMRMG_hang,
                       ctx->AMRMG_normThresh,
                       ctx->AMRMG_verbosity);

    setBottomParameters(ctx->bottom_imax,
                        ctx->bottom_numRestarts,
                        ctx->bottom_eps,
                        ctx->bottom_reps,
                        ctx->bottom_hang,
                        ctx->bottom_small,
                        ctx->bottom_normType,
                        ctx->bottom_verbosity);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
AMRPressureSolver::~AMRPressureSolver ()
{
    undefine();
}


// -----------------------------------------------------------------------------
// Overrides the default AMRMG settings.
// This can only be called _before_ define.
// -----------------------------------------------------------------------------
void AMRPressureSolver::setAMRMGParameters (const int  a_AMRMG_imin,
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
                                            const int  a_AMRMG_verbosity)
{
    // This function cannot be called after define.
    CH_assert(!isDefined());

    m_AMRMG_imin = a_AMRMG_imin;
    m_AMRMG_imax = a_AMRMG_imax;
    m_AMRMG_eps = a_AMRMG_eps;
    m_AMRMG_maxDepth = a_AMRMG_maxDepth;
    m_AMRMG_num_precond_iters = a_AMRMG_num_precond_iters;
    m_AMRMG_num_smooth_down = a_AMRMG_num_smooth_down;
    m_AMRMG_num_smooth_up = a_AMRMG_num_smooth_up;
    m_AMRMG_num_smooth_bottom = a_AMRMG_num_smooth_bottom;
    m_AMRMG_precondMode = a_AMRMG_precondMode;
    m_AMRMG_relaxMode = a_AMRMG_relaxMode;
    m_AMRMG_numMG = a_AMRMG_numMG;
    m_AMRMG_hang = a_AMRMG_hang;
    m_AMRMG_norm_thresh = a_AMRMG_norm_thresh;
    m_AMRMG_verbosity = a_AMRMG_verbosity;
}


// -----------------------------------------------------------------------------
// Overrides the default bottom solver settings.
// This can only be called _before_ define.
// -----------------------------------------------------------------------------
void AMRPressureSolver::setBottomParameters (const Real a_bottom_imax,
                                             const int  a_bottom_numRestarts,
                                             const Real a_bottom_eps,
                                             const Real a_bottom_reps,
                                             const Real a_bottom_hang,
                                             const Real a_bottom_small,
                                             const int  a_bottom_normType,
                                             const int  a_bottom_verbosity)
{
    // This function cannot be called after define.
    CH_assert(!isDefined());

    m_bottom_imax = a_bottom_imax;
    m_bottom_numRestarts = a_bottom_numRestarts;
    m_bottom_eps = a_bottom_eps;
    m_bottom_reps = a_bottom_reps;
    m_bottom_hang = a_bottom_hang;
    m_bottom_small = a_bottom_small;
    m_bottom_normType = a_bottom_normType;
    m_bottom_verbosity = a_bottom_verbosity;
}


// -----------------------------------------------------------------------------
// Allocates everything and leaves object usable.
// This will not erase the solver parameters.
// -----------------------------------------------------------------------------
void AMRPressureSolver::levelDefine (BCMethodHolder           a_bc,
                                     const LevelGeometry&     a_levGeo,
                                     const int                a_numLevels,
                                     const FillJgupInterface* a_customFillJgupPtr)
{
    // Defining this object can be costly.
    // Let's make sure we aren't being wasteful.
    CH_assert(!isDefined());

    // Using this flag is more reliable than checking if lmin == lmax.
    m_isLevelSolver = true;

    // Collect AMR structures. lev corresponds to our vector index,
    // not the true AMR level.
    m_numLevels = a_numLevels;

    Vector<const LevelGeometry*> amrLevGeos(a_numLevels, NULL);
    Vector<DisjointBoxLayout> amrGrids(a_numLevels);
    Vector<IntVect> amrRefRatios(a_numLevels, IntVect::Zero);

    amrLevGeos[a_numLevels-1] = &a_levGeo;
    amrGrids[a_numLevels-1] = a_levGeo.getBoxes();
    amrRefRatios[a_numLevels-1] = IntVect::Zero; // Should not be used.

    for (signed int lev = a_numLevels-2; lev >= 0; --lev) {
        const LevelGeometry* thisLevGeoPtr = amrLevGeos[lev+1]->getCoarserPtr();
        CH_assert(thisLevGeoPtr != NULL);

        amrLevGeos[lev] = thisLevGeoPtr;
        amrGrids[lev] = thisLevGeoPtr->getBoxes();
        amrRefRatios[lev] = thisLevGeoPtr->getFineRefRatio();
    }

    // Create the op factory
    MappedAMRPoissonOpFactory localPoissonOpFactory;
    localPoissonOpFactory.define(amrLevGeos[0]->getDomain(),
                                 amrGrids,
                                 amrRefRatios,
                                 amrLevGeos[0]->getDx(),
                                 a_bc,
                                 m_AMRMG_maxDepth,
                                 m_AMRMG_num_precond_iters,
                                 m_AMRMG_precondMode,
                                 0.0,               // alpha
                                 1.0,               // beta
                                 amrLevGeos,
                                 m_AMRMG_relaxMode,
                                 false, //isHorizontalFactory,
                                 a_customFillJgupPtr);

    // // Compute the lepticity. This will determine how we solve.
    // const RealVect& dx = amrLevGeos[a_numLevels-1]->getDx();
    // const Real H = amrLevGeos[a_numLevels-1]->getDomainLength(SpaceDim-1);
    // Real lepticity = Min(dx[0], dx[SpaceDim-2]) / H;
    // if (lepticity > 2.0) {
    //     m_useLevelLepticSolver = true;
    //     m_useLevelMGSolver = false;
    // } else {
    //     m_useLevelLepticSolver = true;
    //     m_useLevelMGSolver = true;
    // }
    m_useLevelLepticSolver = false;
    m_useLevelMGSolver = true;


    if (m_useLevelLepticSolver) {
        // Allocate a new AMRLepticSolver.
        AMRLepticSolver* lepticPtr = new AMRLepticSolver;
        CH_assert(lepticPtr != NULL);

        // Set solver parameters
        lepticPtr->define(amrLevGeos[0]->getDomain(),
                          localPoissonOpFactory,
                          m_numLevels);

        lepticPtr->m_verbosity = m_AMRMG_verbosity;
        lepticPtr->m_imin = m_AMRMG_imin;
        lepticPtr->setSolverParameters(m_AMRMG_numMG,
                                       m_AMRMG_imax,
                                       m_AMRMG_eps,
                                       m_AMRMG_hang,
                                       m_AMRMG_norm_thresh);

        // Cast down to an AMREllipticSolver.
        m_lepticSolverPtr = lepticPtr;
    }

    if (m_useLevelMGSolver) {
        // Create the bottom solver
        {
            BiCGStabSolver<LevelData<FArrayBox> >* krylovPtr = new BiCGStabSolver<LevelData<FArrayBox> >;
            CH_assert(krylovPtr != NULL);

            krylovPtr->m_eps = m_bottom_eps;
            krylovPtr->m_reps = m_bottom_reps;
            krylovPtr->m_imax = m_bottom_imax;
            krylovPtr->m_numRestarts = m_bottom_numRestarts;
            krylovPtr->m_hang = m_bottom_hang;
            krylovPtr->m_small = m_bottom_small;
            krylovPtr->m_verbosity = m_bottom_verbosity;
            krylovPtr->m_normType = m_bottom_normType;

            m_bottomSolverPtr = krylovPtr;

            // RelaxSolver<LevelData<FArrayBox> >* krylovPtr = new RelaxSolver<LevelData<FArrayBox> >;
            // CH_assert(krylovPtr != NULL);

            // krylovPtr->m_eps = m_bottom_eps;
            // krylovPtr->m_imax = m_bottom_imax;
            // krylovPtr->m_hang = m_bottom_hang;
            // krylovPtr->m_verbosity = m_bottom_verbosity;
            // krylovPtr->m_normType = m_bottom_normType;

            // m_bottomSolverPtr = krylovPtr;
        }

        // Allocate a new AMRMultiGrid solver.
        MappedAMRMultiGrid<LevelData<FArrayBox> >* amrmgPtr
            = new MappedAMRMultiGrid<LevelData<FArrayBox> >;
        CH_assert(amrmgPtr != NULL);

        // Set solver parameters
        amrmgPtr->define(amrLevGeos[0]->getDomain(),
                         localPoissonOpFactory,
                         m_bottomSolverPtr,
                         m_numLevels);

        amrmgPtr->m_verbosity = m_AMRMG_verbosity;
        amrmgPtr->m_imin = m_AMRMG_imin;
        amrmgPtr->setSolverParameters(m_AMRMG_num_smooth_down,
                                      m_AMRMG_num_smooth_up,
                                      m_AMRMG_num_smooth_bottom,
                                      m_AMRMG_numMG,
                                      m_AMRMG_imax,
                                      m_AMRMG_eps,
                                      m_AMRMG_hang,
                                      m_AMRMG_norm_thresh);

        // // Add an inspector
        // RefCountedPtr<MappedAMRMultiGridInspector<LevelData<FArrayBox> > > insPtr(
        //     new OutputMappedAMRMultiGridInspector<LevelData<FArrayBox> >("LevelProj",
        //                                                                  *amrmgPtr)
        // );
        // amrmgPtr->addInspector(insPtr);

        // Cast down to an AMREllipticSolver.
        m_amrmgSolverPtr = amrmgPtr;
    }

    // We are done. This object is now usable.
    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Allocates everything and leaves object usable.
// This will not erase the solver parameters.
// -----------------------------------------------------------------------------
void AMRPressureSolver::define (BCMethodHolder           a_bc,
                                const LevelGeometry&     a_levGeo,
                                const Box&               a_lminDomBox,
                                const int                a_numLevels,
                                const FillJgupInterface* a_customFillJgupPtr)
{
    // Defining this object can be costly.
    // Let's make sure we aren't being wasteful.
    CH_assert(!isDefined());

    // Using this flag is more reliable than checking if lmin == lmax.
    m_isLevelSolver = false;

    // Collect AMR structures. lev corresponds to our vector index,
    // not the true AMR level.
    m_numLevels = a_numLevels;

    // Gather the levgeos we need...
    Vector<const LevelGeometry*> amrLevGeos;
    const int lmin = 0;
    const int lmax = a_numLevels-1;
    gatherAMRLevGeos(amrLevGeos, a_levGeo, a_lminDomBox, lmin, lmax);

    // ...and the corresponding grids and refRatios.
    Vector<DisjointBoxLayout> amrGrids(a_numLevels);
    Vector<IntVect> amrRefRatios(a_numLevels, IntVect::Zero);
    for (int idx = 0; idx <= a_numLevels-1; ++idx) {
        CH_assert(amrLevGeos[idx] != NULL);
        amrGrids[idx] = amrLevGeos[idx]->getBoxes();
        amrRefRatios[idx] = amrLevGeos[idx]->getFineRefRatio();
    }

    // Create the op factory
    MappedAMRPoissonOpFactory localPoissonOpFactory;
    localPoissonOpFactory.define(amrLevGeos[lmin]->getDomain(),
                                 amrGrids,
                                 amrRefRatios,
                                 amrLevGeos[lmin]->getDx(),
                                 a_bc,
                                 m_AMRMG_maxDepth,
                                 m_AMRMG_num_precond_iters,
                                 m_AMRMG_precondMode,
                                 0.0,               // alpha
                                 1.0,               // beta
                                 amrLevGeos,
                                 m_AMRMG_relaxMode,
                                 false, //isHorizontalFactory,
                                 a_customFillJgupPtr);
    // // Compute the lepticity. This will determine how we solve.
    // const RealVect& dx = amrLevGeos[a_numLevels-1]->getDx();
    // const Real H = amrLevGeos[a_numLevels-1]->getDomainLength(SpaceDim-1);
    // Real lepticity = Min(dx[0], dx[SpaceDim-2]) / H;
    // if (lepticity > 2.0) {
    //     m_useLevelLepticSolver = true;
    //     m_useLevelMGSolver = false;
    // } else {
    //     m_useLevelLepticSolver = true;
    //     m_useLevelMGSolver = true;
    // }
    m_useLevelLepticSolver = false;
    m_useLevelMGSolver = true;

    if (s_useAMRLepticSolver) {
        // Allocate a new AMRLepticSolver.
        AMRLepticSolver* lepticPtr = new AMRLepticSolver;
        CH_assert(lepticPtr != NULL);

        // Set solver parameters
        lepticPtr->define(amrLevGeos[0]->getDomain(),
                          localPoissonOpFactory,
                          m_numLevels);

        lepticPtr->m_verbosity = m_AMRMG_verbosity;
        lepticPtr->m_imin = m_AMRMG_imin;
        lepticPtr->setSolverParameters(m_AMRMG_numMG,
                                       m_AMRMG_imax,
                                       m_AMRMG_eps,
                                       m_AMRMG_hang,
                                       m_AMRMG_norm_thresh);

        // Cast down to an AMREllipticSolver.
        m_lepticSolverPtr = lepticPtr;

    }

    if (s_useAMRMGSolver) {
        // Create the bottom solver
        {
            BiCGStabSolver<LevelData<FArrayBox> >* krylovPtr = new BiCGStabSolver<LevelData<FArrayBox> >;
            CH_assert(krylovPtr != NULL);

            krylovPtr->m_eps = m_bottom_eps;
            krylovPtr->m_reps = m_bottom_reps;
            krylovPtr->m_imax = m_bottom_imax;
            krylovPtr->m_numRestarts = m_bottom_numRestarts;
            krylovPtr->m_hang = m_bottom_hang;
            krylovPtr->m_small = m_bottom_small;
            krylovPtr->m_verbosity = m_bottom_verbosity;
            krylovPtr->m_normType = m_bottom_normType;

            m_bottomSolverPtr = krylovPtr;

            // RelaxSolver<LevelData<FArrayBox> >* krylovPtr = new RelaxSolver<LevelData<FArrayBox> >;
            // CH_assert(krylovPtr != NULL);

            // krylovPtr->m_eps = m_bottom_eps;
            // krylovPtr->m_imax = m_bottom_imax;
            // krylovPtr->m_hang = m_bottom_hang;
            // krylovPtr->m_verbosity = m_bottom_verbosity;
            // krylovPtr->m_normType = m_bottom_normType;

            // m_bottomSolverPtr = krylovPtr;
        }

        // Allocate a new AMRMultiGrid solver.
        MappedAMRMultiGrid<LevelData<FArrayBox> >* amrmgPtr
            = new MappedAMRMultiGrid<LevelData<FArrayBox> >;
        CH_assert(amrmgPtr != NULL);

        // Set solver parameters
        amrmgPtr->define(amrLevGeos[0]->getDomain(),
                         localPoissonOpFactory,
                         m_bottomSolverPtr,
                         m_numLevels);

        amrmgPtr->m_verbosity = m_AMRMG_verbosity;
        amrmgPtr->m_imin = m_AMRMG_imin;
        amrmgPtr->setSolverParameters(m_AMRMG_num_smooth_down,
                                      m_AMRMG_num_smooth_up,
                                      m_AMRMG_num_smooth_bottom,
                                      m_AMRMG_numMG,
                                      m_AMRMG_imax,
                                      m_AMRMG_eps,
                                      m_AMRMG_hang,
                                      m_AMRMG_norm_thresh);

        // Cast down to an AMREllipticSolver.
        m_amrmgSolverPtr = amrmgPtr;
    }

    // We are done. This object is now usable.
    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Frees memory and leaves object unusable.
// This will not erase the solver parameters.
// -----------------------------------------------------------------------------
void AMRPressureSolver::undefine ()
{
    // Do we have anything to do?
    if (!isDefined()) return;

    // Free memory
    delete m_lepticSolverPtr;
    m_lepticSolverPtr = NULL;

    delete m_amrmgSolverPtr;
    m_amrmgSolverPtr = NULL;

    delete m_bottomSolverPtr;
    m_bottomSolverPtr = NULL;

    // Reset members
    m_numLevels = -1;
    m_isDefined = false;
}


// -----------------------------------------------------------------------------
// Solves the Poisson problem.
// -----------------------------------------------------------------------------
void AMRPressureSolver::solve (Vector<LevelData<FArrayBox>*>&       a_phi,
                               const Vector<LevelData<FArrayBox>*>& a_rhs,
                               const int                            a_lmin,
                               const int                            a_lmax,
                               const bool                           a_zeroPhi,
                               const bool                           a_forceHomogeneous)
{
    // Sanity checks
    CH_assert(isDefined());
    CH_assert(a_phi.size() == m_numLevels);
    CH_assert(a_phi.size() == a_rhs.size());
    CH_assert(0 <= a_lmin);
    CH_assert(a_lmax < a_phi.size());

    if (a_zeroPhi) {
        setValLevels(a_phi, a_lmin, a_lmax, 0.0);
    }

    if (m_isLevelSolver) {
        // Level solve
        // You can put a loop here to switch back and forth among solvers.
        // I think that a nicer way to go about this would be to hand the
        // leptic solver a pointer to a krylov solver. This way, the leptic
        // method itself can decide how and when to switch methods.
        if (m_useLevelLepticSolver) {
            CH_assert(m_lepticSolverPtr != NULL);
            m_lepticSolverPtr->solve(a_phi,
                                     a_rhs,
                                     a_lmax,
                                     a_lmin,
                                     false, // zero phi?
                                     a_forceHomogeneous);
        }
        if (m_useLevelMGSolver) {
            CH_assert(m_amrmgSolverPtr != NULL);
            m_amrmgSolverPtr->solve(a_phi,
                                    a_rhs,
                                    a_lmax,
                                    a_lmin,
                                    false, // zero phi?
                                    a_forceHomogeneous);
        }
    } else {
        // AMR solve
        // You can put a loop here to switch back and forth among solvers.
        // I think that a nicer way to go about this would be to hand the
        // leptic solver a pointer to a krylov solver. This way, the leptic
        // method itself can decide how and when to switch methods.
        if (s_useAMRLepticSolver) {
            CH_assert(m_lepticSolverPtr != NULL);
            m_lepticSolverPtr->solve(a_phi,
                                     a_rhs,
                                     a_lmax,
                                     a_lmin,
                                     false, // zero phi?
                                     a_forceHomogeneous);
        }
        if (s_useAMRMGSolver) {
            CH_assert(m_amrmgSolverPtr != NULL);
            m_amrmgSolverPtr->solve(a_phi,
                                    a_rhs,
                                    a_lmax,
                                    a_lmin,
                                    false, // zero phi?
                                    a_forceHomogeneous);
        }
    }
}


// -----------------------------------------------------------------------------
// Solves the Poisson problem.
// -----------------------------------------------------------------------------
void AMRPressureSolver::levelSolve (LevelData<FArrayBox>&       a_phi,
                                    const LevelData<FArrayBox>* a_crsePhiPtr,
                                    const LevelData<FArrayBox>& a_rhs,
                                    const bool                  a_zeroPhi,
                                    const bool                  a_forceHomogeneous)
{
    // Sanity checks
    CH_assert(isDefined());
    CH_assert( (a_crsePhiPtr == NULL && m_numLevels == 1) ||
               (a_crsePhiPtr != NULL && m_numLevels == 2) );

    // Package phi and rhs into vectors
    Vector< LevelData<FArrayBox>* > phiVect(m_numLevels, NULL);
    Vector< LevelData<FArrayBox>* > rhsVect(m_numLevels, NULL);
    phiVect[m_numLevels-1] = &a_phi;
    rhsVect[m_numLevels-1] = (LevelData<FArrayBox>*)&a_rhs;
    if (m_numLevels == 2) {
        phiVect[m_numLevels-2] = (LevelData<FArrayBox>*)a_crsePhiPtr;
    }

    // Solve!
    this->solve(phiVect,
                rhsVect,
                m_numLevels-1,
                m_numLevels-1,
                a_zeroPhi,
                a_forceHomogeneous);
}


// -----------------------------------------------------------------------------
// Static utility
// This collects levgeos and puts them into the correct vector index.
// -----------------------------------------------------------------------------
void AMRPressureSolver::gatherAMRLevGeos (Vector<const LevelGeometry*>& a_amrLevGeos,
                                          const LevelGeometry&          a_levGeo,
                                          const Box&                    a_lminDomBox,
                                          const int                     a_lmin,
                                          const int                     a_lmax)
{
    // Gather every levgeo in existence.
    Vector<const LevelGeometry*> allLevGeos = a_levGeo.getAMRLevGeos();

    // Find the index corresponding to lmin.
    int sidx = 0;
    for (; sidx < allLevGeos.size(); ++sidx) {
        if (allLevGeos[sidx]->getDomain().domainBox() == a_lminDomBox) break;
    }

    // Did we find lmin?
    CH_assert(sidx < allLevGeos.size());
    // Is there enough data to take us to lmax?
    CH_assert(allLevGeos.size()-sidx >= a_lmax-a_lmin+1);

    // Populate the output vector
    a_amrLevGeos.resize(a_lmax+1);
    int idx = 0;
    for (; idx < sidx; ++idx) {
        a_amrLevGeos[idx] = NULL;
    }
    CH_assert(idx == sidx);
    for (; idx <= a_lmax; ++idx) {
        a_amrLevGeos[idx] = allLevGeos[sidx];
        ++sidx;
    }
}
