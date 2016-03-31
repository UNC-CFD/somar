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

// This is a rework of code found in the Chombo library. While these functions
// have been reproduced for ease of modification, some of these functions may
// not have been altered at all. To that end, you can find Chombo's copyright
// file at somar/Chombo/Copyright.txt.

#ifndef __MappedMultiGrid_H__INCLUDED__
#define __MappedMultiGrid_H__INCLUDED__

#include "LinearSolver.H"
#include "AnisotropicRefinementTools.H"
#include "ProblemDomain.H"
#include "CH_Timer.H"


// -----------------------------------------------------------------------------
// This class handles the additional tasks of coordinating operations
// between this level and the level one coarser than this 'level'.
// -----------------------------------------------------------------------------
template <typename T>
class MappedMGLevelOp: public LinearOp<T>
{
public:

    // Constructor
    MappedMGLevelOp ()
    : LinearOp<T>(),
      m_mgCrseRefRatio(D_DECL(-1,-1,-1))
    {;}

    // Destructor
    virtual ~MappedMGLevelOp ()
    {;}

    // Create a coarsened  (by two) version of the input data.  This does not include averaging
    // the data.   So if a_fine is over a Box of (0, 0, 0) (63, 63, 63), a_fine should
    // be over a Box (0, 0, 0) (31, 31, 31).
    virtual void createCoarser (T& a_coarse, const T& a_fine, bool ghosted) = 0;

    // Use your relaxtion operator to remove the high frequency wave numbers from
    // the correction so that it may be averaged to a coarser refinement.
    // A point relaxtion scheme, for example takes the form
    // a_correction -= lambda*(L(a_correction) - a_residual).
    virtual void relax (T& a_correction, const T& a_residual, int a_iterations) = 0 ;


    // calculate restricted residual
    // a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
    virtual void restrictResidual (T& a_resCoarse, T& a_phiFine, const T& a_rhsFine) = 0;

    // correct the fine solution based on coarse correction
    // a_phiThisLevel += I[2h->h](a_correctCoarse)
    virtual void prolongIncrement (T& a_phiThisLevel, const T& a_correctCoarse) = 0;

    // Basic accessor for the MG ref ratio
    virtual void setMGCrseRefRatio (const IntVect& a_mgCrseRefRatio)
    {
        // Each dir must be coarsened by 2 or nothing at all.
        D_TERM(CH_assert(a_mgCrseRefRatio[0] == 1 || a_mgCrseRefRatio[0] == 2);,
               CH_assert(a_mgCrseRefRatio[1] == 1 || a_mgCrseRefRatio[1] == 2);,
               CH_assert(a_mgCrseRefRatio[2] == 1 || a_mgCrseRefRatio[2] == 2);)
        // Why use MG if we aren't coarsening at all?
        CH_assert(a_mgCrseRefRatio.product() > 1);
        m_mgCrseRefRatio = a_mgCrseRefRatio;
    }

    // Basic accessor for the MG ref ratio
    virtual IntVect getMGCrseRefRatio () const
    {
        return m_mgCrseRefRatio;
    }

    // Basic accessor for the MG depth. Useful for debugging.
    // If your factory does not set m_mgDepth, nothing tragic will happen.
    virtual int getMGDepth() const
    {
        return m_mgDepth;
    }

protected:
    // Each element must be either 1 or 2. This needs to be set by the factory.
    IntVect m_mgCrseRefRatio;

    // Top level is 0
    int m_mgDepth;

private:
    // Forbidden operations.
    MappedMGLevelOp(const MappedMGLevelOp&);
    MappedMGLevelOp& operator=(const MappedMGLevelOp&);
};



// -----------------------------------------------------------------------------
// Factory class for generating MappedMGLevelOps
// -----------------------------------------------------------------------------
template <class T>
class MappedMGLevelOpFactory
{
public:
    // Base class constructor.
    MappedMGLevelOpFactory ()
    {;}

    //! Destructor.
    virtual ~MappedMGLevelOpFactory ()
    {;}

    // Create an operator at the coarsened index space.
    // Return NULL if no such Multigrid level can be created at this a_depth.
    // If a_homoOnly = true, then only homogeneous boundary conditions will be needed.
    virtual MappedMGLevelOp<T>* MGnewOp(const ProblemDomain&   a_fineIndexSpace,
                                        const int              a_depth,
                                        const bool             a_homoOnly = true,
                                        Vector<IntVect>*       a_allMGRefRatiosPtr = NULL,
                                        const Vector<IntVect>* a_forceAllMGRefRatiosPtr = NULL) = 0;

private:
    // Forbidden operations.
    MappedMGLevelOpFactory(const MappedMGLevelOpFactory&);
    MappedMGLevelOpFactory& operator=(const MappedMGLevelOpFactory&);
};



// -----------------------------------------------------------------------------
// Class to execute geometric multigrid. This class is not meant to be particularly
// user-friendly and a good option for people who want something a tad less raw is to
// use MappedAMRMultiGrid instead.
//
// NOTE: init, cycle, clear, and setBottomSolver are NOT virtual for speed.
// -----------------------------------------------------------------------------
template <class T>
class MappedMultiGrid
{
public:
    // Constructor
    MappedMultiGrid ();

    // Destructor
    virtual ~MappedMultiGrid ();

    // Function to define a MappedMultiGrid object.
    //  a_factory is the factory for generating operators.
    //  a_bottomSolver is called at the bottom of v-cycle.
    //  a_domain is the problem domain at the top of the vcycle.
    //  maxDepth defines the location of the bottom of the v-cycle.
    //  The vcycle will terminate (hit bottom) when the factory returns NULL for a paticular
    //  depth if maxdepth = -1.  Otherwise the vcycle terminates at maxdepth.
    virtual void define (MappedMGLevelOpFactory<T>& a_factory,
                         LinearSolver<T>*           a_bottomSolver,
                         const ProblemDomain&       a_domain,
                         const int                  a_maxDepth = -1,
                         MappedMGLevelOp<T>*        a_finestLevelOp = NULL,
                         const Vector<IntVect>*     a_forceAllMGRefRatiosPtr = NULL);

    // Solve L(a_phi) = a_rhs.   Tolerance is how much you want the norm of the error reduced.
    // verbosity is how chatty you want the function to be.   maxIterations is the maximum number
    // of v-cycles.   This does the whole residual correction switcharoo and calls oneCycle up to
    // maxIterations times, evaluating the residual as it goes.
    virtual void solve (T&         a_phi,
                        const T&   a_rhs,
                        const Real a_tolerance,
                        const int  a_maxIterations);

    // Execute ONE v-cycle of multigrid.
    // If you want the solution to converge, you need to iterate this.
    // See solve() or AMRMultiGrid::solve for a more automatic solve() function.
    // This operates residual-correction form of solution
    // so all boundary conditions are assumed to be homogeneous.
    // L(a_e) = a_res
    virtual void oneCycle (T&       a_e,
                           const T& a_res);

    // Needed by MappedAMRMultiGrid
    void cycle (const int a_depth,
                T&        a_correction,
                const T&  a_residual);

    // used by AMRMultiGrid--not a part of the public API
    void setBottomSolver (LinearSolver<T>* a_bottomSolver);

    // Allocates memory for the residuals and corrections.
    void init (const T& a_correction,
               const T& a_residual);

    // for changing coefficients --- not for the faint of heart.
    virtual Vector<MappedMGLevelOp<T>*> getAllOperators ();

    // Public solver parameters. m_pre and m_post are the ones
    // that usually get set and are the number of relaxations performed
    // before and after multigrid recursion. See AMRMultiGrid for a
    // more user-friendly interface.
    int  m_pre;             // Number of relaxation iters on the way down the V-cycle.
    int  m_post;            // Number of relaxation iters on the way up the V-cycle.
    int  m_bottom;          // Number of relaxation iters before calling the bottom solver.
    int  m_cycle;           // 1 = V-cycle, 2 = W-cycle, etc... (If cycle is < 0, then FMG).
    int  m_depth;           // Stores the current refinement depth. Useful because cycle is recursive.
    int  m_defaultDepth;    // This is where m_depth is when not in a cycle. (max depth + 1).
    int  m_maxForcedDepth;  // This is how deep the coarsening pattern was forced in define.
    bool m_homogeneous;     // Have the BCs already been handled? (Are we working on the residual equation?)
    LinearSolver<T>* m_bottomSolver;
    int m_verbosity;

protected:
    // Returns the ProblemDomain corresponding to a particular depth.
    // If a_depth is zerom you'll just get a_fineIndexSpace back.
    // WARNING: This function does not perform assertions so that it can be used by
    // define. Use at your own risk!
    virtual ProblemDomain getDomainAtDepth (const ProblemDomain& a_fineIndexSpace,
                                            const int            a_depth) const;

    // Undefine this solver and free any used memory.
    void clear ();


    bool                        m_defined;          // Is this solver in a usable state?
    int                         m_bottomCells;      // The number of cells at the max depth.
    ProblemDomain               m_topLevelDomain;   // The uncoarsened domain.
    Vector<MappedMGLevelOp<T>*> m_op;               // The ops at every depth.
    std::vector<bool>           m_ownOp;            // Flags that signal which ops the destructor should delete.
    Vector<T*>                  m_residual;         // The residual at every depth.
    Vector<T*>                  m_correction;       // The correction at every depth.

private:
    MappedMultiGrid(const MappedMultiGrid<T>& a_opin) {
        MayDay::Error("invalid operator");
    }

    void operator=(const MappedMultiGrid<T>& a_opin) {
        MayDay::Error("invalid operator");
    }

    // This is a function defined in LinearOp. In this class and all of its
    // children, a call to this function should throw an error.
    virtual Real dx () const
    {
        MayDay::Error("A mapped op should not call the isotropic dx accessor");
        return 0.0;
    }
};



// ********************* MappedMultiGrid Implementation ************************

// -----------------------------------------------------------------------------
// Returns the ProblemDomain corresponding to a particular depth.
// If a_depth is zerom you'll just get a_fineIndexSpace back.
// WARNING: This function does not perform assertions so that it can be used by
// define. Use at your own risk!
// -----------------------------------------------------------------------------
template <class T>
ProblemDomain MappedMultiGrid<T>::getDomainAtDepth (const ProblemDomain& a_fineIndexSpace,
                                                    const int            a_depth) const
{
    IntVect ref = IntVect::Unit;
    for (int d = 0; d < a_depth; ++d) {
        ref *= m_op[d]->getMGCrseRefRatio();
    }

    ProblemDomain retDom = a_fineIndexSpace;
    coarsen(retDom, a_fineIndexSpace, ref);
    return retDom;
}


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
template <class T>
MappedMultiGrid<T>::MappedMultiGrid ()
: m_pre(3),
  m_post(3),
  m_bottom(0),
  m_cycle(1),
  m_maxForcedDepth(-1),
  m_homogeneous(true),
  m_verbosity(3),
  m_defined(false)
{;}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
template <class T>
MappedMultiGrid<T>::~MappedMultiGrid ()
{
    this->clear();
}


// -----------------------------------------------------------------------------
// Function to define a MappedMultiGrid object.
//  a_factory is the factory for generating operators.
//  a_bottomSolver is called at the bottom of v-cycle.
//  a_domain is the problem domain at the top of the vcycle.
//  maxDepth defines the location of the bottom of the v-cycle.
//  The vycle will terminate (hit bottom) when the factory returns NULL for a
//  particular depth if maxdepth = -1. Otherwise the vycle terminates at maxdepth.
// -----------------------------------------------------------------------------
template <class T>
void MappedMultiGrid<T>::define (MappedMGLevelOpFactory<T>& a_factory,
                                 LinearSolver<T>*           a_bottomSolver,
                                 const ProblemDomain&       a_domain,
                                 const int                  a_maxDepth,
                                 MappedMGLevelOp<T>*        a_finestOp,
                                 const Vector<IntVect>*     a_forceAllMGRefRatiosPtr)
{
    CH_TIME("MappedMultiGrid::define");

    // Reset everything.
    this->clear();
    m_depth = 0;
    m_residual.resize(0);   // the zero'th entry is left null constructed
    m_correction.resize(0); // the zero'th entry is left null constructed
    m_op.resize(0);
    m_ownOp.resize(0);

    // Set what we immediately know.
    m_bottomSolver = a_bottomSolver;
    m_topLevelDomain = a_domain;

    // This is not permanent...we are about to coarsen this.
    m_bottomCells = a_domain.domainBox().numPts();

    // We are about to define operators at every depth. At the top level,
    // we may have been given an operator, so it is a special case.
    MappedMGLevelOp<T>* nextOp;
    bool ownOp;
    if (a_finestOp == NULL) {
        nextOp =  a_factory.MGnewOp(a_domain, 0, true);
        ownOp  = true;
    } else {
        nextOp = a_finestOp;
        ownOp  = false;
    }

    // Now, we have the finest level op. Begin at this level and move on down,
    // defining everything in our path.
    Vector<IntVect> mgRefRatios(0);
    while (nextOp != NULL) {
        // Define everything at this depth.
        m_residual.push_back(new T());
        m_correction.push_back(new T());
        m_op.push_back(nextOp);
        m_ownOp.push_back(ownOp);

        // Move to the next coarser level and get a new op, if necessary.
        m_depth++;
        if ((m_depth < a_maxDepth) || (a_maxDepth < 0)) {
            MappedMGLevelOp<T>* fineOp = nextOp;

            // Create a coarser multigrid operator.
            CH_assert(mgRefRatios.size() == m_depth-1);
            nextOp = a_factory.MGnewOp(a_domain, m_depth, true, &mgRefRatios, a_forceAllMGRefRatiosPtr);
            ownOp = true;

            // Set the finer level's mgCrseRefRatio.
            if (nextOp != NULL) {
                CH_assert(mgRefRatios.size() > 0);
                fineOp->setMGCrseRefRatio(mgRefRatios.back());
            }
        } else {
            nextOp = NULL;
        }
    }

    // We are at the bottom level. Define the bottom solver.
    m_bottomSolver->define(m_op[m_depth - 1], true);
    m_bottomCells = getDomainAtDepth(a_domain, m_depth-1).domainBox().numPts();    // TODO: Test that the indexing is correct.
    m_defaultDepth = m_depth;
    if (a_forceAllMGRefRatiosPtr != NULL) {
        m_maxForcedDepth = a_forceAllMGRefRatiosPtr->size();
    }

    // We are done. This solver is usable.
    m_defined = true;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
void MappedMultiGrid<T>::setBottomSolver (LinearSolver<T>* a_bottomSolver)
{
    m_bottomSolver = a_bottomSolver;
    m_bottomSolver->define(m_op[m_depth - 1], true);
}


// -----------------------------------------------------------------------------
// Allocates memory for the residuals and corrections.
// -----------------------------------------------------------------------------
template <class T>
void MappedMultiGrid<T>::init (const T& a_e,
                               const T& a_residual)
{
    CH_TIME("MutliGrid::init");
    if (m_depth > 1) {
        m_op[0]->createCoarser(*(m_residual[1]), a_residual, false);
        m_op[0]->createCoarser(*(m_correction[1]), a_e, true);
    }

    for (int i = 2; i < m_depth; i++) {
        m_op[i - 1]->createCoarser(*(m_residual[i]), *(m_residual[i - 1]), false);
        m_op[i - 1]->createCoarser(*(m_correction[i]), *(m_correction[i - 1]), true);
    }
}


// -----------------------------------------------------------------------------
// Undefine this solver and free any used memory.
// -----------------------------------------------------------------------------
template <class T>
void MappedMultiGrid<T>::clear()
{
    CH_assert(m_ownOp.size() == m_op.size());
    CH_assert(m_correction.size() == m_op.size());
    CH_assert(m_residual.size() == m_op.size());

    std::vector<bool>::reverse_iterator it = m_ownOp.rbegin();
    for (int i = m_op.size() - 1; i >= 0; --i, ++it) {
        if (*it) delete m_op[i];
        delete m_correction[i];
        delete m_residual[i];
    }
    m_op.resize(0);
    m_residual.resize(0);
    m_correction.resize(0);
    m_defined = false;
}


// -----------------------------------------------------------------------------
// This is the main solve function that should be called externally.
// -----------------------------------------------------------------------------
template <class T>
void MappedMultiGrid<T>::solve (T&         a_phi,
                                const T&   a_rhs,
                                const Real a_tolerance,
                                const int  a_maxIterations)
{
    CH_TIME("MappedMultiGrid::solve");

    // Allocate space for the solve.
    this->init(a_phi, a_rhs);
    T correction, residual;
    m_op[0]->create(correction, a_phi);
    m_op[0]->create(residual, a_rhs);

    // Initialize the fields for the entire solve.
    m_op[0]->setToZero(a_phi);
    m_op[0]->residual(residual, a_phi, a_rhs, false);

    // Initialize the convergence metrics for the entire solve.
    Real errorno = m_op[0]->norm(residual, 0);
    if (m_verbosity > 2) {
        pout() << "multigrid::solve initial residual = " << errorno << std::endl;
    }
    Real compval = a_tolerance * errorno;
    Real epsilon = 1.0e-16;
    compval = Max(compval, epsilon);

    // Time to iterate...
    Real error = errorno;
    int iter = 0;
    while ((error > compval) && (error > a_tolerance * errorno) && (iter < a_maxIterations)) {
        // Initialize the fields for this iteration.
        m_op[0]->setToZero(correction);
        m_op[0]->residual(residual, a_phi, a_rhs, false);   // TODO: I don't think this is needed if iter is 0.

        // Initialize the convergence metrics for this iteration.
        error = m_op[0]->norm(residual, 0);
        if (m_verbosity > 3) {
            pout() << "multigrid::solve iter = " << iter <<  ",  residual = " << error << std::endl;
        }

        // Perform a cycle down to max depth and back.
        this->cycle(0, correction, residual);

        // Update the total correction and move on.
        m_op[0]->incr(a_phi, correction, 1.0);
        iter++;
    }

    // We are finished solving. Inform the user.
    if (m_verbosity > 2) {
        pout() << "multigrid::solve final residual = " << error << std::endl;
    }

    // Free memory that this function allocated.
    m_op[0]->clear(correction);
    m_op[0]->clear(residual);
}


// -----------------------------------------------------------------------------
// Internal function
// This recursively completes one V-cycle of the solver.
// -----------------------------------------------------------------------------
template <class T>
void MappedMultiGrid<T>::oneCycle(T& a_e, const T& a_residual)
{
    CH_TIME("MappedMultiGrid::oneCycle");

    if (m_homogeneous) {
        this->cycle(0, a_e, a_residual);
    } else {
        T correction, residual;
        m_op[0]->create(correction, a_e);
        m_op[0]->create(residual, a_residual);

        m_op[0]->residual(residual, a_e, a_residual);

        m_op[0]->setToZero(correction);
        this->cycle(0, correction, residual);

        m_op[0]->incr(a_e, correction, 1.0);
        m_op[0]->clear(correction);
        m_op[0]->clear(residual);
    }
}


// -----------------------------------------------------------------------------
// Internal function
// -----------------------------------------------------------------------------
template <class T>
void MappedMultiGrid<T>::cycle(int depth, T& correction, const T& residual)
{
    CH_TIME("MappedMultiGrid::cycle");

    if (m_verbosity >= 5) {
        Real rnorm = m_op[depth]->norm(residual, 0);
        pout() << "\tInitial residual norm at depth " << depth << " is " << rnorm << endl;
    }

    // currently I can't drop a timer in this function because it is recursive.  doh
    if (depth == m_depth - 1) {
        if (m_bottomCells == 1) {
            if (m_verbosity >= 5) pout() << "\trelax at depth " << depth << endl;
            m_op[depth  ]->relax(correction, residual, 1);
        } else {
            if (m_verbosity >= 5) pout() << "\trelax at depth " << depth << endl;
            m_op[depth  ]->relax(correction, residual, m_bottom);
            if (m_verbosity >= 5) pout() << "\tsolve at depth " << depth << endl;
            m_bottomSolver->solve(correction, residual);
        }
        //m_op[depth  ]->relax(correction, residual, 1);
    } else {
        int cycles = m_cycle;
        if ( cycles < 0 ) {
            cycles = -cycles;
            // F cycle multigrid
            if (m_verbosity >= 5) {
                pout() << "\trestrict residual from depth " << depth << " to " << depth+1
                       << " with mgRefRatio = " << m_op[depth]->getMGCrseRefRatio() << endl;
            }
            m_op[depth  ]->restrictResidual(*(m_residual[depth + 1]), correction, residual);
            m_op[depth  ]->setToZero(*(m_correction[depth + 1]));

            // recursive call
            cycle(depth + 1, *(m_correction[depth + 1]), *(m_residual[depth + 1]));

            if (m_verbosity >= 5) {
                pout() << "\tProlong increment from depth " << depth+1 << " to " << depth
                       << " with mgRefRatio = " << m_op[depth]->getMGCrseRefRatio() << endl;
            }
            m_op[depth  ]->prolongIncrement(correction, *(m_correction[depth + 1]));

            if (m_verbosity >= 5) pout() << "\trelax at depth " << depth << endl;
            m_op[depth  ]->relax(correction, residual, m_pre);

            for (int img = 0; img < cycles; img++) {
                if (m_verbosity >= 5) {
                    pout() << "\trestrict residual from depth " << depth << " to " << depth+1
                           << " with mgRefRatio = " << m_op[depth]->getMGCrseRefRatio() << endl;
                }
                m_op[depth  ]->restrictResidual(*(m_residual[depth + 1]), correction, residual);
                m_op[depth  ]->setToZero(*(m_correction[depth + 1]));

                m_cycle = 1;          // hack to get a V-cycle
                cycle(depth + 1, *(m_correction[depth + 1]), *(m_residual[depth + 1]));
                m_cycle = -cycles;
                //
                if (m_verbosity >= 5) {
                    pout() << "\tProlong increment from depth " << depth+1 << " to " << depth
                           << " with mgRefRatio = " << m_op[depth]->getMGCrseRefRatio() << endl;
                }
                m_op[depth  ]->prolongIncrement(correction, *(m_correction[depth + 1]));
            }

            m_op[depth  ]->relax(correction, residual, m_post);
        } else {
            if (m_verbosity >= 5) pout() << "\trelax at depth " << depth << endl;
            m_op[depth  ]->relax( correction, residual, m_pre );

            if (m_verbosity >= 5) {
                pout() << "\trestrict residual from depth " << depth << " to " << depth+1
                       << " with mgRefRatio = " << m_op[depth]->getMGCrseRefRatio() << endl;
            }
            m_op[depth  ]->restrictResidual(*(m_residual[depth + 1]), correction, residual);
            m_op[depth  ]->setToZero(*(m_correction[depth + 1]));

            for (int img = 0; img < cycles; img++) {
                cycle(depth + 1, *(m_correction[depth + 1]), *(m_residual[depth + 1]));
            }
            if (m_verbosity >= 5) {
                pout() << "\tProlong increment from depth " << depth+1 << " to " << depth
                       << " with mgRefRatio = " << m_op[depth]->getMGCrseRefRatio() << endl;
            }
            m_op[depth  ]->prolongIncrement(correction, *(m_correction[depth + 1]));

            if (m_verbosity >= 5) pout() << "\trelax at depth " << depth << endl;
            m_op[depth  ]->relax(correction, residual, m_post);
        }
    }

    if (m_verbosity >= 5) {
        T newRes;
        m_op[depth]->create(newRes, residual);
        m_op[depth]->residual(newRes, correction, residual, true);

        Real rnorm = m_op[depth]->norm(newRes, 0);
        pout() << "\tFinal residual norm at depth " << depth << " is " << rnorm << endl;
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
Vector<MappedMGLevelOp<T>*> MappedMultiGrid<T>::getAllOperators ()
{
    Vector<MappedMGLevelOp<T>*> retval = m_op;
    return retval;
}


#endif //!__MappedMultiGrid_H__INCLUDED__
