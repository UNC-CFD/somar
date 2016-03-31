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
#ifndef __MappedAMRMultiGrid_H__INCLUDED__
#define __MappedAMRMultiGrid_H__INCLUDED__

#include "AMREllipticSolver.H"
#include "MappedMultiGrid.H"
#include "AnisotropicRefinementTools.H"
#include "RefCountedPtr.H"
#include "NoOpSolver.H"
#include "AMRIO.H"
#include "Printing.H"

template <class T> class MappedAMRMultiGrid;


// -----------------------------------------------------------------------------
// Operator class for AMR Multigrid
// -----------------------------------------------------------------------------
template <typename T>
class MappedAMRLevelOp : public MappedMGLevelOp<T>
{
public:
    // Constructor
    MappedAMRLevelOp ()
    : MappedMGLevelOp<T>()
    {;}

    // Destructor
    virtual ~MappedAMRLevelOp ()
    {;}

    // return the refinement ratio to next coarser level.
    // return 1 when there are no coarser AMRLevelOp objects.
    virtual IntVect refToCoarser () const = 0;

    // a_residual = a_rhs - L(a_phiFine, a_phi, a_phiCoarse)
    virtual void AMRResidual (T&                   a_residual,
                              const T&             a_phiFine,
                              const T&             a_phi,
                              const T&             a_phiCoarse,
                              const T&             a_rhs,
                              const bool           a_homogeneousDomBC,
                              MappedAMRLevelOp<T>* a_finerOp) = 0;

    // a_residual = a_rhs - L^nf(a_phi, a_phiCoarse)
    // assume no finer AMR level
    virtual void AMRResidualNF (T&         a_residual,
                                const T&   a_phi,
                                const T&   a_phiCoarse,
                                const T&   a_rhs,
                                const bool a_homogeneousBC) = 0;

    // a_residual = a_rhs - L(a_phiFine, a_phi)
    // assume no coarser AMR level
    virtual void AMRResidualNC (T&                   a_residual,
                                const T&             a_phiFine,
                                const T&             a_phi,
                                const T&             a_rhs,
                                const bool           a_homogeneousBC,
                                MappedAMRLevelOp<T>* a_finerOp) = 0;

    // Apply the AMR operator, including coarse-fine matching
    virtual void AMROperator (T&                   a_LofPhi,
                              const T&             a_phiFine,
                              const T&             a_phi,
                              const T&             a_phiCoarse,
                              const bool           a_homogeneousDomBC,
                              MappedAMRLevelOp<T>* a_finerOp) = 0;

    // Apply the AMR operator, including coarse-fine matching.
    // assume no finer AMR level
    virtual void AMROperatorNF (T&         a_LofPhi,
                                const T&   a_phi,
                                const T&   a_phiCoarse,
                                const bool a_homogeneousBC) = 0;

    // Apply the AMR operator, including coarse-fine matching
    // assume no coarser AMR level
    virtual void AMROperatorNC (T&                   a_LofPhi,
                                const T&             a_phiFine,
                                const T&             a_phi,
                                const bool           a_homogeneousBC,
                                MappedAMRLevelOp<T>* a_finerOp) = 0;


    // a_resCoarse = I[h-2h]( a_residual - L(a_correction, a_coarseCorrection))
    virtual void AMRRestrict (T&       a_resCoarse,
                              const T& a_residual,
                              const T& a_correction,
                              const T& a_coarseCorrection) = 0;

    // a_correction += I[2h->h](a_coarseCorrection)
    virtual void AMRProlong (T&       a_correction,
                             const T& a_coarseCorrection) = 0;

    // a_residual = a_residual - L(a_correction, a_coarseCorrection)
    virtual void AMRUpdateResidual (T&       a_residual,
                                    const T& a_correction,
                                    const T& a_coarseCorrection) = 0;

    // Compute norm over all cells on coarse not covered by finer
    virtual Real AMRNorm (const T&       a_coarResid,
                          const T&       a_fineResid,
                          const IntVect& a_refRat,
                          const int      a_ord) {
        return this->norm(a_coarResid, 0);
    }

    //
    virtual void createCoarsened(T&             a_lhs,
                                 const T&       a_rhs,
                                 const IntVect& a_refRat) = 0;


    // Debugging functions...

    //
    virtual void dumpAMR (const Vector<T*>& a_data,
                          const std::string a_name)
    {;}

    //
    virtual void dumpLevel (const T&          a_data,
                            const std::string a_name)
    {;}

    //
    virtual void dumpStuff (const Vector<T*>  a_data,
                            const std::string a_filename)
    {;}

    //
    virtual void outputLevel (const T&          a_rhs,
                              const std::string a_name)
    {;}

    //
    virtual void outputAMR (const Vector<T*>& a_rhs,
                            const std::string a_name)
    {;}


    // -------------------------------------------------------------------------
    // optional optimizations for an AMRLevelOp.  These are not pure virtual
    // functions, since we can build the equivalent algorithmic components from
    // pure virtual functions.  The AMRMultiGrid algorithm actually calls *these*
    // functions, which a smart operator can perform in faster ways.
    // -------------------------------------------------------------------------

    //
    virtual void buildCopier(Copier&  a_copier,
                             const T& a_lhs,
                             const T& a_rhs)
    {;}

    //
    virtual void assignCopier (T&            a_lhs,
                               const T&      a_rhs,
                               const Copier& a_copier)
    {
        this->assign(a_lhs, a_rhs);
    }

    //
    virtual void zeroCovered (T&            a_lhs,
                              T&            a_rhs,
                              const Copier& a_copier)
    {
        this->setToZero(a_rhs);
        this->assignCopier(a_lhs, a_rhs, a_copier);
    }

    //
    virtual Real localMaxNorm (const T& a_phi)
    {
        return this->norm(a_phi, 0);
    }

    /** optimization of AMRProlong that sends in the existing temporary and copier */
    virtual void AMRProlongS (T&            a_correction,
                              const T&      a_coarseCorrection,
                              T&            a_temp,
                              const Copier& a_copier)
    {
        AMRProlong(a_correction, a_coarseCorrection);
    }

    //
    virtual void AMRRestrictS (T&       a_resCoarse,
                               const T& a_residual,
                               const T& a_correction,
                               const T& a_coarseCorrection,
                               T&       a_scratch)
    {
        AMRRestrict(a_resCoarse, a_residual, a_correction, a_coarseCorrection);
    }

    //
    virtual unsigned int orderOfAccuracy () const
    {
        return 2;
    }

    // This routine is for operators with orderOfAccuracy() > 2.
    virtual void enforceCFConsistency (T&       a_coarseCorrection,
                                       const T& a_correction)
    {;}
};



// -----------------------------------------------------------------------------
// Factory to create AMRLevelOps
// -----------------------------------------------------------------------------
template <class T>
class MappedAMRLevelOpFactory : public MappedMGLevelOpFactory<T> {
public:
    virtual ~MappedAMRLevelOpFactory ()
    {;}

    // return a new operator. This is done with a new call.
    // caller is responsible for deletion
    virtual MappedAMRLevelOp<T>* AMRnewOp (const ProblemDomain& a_indexSpace) = 0;

    // return refinement ratio to next finer level.
    virtual IntVect getFineRefRatio (const ProblemDomain& a_indexSpace) const = 0;
};



// -----------------------------------------------------------------------------
// This base class allows one to construct methods for inspecting the
// multigrid algorithm at each of its steps.
// NOTE: I have reproduced this not to add or subtract functionality, but
// to avoid including the file that defines the original AMRMultiGridInspector.
// -----------------------------------------------------------------------------
template <class T>
class MappedAMRMultiGridInspector
{
public:

    // Base class constructor. This must be called by all subclasses.
    MappedAMRMultiGridInspector ()
    {;}

    // Destructor
    virtual ~MappedAMRMultiGridInspector ()
    {;}

    // Override this method to keep track of a multigrid residual computed
    // during a multigrid iteration at the given levels.
    // a_residuals: An array containing the residuals computed by the multigrid
    //              algorithm at each level in the range [a_minLevel, a_maxLevel].
    // a_minLevel: The lowest AMR level at which residuals were computed.
    // a_maxLevel: The highest AMR level at which residuals were computed.
    // a_iter: The multigrid iteration number.
    virtual void recordResiduals (const Vector<T*>& a_residuals,
                                  int               a_minLevel,
                                  int               a_maxLevel,
                                  int               a_iter) = 0;

    // Override this method to keep track of a multigrid correction computed
    // during a V cycle at the given level.
    // a_corrections: An array containing the corrections computed during a V cycle
    //                at each level in the range [a_minLevel, a_maxLevel].
    // a_minLevel: The lowest AMR level at which corrections were computed.
    // a_maxLevel: The highest AMR level at which corrections were computed.
    // a_iter: The multigrid iteration number.
    virtual void recordCorrections (const Vector<T*>& a_corrections,
                                    int               a_minLevel,
                                    int               a_maxLevel,
                                    int               a_iter) = 0;

private:
    // Illegal functions
    MappedAMRMultiGridInspector (const MappedAMRMultiGridInspector&);
    MappedAMRMultiGridInspector& operator= (const MappedAMRMultiGridInspector&);
};


// -----------------------------------------------------------------------------
// This subclass of AMRMultiGridInspector simply writes output files
// containing residuals and corrections at each multigrid iteration.
// -----------------------------------------------------------------------------
template <class T>
class OutputMappedAMRMultiGridInspector: public MappedAMRMultiGridInspector<T>
{
public:
    // Create an AMRMultiGridInspector that inspects the given AMRMultiGrid solver and
    // generates output files containing intermediate data.
    // a_name The name of the inspector. This gets prepended to the output files.
    // a_solver The AMRMultiGrid solver with which this inspector is associated.
    OutputMappedAMRMultiGridInspector (const std::string&     a_name,
                                       MappedAMRMultiGrid<T>& a_solver)
    : MappedAMRMultiGridInspector<T>(),
      m_name(a_name),
      m_solver(&a_solver)
    {;}

    // Destructor.
    ~OutputMappedAMRMultiGridInspector ()
    {;}

    // Overridden methods.
    void recordResiduals (const Vector<T*>& a_residuals,
                          int               a_minLevel,
                          int               a_maxLevel,
                          int               a_iter)
    {
        char name[1024];
        snprintf(name, 1024, "%s.residual.iter.%d.hdf5", m_name.c_str(), a_iter);
        std::string nameStr(name);
        m_solver->outputAMR(const_cast<Vector<T*>&>(a_residuals),
                            nameStr,
                            a_maxLevel,
                            a_minLevel);
    }

    void recordCorrections (const Vector<T*>& a_corrections,
                            int               a_minLevel,
                            int               a_maxLevel,
                            int               a_iter)
    {
        char name[1024];
        snprintf(name, 1024, "%s.correction.iter.%d.hdf5", m_name.c_str(), a_iter);
        std::string nameStr(name);
        m_solver->outputAMR(const_cast<Vector<T*>&>(a_corrections),
                            nameStr,
                            a_maxLevel,
                            a_minLevel);
    }

private:
    std::string m_name;
    MappedAMRMultiGrid<T>* m_solver;
};


// -----------------------------------------------------------------------------
// Class to solve elliptic equations using the Martin and Cartwright algorithm.
// -----------------------------------------------------------------------------
template <class T>
class MappedAMRMultiGrid: public AMREllipticSolver<T> {
public:
    // Constructor
    MappedAMRMultiGrid ();

    // Destructor
    virtual ~MappedAMRMultiGrid ();

    // Define the solver.
    // a_coarseDomain is the index space on the coarsest AMR level.
    // a_factory is the operator factory through which all special information is conveyed.
    // a_bottomSolver is the solver to be used at the termination of multigrid coarsening.
    //  It is the client's responsibility to free up the dynamically-allocated memory.
    // a_numLevels is the number of AMR levels.
    virtual void define (const ProblemDomain&        a_coarseDomain,
                         MappedAMRLevelOpFactory<T>& a_factory,
                         LinearSolver<T>*            a_bottomSolver,
                         const int                   a_numLevels,
                         const int                   a_verbosity = 3);

    // Add an inspector to the list of inspectors maintained by this MappedAMRMultiGrid
    // instance. It will be given the opportunity to record intermediate data.
    virtual void addInspector (RefCountedPtr<MappedAMRMultiGridInspector<T> >& a_inspector)
    {
        CH_assert(!a_inspector.isNull());
        m_inspectors.push_back(a_inspector);
    }

    // Solve L(phi) = rho from l_base to l_max.  To solve over all levels,
    // l_base = 0 and l_max = max_level = numLevels-1.
    virtual void solve (Vector<T*>&       a_phi,
                        const Vector<T*>& a_rhs,
                        const int         a_lMax,
                        const int         a_lBase,
                        const bool        a_zeroPhi = true,
                        const bool        a_forceHomogeneous = false);

    // Same as "solve" except user has taken the reponsibility of having previously
    // called "init" so solver can allocate temporary holders.
    virtual void solveNoInit (Vector<T*>&       a_phi,
                              const Vector<T*>& a_rhs,
                              const int         a_lMax,
                              const int         a_lBase,
                              const bool        a_zeroPhi = true,
                              const bool        a_forceHomogeneous = false);

    // use if you want final residual
    virtual void solveNoInitResid (Vector<T*>&       a_phi,
                                   Vector<T*>&       a_finalResid,
                                   const Vector<T*>& a_rhs,
                                   const int         a_lMax,
                                   const int         a_lBase,
                                   const bool        a_zeroPhi = true,
                                   const bool        a_forceHomogeneous = false);

    //
    virtual void relaxOnlyHomogeneous (Vector<T*>&       a_phi,
                                       const Vector<T*>& a_rhs,
                                       const int         a_lMax,
                                       const int         a_lBase);

    //
    virtual void AMRVCycle(Vector<T*>&       a_correction,
                           const Vector<T*>& a_residual,
                           const int         a_l,
                           const int         a_lMax,
                           const int         a_lBase);

    //
    virtual void setMGCycle (const int a_numMG);

    //
    virtual void init (const Vector<T*>& a_phi,
                       const Vector<T*>& a_rhs,
                       const int         a_lMax,
                       const int         a_lBase);

    // init messes with multigrid depth. This puts it back
    virtual void revert (const Vector<T*>& a_phi,
                         const Vector<T*>& a_rhs,
                         const int         a_lMax,
                         const int         a_lBase);

    //
    virtual MappedAMRLevelOp<T>& levelOp (const int level);

    // resid = L(phi) - rhs
    virtual Real computeAMRResidual (Vector<T*>&       a_resid,
                                     const Vector<T*>& a_phi,
                                     const Vector<T*>& a_rhs,
                                     const int         a_lMax,
                                     const int         a_lBase,
                                     const bool        a_homogeneousBC = false,
                                     const bool        a_computeNorm = true);


    // Just return the normed value of computeAMRResidual. Used for benchmarking.
    virtual Real computeAMRResidual (Vector<T*>&       a_phi,
                                     const Vector<T*>& a_rhs,
                                     const int         a_lMax,
                                     const int         a_lMin);

    // lph = L(phi)
    virtual void computeAMROperator (Vector<T*>&       a_lph,
                                     const Vector<T*>& a_phi,
                                     const int         a_lMax,
                                     const int         a_lBase,
                                     const bool        a_homogeneousBC = false);

    // For changing coefficients.  Use at thy own peril.
    virtual Vector<MappedMGLevelOp<T>*>          getAllOperators ();
    virtual Vector<MappedMGLevelOp<T>*>          getOperatorsOp ();
    virtual Vector<Vector<MappedMGLevelOp<T>*> > getOperatorsMG ();
    virtual Vector<MappedAMRLevelOp<T>*>&        getAMROperators ()
    {
        return m_op;
    }

    // Set parameters of the solve.
    // a_pre is the number of smoothings before averaging.
    // a_post is the number of smoothings after averaging.
    // a_bottom is the number of smoothings at the bottom level.
    // a_numMG = 1 for vcycle, =2 for wcycle (use 1).
    // a_itermax is the max number of v cycles.
    // a_hang is the minimum amount of change per vcycle.
    // a_eps is the solution tolerance.
    // a_normThresh is how close to zero eps*resid is allowed to get.
    virtual void setSolverParameters (const int  a_pre,
                                      const int  a_post,
                                      const int  a_bottom,
                                      const int  a_numMG,
                                      const int  a_iterMax,
                                      const Real a_eps,
                                      const Real a_hang,
                                      const Real a_normThresh);


    // set up bottom solver for internal MG solver
    // This function is normally called by the solve(...) function.
    // However, it must be called if solve will not be called (in particular,
    // if only the V-cycle is being used)
    virtual void setBottomSolver (const int a_lMax,
                                  const int a_lBase);

    //
    virtual void setBottomSolverEpsCushion (const Real a_bottomSolverEpsCushion);

    // Writes a_data to HDF5.
    virtual void outputAMR (const Vector<T*>& a_data,
                            const std::string a_name,
                            const int         a_lmax,
                            const int         a_lbase);

    // Public parameters...
    Real m_eps;
    Real m_hang;
    Real m_normThresh;
    int m_imin;
    int m_iterMax;
    int m_verbosity;
    int m_pre;
    int m_post;
    int m_bottom;
    int m_numMG;

    // max no. of coarsenings -- -1 (default) means coarsen as far as possible
    // If using a value besides the default, need to set it _before_ the define
    // function is called.
    int m_maxDepth;

    // default m_convergenceMetric = 0.:  initial residual will be set to
    // result of computeAMRResidual.
    // if m_convergenceMetric > 0., then initial residual will be set to
    // m_convergenceMetric.
    Real m_convergenceMetric;

    // used to give an additional cushion in the EPS used for bottom solves
    Real m_bottomSolverEpsCushion;

    bool m_solverParamsSet;
    int m_exitStatus;

protected:
    //
    virtual void relax (T&        a_phi,
                        const T&  a_R,
                        const int a_depth,
                        const int a_nRelax = 2);

    //
    virtual void computeAMRResidualLevel (Vector<T*>&       a_resid,
                                          const Vector<T*>& a_phi,
                                          const Vector<T*>& a_rhs,
                                          const int         a_lMax,
                                          const int         a_lBase,
                                          const int         a_ilev,
                                          const bool        a_homogeneousBC);

    // Returns the norm of the residual.
    virtual Real postVCycleOps (Vector<T*>&       a_uberResidual,
                                Vector<T*>&       a_uberCorrection,
                                Vector<T*>&       a_phi,
                                const Vector<T*>& a_rhs,
                                const int         a_lMax,
                                const int         a_lBase,
                                const bool        a_forceHomogeneous);

    //
    virtual void clear ();

    Vector<MappedAMRLevelOp<T>*> m_op;
    Vector<MappedMultiGrid<T>*>  m_mg;
    Vector<T*>                   m_correction;
    Vector<T*>                   m_residual;
    Vector<T*>                   m_resC;
    Vector<Copier>               m_resCopier;
    Vector<Copier>               m_reverseCopier;

    NoOpSolver<T>    m_nosolve;      // Used during initialization before a bottom solver is set.
    LinearSolver<T>* m_bottomSolver;

    Vector<char> m_hasInitBeenCalled;

private:
    // A list of inspectors maintained by this instance.
    Vector<RefCountedPtr<MappedAMRMultiGridInspector<T> > > m_inspectors;

    // Forbidden copiers.
    MappedAMRMultiGrid(const MappedAMRMultiGrid<T>&);
    MappedAMRMultiGrid& operator=(const MappedAMRMultiGrid<T>&);
};



// *********************** AMRMultigrid Implementation *************************


// -----------------------------------------------------------------------------
// Write a_data to HDF5.
// -----------------------------------------------------------------------------
template <class T>
void MappedAMRMultiGrid<T>::outputAMR (const Vector<T*>& a_data,
                                       const std::string a_name,
                                       const int         a_lmax,
                                       const int         a_lbase)
{
    Vector<T*> outputData;
    for (int ilev = a_lbase; ilev <= a_lmax; ilev++) {
        outputData.push_back(a_data[ilev]);
    }
    m_op[a_lbase]->outputAMR(outputData, a_name);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
void MappedAMRMultiGrid<T>::setSolverParameters (const int  a_pre,
                                                 const int  a_post,
                                                 const int  a_bottom,
                                                 const int  a_numMG,
                                                 const int  a_iterMax,
                                                 const Real a_eps,
                                                 const Real a_hang,
                                                 const Real a_normThresh)
{
    m_solverParamsSet = true;
    m_pre        =    a_pre;
    m_post       =    a_post;
    m_bottom     =    a_bottom;
    m_eps        =    a_eps;
    m_hang       =    a_hang;
    m_normThresh =    a_normThresh;
    m_iterMax    =    a_iterMax;
    for (int img = 0; img < m_mg.size(); img++) {
        m_mg[img]->m_pre    = a_pre;
        m_mg[img]->m_post   = a_post;
        m_mg[img]->m_bottom = a_bottom;
    }
    this->setMGCycle(a_numMG);
    m_bottomSolverEpsCushion = 1.0;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
Vector<MappedMGLevelOp<T>*> MappedAMRMultiGrid<T>::getAllOperators ()
{
    Vector<MappedMGLevelOp<T>*> retval;
    for (int iop = 0;  iop < m_op.size(); iop++) {
        MappedMGLevelOp<T>* operPtr = (MappedMGLevelOp<T>*) m_op[iop];
        retval.push_back(operPtr);
    }

    for (int img = 0; img < m_mg.size(); img++) {
        Vector<MappedMGLevelOp<T>*> mgOps = m_mg[img]->getAllOperators();
        retval.append(mgOps);
    }
    return retval;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
Vector<MappedMGLevelOp<T>*> MappedAMRMultiGrid<T>::getOperatorsOp ()
{
    Vector<MappedMGLevelOp<T>*> retval;
    for (int iop = 0;  iop < m_op.size(); iop++) {
        MappedMGLevelOp<T>* operPtr = (MappedMGLevelOp<T>*) m_op[iop];
        retval.push_back(operPtr);
    }
    return retval;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
Vector<Vector<MappedMGLevelOp<T>*> > MappedAMRMultiGrid<T>::getOperatorsMG()
{
    Vector<Vector<MappedMGLevelOp<T>*> > retval(m_mg.size());

    for (int img = 0; img < m_mg.size(); img++) {
        retval[img] = m_mg[img]->getAllOperators();
    }
    return retval;
}


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
template <class T>
MappedAMRMultiGrid<T>::MappedAMRMultiGrid ()
: m_eps(1E-6),
  m_hang(1E-15),
  m_normThresh(1E-30),
  m_imin(5),
  m_iterMax(20),
  m_verbosity(3),
  m_pre(2),
  m_post(2),
  m_bottom(2),
  m_numMG(1),
  m_maxDepth(-1),
  m_convergenceMetric(0.),
  m_bottomSolverEpsCushion(1.0),
  m_solverParamsSet(false),
  m_bottomSolver(NULL),
  m_inspectors()
{;}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
void MappedAMRMultiGrid<T>::setMGCycle (const int a_numMG)
{
    for (int ilev = 0; ilev < m_op.size(); ilev++) {
        // m_mg[ilev]->m_numMG = a_numMG;
        m_mg[ilev]->m_cycle = a_numMG;
    }
    m_numMG = a_numMG;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
void MappedAMRMultiGrid<T>::relax (T&        a_correction,
                                   const T&  a_residual,
                                   const int a_depth,   // this is actually ilev
                                   const int a_numSmooth)
{
    CH_TIME("MappedAMRMultiGrid::relax");

    // Get the AMR ref ratio to the coarser level.
    const IntVect& ref = m_op[a_depth]->refToCoarser();

    if (D_TERM(ref[0] > 2, || ref[1] > 2, || ref[2] > 2)) {
        // intermediate multigrid levels exist between this level and the
        // next coarser AMR level, so do a mini v-cyle to smooth on those
        CH_assert(m_mg[a_depth]->m_maxForcedDepth > 0);

        const int tmp = m_mg[a_depth]->m_depth;
        m_mg[a_depth]->m_depth = m_mg[a_depth]->m_maxForcedDepth + 1;
        if (m_verbosity >= 5) {
            pout() << "cycle at level " << a_depth << " with mini V-cycle max depth = "
                   << m_mg[a_depth]->m_maxForcedDepth << endl;
        }
        m_mg[a_depth]->cycle(0, a_correction, a_residual);
        m_mg[a_depth]->m_depth = tmp;

    } else {
        // No intermediate multigrid levels between AMR levels, so
        // just call relax on this level and be done with it.
        if (m_verbosity >= 5) pout() << "relax at level " << a_depth << endl;
        m_op[a_depth]->relax(a_correction, a_residual, a_numSmooth);  //numSmoothDown
    }
}



// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
template <class T>
MappedAMRMultiGrid<T>::~MappedAMRMultiGrid()
{
    CH_TIME("~MappedAMRMultiGrid");
    this->clear();
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
MappedAMRLevelOp<T>&  MappedAMRMultiGrid<T>::levelOp (const int a_level)
{
    return *(m_op[a_level]);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
Real MappedAMRMultiGrid<T>::computeAMRResidual (Vector<T*>&       a_resid,
                                                const Vector<T*>& a_phi,
                                                const Vector<T*>& a_rhs,
                                                const int         l_max,
                                                const int         l_base,
                                                const bool        a_homogeneousBC,
                                                const bool        a_computeNorm)
{
    CH_TIME("MappedAMRMultiGrid::computeAMRResidual");

    Real rnorm = 0;
    Real localNorm = 0;
    for (int ilev = l_base; ilev <= l_max; ilev++) {
        //always used at top level where bcs are inhomogeneous
        computeAMRResidualLevel(a_resid,
                                a_phi,
                                a_rhs,
                                l_max, l_base, ilev, a_homogeneousBC);
        if (a_computeNorm) {
            if (ilev == l_max) {
                localNorm = m_op[ilev]->localMaxNorm(*a_resid[ilev]);
            } else {
                m_op[ilev]->zeroCovered(*a_resid[ilev], *m_resC[ilev + 1], m_resCopier[ilev + 1]);
                localNorm = m_op[ilev]->localMaxNorm(*a_resid[ilev]);
            }
            rnorm = Max(localNorm, rnorm);
        }
    }
#ifdef CH_MPI
    if (a_computeNorm) {
        CH_TIME("MPI_Allreduce");
        Real recv;
        int result = MPI_Allreduce(&rnorm, &recv, 1, MPI_CH_REAL,
                                   MPI_MAX, Chombo_MPI::comm);
        if (result != MPI_SUCCESS) {
            //bark!!!
            MayDay::Error("sorry, but I had a communcation error on norm");
        }
        rnorm = recv;
    }
#endif

    return rnorm; // if a_computeNorm is false, then this just returns zero.
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
Real MappedAMRMultiGrid<T>::computeAMRResidual (Vector<T*>&       a_phi,
                                                const Vector<T*>& a_rhs,
                                                const int         l_max,
                                                const int         l_min)
{
    return computeAMRResidual(m_residual,
                              a_phi,
                              a_rhs,
                              l_max,
                              l_min,
                              false,
                              true);


}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
void MappedAMRMultiGrid<T>::computeAMROperator(Vector<T*>&       a_lph,
                                               const Vector<T*>& a_phi,
                                               const int         l_max,
                                               const int         l_base,
                                               const bool        a_homogeneousBC)
{
    for (int ilev = l_base; ilev <= l_max; ilev++) {
        m_op[ilev]->create(*m_residual[ilev], *a_lph[ilev]);
        m_op[ilev]->setToZero(*m_residual[ilev]);
    }
    computeAMRResidual(a_lph, a_phi, m_residual, l_max, l_base, a_homogeneousBC, false);
    // fixing the bug of not negating the result --- Qinghai Zhang 12/10/2009
    // In loving memory of the two weeks I lost over this bug!!!
    for (int ilev = l_base; ilev <= l_max; ilev++) {
        m_op[ilev]->scale(*a_lph[ilev], -1.);
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
void MappedAMRMultiGrid<T>::computeAMRResidualLevel (Vector<T*>&       a_resid,
                                                     const Vector<T*>& a_phi,
                                                     const Vector<T*>& a_rhs,
                                                     const int         l_max,
                                                     const int         l_base,
                                                     const int         ilev,
                                                     const bool        a_homogeneousBC)
{
    CH_TIME("MappedAMRMultiGrid<T>::computeAMRResidualLevel");
    CH_assert(m_hasInitBeenCalled[ilev] == 't');

    //m_op[ilev]->setToZero(*(a_resid[l_max]));
    if (l_max != l_base) {
        if (ilev == l_max) {
            m_op[l_max]->AMRResidualNF(*(a_resid[l_max]), *(a_phi[l_max]),
                                       *(a_phi[l_max - 1]), *(a_rhs[l_max]),
                                       a_homogeneousBC);
        } else if (ilev == l_base) {
            if (l_base == 0) {
                m_op[l_base]->AMRResidualNC(*(a_resid[l_base]), *(a_phi[l_base + 1]),
                                            *(a_phi[l_base]),  *(a_rhs[l_base]),
                                            a_homogeneousBC, m_op[l_base + 1]);
            } else {
                m_op[l_base]->AMRResidual(*a_resid[l_base], *a_phi[l_base + 1], *a_phi[l_base],
                                          *a_phi[l_base - 1], *a_rhs[l_base],
                                          a_homogeneousBC, m_op[l_base + 1]);
            }
        } else {
            m_op[ilev]->AMRResidual(*a_resid[ilev], *a_phi[ilev + 1], *a_phi[ilev],
                                    *a_phi[ilev - 1], *a_rhs[ilev],
                                    a_homogeneousBC, m_op[ilev + 1]);
        }
    } else {
        CH_assert(ilev == l_base);
        if (l_base == 0) {
            m_op[l_max]->residual(*a_resid[l_max], *a_phi[l_max], *a_rhs[l_max], a_homogeneousBC);
        } else {
            m_op[l_max]->AMRResidualNF(*(a_resid[l_max]), *(a_phi[l_max]),
                                       *(a_phi[l_max - 1]), *(a_rhs[l_max]),
                                       a_homogeneousBC);
        }
    }

}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template<class T>
void MappedAMRMultiGrid<T>::solve (Vector<T*>&       a_phi,
                                   const Vector<T*>& a_rhs,
                                   const int         l_max,
                                   const int         l_base,
                                   const bool        a_zeroPhi,
                                   const bool        a_forceHomogeneous)
{
    CH_TIME("MappedAMRMultiGrid::solve");

    init(a_phi, a_rhs, l_max, l_base);
    solveNoInit(a_phi, a_rhs, l_max, l_base, a_zeroPhi, a_forceHomogeneous);
    //puts multigrid depth back where it started so the solver can be reused
    revert(a_phi, a_rhs, l_max, l_base);
}


template<class T>
void MappedAMRMultiGrid<T>::solveNoInit (Vector<T*>&       a_phi,
                                         const Vector<T*>& a_rhs,
                                         const int         l_max,
                                         const int         l_base,
                                         const bool        a_zeroPhi,
                                         const bool        a_forceHomogeneous)
{
    Vector<T*> uberResidual(a_rhs.size());
    int lowlim = l_base;
    if (l_base > 0)  // we need an ubercorrection one level lower than l_base
        --lowlim;
    // for (int ilev = l_base; ilev <= l_max; ilev++)
    for (int ilev = lowlim; ilev <= l_max; ilev++) {
        uberResidual[ilev] = new T();
        if (ilev >= l_base) {
            m_op[ilev]->create(*uberResidual[ilev], *a_rhs[ilev]);
        }
    }
    solveNoInitResid(a_phi, uberResidual, a_rhs, l_max, l_base, a_zeroPhi, a_forceHomogeneous);
    for (int i = lowlim; i <= l_max; i++) {
        m_op[i]->clear(*uberResidual[i]);
        delete uberResidual[i];
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template<class T>
void MappedAMRMultiGrid<T>::solveNoInitResid (Vector<T*>&       a_phi,
                                              Vector<T*>&       uberResidual,
                                              const Vector<T*>& a_rhs,
                                              const int         l_max,
                                              const int         l_base,
                                              const bool        a_zeroPhi,
                                              const bool        a_forceHomogeneous)
{
    CH_TIMERS("MappedAMRMultiGrid::solveNoInit");
    CH_TIMER("MappedAMRMultiGrid::AMRVcycle", vtimer);
    setBottomSolver(l_max, l_base);

    CH_assert(l_base <= l_max);
    CH_assert(a_rhs.size() == a_phi.size());

    //these correspond to the residual and correction
    //that live in AMRSolver
    Vector<T*> uberCorrection(a_rhs.size());

    // This stores the best solution
    Vector<T*> bestPhi(a_phi.size());

    int lowlim = l_base;
    if (l_base > 0)  // we need an ubercorrection one level lower than l_base
        lowlim--;
    bool outputIntermediates = false;

    for (int ilev = lowlim; ilev <= l_max; ilev++) {
        uberCorrection[ilev] = new T();
        m_op[ilev]->create(*uberCorrection[ilev], *a_phi[ilev]);
        if (ilev >= l_base) {
            m_op[ilev]->create(*uberResidual[ilev], *a_rhs[ilev]);
            m_op[ilev]->setToZero(*(uberResidual[ilev]));
        }
        m_op[ilev]->setToZero(*(uberCorrection[ilev]));

        // Allocate bestPhi
        bestPhi[ilev] = new T();
        m_op[ilev]->create(*bestPhi[ilev], *a_phi[ilev]);
    }
    if (a_zeroPhi)
        for (int ilev = l_base; ilev <= l_max; ++ilev) {
            m_op[ilev]->setToZero(*(a_phi[ilev]));
        }
    //compute initial residual and initialize internal residual to it

    // Set bestPhi equal to initial guess
    for (int ilev = lowlim; ilev <= l_max; ilev++) {
        m_op[ilev]->assign(*bestPhi[ilev], *a_phi[ilev]);
    }

    Real initial_rnorm = 0;
    {
        CH_TIME("Initial AMR Residual");
        initial_rnorm = computeAMRResidual(uberResidual, a_phi, a_rhs, l_max, l_base, a_forceHomogeneous, true);
    }

    if (m_convergenceMetric != 0.) {
        initial_rnorm = m_convergenceMetric;
    }

    Real rnorm = initial_rnorm;
    Real norm_last = 2 * initial_rnorm;
    Real best_rnorm = rnorm;
    bool useBestPhi = false;
    bool somethingConverged = false;  // This remains false if no iters converge to a sol'n better than what was given to us.

    /// set bottom solver convergence norm and solver tolerance
    m_bottomSolver->setConvergenceMetrics(initial_rnorm, m_bottomSolverEpsCushion * m_eps);

    int iter = 0;
    if (m_verbosity == 2) {
        pout() << "    MappedAMRMultiGrid:: " << std::flush;
    } else if (m_verbosity > 2) {
        pout() << "    MappedAMRMultiGrid:: iteration = " << std::fixed << iter
               << ", residual norm = " << std::scientific << rnorm << std::endl;
    }

    bool goNorm = rnorm > m_normThresh;                        //iterate if norm is not small enough
    bool goRedu = rnorm > m_eps * initial_rnorm;               //iterate if initial norm is not reduced enough
    bool goIter = iter < m_iterMax;                            //iterate if iter < max iteration count
    bool goHang = iter < m_imin || rnorm < (1 - m_hang) * norm_last; //iterate if we didn't hang

    while (goIter && goRedu && goHang && goNorm) {
        // Report the residuals to the inspectors.
        for (int i = 0; i < m_inspectors.size(); ++i)
            m_inspectors[i]->recordResiduals(uberResidual, l_base, l_max, iter);

        if (outputIntermediates) {
            char strresname[100];
            sprintf(strresname, "amrmg.res.iter.%03d.hdf5", iter);
            string nameres(strresname);
            outputAMR(uberResidual, nameres, l_max, l_base);
        }

        norm_last = rnorm;

        //this generates a correction from the current residual
        CH_START(vtimer);
        AMRVCycle(uberCorrection, uberResidual, l_max, l_max, l_base);
        CH_STOP(vtimer);


        // Report the corrections to the inspectors.
        for (int i = 0; i < m_inspectors.size(); ++i)
            m_inspectors[i]->recordCorrections(uberCorrection, l_base, l_max, iter);

        // Do post VCycle stuff
        rnorm = postVCycleOps(uberResidual, uberCorrection, a_phi, a_rhs, l_max, l_base, a_forceHomogeneous);
        iter++;

        // Check if residual has dropped.
        if (rnorm <= best_rnorm) {
            // OK. Save this as the best solution and move on.
            best_rnorm = rnorm;
            for (int ilev = l_base; ilev <= l_max; ++ilev) {
                m_op[ilev]->assign(*bestPhi[ilev], *a_phi[ilev]);
            }
            useBestPhi = false;
            somethingConverged = true;
        } else {
            // Bad. VCycle did not produce a better solution.
            useBestPhi = true;
        }

        // Let the user know what happened
        if (m_verbosity >= 3) { ////
            if (useBestPhi) {
                pout() << "   [D] ";   // Mark iterations that have diverged.
            } else {
                pout() << "       ";
            }
            pout() << "MappedAMRMultiGrid:: iteration = " << std::fixed << iter
                   << ", residual norm = " << std::scientific << rnorm;
            if (rnorm > 0.0) {
                pout() << ", rate = " << std::fixed << norm_last / rnorm;
            }
            pout() << std::endl;
        }

        goNorm = rnorm > m_normThresh;                        //keep iterating if norm is not small enough
        goRedu = rnorm > m_eps * initial_rnorm;               //keep iterating if initial norm is not reduced enough
        goIter = iter < m_iterMax;                            //keep iterating if iter < max iteration count
        goHang = iter < m_imin || rnorm < (1 - m_hang) * norm_last; //keep iterating if we didn't hang
    }

    // Use the best solution available
    if (useBestPhi) {
        rnorm = best_rnorm;
        for (int ilev = l_base; ilev <= l_max; ++ilev) {
            m_op[ilev]->assign(*a_phi[ilev], *bestPhi[ilev]);
        }
    }

    // Did we blow up?
    if ((rnorm > 10.*initial_rnorm) && (rnorm > 10.*m_eps)) {
        pout() << "solver seems to have blown up" << endl;
        MayDay::Error("kaboom");

        // TODO: Maybe we shouldn't just give up.
    }

    if (!somethingConverged && rnorm >= initial_rnorm && rnorm >= m_eps) {
        pout() << "MappedAMRMultiGrid solver blew up." << endl;
        // std::cout << "MappedAMRMultiGrid solver blew up." << endl;
        MayDay::Error("MappedAMRMultiGrid solver blew up");
    }

    // The solver has finished. Figure out the final state of the solution.
    m_exitStatus = int(!goRedu) + int(!goIter) * 2 + int(!goHang) * 4 + int(!goNorm) * 8;

    if (m_verbosity == 2) {
        if (initial_rnorm != 0.0) {
            pout() << std::fixed << iter << " iters\t rel res = "
                   << std::scientific << rnorm/initial_rnorm << std::endl;
        } else {
            pout() << std::fixed << iter << " iters\t rel res = "
                   << std::scientific << rnorm << " / " << initial_rnorm << std::endl;
        }
    } else if (m_verbosity > 2) {
        pout() << "    MappedAMRMultiGrid:: iteration = " << std::fixed << iter
               << ", residual norm = " << std::scientific << rnorm << std::endl;
    }

    if (m_verbosity > 1) {
        if (!goIter && goRedu && goNorm) { // goRedu=T, goIter=F, goHang=?, goNorm=T
            // m_exitStatus == 0 + 2 + 0|4 + 0 = 2|6
            pout() << "    MappedAMRMultiGrid:: WARNING: Exit because max iteration count exceeded" << std::endl;
        }
        if (!goHang && goRedu && goNorm) { // goRedu=T, goIter=?, goHang=F, goNorm=T
            // m_exitStatus == 0 + 0|2 + 4 + 0 = 4|6
            pout() << "    MappedAMRMultiGrid:: WARNING: Exit because of solver hang" << std::endl;
        }
        if (m_verbosity > 4) {
            pout() << "    MappedAMRMultiGrid:: exitStatus = " << std::fixed << m_exitStatus << std::endl;
        }
    }

    // Clean up after ourselves
    for (int i = lowlim; i <= l_max; i++) {
        m_op[i]->clear(*uberCorrection[i]);
        delete uberCorrection[i];
        delete bestPhi[i];
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template<class T>
Real MappedAMRMultiGrid<T>::postVCycleOps (Vector<T*>&       a_uberResidual,
                                           Vector<T*>&       a_uberCorrection,
                                           Vector<T*>&       a_phi,
                                           const Vector<T*>& a_rhs,
                                           const int         l_max,
                                           const int         l_base,
                                           const bool        a_forceHomogeneous)
{
    // Increment phi by correction and reset correction to zero
    for (int ilev = l_base; ilev <= l_max; ilev++) {
        m_op[ilev]->incr(*(a_phi[ilev]), *(a_uberCorrection[ilev]), 1.0);
        m_op[ilev]->setToZero(*(a_uberCorrection[ilev]));
    }

    // For solvers with accuracy higher than 2nd order
    //  consistency between levels has to be explicitly enforced.
    // Qinghai Zhang
    if (m_op[0]->orderOfAccuracy() > 2) {
        for (int ilev = l_max; ilev > l_base; ilev--) {
            m_op[ilev]->enforceCFConsistency(*a_phi[ilev - 1], *a_phi[ilev]);
        }
    }
    //------end enforcing consistency.

    // recompute residual
    return computeAMRResidual(a_uberResidual, a_phi, a_rhs, l_max, l_base, a_forceHomogeneous, true);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template<class T>
void MappedAMRMultiGrid<T>::relaxOnlyHomogeneous (Vector<T*>&       a_phi,
                                                  const Vector<T*>& a_rhs,
                                                  const int         l_max,
                                                  const int         l_base)
{
    CH_TIME("MappedAMRMultiGrid::relaxOnly");
    CH_assert(l_max == l_base);
    CH_assert(l_max == 0);
    init(a_phi, a_rhs, l_max, l_base);

    CH_assert(a_rhs.size() == a_phi.size());

    Vector<T*> uberResidual(a_rhs.size());

    for (int ilev = 0; ilev < m_op.size(); ilev++) {
        uberResidual[ilev] = new T();
        m_op[ilev]->create(*uberResidual[ilev], *a_rhs[ilev]);
    }

    //compute initial residual and initialize internal residual to it
    m_op[0]->residual(*uberResidual[0], *a_phi[0], *a_rhs[0], true);
    Real initial_rnorm = m_op[0]->norm(*uberResidual[0], 0);
    Real rnorm = initial_rnorm;
    Real norm_last = 2 * initial_rnorm;


    int iter = 0;
    if (m_verbosity >= 3) {
        pout() << "    MappedAMRMultiGrid::relaxOnly iteration = "
               << std::fixed << iter
               << ", residual norm = " << std::scientific << rnorm
               << std::endl;
    }
    bool goNorm = rnorm > m_normThresh;                        //iterate if norm is not small enough
    bool goRedu = rnorm > m_eps * initial_rnorm;               //iterate if initial norm is not reduced enough
    bool goIter = iter < m_iterMax;                            //iterate if iter < max iteration count
    bool goHang = iter < m_imin || rnorm < (1 - m_hang) * norm_last; //iterate if we didn't hang
    while (goIter && goRedu && goHang && goNorm) {
        norm_last = rnorm;
        m_op[0]->relax(*a_phi[0], *a_rhs[0], 1);

        iter++;
        //recompute residual
        m_op[0]->residual(*uberResidual[0], *a_phi[0], *a_rhs[0], true);
        rnorm = m_op[0]->norm(*uberResidual[0], 2);
        if (m_verbosity >= 4) {
            pout() << "    MappedAMRMultiGrid::relaxOnly iteration = "
                   << std::fixed << iter
                   << ", residual norm = " << std::scientific << rnorm
                   << ", rate = " << std::fixed << norm_last / rnorm
                   << std::endl;
        }
        goNorm = rnorm > m_normThresh;                        //keep iterating if norm is not small enough
        goRedu = rnorm > m_eps * initial_rnorm;               //keep iterating if initial norm is not reduced enough
        goIter = iter < m_iterMax;                            //keep iterating if iter < max iteration count
        goHang = iter < m_imin || rnorm < (1 - m_hang) * norm_last; //keep iterating if we didn't hang
    }
    m_exitStatus = 0;
    for (int i = 0; i < m_op.size(); i++) {
        m_op[i]->clear(*uberResidual[i]);
        delete uberResidual[i];
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
void MappedAMRMultiGrid<T>::clear ()
{
    CH_assert(m_correction.size() == m_op.size());
    CH_assert(m_residual.size() == m_op.size());
    CH_assert(m_resC.size() == m_op.size());
    CH_assert(m_mg.size() == m_op.size());

    for (int i = 0; i < m_op.size(); i++) {
        m_op[i]->clear(*m_correction[i]);
        m_op[i]->clear(*m_residual[i]);
        m_op[i]->clear(*m_resC[i]);

        delete m_correction[i];
        delete m_residual[i];
        delete m_resC[i];
        delete m_op[i];
        delete m_mg[i];

        m_correction[i] = NULL;
        m_residual[i] = NULL;
        m_resC[i] = NULL;
        m_op[i] = NULL;
        m_mg[i] = NULL;
    }
    m_hasInitBeenCalled.resize(0);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
void MappedAMRMultiGrid<T>::init (const Vector<T*>& a_phi,
                                  const Vector<T*>& a_rhs,
                                  const int         l_max,
                                  const int         l_base)
{
    CH_TIME("MappedAMRMultiGrid::init");

    for (int i = l_base; i <= l_max; ++i) {
        m_hasInitBeenCalled[i] = 't';
        MappedAMRLevelOp<T>& op = *(m_op[i]);

        op.create(*m_correction[i], *a_phi[i]);
        op.create(*m_residual[i],   *a_rhs[i]);

        // This triggers the mini V-cycles when the AMR ref ratio > 2.
        if (i != l_base) {
            if (m_mg[i]->m_maxForcedDepth > 0) {
                m_mg[i]->m_depth = m_mg[i]->m_maxForcedDepth + 1;
            } else {
                m_mg[i]->m_depth = 0;
            }
        }

        m_mg[i]->init(*a_phi[i], *a_rhs[i]);

        if (i != l_base) {
            const IntVect r = op.refToCoarser();
            op.createCoarsened(*m_resC[i], *a_rhs[i], r);
            op.buildCopier(m_resCopier[i], *a_rhs[i-1], *m_resC[i]);
            m_reverseCopier[i] = m_resCopier[i];
            m_reverseCopier[i].reverse();
        }
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
void MappedAMRMultiGrid<T>::revert (const Vector<T*>& a_phi,
                                    const Vector<T*>& a_rhs,
                                    const int         l_max,
                                    const int         l_base)
{
    CH_TIME("MappedAMRMultiGrid::revert");
    for (int i = l_base; i <= l_max; i++) {
        if (i != l_base) {
            m_mg[i]->m_depth = m_mg[i]->m_defaultDepth;
        }
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
void MappedAMRMultiGrid<T>::setBottomSolver (const int l_max,
                                             const int l_base)
{
    for (int ilev = l_base; ilev <= l_max; ilev++) {
        // only do this if we haven't already set these
        if (!m_solverParamsSet) {
            m_mg[ilev]->m_pre = m_pre;
            m_mg[ilev]->m_post = m_post;
            m_mg[ilev]->m_bottom = m_bottom;
        }
        m_mg[ilev]->m_bottomSolver = &m_nosolve;
    }

    m_mg[l_base]->setBottomSolver(m_bottomSolver);
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
void MappedAMRMultiGrid<T>::setBottomSolverEpsCushion (Real a_bottomSolverEpsCushion)
{
    CH_assert((a_bottomSolverEpsCushion >  0.) &&
              (a_bottomSolverEpsCushion <= 1.));
    m_bottomSolverEpsCushion = a_bottomSolverEpsCushion;
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template <class T>
void MappedAMRMultiGrid<T>::define (const ProblemDomain&        a_coarseDomain,
                                    MappedAMRLevelOpFactory<T>& a_factory,
                                    LinearSolver<T>*            a_bottomSolver,
                                    const int                   a_maxAMRLevels,
                                    const int                   a_verbosity)
{
    CH_TIME("MappedAMRMultiGrid::define");

    m_verbosity = a_verbosity;

    this->clear();
    m_op.resize( a_maxAMRLevels, NULL);
    m_mg.resize( a_maxAMRLevels, NULL);
    m_hasInitBeenCalled.resize(a_maxAMRLevels, 'f');

    m_correction.resize( a_maxAMRLevels, NULL);
    m_residual.  resize( a_maxAMRLevels, NULL);
    m_resC.      resize( a_maxAMRLevels, NULL);
    m_resCopier. resize( a_maxAMRLevels);
    m_reverseCopier.resize(  a_maxAMRLevels);
    m_bottomSolver = a_bottomSolver;

    ProblemDomain current = a_coarseDomain;
    for (int i = 0; i < a_maxAMRLevels; i++) {
        CH_assert(m_correction[i] == NULL);
        CH_assert(m_residual[i] == NULL);
        CH_assert(m_resC[i] == NULL);

        m_correction[i] = new T();
        m_residual[i]   = new T();
        m_resC[i]       = new T();

        if (m_verbosity >= 5) {
            pout() << "calling a_factory.AMRnewOp(current) with lev = "
                   << std::fixed << i
                   << " and index space = " << current.domainBox() << endl;
        }

        CH_assert(m_op[i] == NULL);
        m_op[i] = a_factory.AMRnewOp(current);
        CH_assert(m_op[i] != NULL);

        m_mg[i] = new MappedMultiGrid<T>();
        if (i == 0) {
            // Let the MG define function decide its own coarsening pattern.
            m_mg[i]->define(a_factory, &m_nosolve, current, m_maxDepth, m_op[i]);
            m_mg[i]->m_verbosity = m_verbosity;
        } else {
            // Mini V-cycle will need a special coarsening pattern.

            // Determine how much refinement the mini V-cycle will need.
            Vector<IntVect> allMGRefRatios(0);
            int intermediateDepth = 0;
            IntVect r = m_op[i]->refToCoarser();

            while (D_TERM(r[0] > 2, || r[1] > 2, || r[2] > 2)) {
                IntVect thisMGref = IntVect::Unit;
                D_TERM(if (r[0] > 2) {r[0] /= 2; thisMGref[0] = 2;},
                       if (r[1] > 2) {r[1] /= 2; thisMGref[1] = 2;},
                       if (r[2] > 2) {r[2] /= 2; thisMGref[2] = 2;})
                if (thisMGref.product() > 1) {
                    allMGRefRatios.push_back(thisMGref);
                    intermediateDepth++;
                }
            }

            // Flip the ref ratios so anisotropic refinement occurs first.
            Vector<IntVect> allMGRefRatiosReverse(intermediateDepth);
            for (int idx = 0; idx < intermediateDepth; ++idx) {
                allMGRefRatiosReverse[intermediateDepth - idx - 1] = allMGRefRatios[idx];
            }
            allMGRefRatios.resize(0);

            // This also sets m_mg[i]->m_maxForcedDepth.
            m_mg[i]->define(a_factory, &m_nosolve, current, m_maxDepth, m_op[i], &allMGRefRatiosReverse);
            m_mg[i]->m_verbosity = m_verbosity;
        }

        // Only do this if it will be used (avoiding a reference to invalid
        // and/or unavailable refinement ratios)
        if (i < a_maxAMRLevels - 1) {
            ProblemDomain notSoCurrent = current;
            refine(current, notSoCurrent, a_factory.getFineRefRatio(current));
        }
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
template<class T>
void MappedAMRMultiGrid<T>::AMRVCycle (Vector<T*>&       a_uberCorrection,
                                       const Vector<T*>& a_uberResidual,
                                       const int         ilev,
                                       const int         l_max,
                                       const int         l_base)
{
    if (ilev == l_max) {
        for (int level = l_base; level <= l_max; level++) {
            m_op[level]->assignLocal(*m_residual[level], *a_uberResidual[level]);
            m_op[level]->setToZero(*m_correction[level]);
        }
    }

    if (l_max == l_base) {
        CH_assert(ilev == l_base);
        if (m_verbosity >= 5) pout() << "oneCycle " << std::fixed << ilev << "\n";
        m_mg[l_base]->oneCycle(*(a_uberCorrection[ilev]), *(a_uberResidual[ilev]));
    } else if (ilev == l_base) {
        if (m_verbosity >= 5) pout() << "oneCycle " << std::fixed << ilev << "\n";
        m_mg[l_base]->oneCycle(*(m_correction[ilev]), *(m_residual[ilev]));

        m_op[ilev]->incr(*(a_uberCorrection[ilev]), *(m_correction[ilev]), 1.0);
    } else {
        //============= Downsweep ========================

        if (m_verbosity >= 5) pout() << "Calling relax on level " << std::fixed << ilev << "...\n";
        this->relax(*(m_correction[ilev]), *(m_residual[ilev]), ilev, m_pre);
        m_op[ilev]->incr(*(a_uberCorrection[ilev]), *(m_correction[ilev]), 1.0);

        // Set next coarser level correction to zero
        m_op[ilev - 1]->setToZero(*(m_correction[ilev - 1]));

        // Recompute residual on next coarser level
        //  for the valid region NOT covered by this level.
        computeAMRResidualLevel(m_residual,
                                a_uberCorrection,
                                a_uberResidual,
                                l_max, l_base, ilev - 1,
                                true);

        // Compute the restriction of the residual to the coarser level resC.
        if (m_verbosity >= 5) pout() << "restrict residual " << std::fixed << ilev << " to " << ilev - 1 << "\n";
        m_op[ilev]->AMRRestrictS(*(m_resC[ilev]),
                                 *(m_residual[ilev]),
                                 *(m_correction[ilev]),
                                 *(m_correction[ilev - 1]),
                                 *(a_uberCorrection[ilev]));

        // Overwrite residual on the valid region of the next coarser level
        //  with coarsened residual from this level
        m_op[ilev - 1]->assignCopier(*m_residual[ilev - 1], *(m_resC[ilev]), m_resCopier[ilev]);

        //============finish Compute residual for the next coarser level======

        for (int img = 0; img < m_numMG; img++) {
            AMRVCycle(a_uberCorrection, a_uberResidual, ilev - 1, l_max, l_base);
        }

        //================= Upsweep ======================
        //increment the correction with coarser version
        if (m_verbosity >= 5) pout() << "Prolong correction " << std::fixed << ilev - 1 << " to " << ilev << "\n";
        m_op[ilev]->AMRProlongS(*(m_correction[ilev]), *(m_correction[ilev - 1]),
                                *m_resC[ilev], m_reverseCopier[ilev]);

        //recompute residual
        m_op[ilev]->AMRUpdateResidual(*(m_residual[ilev]), *(m_correction[ilev]), *(m_correction[ilev - 1]));

        //compute correction to the correction
        if (m_verbosity >= 5) pout() << "Calling relax on level " << std::fixed << ilev << "...\n";
        T& dCorr = *(a_uberCorrection[ilev]); // user uberCorrection as holder for correction to correction
        m_op[ilev]->setToZero(dCorr);
        this->relax(dCorr, *(m_residual[ilev]), ilev, m_post);

        //correct the correction with the correction to the correction
        if (m_verbosity >= 5) {
            T newRes;
            m_op[ilev]->create(newRes, *(m_correction[ilev]));
            m_op[ilev]->residual(newRes, *(m_correction[ilev]), *(m_residual[ilev]), true);

            T fineRes;
            IntVect ref = IntVect::Unit;
            Real rnorm = m_op[ilev]->AMRNorm(newRes, fineRes, ref, 0);

            pout() << "Pre residual max norm = " << std::scientific << rnorm << endl;
        }
        m_op[ilev]->incr(*(m_correction[ilev]), dCorr, 1.0);
        if (m_verbosity >= 5) {
            T newRes;
            m_op[ilev]->create(newRes, *(m_correction[ilev]));
            m_op[ilev]->residual(newRes, *(m_correction[ilev]), *(m_residual[ilev]), true);

            T fineRes;
            IntVect ref = IntVect::Unit;
            Real rnorm = m_op[ilev]->AMRNorm(newRes, fineRes, ref, 0);

            pout() << "Post residual max norm = " << std::scientific << rnorm << endl;
        }

        m_op[ilev]->assignLocal(*(a_uberCorrection[ilev]), *(m_correction[ilev]));
    }
}


#endif //!__MappedAMRMultiGrid_H__INCLUDED__
