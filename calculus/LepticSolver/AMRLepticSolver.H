/*******************************************************************************
 *    FILE: AMRLepticSolver.H
 *    DESC: The AMR version of the leptic solver, as the file name suggests.
 *    DATE: Sun 18 May 2014 12:49:21 PM EDT  (UTC-0400)
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef AMRLepticSolver_H__INCLUDED__
#define AMRLepticSolver_H__INCLUDED__

#include "AMREllipticSolver.H"
#include "LevelLepticSolver.H"


class AMRLepticSolver: public AMREllipticSolver<LevelData<FArrayBox> >
{
public:
    // Constructor
    AMRLepticSolver ();

    // Destructor
    virtual ~AMRLepticSolver ();

    // TODO: Change params.
    void setSolverParameters (const int  a_numMG,
                              const int  a_iterMax,
                              const Real a_eps,
                              const Real a_hang,
                              const Real a_normThresh);

    void define (const ProblemDomain&                            a_coarseDomain,
                 MappedAMRLevelOpFactory<LevelData<FArrayBox> >& a_factory,
                 const int                                       a_maxAMRLevels,
                 const int                                       a_verbosity = 3);

    virtual void clear ();

    virtual void init (const Vector<LevelData<FArrayBox>*>& a_phi,
                       const Vector<LevelData<FArrayBox>*>& a_rhs,
                       const int                            l_max,
                       const int                            l_base);

    virtual void solve (Vector<LevelData<FArrayBox>*>&       a_phi,
                        const Vector<LevelData<FArrayBox>*>& a_rhs,
                        const int                            l_max,
                        const int                            l_base,
                        const bool                           a_zeroPhi,
                        const bool                           a_forceHomogeneous);

    virtual void solveNoInit (Vector<LevelData<FArrayBox>*>&       a_phi,
                              const Vector<LevelData<FArrayBox>*>& a_rhs,
                              const int                            l_max,
                              const int                            l_base,
                              const bool                           a_zeroPhi,
                              const bool                           a_forceHomogeneous);

    virtual void solveNoInitResid (Vector<LevelData<FArrayBox>*>&       a_phi,
                                   Vector<LevelData<FArrayBox>*>&       uberResidual,
                                   const Vector<LevelData<FArrayBox>*>& a_rhs,
                                   const int                            l_max,
                                   const int                            l_base,
                                   const bool                           a_zeroPhi,
                                   const bool                           a_forceHomogeneous);

    virtual void AMRVCycle (Vector<LevelData<FArrayBox>*>&       a_correction,
                            const Vector<LevelData<FArrayBox>*>& a_residual,
                            const int                            a_l,
                            const int                            a_lMax,
                            const int                            a_lBase);

    virtual Real postVCycleOps (Vector<LevelData<FArrayBox>*>&       a_uberResidual,
                                Vector<LevelData<FArrayBox>*>&       a_uberCorrection,
                                Vector<LevelData<FArrayBox>*>&       a_phi,
                                const Vector<LevelData<FArrayBox>*>& a_rhs,
                                const int                            l_max,
                                const int                            l_base,
                                const bool                           a_forceHomogeneous);

    virtual void computeAMRResidualLevel (Vector<LevelData<FArrayBox>*>&       a_resid,
                                          const Vector<LevelData<FArrayBox>*>& a_phi,
                                          const Vector<LevelData<FArrayBox>*>& a_rhs,
                                          const int                            l_max,
                                          const int                            l_base,
                                          const int                            ilev,
                                          const bool                           a_homogeneousBC);

    virtual Real computeAMRResidual (Vector<LevelData<FArrayBox>*>&       a_resid,
                                     const Vector<LevelData<FArrayBox>*>& a_phi,
                                     const Vector<LevelData<FArrayBox>*>& a_rhs,
                                     const int                            l_max,
                                     const int                            l_base,
                                     const bool                           a_homogeneousBC,
                                     const bool                           a_computeNorm);

    // Write a_data to HDF5.
    virtual void outputAMR (const Vector<LevelData<FArrayBox>*>& a_data,
                            const std::string                    a_name,
                            const int                            a_lmax,
                            const int                            a_lbase);

    // Public parameters...
    Real m_eps;
    Real m_hang;
    Real m_normThresh;
    int m_imin;
    int m_iterMax;
    int m_exitStatus;
    int m_numMG;
    int m_verbosity;

    // default m_convergenceMetric = 0.:  initial residual will be set to
    // result of computeAMRResidual.
    // if m_convergenceMetric > 0., then initial residual will be set to
    // m_convergenceMetric.
    Real m_convergenceMetric;

protected:
    Vector<LevelLepticSolver*>                       m_amrLepticSolver;
    Vector<MappedAMRLevelOp<LevelData<FArrayBox> >*> m_op;
    Vector<LevelData<FArrayBox>*>                    m_correction;
    Vector<LevelData<FArrayBox>*>                    m_residual;
    Vector<LevelData<FArrayBox>*>                    m_resC;
    Vector<Copier>                                   m_resCopier;
    Vector<Copier>                                   m_reverseCopier;

private:
    // Forbidden copiers.
    AMRLepticSolver(const AMRLepticSolver&);
    AMRLepticSolver& operator=(const AMRLepticSolver&);
};



#endif //!AMRLepticSolver_H__INCLUDED__
