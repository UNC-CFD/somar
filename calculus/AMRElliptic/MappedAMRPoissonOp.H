/*******************************************************************************
 *  FILE: MappedAMRPoissonOp.H
 *  DESC: Operator for solving (alpha I + beta*Laplacian)(phi) = rho over an
 *        AMR hierarchy in mapped coordiantes.
 *  DATE: Wed 02 May 2012 09:05:58 PM EDT
 *  MAIL: santilli@physics.unc.edu
 ******************************************************************************/
#ifndef __MAPPEDAMRPOISSONOP_H__INCLUDED__
#define __MAPPEDAMRPOISSONOP_H__INCLUDED__

#include "RestrictionStrategy.H"
#include "ProlongationStrategy.H"
#include "RelaxationMethod.H"
#include "LevelData.H"
#include "LevelDataOps.H"
#include "MappedAMRTGA.H"
#include "LepticOperator.H"
#include "EllipticBCInterface.H"
#include "MappedQuadCFInterp.H"
#include "MappedLevelFluxRegister.H"
#include "CFRegion.H"
#include "FluxBox.H"
#include "GNUC_Extensions.H"
class MappedAMRPoissonOpFactory;
class RealVect;


// -----------------------------------------------------------------------------
// Operator for solving (alpha I + beta*Laplacian)(phi) = rho
// over an AMR hierarchy in mapped coordinates.
// -----------------------------------------------------------------------------
class MappedAMRPoissonOp : public MappedLevelTGAHelmOp<LevelData<FArrayBox>, FluxBox>,
                           public LepticOperator
{
    friend class MappedAMRPoissonOpFactory;
    friend void test_LinearOp ();
    friend void test_HorizontalLinearOp ();

public:
    // Default constructor
    MappedAMRPoissonOp ()
    : m_isDiagonal(false),   // Factory sets this. I just wanted to shut valgrind up.
      m_horizontalOp(false),
      m_restrictPtr(NULL),
      m_prolongPtr(NULL),
      m_relaxPtr(NULL),
      m_precondRelaxPtr(NULL)
    {;}

    // Destructor
    virtual ~MappedAMRPoissonOp ()
    {
        delete m_restrictPtr;
        m_restrictPtr = NULL;

        delete m_prolongPtr;
        m_prolongPtr = NULL;

        if (m_relaxPtr == m_precondRelaxPtr) {
            delete m_relaxPtr;
            m_relaxPtr = NULL;
            m_precondRelaxPtr = NULL;
        } else {
            delete m_relaxPtr;
            m_relaxPtr = NULL;

            delete m_precondRelaxPtr;
            m_precondRelaxPtr = NULL;
        }
    }

    // Defines this op to simply calculate Laplacians and residuals.
    void define(BCMethodHolder&                             a_bc,
                const LevelGeometry&                        a_levGeo,
                const Real                                  a_alpha,
                const Real                                  a_beta,
                const RefCountedPtr<LevelData<FArrayBox> >& a_lapDiagsPtr = RefCountedPtr<LevelData<FArrayBox> >(NULL));

    // Creates a simple LinearOp for a LinearSolver.
    // If you have a coarser level, this will only handle homogeneous CFBCs.
    // Take care of your CFBCs, then solve the residual equation.
    // If a_JinvPtr.isNull(), this op will just use Jinv=1.
    void linearOpDefine (BCMethodHolder&                             a_bc,
                         const RealVect&                             a_dx,
                         const RealVect&                             a_dxCrse,
                         const Real                                  a_alpha,
                         const Real                                  a_beta,
                         const bool                                  a_isDiagonal,
                         const bool                                  a_horizontalOp,
                         const RefCountedPtr<LevelData<FluxBox> >&   a_JgupPtr,
                         const RefCountedPtr<LevelData<FArrayBox> >& a_JinvPtr = RefCountedPtr<LevelData<FArrayBox> >(NULL),
                         const int                                   a_precondMode = ProblemContext::PrecondMode::DiagRelax);

protected:
    // Full define function for operator with both coarser and finer levels
    void define(const DisjointBoxLayout&   a_grids,
                const DisjointBoxLayout&   a_gridsFiner,
                const DisjointBoxLayout&   a_gridsCoarser,
                const RealVect&            a_dxLevel,
                const IntVect&             a_refRatio,
                const IntVect&             a_refRatioFiner,
                const ProblemDomain&       a_domain,
                BCMethodHolder&            a_bc,
                const Copier&              a_exchange,
                const CFRegion&            a_cfregion);

    // Full define function for operator with finer levels, but no coarser
    void define(const DisjointBoxLayout&   a_grids,
                const DisjointBoxLayout&   a_gridsFiner,
                const RealVect&            a_dxLevel,
                const IntVect&             a_refRatio, // dummy arg, send in 1
                const IntVect&             a_refRatioFiner,
                const ProblemDomain&       a_domain,
                BCMethodHolder&            a_bc,
                const Copier&              a_exchange,
                const CFRegion&            a_cfregion);

    // Define function for operator which has no finer AMR level
    void define(const DisjointBoxLayout&   a_grids,
                const DisjointBoxLayout&   a_baseBAPtr,
                const RealVect&            a_dxLevel,
                const IntVect&             a_refRatio,
                const ProblemDomain&       a_domain,
                BCMethodHolder&            a_bc,
                const Copier&              a_exchange,
                const CFRegion&            a_cfregion);

    // Define function for operator which has no finer or coarser AMR level
    void define(const DisjointBoxLayout&   a_grids,
                const RealVect&            a_dx,
                const ProblemDomain&       a_domain,
                BCMethodHolder&            a_bc,
                const Copier&              a_exchange,
                const CFRegion&            a_cfregion);

    // Define function for operator which has no finer or coarser AMR level
    // This version is not used by the factory.
    void define(const DisjointBoxLayout&   a_grids,
                const RealVect&            a_dx,
                const ProblemDomain&       a_domain,
                BCMethodHolder&            a_bc) __deprecated;

    // Full define function that mimics the old PoissonOp.
    // Deprecated, throws an error.
    void define(const DisjointBoxLayout& a_grids,
                const DisjointBoxLayout* a_baseBAPtr,
                const RealVect&          a_dxLevel,
                const IntVect&           a_refRatio,
                const ProblemDomain&     a_domain,
                BCMethodHolder&          a_bc) __deprecated;
public:

    // Installs a new metric
    virtual void setJgup(RefCountedPtr<LevelData<FluxBox> >& a_JgupPtr);

    virtual void setJgup(const FillJgupInterface* a_customFillJgupPtr);

    // Sets the time used to evaluate BCs.
    virtual void setTime(Real a_time);

    // Returns the time used to evaluate BCs.
    virtual Real getTime() const;

    // For now, returns (1,1,1) for a full op or (1,1,0) for a horizontal op.
    inline virtual IntVect getActiveDirs () const;


    // LinearOp functions ------------------------------------------------------

    virtual void residual(LevelData<FArrayBox>&       a_lhs,
                          const LevelData<FArrayBox>& a_phi,
                          const LevelData<FArrayBox>& a_rhs,
                          bool                        a_homogeneous = false);

    /// despite what you might think, the "I" here means "Ignore the coarse-fine boundary"
    virtual void residualI(LevelData<FArrayBox>&       a_lhs,
                           const LevelData<FArrayBox>& a_phi,
                           const LevelData<FArrayBox>& a_rhs,
                           bool                        a_homogeneous = false);

    virtual void preCond(LevelData<FArrayBox>&       a_correction,
                         const LevelData<FArrayBox>& a_residual);

    virtual void applyOp(LevelData<FArrayBox>&       a_lhs,
                         const LevelData<FArrayBox>& a_phi,
                         bool                        a_homogeneous = false);

    /// despite what you might think, the "I" here means "Ignore the coarse-fine boundary"
    virtual void applyOpI(LevelData<FArrayBox>&       a_lhs,
                          const LevelData<FArrayBox>& a_phi,
                          bool                        a_homogeneous = false) __hot;

    virtual void applyOpNoBoundary(LevelData<FArrayBox>&       a_lhs,
                                   const LevelData<FArrayBox>& a_rhs,
                                   const bool                  a_homogeneous);

    inline virtual void create(LevelData<FArrayBox>&       a_lhs,
                               const LevelData<FArrayBox>& a_rhs);

    virtual void createCoarsened(LevelData<FArrayBox>&       a_lhs,
                                 const LevelData<FArrayBox>& a_rhs,
                                 const IntVect&              a_refRat);

    inline virtual void assign(LevelData<FArrayBox>&       a_lhs,
                               const LevelData<FArrayBox>& a_rhs);

    inline virtual void assignLocal(LevelData<FArrayBox>&       a_lhs,
                                    const LevelData<FArrayBox>& a_rhs);

    inline virtual void buildCopier(Copier& a_copier,
                                    const LevelData<FArrayBox>& a_lhs,
                                    const LevelData<FArrayBox>& a_rhs);

    inline virtual void assignCopier(LevelData<FArrayBox>&       a_lhs,
                                     const LevelData<FArrayBox>& a_rhs,
                                     const Copier&               a_copier);

    inline virtual void zeroCovered(LevelData<FArrayBox>& a_lhs,
                                    LevelData<FArrayBox>& a_rhs,
                                    const Copier&         a_copier);

    inline virtual Real dotProduct(const LevelData<FArrayBox>& a_1,
                                   const LevelData<FArrayBox>& a_2);

    /* multiple dot products (for GMRES) */
    inline virtual void mDotProduct(const LevelData<FArrayBox>& a_1,
                                    const int a_sz,
                                    const LevelData<FArrayBox> a_2[],
                                    Real a_mdots[]);

    inline virtual void incr(LevelData<FArrayBox>&       a_lhs,
                             const LevelData<FArrayBox>& a_x,
                             Real                        a_scale);

    inline virtual void axby(LevelData<FArrayBox>&       a_lhs,
                             const LevelData<FArrayBox>& a_x,
                             const LevelData<FArrayBox>& a_y,
                             Real                        a_a,
                             Real                        a_b);

    inline virtual void scale(LevelData<FArrayBox>& a_lhs,
                              const Real&           a_scale);

    inline virtual Real norm(const LevelData<FArrayBox>& a_x,
                             int                         a_ord);

    virtual Real localMaxNorm(const LevelData<FArrayBox>& a_x);

    inline virtual void setToZero(LevelData<FArrayBox>& a_x);


    // MGLevelOp functions -----------------------------------------------------

    virtual void relax(LevelData<FArrayBox>&       a_e,
                       const LevelData<FArrayBox>& a_residual,
                       int                         a_iterations);

    // MG function to define (but not fill) the coarser level data holder.
    virtual void createCoarser(LevelData<FArrayBox>&       a_coarse,
                               const LevelData<FArrayBox>& a_fine,
                               bool                        a_ghosted);

    // calculate restricted residual
    // a_resCoarse[2h] = I[h->2h] (rhsFine[h] - L[h](phiFine[h])
    virtual void restrictResidual(LevelData<FArrayBox>&       a_resCoarse,
                                  LevelData<FArrayBox>&       a_phiFine,
                                  const LevelData<FArrayBox>& a_rhsFine);

    // correct the fine solution based on coarse correction
    // a_phiThisLevel += I[2h->h](a_correctCoarse)
    inline virtual void prolongIncrement(LevelData<FArrayBox>&       a_phiThisLevel,
                                         const LevelData<FArrayBox>& a_correctCoarse);


    // AMRLevelOp functions ----------------------------------------------------

    // a_residual = a_rhs - L(a_phi, a_phiFine, a_phiCoarse)
    virtual void AMRResidual(LevelData<FArrayBox>&                    a_residual,
                             const LevelData<FArrayBox>&              a_phiFine,
                             const LevelData<FArrayBox>&              a_phi,
                             const LevelData<FArrayBox>&              a_phiCoarse,
                             const LevelData<FArrayBox>&              a_rhs,
                             const bool                               a_homogeneousPhysBC,
                             MappedAMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

    // residual assuming no more coarser AMR levels
    virtual void AMRResidualNC(LevelData<FArrayBox>&                    a_residual,
                               const LevelData<FArrayBox>&              a_phiFine,
                               const LevelData<FArrayBox>&              a_phi,
                               const LevelData<FArrayBox>&              a_rhs,
                               const bool                               a_homogeneousPhysBC,
                               MappedAMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

    // a_residual = a_rhs - L(a_phi, a_phiCoarse)
    virtual void AMRResidualNF(LevelData<FArrayBox>&       a_residual,
                               const LevelData<FArrayBox>& a_phi,
                               const LevelData<FArrayBox>& a_phiCoarse,
                               const LevelData<FArrayBox>& a_rhs,
                               const bool                  a_homogeneousPhysBC);

    // Apply the AMR operator, including coarse-fine matching
    virtual void AMROperator(LevelData<FArrayBox>&                    a_LofPhi,
                             const LevelData<FArrayBox>&              a_phiFine,
                             const LevelData<FArrayBox>&              a_phi,
                             const LevelData<FArrayBox>&              a_phiCoarse,
                             const bool                               a_homogeneousDomBC,
                             MappedAMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

    //  Apply the AMR operator, including coarse-fine matching
    //  assume no coarser AMR level
    virtual void AMROperatorNC(LevelData<FArrayBox>&                    a_LofPhi,
                               const LevelData<FArrayBox>&              a_phiFine,
                               const LevelData<FArrayBox>&              a_phi,
                               const bool                               a_homogeneousBC,
                               MappedAMRLevelOp<LevelData<FArrayBox> >* a_finerOp);


    //  Apply the AMR operator, including coarse-fine matching.
    //  assume no finer AMR level
    virtual void AMROperatorNF(LevelData<FArrayBox>&       a_LofPhi,
                               const LevelData<FArrayBox>& a_phi,
                               const LevelData<FArrayBox>& a_phiCoarse,
                               const bool                  a_homogeneousBC);

    //  a_resCoarse = I[h->2h]( a_residual - L(a_correction, a_coarseCorrection) )
    //  it is assumed that a_resCoarse has already been filled in with the coarse
    //  version of AMRResidualNF and that this operation is free to overwrite
    //  in the overlap regions.
    virtual void AMRRestrict(LevelData<FArrayBox>&       a_resCoarse,
                             const LevelData<FArrayBox>& a_residual,
                             const LevelData<FArrayBox>& a_correction,
                             const LevelData<FArrayBox>& a_coarseCorrection);

    virtual void AMRRestrictS(LevelData<FArrayBox>&       a_resCoarse,
                              const LevelData<FArrayBox>& a_residual,
                              const LevelData<FArrayBox>& a_correction,
                              const LevelData<FArrayBox>& a_coarseCorrection,
                              LevelData<FArrayBox>&       a_scratch);

    //  a_correction += I[2h->h](a_coarseCorrection)
    virtual void AMRProlong(LevelData<FArrayBox>&       a_correction,
                            const LevelData<FArrayBox>& a_coarseCorrection);

    //  optimization of AMRProlong that sends in the existing temporary and copier
    virtual void AMRProlongS(LevelData<FArrayBox>&       a_correction,
                             const LevelData<FArrayBox>& a_coarseCorrection,
                             LevelData<FArrayBox>&       a_temp,
                             const Copier&               a_copier);

    //  a_residual = a_residual - L(a_correction, a_coarseCorrection)
    virtual void AMRUpdateResidual(LevelData<FArrayBox>&       a_residual,
                                   const LevelData<FArrayBox>& a_correction,
                                   const LevelData<FArrayBox>& a_coarseCorrection);

    // Compute norm over all cells on coarse not covered by finer
    virtual Real AMRNorm(const LevelData<FArrayBox>& a_coarseResid,
                         const LevelData<FArrayBox>& a_fineResid,
                         const IntVect&              a_refRat,
                         const int                   a_ord);

    // For MappedAMRMultiGrid
    virtual IntVect refToCoarser () const {
        return m_refToCoarser;
    }

    // For TGA to reset parameters
    virtual void setAlphaAndBeta(const Real& a_alpha,
                                 const Real& a_beta);

    // Change boundary conditions
    inline virtual void setBC(const BCMethodHolder& a_bc);

    // For TGA stuff---in this case a noop
    inline virtual void diagonalScale(LevelData<FArrayBox>& a_rhs,
                                      bool                  a_kappaWeighted);

    // For TGA stuff---in this case a noop
    inline virtual void divideByIdentityCoef(LevelData<FArrayBox>& a_rhs);

    // For TGA stuff -- in this case, a no-op.
    // This could be used to compute the flux to be saved and retrieved by
    // getFlux, but that's only necessary when the flux needs to be defined
    // over the entire level in one shot.
    inline virtual void fillGrad (const LevelData<FArrayBox>& a_phi);

    // For TGA stuff
    // Since fillGrad is a noop, this function must perform the flux calc.
    virtual void getFlux(FluxBox&                    a_flux,
                         const LevelData<FArrayBox>& a_data,
                         const Box&                  a_grid,
                         const DataIndex&            a_dit,
                         Real                        a_scale);

    virtual void reflux(const LevelData<FArrayBox>&              a_phiFine,
                        const LevelData<FArrayBox>&              a_phi,
                        LevelData<FArrayBox>&                    a_residual,
                        MappedAMRLevelOp<LevelData<FArrayBox> >* a_finerOp);

    virtual void write(const LevelData<FArrayBox>* a_data,
                       const char*                 a_filename);

    virtual void outputAMR (const Vector<LevelData<FArrayBox>*>& a_rhs,
                            const std::string a_name);

    // For LepticSolver...
    inline virtual const RealVect& getDx () const;
    inline virtual const RealVect& getDxCrse () const;
    inline virtual void setDxCrse (const RealVect& a_dxCrse);
    inline virtual const RefCountedPtr<LevelData<FluxBox> > getFCJgup () const;
    inline virtual const RefCountedPtr<LevelData<FArrayBox> > getCCJinv () const;
    inline virtual BCMethodHolder getBCs () const;

    // This is this operator's smallest allowable grid size.
    // (Just an accessor for s_maxCoarse.)
    inline virtual int smallestGridSize () const;

    // This is just for debugging purposes
    virtual const ProblemDomain& getDomain () const {return m_domain;}

    // Public constants
    Real m_alpha, m_beta, m_aCoef, m_bCoef;

    int m_preCondSmoothIters;

    RefCountedPtr<LevelData<FluxBox> >   m_FCJgup;
    RefCountedPtr<LevelData<FArrayBox> > m_CCJinv;
    RefCountedPtr<LevelData<FArrayBox> > m_lapDiag;

    // Needed for homogeneous interpolation
    // Set by the factory
    RealVect  m_dxCrse;

    static int s_maxCoarse;

protected:
    std::string             m_levid;

    Real                    m_time;

    RealVect                m_dx;
    ProblemDomain           m_domain;
    Box                     m_validDomain;

    LevelDataOps<FArrayBox> m_levelOps;

    BCMethodHolder          m_bc;

    CFRegion                m_cfregion;

    Copier                  m_exchangeCopier;

    MappedQuadCFInterp      m_interpWithCoarser;

    MappedLevelFluxRegister m_levfluxreg;

    DisjointBoxLayout       m_coarsenedMGrids;

    IntVect                 m_refToCoarser; // AMR, not MG
    IntVect                 m_refToFiner;   // AMR, not MG

    bool                    m_isDiagonal;

    // Set to false for normal solves, and true for horizontal leptic solves.
    bool                    m_horizontalOp;

    // Multigrid methods and strategies.
    RestrictionStrategy*    m_restrictPtr;
    ProlongationStrategy*   m_prolongPtr;
    RelaxationMethod*       m_relaxPtr;
    RelaxationMethod*       m_precondRelaxPtr;
    int                     m_precondMode;



    virtual void getFlux(FArrayBox&       a_flux,
                         const FArrayBox& a_data,
                         const DataIndex& a_di,
                         int              a_dir,
                         int              a_ref = 1) const;

    virtual void getFlux(FArrayBox&       a_flux,
                         const FArrayBox& a_data,
                         const Box&       a_edgebox,
                         const DataIndex& a_di,
                         int              a_dir,
                         int              a_ref = 1) const;

    virtual void getFluxDuringExchange(FArrayBox&       a_flux,
                                       const FArrayBox& a_data,
                                       const FArrayBox& a_extrap,
                                       const Box&       a_edgebox,
                                       const DataIndex& a_di,
                                       int              a_dir,
                                       int              a_ref = 1) const;

    virtual void getFluxAfterExchange(FArrayBox&       a_flux,
                                      const FArrayBox& a_data,
                                      const FArrayBox& a_extrap,
                                      const Box&       a_edgebox,
                                      const DataIndex& a_di,
                                      int              a_dir,
                                      int              a_ref = 1) const;

    virtual void getFluxComplete (FArrayBox&       a_flux,
                                  const FArrayBox& a_data,
                                  const FArrayBox& a_extrap,
                                  const Box&       a_edgebox,
                                  const DataIndex& a_di,
                                  int              a_dir,
                                  int              a_ref = 1) const;

    // Calls MappedQuadCFInterp or homogeneousCFInterp. Also fills corner ghosts.
    // This operator is the finer level.
    virtual void interpCFGhosts (LevelData<FArrayBox>&       a_phi,
                                 const LevelData<FArrayBox>* a_phiCoarsePtr,
                                 const bool                  a_isHomogeneous) const;

    // Exchanges all necessary ghosts.
    virtual void exchangeComplete (LevelData<FArrayBox>& a_phi) const;

    // Copies a_phi and extrapolates ghosts at physical boundaries.
    virtual void fillExtrap (FArrayBox&       a_extrap,
                             const FArrayBox& a_phi,
                             const int        a_order) const;

private:
    virtual void createCoarsened(LevelData<FArrayBox>&       a_lhs,
                                 const LevelData<FArrayBox>& a_rhs,
                                 const int&                  a_refRat)
    {
        MayDay::Error("Illegal function");
    }

    virtual Real dx() const
    {
        MayDay::Error("Illegal function");
        return 0.0;
    }
};


// **************************** Inline functions *******************************

// -----------------------------------------------------------------------------
// Creates a data holder a_lhs that mirrors a_rhs. Does not copy the data of
// a_rhs, just makes a holder the same size.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::create(LevelData<FArrayBox>&       a_lhs,
                                const LevelData<FArrayBox>& a_rhs)
{
    m_levelOps.create(a_lhs, a_rhs);
}


// -----------------------------------------------------------------------------
// Sets a_lhs equal to a_rhs. This function will do exchanges if the data
// holders do not share the same grids.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::assign(LevelData<FArrayBox>&       a_lhs,
                                const LevelData<FArrayBox>& a_rhs)
{
    m_levelOps.assign(a_lhs, a_rhs);
}


// -----------------------------------------------------------------------------
// Set a_lhs equal to a_rhs assuming a_lhs and a_rhs share the same grids.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::assignLocal(LevelData<FArrayBox>&        a_lhs,
                                     const  LevelData<FArrayBox>& a_rhs)
{
    DataIterator dit = a_lhs.dataIterator();
    for (dit.reset(); dit.ok(); ++dit)
        a_lhs[dit].copy(a_rhs[dit]);
}


// -----------------------------------------------------------------------------
// Creates a copier that sends data from a_rhs to a_lhs.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::buildCopier(Copier&                     a_copier,
                                     const LevelData<FArrayBox>& a_lhs,
                                     const LevelData<FArrayBox>& a_rhs)
{
    a_copier.define(a_rhs.disjointBoxLayout(),
                    a_lhs.disjointBoxLayout(),
                    IntVect::Zero);
}


// -----------------------------------------------------------------------------
// Copies data from a_rhs to a_lhs using a_copier.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::assignCopier(LevelData<FArrayBox>&       a_lhs,
                                      const LevelData<FArrayBox>& a_rhs,
                                      const Copier&               a_copier)
{
    a_rhs.copyTo(a_rhs.interval(), a_lhs, a_lhs.interval(), a_copier);
}

// -----------------------------------------------------------------------------
// Sets a_rhs to zero, then copies the null data to a_lhs using a_copier.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::zeroCovered(LevelData<FArrayBox>& a_lhs,
                                     LevelData<FArrayBox>& a_rhs,
                                     const Copier&         a_copier)
{
    m_levelOps.copyToZero(a_lhs, a_copier);
}


// -----------------------------------------------------------------------------
// Calculates a_1 . a_2 without using any notion of a metric.
// -----------------------------------------------------------------------------
Real MappedAMRPoissonOp::dotProduct(const LevelData<FArrayBox>& a_1,
                                    const LevelData<FArrayBox>& a_2)
{
    return m_levelOps.dotProduct(a_1, a_2);
}


// -----------------------------------------------------------------------------
// Same as dotProduct, but performs the calculation over a_sz different a_2s.
// All a_sz results are placed in the entries of a_mdots.
// a_mdots needs to be properly allocated prior to call.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::mDotProduct(const LevelData<FArrayBox>& a_1,
                                     const int                   a_sz,
                                     const LevelData<FArrayBox>  a_2[],
                                     Real                        a_mdots[])
{
    m_levelOps.mDotProduct(a_1, a_sz, a_2, a_mdots);
}


// -----------------------------------------------------------------------------
// Increment by scaled amount (a_lhs += a_scale*a_x).
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::incr(LevelData<FArrayBox>&       a_lhs,
                              const LevelData<FArrayBox>& a_x,
                              Real                        a_scale)
{
    m_levelOps.incr(a_lhs, a_x, a_scale);
}


// -----------------------------------------------------------------------------
// Set input to a scaled sum (a_lhs = a_a*a_x + a_b*a_y).
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::axby(LevelData<FArrayBox>&       a_lhs,
                              const LevelData<FArrayBox>& a_x,
                              const LevelData<FArrayBox>& a_y,
                              Real                        a_a,
                              Real                        a_b)
{
    m_levelOps.axby(a_lhs, a_x, a_y, a_a, a_b);
}


// -----------------------------------------------------------------------------
// Multiply the input by a given scale (a_lhs *= a_scale).
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::scale(LevelData<FArrayBox>& a_lhs,
                               const Real&           a_scale)
{
    m_levelOps.scale(a_lhs, a_scale);
}


// -----------------------------------------------------------------------------
// Return the norm of a_rhs.
// a_ord == 0 max norm, a_ord == 1 sum(abs(a_rhs)), else, L(a_ord) norm.
// -----------------------------------------------------------------------------
Real MappedAMRPoissonOp::norm(const LevelData<FArrayBox>& a_x,
                              int                         a_ord)
{
    return CH_XD::norm(a_x, a_x.interval(), a_ord);
}


// -----------------------------------------------------------------------------
// Set a_lhs to zero.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::setToZero(LevelData<FArrayBox>& a_lhs)
{
    m_levelOps.setToZero(a_lhs);
}


// -----------------------------------------------------------------------------
// Accessor for the LepticSolver.
// -----------------------------------------------------------------------------
const RealVect& MappedAMRPoissonOp::getDx () const
{
    return m_dx;
}


// -----------------------------------------------------------------------------
// Accessor for the LepticSolver.
// -----------------------------------------------------------------------------
const RealVect& MappedAMRPoissonOp::getDxCrse () const
{
    return m_dxCrse; // AMR, not MG
}


// -----------------------------------------------------------------------------
// Accessor for the LepticSolver.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::setDxCrse (const RealVect& a_dxCrse)
{
    m_dxCrse = a_dxCrse; // AMR, not MG
}


// -----------------------------------------------------------------------------
// Accessor for the LepticSolver.
// -----------------------------------------------------------------------------
const RefCountedPtr<LevelData<FluxBox> > MappedAMRPoissonOp::getFCJgup () const
{
    return m_FCJgup;
}


// -----------------------------------------------------------------------------
// Accessor for the LepticSolver.
// -----------------------------------------------------------------------------
const RefCountedPtr<LevelData<FArrayBox> > MappedAMRPoissonOp::getCCJinv () const
{
    return m_CCJinv;
}


// -----------------------------------------------------------------------------
// This is this operator's smallest allowable grid size.
// (Just an accessor for s_maxCoarse.)
// -----------------------------------------------------------------------------
int MappedAMRPoissonOp::smallestGridSize () const
{
    return s_maxCoarse;
}


// -----------------------------------------------------------------------------
// Accessor for the LepticSolver.
// -----------------------------------------------------------------------------
BCMethodHolder MappedAMRPoissonOp::getBCs () const
{
    return m_bc;
}


// -----------------------------------------------------------------------------
// Access function.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::setBC(const BCMethodHolder& a_bc)
{
    m_bc = a_bc;
}


// -----------------------------------------------------------------------------
// Computes prolongation of correction to finer level
// phi[h] += I[mh->h](coarse[mh])
// This op is on the finer grid.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::prolongIncrement (LevelData<FArrayBox>&       a_phiThisLevel,
                                           const LevelData<FArrayBox>& a_correctCoarse)
{
    CH_assert(m_prolongPtr != NULL);
    m_prolongPtr->prolongIncrement(a_phiThisLevel, a_correctCoarse);
}


// -----------------------------------------------------------------------------
// For now, returns (1,1,1) for a full op or (1,1,0) for a horizontal op.
// -----------------------------------------------------------------------------
IntVect MappedAMRPoissonOp::getActiveDirs () const
{
    if (m_horizontalOp) return IntVect::Unit - BASISV(SpaceDim-1);
    return IntVect::Unit;
}


// -----------------------------------------------------------------------------
// For TGA stuff---in this case a noop
// Used when solving f(x) * dphi/dt = L[phi].
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::diagonalScale(LevelData<FArrayBox>& a_rhs,
                                       bool                  a_kappaWeighted)
{;}


// -----------------------------------------------------------------------------
// For TGA stuff---in this case a noop
// I think this is similar to diagonalScale.
// Unfortunately, this function is undocumented.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::divideByIdentityCoef(LevelData<FArrayBox>& a_rhs)
{;}


// -----------------------------------------------------------------------------
// For TGA stuff -- in this case, a no-op.
// This could be used to compute the flux to be saved and retrieved by
// getFlux, but that's only necessary when the flux needs to be defined
// over the entire level in one shot.
// -----------------------------------------------------------------------------
void MappedAMRPoissonOp::fillGrad (const LevelData<FArrayBox>& a_phi)
{;}



#endif //!__MAPPEDAMRPOISSONOP_H__INCLUDED__

