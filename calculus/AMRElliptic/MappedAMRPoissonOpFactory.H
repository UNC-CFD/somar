/*******************************************************************************
 *  FILE: MappedAMRPoissonOpFactory.hpp
 *  DESC: Factory class to create MappedAMRPoissonOps
 *  DATE: Wed 02 May 2012 10:20:37 PM EDT
 *  MAIL: santilli@physics.unc.edu
 ******************************************************************************/
#ifndef __MAPPEDAMRPOISSONOPFACTORY_HPP__INCLUDED__
#define __MAPPEDAMRPOISSONOPFACTORY_HPP__INCLUDED__

#include "MappedAMRPoissonOp.H"
#include "LevelGeometry.H"


// -----------------------------------------------------------------------------
// Factory to create MappedAMRPoissonOps
// -----------------------------------------------------------------------------
class MappedAMRPoissonOpFactory: public MappedAMRLevelOpFactory<LevelData<FArrayBox> >
{
public:
    // The default destructor
    MappedAMRPoissonOpFactory()
    : m_bypassValidation(false),
      m_customFillJgupPtr(NULL)
    {;}

    // The default destructor
    virtual ~MappedAMRPoissonOpFactory();

    // Full AMR factory define
    // Assembles the hierarchy of data via a_baseLevGeoPtr->getAMRLevGeos().
    // Really, you can send in any levGeoPtr and this function will find the base.
    void define(const LevelGeometry*     a_baseLevGeoPtr,
                Real                     a_alpha,
                Real                     a_beta,
                BCMethodHolder&          a_bc,
                int                      a_maxDepth,
                int                      a_preCondSmoothIters,
                int                      a_precondMode,
                int                      a_relaxMode,
                bool                     a_horizontalFactory = false,
                const FillJgupInterface* a_customFillJgupPtr = NULL);

    // Full AMR factory define
    // a_coarseDomain: the domain at the coarsest level.
    // a_grids: the AMR hierarchy.
    // a_refRatios: the refinement ratios between levels. The ratio lives with
    //              the coarser level so a_refRatios[ilev] is the ratio between
    //              ilev and ilev+1
    // a_coarseDx: the grid spacing at the coarsest level.
    // a_bc: holds the boundary conditions.
    // a_alpha: the identity coefficient
    // a_beta: the laplacian coefficient.
    void define(const ProblemDomain&                a_coarseDomain,
                const Vector<DisjointBoxLayout>&    a_grids,
                const Vector<IntVect>&              a_refRatios,
                const RealVect&                     a_coarsedx,
                BCMethodHolder&                     a_bc,
                int                                 a_maxDepth,
                int                                 a_preCondSmoothIters,
                int                                 a_precondMode,
                Real                                a_alpha,
                Real                                a_beta,
                const Vector<const LevelGeometry*>& a_vlevGeoPtr,
                int                                 a_relaxMode,
                bool                                a_horizontalFactory = false,
                const FillJgupInterface*            a_customFillJgupPtr = NULL);

    // Single level AMR factory define -- used by the Leptic solver.
    void define(const RefCountedPtr<LevelData<FluxBox> >&   a_JgupPtr,
                const RefCountedPtr<LevelData<FArrayBox> >& a_JinvPtr,
                const Copier&                               a_copier,
                const CFRegion&                             a_cfregion,
                const RealVect&                             a_dx,
                Real                                        a_alpha,
                Real                                        a_beta,
                BCMethodHolder&                             a_bc,
                int                                         a_maxDepth,
                int                                         a_preCondSmoothIters,
                int                                         a_precondMode,
                int                                         a_relaxMode,
                bool                                        a_horizontalFactory = false);

    // Clears all memory pools and undefines this object.
    void clear ();

    // Create an operator at an index space = coarsen(fineIndexSpace, 2^depth).
    // Return NULL if no such Multigrid level can be created at this a_depth.
    // If a_homoOnly = true, then only homogeneous BCs will be needed.
    virtual MappedMGLevelOp<LevelData<FArrayBox> >* MGnewOp(const ProblemDomain&   a_FineindexSpace,
                                                            const int              a_depth,
                                                            const bool             a_homoOnly = true,
                                                            Vector<IntVect>*       a_allMGRefRatiosPtr = NULL,
                                                            const Vector<IntVect>* a_forceAllMGRefRatiosPtr = NULL);

    // Return a new operator. This is done with a new call.
    // Caller is responsible for deletion.
    virtual MappedAMRLevelOp<LevelData<FArrayBox> >* AMRnewOp(const ProblemDomain& a_indexSpace);

    // Returns the refinement ratio to the next finer level
    virtual IntVect getFineRefRatio (const ProblemDomain& a_indexSpace) const;

    // When setting up horizontal solves, this level's LevGeo will not be tied
    // to a coarse-level LevGeo. Therefore, we need to be told what the coarse
    // level dx is.
    virtual void forceDxCrse(const RealVect& a_dxCrse);

protected:
    // Upon completion, the metric field ptrs will be defined properly.
    virtual void validateMetricPtrs(const int                a_AMRlevel,
                                    const int                a_MGdepth,
                                    const IntVect&           a_mgRefRatio,
                                    const IntVect&           a_coarsening,
                                    const DisjointBoxLayout& a_layout);

    // Fill the metric fields with data.
    // This function was made special for filling MG levels.
    virtual void fill_MGfields(const int      a_AMRlevel,
                               const int      a_MGdepth,
                               const IntVect& a_mgRefRatio,
                               const IntVect& a_coarsening);

    // The maximum MG depth (-1 for as deep as possible)
    int m_maxDepth;
    int m_preCondSmoothIters;
    int m_precondMode;
    int m_relaxMode;

    Vector<ProblemDomain>     m_domains;
    Vector<DisjointBoxLayout> m_boxes;
    Vector<Copier>            m_exchangeCopiers;
    Vector<CFRegion>          m_cfregion;

    Vector<RealVect> m_dx;
    Vector<IntVect>  m_refRatios; // refinement to next finer level

    BCMethodHolder m_bc;

    Real m_alpha;
    Real m_beta;

    bool m_horizontalFactory;
    IntVect m_maskedMaxCoarse;

    // The problem's geometry object. This class does not own the pointer.
    Vector<const LevelGeometry*> m_vlevGeoPtr;
    const FillJgupInterface*     m_customFillJgupPtr;

    // TODO: Would be more memory efficient as a std::map<BoxLayout, RCP>
    // These are the AMRMG metric fields. The format is
    // m_fieldPtr[AMR_level][i] where MG_coarsening_factor = 2^i
    Vector< Vector< RefCountedPtr< LevelData<FluxBox> > > >   m_vvJgup;
    Vector< Vector< RefCountedPtr< LevelData<FArrayBox> > > > m_vvJinv;
    Vector< Vector< RefCountedPtr< LevelData<FArrayBox> > > > m_vvlapDiag;

    // Be careful with this.
    bool m_bypassValidation;

    // When setting up horizontal solves, this level's LevGeo will not be tied
    // to a coarse-level LevGeo. Therefore, we need to be told what the coarse
    // level dx is.
    bool     m_useForceDxCrse;
    RealVect m_forceDxCrse;
};


#endif //!__MAPPEDAMRPOISSONOPFACTORY_HPP__INCLUDED__

