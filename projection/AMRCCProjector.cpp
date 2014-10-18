#include "AMRCCProjector.H"
#include "Debug.H"
#include "Copier.H"
#include "CornerCopier.H"


// -----------------------------------------------------------------------------
// Default constructor
// This sets the solver parameters, but leaves object unusable.
// -----------------------------------------------------------------------------
AMRCCProjector::AMRCCProjector ()
: m_isDefined(false)
{
    const ProblemContext* ctx = ProblemContext::getInstance();

    m_solver.setAMRMGParameters(ctx->syncProjection_AMRMG_imin,
                                ctx->syncProjection_AMRMG_imax,
                                ctx->syncProjection_AMRMG_eps,
                                ctx->syncProjection_AMRMG_maxDepth,
                                ctx->syncProjection_AMRMG_num_smooth_precond,
                                ctx->syncProjection_AMRMG_num_smooth_down,
                                ctx->syncProjection_AMRMG_num_smooth_up,
                                ctx->syncProjection_AMRMG_num_smooth_bottom,
                                ctx->syncProjection_AMRMG_precondMode,
                                ctx->syncProjection_AMRMG_relaxMode,
                                ctx->syncProjection_AMRMG_numMG,
                                ctx->syncProjection_AMRMG_hang,
                                ctx->syncProjection_AMRMG_normThresh,
                                ctx->syncProjection_AMRMG_verbosity);

    m_solver.setBottomParameters(ctx->syncProjection_bottom_imax,
                                 ctx->syncProjection_bottom_numRestarts,
                                 ctx->syncProjection_bottom_eps,
                                 ctx->syncProjection_bottom_reps,
                                 ctx->syncProjection_bottom_hang,
                                 ctx->syncProjection_bottom_small,
                                 ctx->syncProjection_bottom_normType,
                                 ctx->syncProjection_bottom_verbosity);

    m_applySyncProjection = ctx->applySyncCorrection;
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
AMRCCProjector::~AMRCCProjector ()
{
    undefine();
}


// -----------------------------------------------------------------------------
// Allocates memory and leaves object useable.
// -----------------------------------------------------------------------------
void AMRCCProjector::define (Vector<LevelData<FArrayBox>*>& a_amrPressure,
                             const PhysBCUtil&              a_physBCUtil,
                             const LevelGeometry&           a_levGeo,
                             const FillJgupInterface*       a_customFillJgupPtr)
{
    // Clear the old data.
    undefine();

    // Collect composite structures
    const int numLevels = a_amrPressure.size();
    m_amrLevGeo = a_levGeo.getAMRLevGeos();
    m_pressure = a_amrPressure;

#   ifndef NDEBUG
        // Sanity checks...

        // Find the lowest defined level.
        int lbase = 0;
        for (; lbase < numLevels; ++lbase) {
            if (m_pressure[lbase] != NULL) break;
        }

        // All levels above lbase MUST be well-defined.
        int lev = lbase;
        for (; lev < numLevels; ++lev) {
            // Is pressure allocated?
            if (m_pressure[lev] == NULL) {
                ostringstream msg;
                msg << "AMRCCProjector::define found lbase = " << lbase
                    << " but pressure is NULL at level = " << lev;
                MayDay::Error(msg.str().c_str());
            }

            // Is pressure defined?
            if (!(m_pressure[lev]->isDefined())) {
                ostringstream msg;
                msg << "AMRCCProjector::define found lbase = " << lbase
                    << " but pressure is undefined at level = " << lev;
                MayDay::Error(msg.str().c_str());
            }

            // Do the pressure and levGeo grids match?
            if (!(m_pressure[lev]->getBoxes() == m_amrLevGeo[lev]->getBoxes())) {
                pout() << "m_pressure [" << lev << "]->getBoxes() = " << m_pressure[lev]->getBoxes();
                pout() << "m_amrLevGeo[" << lev << "]->getBoxes() = " << m_amrLevGeo[lev]->getBoxes();
                pout() << std::flush;

                ostringstream msg;
                msg << "AMRCCProjector::define found lbase = " << lbase
                    << " but pressure and levGeo at level = " << lev
                    << " do not have matching grids";
                MayDay::Error(msg.str().c_str());
            }
        }
#   endif

    // Collect BCs.
    const ProblemContext* ctx = ProblemContext::getInstance();
    const bool isViscous = (ctx->nu > 0.0);

    m_solverBC = a_physBCUtil.SyncProjFuncBC();
    m_divBC = a_physBCUtil.uStarFuncBC(isViscous);
    m_gradBC = a_physBCUtil.gradESyncFuncBC();

    // Define CF-BC interpolators.
    m_pressureCFInterp.resize(numLevels);
    m_velCFInterp.resize(numLevels);

    for (int lev = 0; lev < numLevels; ++lev) {
        m_pressureCFInterp[lev] = RefCountedPtr<MappedQuadCFInterp>(new MappedQuadCFInterp);
        m_velCFInterp[lev] = RefCountedPtr<MappedQuadCFInterp>(new MappedQuadCFInterp);
    }

    for (int lev = 1; lev < numLevels; ++lev) {
        const LevelGeometry& levGeo = *m_amrLevGeo[lev];

        m_pressureCFInterp[lev]->define(levGeo.getBoxes(),
                                        &(m_amrLevGeo[lev-1]->getBoxes()),
                                        levGeo.getDx(),
                                        levGeo.getCrseRefRatio(),
                                        1, // num comps
                                        levGeo.getDomain());

        m_velCFInterp[lev]->define(levGeo.getBoxes(),
                                   &(m_amrLevGeo[lev-1]->getBoxes()),
                                   levGeo.getDx(),
                                   levGeo.getCrseRefRatio(),
                                   SpaceDim, // num comps
                                   levGeo.getDomain());
    }

    // Define the pressure solver.
    const Box lminDomBox = m_amrLevGeo[0]->getDomain().domainBox();
    m_solver.define(m_solverBC, a_levGeo, lminDomBox, numLevels, a_customFillJgupPtr);

    // This object is ready to be used.
    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Frees memory and leaves object unuseable.
// -----------------------------------------------------------------------------
void AMRCCProjector::undefine ()
{
    if (m_isDefined) {
        // Base members
        m_time = BOGUS_TIME;
        m_pressure.clear();
        m_solver.undefine();

        // Our members
        m_amrLevGeo.clear();
        m_pressureCFInterp.clear();
        m_velCFInterp.clear();
        m_isDefined = false;
    }
}


// -----------------------------------------------------------------------------
// From BaseProjector:
//  Computes J^{-1}\partial_i(J u^i) over an AMR hierarchy.
//  This must be overriden or an error will be thrown.
// -----------------------------------------------------------------------------
void AMRCCProjector::computeDiv (Vector<LevelData<FArrayBox>*>&       a_div,
                                 const Vector<LevelData<FArrayBox>*>& a_flux,
                                 const int                            a_lmin,
                                 const int                            a_lmax) const
{
    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(0 <= a_lmin);
    CH_assert(a_lmin <= a_lmax);
    CH_assert(a_div .size() >= a_lmax);
    CH_assert(a_flux.size() >= a_lmax);

    // Loop over levels and compute Div(Flux)
    for (int lev = a_lmin; lev <= a_lmax; ++lev) {
        // More sanity checks
        CH_assert(a_div[lev] != NULL);
        CH_assert(a_div[lev]->getBoxes().physDomain() == m_amrLevGeo[lev]->getBoxes().physDomain());
        CH_assert(a_div[lev]->getBoxes() == m_amrLevGeo[lev]->getBoxes());

        CH_assert(a_flux[lev] != NULL);
        CH_assert(a_flux[lev]->getBoxes().physDomain() == m_amrLevGeo[lev]->getBoxes().physDomain());
        CH_assert(a_flux[lev]->getBoxes() == m_amrLevGeo[lev]->getBoxes());

        // Grab the coarse velocity pointer
        LevelData<FArrayBox>* crseFluxPtr = NULL;
        if (lev > 0) {
            crseFluxPtr = a_flux[lev-1];
        }

        // Grab the fine velocity pointer
        LevelData<FArrayBox>* fineFluxPtr = NULL;
        if (lev < a_lmax) {
            fineFluxPtr = a_flux[lev+1];
        }

        // Just in case...
        Interval velComps(0,SpaceDim-1);
        a_flux[lev]->exchange(velComps);      // TODO: I don't think this is needed!

        // Compute Div[flux]
        Divergence::compDivergenceCC(*a_div[lev],
                                     *a_flux[lev],
                                     crseFluxPtr,
                                     fineFluxPtr,
                                     true,  // doQuadInterp
                                     *m_amrLevGeo[lev],
                                     m_time,
                                     &m_divBC);
    }
}


// -----------------------------------------------------------------------------
// From BaseProjector:
//  Computes Jg^{i,j}\partial_j(phi) over an AMR hierarchy.
//  This must be overriden or an error will be thrown.
// -----------------------------------------------------------------------------
void AMRCCProjector::computeGrad (Vector<LevelData<FArrayBox>*>&       a_flux,
                                  const Vector<LevelData<FArrayBox>*>& a_phi,
                                  const int                            a_lmin,
                                  const int                            a_lmax) const
{
    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(0 <= a_lmin);
    CH_assert(a_lmin <= a_lmax);
    CH_assert(a_flux.size() >= a_lmax);
    CH_assert(a_phi .size() >= a_lmax);

    // Loop over levels and compute Div(Flux)
    for (int lev = a_lmin; lev <= a_lmax; ++lev) {
        // More sanity checks
        CH_assert(a_phi[lev] != NULL);
        CH_assert(a_phi[lev]->getBoxes().physDomain() == m_amrLevGeo[lev]->getBoxes().physDomain());
        CH_assert(a_phi[lev]->getBoxes() == m_amrLevGeo[lev]->getBoxes());

        CH_assert(a_flux[lev] != NULL);
        CH_assert(a_flux[lev]->getBoxes().physDomain() == m_amrLevGeo[lev]->getBoxes().physDomain());
        CH_assert(a_flux[lev]->getBoxes() == m_amrLevGeo[lev]->getBoxes());

        // Grab the coarse phi pointer
        LevelData<FArrayBox>* crsePhiPtr = NULL;
        if (lev > 0) {
            crsePhiPtr = a_phi[lev-1];
        }

        // Grab the fine phi pointer
        LevelData<FArrayBox>* finePhiPtr = NULL;
        if (lev < a_lmax) {
            finePhiPtr = a_phi[lev+1];
        }

        // Just in case...
        {
            const DisjointBoxLayout& grids = a_phi[lev]->getBoxes();
            const ProblemDomain& domain = grids.physDomain();
            const IntVect& ghostVect = a_phi[lev]->ghostVect();

            Copier excp(grids, grids, domain, ghostVect, true);
            a_phi[lev]->exchange(excp);      // TODO: Is this needed?

            CornerCopier exccp(grids, grids, domain, ghostVect, true);
            a_phi[lev]->exchange(exccp);     // TODO: Is this needed?
        }

        // Compute Grad[phi]
        Gradient::compGradientCC(*a_flux[lev],
                                 *a_phi[lev],
                                 crsePhiPtr,
                                 finePhiPtr,
                                 (MappedQuadCFInterp&)(*m_pressureCFInterp[lev]),
                                 *m_amrLevGeo[lev],
                                 m_time,
                                 &m_gradBC);
    }
}


// -----------------------------------------------------------------------------
// Applies vel = vel - corr. Dt is built into corr.
// -----------------------------------------------------------------------------
void AMRCCProjector::applyCorrection (Vector<LevelData<FArrayBox>*>&       a_amrJVel,
                                      const Vector<LevelData<FArrayBox>*>& a_amrCorr,
                                      const Real                           a_dt,
                                      const int                            a_lmin,
                                      const int                            a_lmax) const
{
    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(0 <= a_lmin);
    CH_assert(a_lmin <= a_lmax);
    CH_assert(a_amrJVel.size() >= a_lmax);
    CH_assert(a_amrCorr.size() >= a_lmax);

    // Do we want to apply the projection?
    if (!m_applySyncProjection) return;

    const Real dtScale = (a_dt == 0.0)? -1.0: -a_dt;

    for (signed int lev = a_lmax; lev >= a_lmin; --lev) {
        // More sanity checks
        CH_assert(a_amrJVel[lev] != NULL);
        CH_assert(a_amrJVel[lev]->getBoxes().physDomain() == m_amrLevGeo[lev]->getBoxes().physDomain());
        CH_assert(a_amrJVel[lev]->getBoxes() == m_amrLevGeo[lev]->getBoxes());

        CH_assert(a_amrCorr[lev] != NULL);
        CH_assert(a_amrCorr[lev]->getBoxes().physDomain() == m_amrLevGeo[lev]->getBoxes().physDomain());
        CH_assert(a_amrCorr[lev]->getBoxes() == m_amrLevGeo[lev]->getBoxes());

        // Apply correction
        DataIterator dit = a_amrJVel[lev]->dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& JVelFAB = (*a_amrJVel[lev])[dit];
            const FArrayBox& corrFAB = (*a_amrCorr[lev])[dit];

            JVelFAB.plus(corrFAB, dtScale);
        }

        // Average valid velocity down
        if (lev < a_lmax) {
            const LevelGeometry* fineLevGeoPtr = m_amrLevGeo[lev+1];
            const bool considerCellSizes = false; // J is built into a_amrJVel.

            MappedCoarseAverage avgDownObj(fineLevGeoPtr->getBoxes(),
                                           SpaceDim,    // comps
                                           fineLevGeoPtr->getCrseRefRatio());

            avgDownObj.averageToCoarse(*a_amrJVel[lev],
                                       *a_amrJVel[lev+1],
                                       fineLevGeoPtr,
                                       considerCellSizes);
        }
    }
}
