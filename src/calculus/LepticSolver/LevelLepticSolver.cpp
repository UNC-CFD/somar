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
#include "LevelLepticSolver.H"
#include "LevelLepticSolverF_F.H"
#include "LepticOperator.H"
#include "LepticBoxUtils.H"
#include "LepticMeshRefine.H"
#include "EllipticBCUtils.H"
#include "ExtrapolationUtils.H"
#include "HomogeneousCFInterp.H"
#include "TridiagUtilsF_F.H"
#include "DivCurlGradF_F.H"
#include "computeMappedSum.H"
#include "Subspace.H"
#include "SubspaceF_F.H"
#include "SetValLevel.H"
#include "MiscUtils.H"
#include "Debug.H"
#include "Constants.H"
#include "AMRLESMeta.H"
#include "ProblemContext.H"


// LevelLepticSolver static members
IntVect LevelLepticSolver::s_vmask = BASISV(CH_SPACEDIM-1);
IntVect LevelLepticSolver::s_hmask = IntVect::Unit - LevelLepticSolver::s_vmask;


// -----------------------------------------------------------------------------
// Default constructor -- leaves object unusable.
// -----------------------------------------------------------------------------
LevelLepticSolver::LevelLepticSolver ()
: m_isDefined(false)
, m_origOpPtr(NULL)             // We don't own this.
, m_crsePhiPtr(NULL)            // We don't own this either.
, m_horizOpFactoryPtr(NULL)     // We own this.
, m_horizSolverPtr(NULL)        // We own this.
, m_horizBottomSolverPtr(NULL)  // That's ours too.
{
    this->undefine();
    this->setDefaultParameters();
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
LevelLepticSolver::~LevelLepticSolver ()
{
    this->undefine();
}


// -----------------------------------------------------------------------------
// Deallocates all memory and brings this solver back to an undefined state.
// -----------------------------------------------------------------------------
void LevelLepticSolver::undefine ()
{
    m_homogeneous = false;
    m_resNorms.clear();
    m_exitStatus = ExitStatus::NONE;

    // LepticOperator stuff
    m_origOpPtr = NULL;
    m_crsePhiPtr = NULL;
    m_dx = RealVect::Zero;
    m_dxCrse = RealVect::Zero;
    m_origJgupPtr = RefCountedPtr<LevelData<FluxBox> >(NULL);
    m_origJinvPtr = RefCountedPtr<LevelData<FArrayBox> >(NULL);

    // Vertical grid stuff
    m_domain = ProblemDomain();
    m_grids = DisjointBoxLayout();
    m_exCopier.clear();
    m_origToVertCopier.clear();
    m_vertToOrigCopier.clear();
    m_CFRegion = CFRegion();
    m_opPtr = RefCountedPtr<MappedAMRLevelOp<LevelData<FArrayBox> > >(NULL);
    m_JgupPtr = RefCountedPtr<LevelData<FluxBox> >(NULL);
    m_JinvPtr = RefCountedPtr<LevelData<FArrayBox> >(NULL);
    m_flatGrids = BoxLayout();
    m_flatDI.clear();
    m_flatDIComplement.clear();
    // m_vertBCTypes.clear(); // LayoutData has no clear method.

    // Horizontal grid stuff
    m_doHorizSolve = true;
    m_horizDomain = ProblemDomain();
    m_horizGrids = DisjointBoxLayout();
    m_flatToHorizCopier.clear();
    m_horizToFlatCopier.clear();
    m_horizRemoveAvg = false;

    delete m_horizOpFactoryPtr;
    m_horizOpFactoryPtr = NULL;

    delete m_horizSolverPtr;
    m_horizSolverPtr = NULL;

    delete m_horizBottomSolverPtr;
    m_horizBottomSolverPtr = NULL;

    // This object is no longer defined
    m_isDefined = false;
}


// -----------------------------------------------------------------------------
// Full define -- leave object in a useable state.
// This is an override of the pure virtual LinearSolver function.
//
// The operator is needed to calculate residuals, to provide metric data,
// to provide BC types, etc. This solver asks more of the op than the
// LinearOp interface can provide. Therefore, the op MUST be able to be
// cast into a LepticOperator. MappedAMRPoissonOp can do this. This function
// asks for a LinearOp and not a LepticOperator so that it can be used as a
// bottom solver.
// -----------------------------------------------------------------------------
void LevelLepticSolver::define (LinearOp<LevelData<FArrayBox> >* a_operator,
                                bool                             a_homogeneous)
{
    if (m_isDefined) {
        this->undefine();
    }

    // Will this be a homogeneous solve? If so, the CF-BCs will be simpler.
    this->setHomogeneous(a_homogeneous);

    // Make sure we were handed an appropriate operator.
    m_origOpPtr = a_operator;
    const LepticOperator* lepticOpPtr = dynamic_cast<LepticOperator*>(a_operator);
    CH_assert(lepticOpPtr != NULL);

    // Extract what we need from the operator.
    m_dx = lepticOpPtr->getDx();
    m_dxCrse = lepticOpPtr->getDxCrse();
    m_origJgupPtr = lepticOpPtr->getFCJgup();
    m_origJinvPtr = lepticOpPtr->getCCJinv();
    int blockFactor = lepticOpPtr->smallestGridSize();
    BCMethodHolder bcHolder = lepticOpPtr->getBCs();


    // Gather info about the original grid structure, domain, periodicity, etc.
    // This will also set up arrays boxes for use in the horizontal solver.

    // Estimate maxBoxSize from the original grids.
    const DisjointBoxLayout& origGrids = m_origJgupPtr->getBoxes();
    const Vector<Box> origBoxArray = origGrids.boxArray();
    const IntVect maxBoxSize = LepticBoxUtils::getMaxBoxSize(origBoxArray);

    // Get domain and pack periodicity info into a C array.
    m_domain = origGrids.physDomain();
    bool isPeriodic[CH_SPACEDIM] = {D_DECL(
        m_domain.isPeriodic(0),
        m_domain.isPeriodic(1),
        m_domain.isPeriodic(2)
    )};

    // Vertical grid stuff...
    {
        // Create boxes suitable for the vertical solver.
        Vector<Box> vertBoxArray;
        LepticBoxUtils::createVerticalSolverGrids(vertBoxArray,
                                                  origBoxArray,
                                                  m_domain.domainBox());
        m_grids.defineAndLoadBalance(vertBoxArray, NULL, m_domain);
        // DataIterator dit = m_grids.dataIterator();

        // Copy metric data. Note that we only copy Jgup and set Jinv to 1.
        // We do this because solve() sets up the residual equation scaled by J to simplify things.
        m_JgupPtr = RefCountedPtr<LevelData<FluxBox> >(new LevelData<FluxBox>(m_grids, SpaceDim));
        CH_assert(!m_JgupPtr.isNull());
        debugInitLevel(*m_JgupPtr);
        m_origJgupPtr->copyTo(*m_JgupPtr);

        m_JinvPtr = RefCountedPtr<LevelData<FArrayBox> >(new LevelData<FArrayBox>(m_grids, 1));
        setValLevel(*m_JinvPtr, 1.0);

        // This helps copy the rhs.
        m_origToVertCopier.define(origGrids, m_grids, m_domain, IntVect::Zero, false);
        m_vertToOrigCopier = m_origToVertCopier;
        m_vertToOrigCopier.reverse();

        // Create exchange copiers.
        m_exCopier.exchangeDefine(m_grids, IntVect::Unit);
        m_exCornerCopier.define(m_grids, m_grids, m_domain, IntVect::Unit, true);

        // Create CFRegion.
        m_CFRegion.define(m_grids, m_domain);

        // Create an operator that can calulate residuals and such.
        {
            MappedAMRPoissonOpFactory opFact;
            opFact.define(m_JgupPtr,
                          m_JinvPtr,
                          m_exCopier,
                          m_CFRegion,
                          m_dx,
                          0.0,      // alpha
                          1.0,      // beta
                          bcHolder,
                          0,
                          0,        // preCondSmoothIters
                          0,        // precondMode
                          ProblemContext::RelaxMode::NORELAX,
                          false);   // horizontal factory?

            m_opPtr = RefCountedPtr<MappedAMRLevelOp<LevelData<FArrayBox> > >(opFact.AMRnewOp(m_domain));
            LepticOperator* lepticOpPtr = dynamic_cast<LepticOperator*>(&*m_opPtr);
            CH_assert(lepticOpPtr != NULL);
            lepticOpPtr->setDxCrse(m_dxCrse);
        }

        // Gather vertical BC types.
        m_vertBCTypes.define(m_grids);
        LevelLepticSolver::gatherVerticalBCTypes(m_vertBCTypes,
                                                 m_doHorizSolve,
                                                 m_domain,
                                                 m_CFRegion,
                                                 bcHolder);
        // pout() << "doHorizSolve = " << m_doHorizSolve << endl;
    } // end of vertical grid stuff.

    // Set up horizontal stuctures, if needed.
    if (m_doHorizSolve) {
        const Box& domBox = m_domain.domainBox();
        const int loIdx = domBox.smallEnd(SpaceDim-1);

        // 1. We need to create a flat set of grids that are compatible with
        // m_grids. The result will NOT be disjoint.
        m_flatGrids.deepCopy(m_grids);
        LepticBoxUtils::FlattenTransform stampy(loIdx);
        m_flatGrids.transform(stampy);
        m_flatGrids.close();
        CH_assert(m_grids.compatible(m_flatGrids));

        // 2. Create an array of data indices that point to grids that vertically
        // span the domain. This subset of m_flatGrids should be disjoint.
        m_flatDI.clear();
        m_flatDIComplement.clear();

        DataIterator dit = m_grids.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            const DataIndex& di = dit();
            const Box& valid = m_grids[di];

            if (LepticBoxUtils::vertSpanCheck(valid, domBox)) {
                m_flatDI.push_back(di);
            } else {
                m_flatDIComplement.push_back(di);
            }
        }

        // 3. We need to make a version of m_flatGrids that are suitable for
        // horizontal solves. That is, we want to remove the unused grids and
        // perform load balancing.

        // Calculate horizontal domain.
        Box flatDomBox = flattenBox(domBox, SpaceDim-1);
        m_horizDomain.define(flatDomBox, isPeriodic);

        // Create the load-balanced horizontal grids.
        Vector<Box> horizBoxArray;
        LepticBoxUtils::createHorizontalSolverGrids(horizBoxArray,
                                                    m_grids.boxArray(),
                                                    domBox,
                                                    blockFactor);
        CH_assert(horizBoxArray.size() > 0);
        m_horizGrids.defineAndLoadBalance(horizBoxArray, NULL, m_horizDomain);

        // Create the copiers...
        m_flatToHorizCopier.define(m_flatGrids, m_horizGrids, m_horizDomain, IntVect::Zero, false);
        m_horizToFlatCopier = m_flatToHorizCopier;
        m_horizToFlatCopier.reverse();

        // Create the LevelGeometry used by the horizontal solvers.
        LevelGeometry horizLevGeo;
        horizLevGeo.define(m_dx);
        horizLevGeo.regridVertAvg(m_horizGrids, domBox);

        // Solver stuff...

        // Define the BCs.
        RefCountedPtr<BCFluxClass> bcFluxPtr(new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                                                              RealVect::Zero,
                                                                              s_hmask,
                                                                              s_hmask));
        BCMethodHolder horizBCHolder;
        horizBCHolder.addBCMethod(bcFluxPtr);

        // Define the AMR op factory.
        m_horizOpFactoryPtr = new MappedAMRPoissonOpFactory;
        m_horizOpFactoryPtr->define(&horizLevGeo,
                                    0.0,     // alpha
                                    1.0,     // beta
                                    horizBCHolder,
                                    m_horiz_maxDepth,
                                    m_horiz_numSmoothPreCond,
                                    0,        // precondMode
                                    m_horiz_relaxType,
                                    true);   // horizontal factory?
        m_horizOpFactoryPtr->forceDxCrse(m_dxCrse);

        // Define the bottom solver.
        m_horizBottomSolverPtr = new BiCGStabSolver<LevelData<FArrayBox> >;
        m_horizBottomSolverPtr->m_numRestarts = m_horizBottom_numRestarts;
        m_horizBottomSolverPtr->m_imax = m_horizBottom_imax;
        m_horizBottomSolverPtr->m_normType = m_normType;
        m_horizBottomSolverPtr->m_verbosity = m_horizBottom_verbosity;
        m_horizBottomSolverPtr->m_eps = m_horizBottom_eps;
        m_horizBottomSolverPtr->m_hang = m_horizBottom_hang;

        // Define the AMR solver.
        m_horizSolverPtr = new MappedAMRMultiGrid<LevelData<FArrayBox> >;
        m_horizSolverPtr->define(m_horizDomain,
                                 *m_horizOpFactoryPtr,
                                 m_horizBottomSolverPtr,
                                 1); // numLevels
        m_horizSolverPtr->m_verbosity = m_horiz_verbosity;
        m_horizSolverPtr->m_imin = m_horiz_imin;
        m_horizSolverPtr->setSolverParameters(m_horiz_numSmoothDown,
                                              m_horiz_numSmoothUp,
                                              m_horiz_numSmoothBottom,
                                              m_horiz_numMG,
                                              m_horiz_imax,
                                              m_horiz_eps,
                                              m_horiz_hang,
                                              m_horiz_normThresh);

        // Do we need to remove the average from the horizontal solution?
        do {
            // I figure that we will need to remove the average if the
            // horizontal domain is completely spanned by m_horizGrids.
            // To check, we can just count the number of cells in each.
            // NOTE: This assumes Neumann BCs all around the domain!

            // Get the number of cells in the horizontal domain.
            const Box& horizDomBox = m_horizDomain.domainBox();
            if (!horizDomBox.numPtsOK()) {
                MayDay::Warning("LevelLepticSolver: horizDomBox.numPtsOK() failed");
                break;
            }
            const long domNumPts = horizDomBox.numPts();

            // Add up the number of cells in each grid.
            long numPts = 0;
            LayoutIterator lit = m_horizGrids.layoutIterator();
            for (lit.reset(); lit.ok(); ++lit) {
                numPts += m_horizGrids[lit].numPts();
            }

            // Compare.
            m_horizRemoveAvg = (numPts == domNumPts);
        } while(0);
    }

    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Reset whether the solver is homogeneous.
// This is an override of the pure virtual LinearSolver function.
// -----------------------------------------------------------------------------
void LevelLepticSolver::setHomogeneous (bool a_homogeneous)
{
    m_homogeneous = a_homogeneous;
}



// -----------------------------------------------------------------------------
// Used for setting CFBCs. Set to NULL for homogeneous CFBCs.
// -----------------------------------------------------------------------------
void LevelLepticSolver::setCrsePhiPtr (const LevelData<FArrayBox>* a_crsePhiPtr)
{
    m_crsePhiPtr = a_crsePhiPtr;
}


// -----------------------------------------------------------------------------
// Initializes all parameters with default values.
// -----------------------------------------------------------------------------
void LevelLepticSolver::setDefaultParameters ()
{
    setParameters(4,        // a_maxOrder
                  1.0e-15,  // a_eps
                  1.0e-15,  // a_hang
                  0,        // a_normType
                  0,        // a_verbosity
                  1.0e-14); // a_horizRhsTol

    setHorizMGParameters(5,         // a_imin
                         20,        // a_imax
                         4,         // a_numSmoothDown
                         4,         // a_numSmoothBottom
                         4,         // a_numSmoothUp
                         2,         // a_numSmoothPreCond
                         ProblemContext::RelaxMode::LEVEL_GSRB,
                         -1,        // a_maxDepth
                         1.0e-12,   // a_eps
                         1.0e-15,   // a_hang
                         1.0e-30,   // a_normThresh
                         0);        // a_verbosity

    setHorizBottomParameters(80,        // a_imax
                             1.0e-12,   // a_eps
                             5,         // a_numRestarts
                             1.0e-15,   // a_hang
                             0);        // a_verbosity
}


// -----------------------------------------------------------------------------
// Sets the solver's overall parameters.
// a_horizRhsTol is set in this function and not in one of the horiz param
// functions because it is a statement that controls the overall performance
// of the leptic solver, not the horizontal solver. If |horizontal rhs| is less
// than a_horizRhsTol, then a horizontal solve will not be performed.
// -----------------------------------------------------------------------------
void LevelLepticSolver::setParameters (const int  a_maxOrder,
                                       const Real a_eps,
                                       const Real a_hang,
                                       const int  a_normType,
                                       const int  a_verbosity,
                                       const Real a_horizRhsTol)
{
    m_maxOrder = a_maxOrder;
    m_eps = a_eps;
    m_hang = a_hang;
    m_normType = a_normType;
    m_verbosity = a_verbosity;
    m_horiz_rhsTol = a_horizRhsTol;
}


// -----------------------------------------------------------------------------
// Sets the horizontal MG solver's parameters.
// -----------------------------------------------------------------------------
void LevelLepticSolver::setHorizMGParameters (const int  a_imin,
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
                                              const int  a_verbosity)
{
    m_horiz_imin = a_imin;
    m_horiz_imax = a_imax;
    m_horiz_numSmoothDown = a_numSmoothDown;
    m_horiz_numSmoothUp = a_numSmoothUp;
    m_horiz_numSmoothBottom = a_numSmoothBottom;
    m_horiz_numSmoothPreCond = a_numSmoothPreCond;
    m_horiz_relaxType = a_relaxType;
    m_horiz_numMG = 1; // I'm hardcoding this because 1 is all that seems to work.
    m_horiz_maxDepth = a_maxDepth;
    m_horiz_eps = a_eps;
    m_horiz_hang = a_hang;
    m_horiz_normThresh = a_normThresh;
    m_horiz_verbosity = a_verbosity;
}


// -----------------------------------------------------------------------------
// Sets the horizontal bottom solver's parameters.
// -----------------------------------------------------------------------------
void LevelLepticSolver::setHorizBottomParameters (const int  a_imax,
                                                  const Real a_eps,
                                                  const int  a_numRestarts,
                                                  const Real a_hang,
                                                  const int  a_verbosity)
{
    m_horizBottom_imax = a_imax;
    m_horizBottom_eps = a_eps;
    m_horizBottom_numRestarts = a_numRestarts;
    m_horizBottom_hang = a_hang;
    m_horizBottom_verbosity = a_verbosity;
}


// -----------------------------------------------------------------------------
// Solve L[phi] = rhs.
// This is an override of the pure virtual LinearSolver function.
// -----------------------------------------------------------------------------
void LevelLepticSolver::solve (LevelData<FArrayBox>&       a_phi,
                               const LevelData<FArrayBox>& a_rhs)
{
    CH_TIME("LevelLepticSolver::solve");

    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(a_phi.nComp() == 1);
    CH_assert(a_rhs.nComp() == 1);
    CH_assert(a_phi.getBoxes().physDomain() == m_domain);
    CH_assert(a_phi.getBoxes().compatible(m_origJgupPtr->getBoxes()));
    CH_assert(a_rhs.getBoxes().compatible(m_origJgupPtr->getBoxes()));

    int maxOrder = m_maxOrder;
    const DisjointBoxLayout& origGrids = a_phi.getBoxes();
    const Real H = LevelGeometry::getDomainLength(SpaceDim-1);
    CH_assert(m_dx[SpaceDim-1] * m_domain.size(SpaceDim-1) == H);

    const bool useKrylovSolver = true;
    BiCGStabSolver<LevelData<FArrayBox> > krylovSolver;
    krylovSolver.define(&*m_opPtr, true);
    krylovSolver.m_imax = 80;
    krylovSolver.m_verbosity = 0;
    krylovSolver.m_numRestarts = 0;
    krylovSolver.m_normType = m_normType;


    // Initialization ----------------------------------------------------------

    // Allocation.
    LevelData<FArrayBox> phiTotal(m_grids, 1, IntVect::Unit);   // The cummulative solution.
    LevelData<FArrayBox> vertPhi(m_grids, 1, IntVect::Unit);    // The current vertical solution.
    RefCountedPtr<LevelData<FArrayBox> > rhsPtr(new LevelData<FArrayBox>(m_grids, 1)); // The current vertical rhs.
    RefCountedPtr<LevelData<FArrayBox> > tmpRhsPtr(new LevelData<FArrayBox>(m_grids, 1)); // Temp storage

    BoxLayoutData<FArrayBox> excess;        // The excess function.
    BoxLayoutData<FArrayBox> flatRhs;       // Used to compute the horizontal rhs.
    LevelData<FArrayBox>     horizPhi;      // The horizontal problem's solution.
    LevelData<FArrayBox>     horizRhs;      // The horizontal problem's rhs.

    if (m_doHorizSolve) {
        excess  .define(m_flatGrids, 1);             // The excess function.
        flatRhs .define(m_flatGrids, 1);             // Used to compute the horizontal rhs.
        horizPhi.define(m_horizGrids, 1, s_hmask);   // The horizontal problem's solution.
        horizRhs.define(m_horizGrids, 1);            // The horizontal problem's rhs.
    }

    // These avoid the LevelData::clear() issues.
    bool useExcess = m_doHorizSolve;
    bool useHorizPhi = m_doHorizSolve;

    // Send J*residual to vertical grids.
    {
        LevelData<FArrayBox> res(origGrids, 1);
        if (m_crsePhiPtr == NULL) {
            m_origOpPtr->residual(res, a_phi, a_rhs, m_homogeneous);
        } else {
            typedef MappedAMRLevelOp<LevelData<FArrayBox> > AMROpType;
            AMROpType* castOp = dynamic_cast<AMROpType*>(m_origOpPtr);
            CH_assert(castOp != NULL);
            castOp->AMRResidualNF(res, a_phi, *m_crsePhiPtr, a_rhs, m_homogeneous);
        }

        DataIterator dito = origGrids.dataIterator();
        for (dito.reset(); dito.ok(); ++dito) {
            res[dito].divide((*m_origJinvPtr)[dito]);
        }

        res.copyTo(*rhsPtr, m_origToVertCopier);
    }

    // Initialize convergence metrics
    m_resNorms.clear();
    m_resNorms.reserve(maxOrder+1);
    Real resNorm = m_opPtr->norm(*rhsPtr, m_normType);
    m_resNorms.push_back(resNorm);
    if (m_verbosity >= 4) {
        pout() << endl;
        pout() << " Absolute initial residual norm = "
               << std::scientific << resNorm << "\n";
        pout() << " Relative residual norm = 1.0" << endl;
    }

    // Initialize BC values.
    BoundaryData<Real> bdryData(m_grids, m_domain);
    CH_assert(!bdryData.isFlat());
    bdryData.setVal(0.0);
#ifndef NDEBUG
    // Check for consistency. Integral[rhs]-Integral[fluxes] should be zero.
    // This only needs to be done if we have Neum or periodic BCs all around.
    if (m_doHorizSolve && m_horizRemoveAvg && m_verbosity >= 4) {
        if (m_verbosity >= 4) {
            pout() << "\t" << SpaceDim << "D consistency check (should be zero) = "
                   << bdryData.consistencyCheck(*rhsPtr, m_dx) << endl;
        }
    }
#endif

    // Initialize solution. Remember, this will accumulate, so start with zero.
    setValLevel(phiTotal, 0.0);


    // Solve -------------------------------------------------------------------
    for (int order = 0; order <= maxOrder; ++order) {
        if (m_verbosity >= 4) {
            pout() << "O(eps^" << order << "):" << endl;
        }

        // Set up the BCs for all-Neum problem...
        // Do this even if we have abandoned the horizontal solves!
        if (m_doHorizSolve) {
            // Compute -Grad[phi_{order-1}]^z or leave BCs zero.
            // NOTE: This clobbers vertPhi!
            if (order >= 1) {
                this->levelVertHorizGradient(bdryData, vertPhi, -1.0);
            }

            // Add previous excess to upper BC
            if (order >= 1 && useExcess) {
                bdryData.vertPlus(excess, 1.0, Side::Hi);
            }

            // Compute new excess or delete if no longer needed.
            if (useExcess) {
                this->computeVerticalExcess(excess, *rhsPtr, bdryData);

                // Report size of excess
                const Real excessNorm = norm(excess, excess.interval(), m_normType);
                if (m_verbosity >= 5) {
                    pout() << "\tAbsolute excess norm = " << excessNorm << endl;
                }

                // The algo does not require the excess function after O(1).
                if (order == 1) {
                    useExcess = false;
                }
            }

            // Remove excess from upper vertical BC.
            if (order == 0 && useExcess) {
                bdryData.vertPlus(excess, -1.0, Side::Hi);
            }
            // At this point, BCs and rhs should jive. The vertical line solver
            // will check for us.

        } // end if Neum-Neum BCs

        // Solve the vertical problem and accumulate correction...
        this->verticalLineSolver(vertPhi, *rhsPtr, m_vertBCTypes, bdryData);

        // Set up the horizontal problem if needed...
        if (useHorizPhi) {
            // Compute \partial_m \bar{Jg^{mj} \partial_j \phi_p^v}. This exchanges vertPhi.
            this->computeHorizRHS(flatRhs, vertPhi, bdryData);

            // Subtract excess/H from horizontal rhs if needed.
            if (useExcess) {
                for (int idx = 0; idx < m_flatDI.size(); ++idx) {
                    const DataIndex& di = m_flatDI[idx];
                    flatRhs[di].plus(excess[di], -1.0/H);
                }
            }

            // Send flatRhs to m_horizGrids.
            setValLevel(horizRhs, 0.0);
            flatRhs.addTo(flatRhs.interval(), horizRhs, horizRhs.interval(), m_horizDomain, m_flatToHorizCopier);

            // TODO: Find a better scaling for the horizontal norm.
            Real horizRhsNorm = norm(horizRhs, horizRhs.interval(), m_normType);
            if (m_verbosity >= 5) {
                pout() << "\tHorizontal relative rhs norm = " << horizRhsNorm / m_resNorms[0] << endl;
            }

            if (m_horiz_rhsTol * m_resNorms[0] > horizRhsNorm) {
                if (m_verbosity >= 4) {
                    pout() << "\tAbandoning horizontal solves." << endl;
                }
                useHorizPhi = false;
            }
        }

        // Solve the horizontal problem and accumulate the correction...
        if (useHorizPhi) {
            this->horizontalSolver(horizPhi, horizRhs);
            this->addHorizontalCorrection(vertPhi, horizPhi);
            // At this point, vertPhi holds phi_p^v + phi_p^h.
        }


        // Finalize order...
        {
            // Check norms. This exchanges vertPhi.
            m_opPtr->residual(*tmpRhsPtr, vertPhi, *rhsPtr, true);
            resNorm = m_opPtr->norm(*tmpRhsPtr, m_normType);

            Real relResNorm = resNorm / m_resNorms[0];
            Real prevRelResNorm = m_resNorms.back() / m_resNorms[0];
            Real redu = prevRelResNorm - relResNorm;

            if (m_verbosity >= 4) {
                pout() << " Relative residual norm = " << relResNorm << endl;
            }

            // If we are not converging, use a Krylov solver.
            if (redu <= m_hang && useKrylovSolver) {
                // Solve with initial guess = 0.
                m_opPtr->setToZero(vertPhi);
                krylovSolver.solve(vertPhi, *rhsPtr);

                // Apply the Krylov solver's correction.
                m_opPtr->residual(*tmpRhsPtr, vertPhi, *rhsPtr, true);
                resNorm = m_opPtr->norm(*tmpRhsPtr, m_normType);
                relResNorm = resNorm / m_resNorms[0];

                if (m_verbosity >= 4) {
                    pout() << " Relative residual norm (krylov solver) = " << relResNorm << endl;
                }
            }

            // Set up RHS.
            std::swap(rhsPtr, tmpRhsPtr);
            m_resNorms.push_back(resNorm);

            // Check convergence status.
            redu = prevRelResNorm - relResNorm;
            if (redu > m_hang) {
                // We are converging. Use correction.
                this->addVerticalCorrection(phiTotal, vertPhi);

                if (order < maxOrder-1) {
                    // We are converging and still have more iteratin' to do.
                    m_exitStatus = ExitStatus::CONVERGE;
                } else {
                    // We are converging, but we reached max iters.
                    m_exitStatus = ExitStatus::ITER;
                }

            } else if (-redu > m_hang) {
                // We are diverging.
                if (m_verbosity >= 4) {
                    pout() << " We are diverging. Abandoning solve." << endl;
                }
                if (order == 0) {
                    m_exitStatus = ExitStatus::KABOOM;
                } else {
                    m_exitStatus = ExitStatus::DIVERGE;
                }
                break;

            } else {
                // We are hanging.
                if (m_verbosity >= 4) {
                    pout() << " We are hanging. Abandoning solve." << endl;
                }
                if (order == 0) {
                    m_exitStatus = ExitStatus::KABOOM;
                } else {
                    m_exitStatus = ExitStatus::HANG;
                }
                break;
            }

            // Should we turn off the horizontal solver?
            if (LevelGeometry::isDiagonal()) {
                useHorizPhi = false;
            }
        }
    } // end loop over orders

    // Finale ------------------------------------------------------------------

    if (m_exitStatus != ExitStatus::KABOOM) {
        // Send result back to original holder.
        LevelData<FArrayBox> corr(origGrids, 1);
        phiTotal.copyTo(corr, m_vertToOrigCopier);

        for (DataIterator dit(origGrids); dit.ok(); ++dit) {
            a_phi[dit].plus(corr[dit], 1.0);
        }

        // Inform user of the result if we didn't already do so.
        if (m_verbosity == 3) {
            pout() << " Final relative residual norm = "
                   << m_resNorms.back() / m_resNorms[0] << endl;
        }

    } else {
        // Report the blow up.
        if (m_verbosity == 3) {
            pout() << " Leptic solver blew up." << endl;
        }
    }
}


// -----------------------------------------------------------------------------
// Inform the user why we stopped solving.
// This returns a human-readable status.
// -----------------------------------------------------------------------------
char* LevelLepticSolver::exitStatusStr () const
{
    if (m_exitStatus == ExitStatus::NONE) return "NONE";
    else if (m_exitStatus == ExitStatus::CONVERGE) return "CONVERGE";
    else if (m_exitStatus == ExitStatus::ITER) return "ITER";
    else if (m_exitStatus == ExitStatus::HANG) return "HANG";
    else if (m_exitStatus == ExitStatus::DIVERGE) return "DIVERGE";
    else if (m_exitStatus == ExitStatus::KABOOM) return "KABOOM";

    return "UNKNOWN!";
}


// -----------------------------------------------------------------------------
// Compute -\partial_m \bar{Jgup^{mj} \partial_j \phi}.
// This does set BCs and exchanges phi. This does not subtract excess/H.
// Careful! a_rhs will be defined over m_flatGrids, not m_horizGrids!
// -----------------------------------------------------------------------------
void LevelLepticSolver::computeHorizRHS (BoxLayoutData<FArrayBox>& a_rhs,
                                         LevelData<FArrayBox>&     a_phi,
                                         const BoundaryData<Real>& a_bdryData) const
{
    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(!a_bdryData.isFlat());
    CH_assert(a_phi.getBoxes() == m_grids);
    CH_assert(a_rhs.boxLayout() == m_flatGrids);
    CH_assert(a_phi.nComp() == 1);
    CH_assert(a_rhs.nComp() == 1);

    // Gather geometric info
    const Box& domBox = m_domain.domainBox();
    const int Nz = domBox.size(SpaceDim-1);

    Vector<Box> FCDomInterior(SpaceDim-1, domBox);
    for (int FCdir = 0; FCdir < SpaceDim-1; ++FCdir) {
        FCDomInterior[FCdir].surroundingNodes(FCdir);
        FCDomInterior[FCdir].grow(FCdir, -1);
    }

    // Fill all of phi's ghosts
    extrapAllGhosts(a_phi, 2); // TODO: Should we use a_bdryData?
    homogeneousCFInterp(a_phi, m_dx, m_dxCrse, m_CFRegion);
    a_phi.exchange(m_exCopier);

    if (!LevelGeometry::isDiagonal()) {
        ExtrapolateCFEV(a_phi, m_CFRegion, 2);
        // for (dit.reset(); dit.ok(); ++dit) {
        //     ExtrapolateBCNoEV(a_phi[dit], m_grids[dit], m_domain, 2);
        // }
        a_phi.exchange(m_exCornerCopier);
    }

    // Initialize rhs to zero.
    setValLevel(a_rhs, 0.0);

    // Loop over grids
    for (int idx = 0; idx < m_flatDI.size(); ++idx) {
        const DataIndex& di = m_flatDI[idx];
        FArrayBox& rhsFAB = a_rhs[di];
        const FArrayBox& phiFAB = a_phi[di];
        const Box& valid = m_grids[di];

        // Loop over dirs accumulating -div[fluxes].
        for (int FCdir = 0; FCdir < SpaceDim-1; ++FCdir) {
            const FArrayBox& JgupFAB = (*m_JgupPtr)[di][FCdir];
            const Box faceBox = surroundingNodes(valid, FCdir);
            const Box interiorFaceBox = faceBox & FCDomInterior[FCdir];
            const Box flatFaceBox = surroundingNodes(m_flatGrids[di], FCdir);

            // Compute gradient...
            FArrayBox gradPhiFAB(faceBox, 1);
            debugInit(gradPhiFAB);

            if (LevelGeometry::isDiagonal()) {
                // Use simpler orthogonal version
                const Real dxDir = m_dx[FCdir];
                FORT_MAPPEDMACGRADORTHO (
                    CHF_FRA1(gradPhiFAB,0),
                    CHF_CONST_FRA1(phiFAB,0),
                    CHF_CONST_FRA1(phiFAB,0),
                    CHF_CONST_FRA(JgupFAB),
                    CHF_BOX(interiorFaceBox),
                    CHF_CONST_REAL(dxDir),
                    CHF_INT(FCdir),
                    CHF_INT(FCdir));

            } else {
                // Use full non-orthogonal version
                FORT_MAPPEDMACGRAD (
                    CHF_FRA1(gradPhiFAB,0),
                    CHF_CONST_FRA1(phiFAB,0),
                    CHF_CONST_FRA1(phiFAB,0),
                    CHF_CONST_FRA(JgupFAB),
                    CHF_BOX(interiorFaceBox),
                    CHF_CONST_REALVECT(m_dx),
                    CHF_INT(FCdir),
                    CHF_INT(FCdir));
            }

            // Set gradient at boundaries
            for (SideIterator sit; sit.ok(); ++sit) {
                const FArrayBox& bcFAB = a_bdryData.getData(di, FCdir, sit());
                if (!bcFAB.box().isEmpty()) {
                    gradPhiFAB.copy(bcFAB);
                }
            }

            // Compute vertical average of gradient...
            FArrayBox avgGradPhiFAB(flatFaceBox, 1);
            avgGradPhiFAB.setVal(0.0); // Integral will accumulate.

            const IntVect flatShift = flatFaceBox.smallEnd() * s_vmask;
            const Real dzScale = 1.0 / Real(Nz);

            FORT_UNMAPPEDVERTINTEGRAL(
                CHF_FRA1_SHIFT(avgGradPhiFAB, 0, flatShift),
                CHF_CONST_FRA1(gradPhiFAB,0),
                CHF_BOX(faceBox),
                CHF_CONST_REAL(dzScale));

            // Accumulate -divergence...
            const Real dxScale = -1.0 / m_dx[FCdir];
            const Box& CCregion = rhsFAB.box();

            FORT_LEPTICACCUMDIV (
                CHF_FRA1(rhsFAB,0),
                CHF_CONST_FRA1(avgGradPhiFAB,0),
                CHF_BOX(CCregion),
                CHF_CONST_REAL(dxScale),
                CHF_CONST_INT(FCdir));

        } // end loop over directions (FCdir)
    } // end loop over grids (idx)
}


// -----------------------------------------------------------------------------
// Computes scale*Jgup^{SpaceDim-1,m}*D[phi]/dx^[m], m=[0,SpaceDim-1).
// Used to generate BCs for the vertical problems.
// - No ghosts need to be set. All ghosts will be extrapolated.
// - It only makes sense to use this function on Neum-Neum problems.
// - It will be up to you to add the excess if needed.
// -----------------------------------------------------------------------------
void LevelLepticSolver::levelVertHorizGradient (BoundaryData<Real>&   a_vertBdryData,
                                                LevelData<FArrayBox>& a_phi,
                                                const Real            a_scale) const
{
    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(a_phi.getBoxes() == m_grids);

    // // This is the correct way to do things, but the extrap block in the loop
    // // is non-blocking, fast, and seems to work well.
    // extrapAllGhosts(a_phi, 2);
    // a_phi.exchange(m_exCopier);
    // a_phi.exchange(m_exCornerCopier);

    // Initialize to zero.
    a_vertBdryData.setVal(0.0);

    // If the metric is orthogonal, just do the easy thing.
    if (LevelGeometry::isDiagonal()) return;

    SideIterator sit;
    for (sit.reset(); sit.ok(); ++sit) {
        const Side::LoHiSide& iside = sit();
        const int sideComp = int(iside);
        const int isign = sign(iside);

        const Box& domBox = m_domain.domainBox();
        const Box domBdry = bdryBox(domBox, SpaceDim-1, iside, 1);

        DataIterator dit = a_phi.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            // Create references
            FArrayBox& bcValsFAB = a_vertBdryData.getData(dit(), SpaceDim-1, iside);
            FArrayBox& phiFAB = a_phi[dit];
            const FArrayBox& JgupFAB = (*m_JgupPtr)[dit][SpaceDim-1];
            const Box destFCBox = bcValsFAB.box();
            CH_assert(destFCBox == bdryBox(m_grids[dit], SpaceDim-1, iside, 1));

            // Extrapolate ghosts. Just do this in-place since we
            // won't need phi from the previous order anymore.
            const int extrapOrder = 2;
            Box domValid = phiFAB.box();
            domValid &= domBox;
            ExtrapolateFaceAndCopy(phiFAB, phiFAB, domValid, SpaceDim-1, iside, extrapOrder);  // TODO: Should we extrap edges?

            // Shift bcValsFAB
            VertShifter<Real> shifty(bcValsFAB, destFCBox);

            // Make sure centerings are correct
            CH_assert(bcValsFAB.box().contains(destFCBox));
            CH_assert(JgupFAB.box().contains(destFCBox));

            // Compute!
            FORT_LEPTICVERTHORIZGRAD (
                CHF_FRA1(bcValsFAB, 0),
                CHF_CONST_FRA1(phiFAB, 0),
                CHF_CONST_FRA(JgupFAB),
                CHF_BOX(destFCBox),
                CHF_CONST_INT(isign),
                CHF_CONST_REALVECT(m_dx),
                CHF_CONST_REAL(a_scale));

            // Restore bcValsFAB.
            shifty.restore();
        }
    }

    // Remember how I said we wouldn't need the contents phi anymore?
    // Well, let's make sure I am right!
    debugInitLevel(a_phi);
}


// -----------------------------------------------------------------------------
// Computes the excess function, excess = hiNeumBC - loNeumBC - Integral[rhs].
// -----------------------------------------------------------------------------
void LevelLepticSolver::computeVerticalExcess (BoxLayoutData<FArrayBox>&   a_excess,
                                               const LevelData<FArrayBox>& a_rhs,
                                               const BoundaryData<Real>&   a_bdryData) const
{
    // Sanity checks
    CH_assert(m_isDefined);

    CH_assert(m_doHorizSolve);
    CH_assert(a_excess.isDefined());

    CH_assert(a_excess.boxLayout() == m_flatGrids);
    CH_assert(a_rhs.getBoxes() == m_grids);

    CH_assert(a_excess.nComp() == 1);
    CH_assert(a_rhs   .nComp() == 1);

    // Initialize phi with a bogus value
    debugInitLevel(a_excess);

    // Set the excess to zero in regions that don't need it.
    setToZero(a_excess, m_flatDIComplement);

    // Loop over grids that require an excess
    for (int idx = 0; idx < m_flatDI.size(); ++idx) {
        const DataIndex& di = m_flatDI[idx];
        FArrayBox& excessFAB = a_excess[di];

        { // Copy hi BC
            const FArrayBox& hiBCFAB = a_bdryData.getData(di, SpaceDim-1, Side::Hi);
            CH_assert(!hiBCFAB.box().isEmpty());

            VertShifter<Real> shifty(excessFAB, hiBCFAB);
            excessFAB.copy(hiBCFAB);
            shifty.restore();
        }
        { // Subtract lo BC
            const FArrayBox& loBCFAB = a_bdryData.getData(di, SpaceDim-1, Side::Lo);
            CH_assert(!loBCFAB.box().isEmpty());

            VertShifter<Real> shifty(excessFAB, loBCFAB);
            excessFAB.plus(loBCFAB, -1.0);
            shifty.restore();
        }
        { // Subtract integral
            const FArrayBox& rhsFAB = a_rhs[di];
            const Box& valid = m_grids[di];
            const Real dzScale = -1.0 * m_dx[SpaceDim-1];
            const IntVect shift = excessFAB.box().smallEnd() * s_vmask;

            FORT_UNMAPPEDVERTINTEGRAL(
                CHF_FRA1_SHIFT(excessFAB,0,shift),
                CHF_CONST_FRA1(rhsFAB,0),
                CHF_BOX(valid),
                CHF_CONST_REAL(dzScale));
        }
    }
}


// -----------------------------------------------------------------------------
// Computes the vertical solutions.
// vertRhs will be temporarily altered, but then restored.
// No BCs will be set on vertPhi and this is a non-blocking function.
// -----------------------------------------------------------------------------
void LevelLepticSolver::verticalLineSolver (LevelData<FArrayBox>&            a_vertPhi,
                                            LevelData<FArrayBox>&            a_vertRhs,
                                            const LayoutData<Tuple<int,2> >& a_vertBCTypes,
                                            const BoundaryData<Real>&        a_bdryData) const
{
    CH_TIME("LevelLepticSolver::verticalLineSolver");

    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(!m_JgupPtr.isNull());

    CH_assert(a_vertPhi.ghostVect()[CH_SPACEDIM-1] >= 1);

    CH_assert(a_vertPhi    .getBoxes()  == m_grids);
    CH_assert(a_vertRhs    .getBoxes()  == m_grids);
    CH_assert(a_vertBCTypes.boxLayout() == m_grids);
    CH_assert(m_JgupPtr->getBoxes().compatible(m_grids));

    CH_assert(a_vertPhi    .nComp() == 1);
    CH_assert(a_vertRhs    .nComp() == 1);

    // Initialize phi with a bogus value
    debugInitLevel(a_vertPhi);

    // Gather geometric info
    const Real dz = m_dx[SpaceDim-1];
    const Real dzCrse = m_dxCrse[SpaceDim-1];
    DataIterator dit = m_grids.dataIterator();
    SideIterator sit;

#ifndef NDEBUG
    // If we have Neum-Neum BCs, check that each Poisson problem is consistent.
    if (m_doHorizSolve && m_verbosity >= 5) {
        BoxLayoutData<FArrayBox> consistency(m_flatGrids, 1);
        this->computeVerticalExcess(consistency, a_vertRhs, a_bdryData);
        const Real consistencyNorm = norm(consistency, consistency.interval(), m_normType);
        pout() << "\tVertical problem consistency check (should be zero) = " << consistencyNorm << endl;
    }
#endif

    // Loop over grids and solve.
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& phiFAB = a_vertPhi[dit];
        FArrayBox& rhsFAB = a_vertRhs[dit];
        const Tuple<int,2>& bcType = a_vertBCTypes[dit];
        const Box& valid = m_grids[dit];
        const int Nz = valid.size(SpaceDim-1);

        const Interval vertInt(SpaceDim-1, SpaceDim-1);
        const FluxBox& JgupFlub = (*m_JgupPtr)[dit];
        const FArrayBox JgzzFAB(vertInt, (FArrayBox&)JgupFlub[SpaceDim-1]);

        // Sanity checks
        CH_assert(rhsFAB.box().size(SpaceDim-1) == Nz);
        CH_assert(phiFAB.box().size(SpaceDim-1) == Nz + 2);
        CH_assert(phiFAB.box().contains(valid));
        CH_assert(rhsFAB.box().contains(valid));
        CH_assert(enclosedCells(JgzzFAB.box()).contains(valid));

        // A few points are worth mentioning.
        // 1. The BCTypes can not change from vertical line to vertical line throughout valid.
        //    If they do, then m_grids are ill-formed. At this point, we won't check.
        // 2. If the BCType is Neum, then we roll the BC values at the upper boundary into
        //    the rhsFAB. For all other BCs, we can assume they are homog and do nothing.
        // 3. For periodic BCs, we will throw an error for now. It may be better to use
        //    a 1D spectral solver in that case.
        // 4. For all other types of BC, specifically BCType::Undefined, we again throw
        //    an error.

        // Gather BCType info for this box and if needed, roll BCs in using
        // rhs = rhs -/+ NeumBCVal/dz.
        CH_assert(!m_domain.isPeriodic(SpaceDim-1));

        Vector<RefCountedPtr<FArrayBox> > rollInFAB(2);

        for (sit.reset(); sit.ok(); ++sit) {
            // Create references for convenience
            const Side::LoHiSide iside = sit();
            const int sideComp = int(iside);
            const int isign = sign(iside);

            if (bcType[sideComp] == BCType::Neum) {
                // Roll in Neumann BCs...

                // This is where the roll-in will happen.
                Box destBox = adjCellBox(valid, SpaceDim-1, iside, 1);
                destBox.shift(SpaceDim-1, -isign);
                destBox.shiftHalf(SpaceDim-1, isign);

                rollInFAB[sideComp] = RefCountedPtr<FArrayBox>(new FArrayBox);
                rollInFAB[sideComp]->define(destBox, 1);
                rollInFAB[sideComp]->setVal(0.0);

                // This is where the BC values live.
                const FArrayBox& bcValsFAB = a_bdryData.getData(dit(), SpaceDim-1, iside);
                const Box& srcBox = bcValsFAB.box();
                CH_assert(srcBox.size() == destBox.size());

                // Calculate roll-in
                const Real scale = -Real(isign) / dz;
                rollInFAB[sideComp]->plus(bcValsFAB, srcBox, destBox, scale, 0, 0, 1);

                // Set to rhsFAB centering
                rollInFAB[sideComp]->shiftHalf(SpaceDim-1, -isign);
            }
        } // end loop over sides (sit)

        // Roll in BCs if necessary
        if (!rollInFAB[0].isNull()) rhsFAB.plus(*rollInFAB[0], 1.0);
        if (!rollInFAB[1].isNull()) rhsFAB.plus(*rollInFAB[1], 1.0);

        // Time to solve...
        if (bcType[0] == BCType::Neum && bcType[1] == BCType::Neum) {
            // Use the modified homogeneous Neumann tridiagonal solver.
            const int vertDir = SpaceDim-1;
            Box bottomBox = adjCellLo(valid, vertDir, 1);
            bottomBox.shift(vertDir, 1);

            FORT_TRIDIAGPOISSONNN1DFAB (
                CHF_FRA(phiFAB),
                CHF_CONST_FRA(rhsFAB),
                CHF_CONST_FRA1(JgzzFAB,0),
                CHF_BOX(bottomBox),
                CHF_CONST_INT(Nz),
                CHF_CONST_REAL(dz),
                CHF_CONST_INT(vertDir));

        } else {
            // Just use regular homogeneous tridiagonal solver.

            // Create workspace for Fortran
            Box workspace(IntVect::Unit, IntVect(D_DECL(Nz,1,1)));
            BaseFab<Real> D(workspace, 1);
            BaseFab<Real> B(workspace, 1);

            workspace.setSmall(0, 2);
            BaseFab<Real> DU(workspace, 1);

            workspace.setSmall(0, 1);
            workspace.setBig(0, Nz-1);
            BaseFab<Real> DL(workspace, 1);

            const IntVect validShift = valid.smallEnd();

            // Solve!
            FORT_LEPTICLAPACKVERTICALSOLVER(
                CHF_FRA1_SHIFT(phiFAB, 0, validShift),
                CHF_CONST_FRA1_SHIFT(rhsFAB, 0, validShift),
                CHF_CONST_FRA1_SHIFT(JgzzFAB, 0, validShift),
                CHF_BOX_SHIFT(valid, validShift),
                CHF_CONST_REAL(dz),
                CHF_CONST_REAL(dzCrse),
                CHF_FRA1(DU,0),
                CHF_FRA1(D,0),
                CHF_FRA1(DL,0),
                CHF_FRA1(B,0),
                CHF_CONST_INT(bcType[0]),
                CHF_CONST_INT(bcType[1]));

            // Free memory
            B.clear();
            D.clear();
            DU.clear();
            DL.clear();
        }

        // Unroll BCs if necessary
        if (!rollInFAB[0].isNull()) rhsFAB.plus(*rollInFAB[0], -1.0);
        if (!rollInFAB[1].isNull()) rhsFAB.plus(*rollInFAB[1], -1.0);

    } // end loop over grids (dit)

    // writeLevelHDF5(a_vertRhs, 0.0, false);
    // writeLevelHDF5(a_vertPhi, 0.0, false);
}


// -----------------------------------------------------------------------------
// Computes the horizontal solutions.
// -----------------------------------------------------------------------------
void LevelLepticSolver::horizontalSolver (LevelData<FArrayBox>&       a_phi,
                                          const LevelData<FArrayBox>& a_rhs)
{
    CH_TIME("LevelLepticSolver::horizontalSolver");

    // Sanity checks
    CH_assert(m_isDefined);
    CH_assert(m_doHorizSolve);
    CH_assert(a_phi.getBoxes() == m_horizGrids);
    CH_assert(a_rhs.getBoxes() == m_horizGrids);

#ifndef NDEBUG
    if (m_horizRemoveAvg && m_verbosity >= 5) {
        // Check for consistency. Integral[rhs]-Integral[fluxes] should be zero.
        BoundaryData<Real> horizBdryData(m_horizGrids, m_horizDomain);
        horizBdryData.setVal(0.0);
        Real consistency = horizBdryData.consistencyCheck(a_rhs, m_dx);
        pout() << "\tFlat, " << SpaceDim-1 << "D consistency check (should be zero) = "
               << consistency << endl;
    }
#endif

    Vector<LevelData<FArrayBox>*> horizPhi(1, &a_phi);
    Vector<LevelData<FArrayBox>*> horizRhs(1, const_cast<LevelData<FArrayBox>*>(&a_rhs));

    m_horizSolverPtr->solve(horizPhi,
                            horizRhs,
                            0,      // lMax
                            0,      // lBase
                            true,   // zero phi
                            true);  // force homog

    if (m_horizRemoveAvg) {
        this->setZeroAvg(a_phi);
    }
}


// -----------------------------------------------------------------------------
// Adds a vertical correction to the vertical solution.
// -----------------------------------------------------------------------------
void LevelLepticSolver::addVerticalCorrection (LevelData<FArrayBox>&       a_phiTotal,
                                               const LevelData<FArrayBox>& a_vertCor) const
{
    CH_assert(m_isDefined);
    CH_assert(a_phiTotal.getBoxes() == m_grids);
    CH_assert(a_vertCor .getBoxes() == m_grids);

    DataIterator dit = m_grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& totalFAB = a_phiTotal[dit];
        const FArrayBox& corFAB = a_vertCor[dit];
        const Box& valid = m_grids[dit];

        totalFAB.plus(corFAB, valid, 0, 0, 1);
    }
}


// -----------------------------------------------------------------------------
// Adds a horizontal correciton to the vertical solution.
// -----------------------------------------------------------------------------
void LevelLepticSolver::addHorizontalCorrection (LevelData<FArrayBox>&       a_vertPhi,
                                                 const LevelData<FArrayBox>& a_horizCor) const
{
    CH_assert(m_isDefined);
    CH_assert(a_vertPhi .nComp() == 1);
    CH_assert(a_horizCor.nComp() == 1);
    CH_assert(a_vertPhi. getBoxes() == m_grids);
    CH_assert(a_horizCor.getBoxes() == m_horizGrids);

    // Send correction to m_flatGrids.
    BoxLayoutData<FArrayBox> flatData(m_flatGrids, 1);
    debugInitLevel(flatData);
    a_horizCor.copyTo(flatData, m_horizToFlatCopier);

    // Add the extrusion to the vertical solution.
    for (int idx = 0; idx < m_flatDI.size(); ++idx) {
        const DataIndex& di = m_flatDI[idx];
        FArrayBox& phiFAB = a_vertPhi[di];
        const FArrayBox& flatFAB = flatData[di];
        const Box& fullRegion = m_grids[di];   // This was phiFAB.box() and caused an invalid memory access error.

        FORT_ADDEXTRUSION(
            CHF_FRA1(phiFAB, 0),
            CHF_CONST_FRA1(flatFAB, 0),
            CHF_BOX(fullRegion));
    }
}


// -----------------------------------------------------------------------------
// Collects each boxes vertical BC type. If all BCTypes are Neum or periodic,
// a_doHorizSolve will be set to true.
// -----------------------------------------------------------------------------
void LevelLepticSolver::gatherVerticalBCTypes (LayoutData<Tuple<int,2> >& a_BCTypes,
                                               bool&                      a_doHorizSolve,
                                               const ProblemDomain&       a_domain,
                                               CFRegion&                  a_CFRegion,
                                               const BCMethodHolder&      a_bc)
{
    // Gather grid info.
    const BoxLayout& grids = a_BCTypes.boxLayout();
    DataIterator dit = grids.dataIterator();
    const Box& domBox = a_domain.domainBox();
    const int isPeriodic = a_domain.isPeriodic(SpaceDim-1);

    // Gather BCTypes from the BCDescriptors.
    // NOTE: Is is assumed that if we do not have Neum BCs, then the fluxBCType is irrelevant.
    Tuple<int, 2> physBCType;
    {
        Tuple<int, 2> ghostBCType;
        ghostBCType[Side::Lo] = a_bc.getGhostDescriptor()[SpaceDim-1][Side::Lo];
        ghostBCType[Side::Hi] = a_bc.getGhostDescriptor()[SpaceDim-1][Side::Hi];

        Tuple<int, 2> fluxBCType;
        fluxBCType[Side::Lo] = a_bc.getFluxDescriptor()[SpaceDim-1][Side::Lo];
        fluxBCType[Side::Hi] = a_bc.getFluxDescriptor()[SpaceDim-1][Side::Hi];

        physBCType = ghostBCType;
        if (ghostBCType[Side::Lo] == BCType::Neum || fluxBCType[Side::Lo] == BCType::Neum) {
            physBCType[Side::Lo] = BCType::Neum;
        }
        if (ghostBCType[Side::Hi] == BCType::Neum || fluxBCType[Side::Hi] == BCType::Neum) {
            physBCType[Side::Hi] = BCType::Neum;
        }
    }

    for (SideIterator sit; sit.ok(); ++sit) {
        const Side::LoHiSide iside = sit();
        CH_assert(iside == 0 || iside == 1);
        const int domSideIdx = domBox.sideEnd(iside)[SpaceDim-1];

        for (dit.reset(); dit.ok(); ++dit) {
            int& thisBCType = a_BCTypes[dit][iside];
            const Box& valid = grids[dit];
            const int validSideIdx = valid.sideEnd(iside)[SpaceDim-1];

            // Initialize to error BCs.
            thisBCType = BCType::Undefined;

            if (validSideIdx == domSideIdx) {
                // We are at a physical boundary...

                if (isPeriodic) {
                    thisBCType = BCType::Periodic;
                } else {
                    thisBCType = physBCType[iside];
                    if (thisBCType != BCType::Neum) {
                        a_doHorizSolve = false;
                    }
                }

            } else {
                // We are NOT at a physical boundary...

                const Box ghostBox = adjCellBox(valid, SpaceDim-1, iside, 1);
                const CFIVS& cfivs = (iside == Side::Lo?
                                      a_CFRegion.loCFIVS(dit(), SpaceDim-1):
                                      a_CFRegion.hiCFIVS(dit(), SpaceDim-1));

                // If the cfivs is empty, we are not at a CF interface.
                // We throw an error because this implies decomposition in the vertical.
                if (cfivs.isEmpty()) {
                    MayDay::Error("Vertical grids are ill-formed");

                } else  if (cfivs.isPacked()) {
                    const Box& cfBox = cfivs.packedBox();
                    if (cfBox.contains(ghostBox)) {
                        thisBCType = BCType::CF;
                        a_doHorizSolve = false;
                    } else {
                        MayDay::Error("Vertical grids are ill-formed");
                    }

                } else {
                    const IntVectSet& ivs = cfivs.getIVS();
                    if (ivs.contains(ghostBox)) {
                        thisBCType = BCType::CF;
                        a_doHorizSolve = false;
                    } else {
                        MayDay::Error("Vertical grids are ill-formed");
                    }
                }
            } // end if region abuts physical boundary or not
        } // end loop over grids (dit)
    } // end loop over sides (sit)

    // Search for Neum-Neum BCs, which demand horizontal solves.
    a_doHorizSolve = false;
    for (dit.reset(); dit.ok(); ++dit) {
        const int loBCType = a_BCTypes[dit][0];
        const int hiBCType = a_BCTypes[dit][1];

        const bool loNotDiri = (loBCType == BCType::Neum || loBCType == BCType::Periodic);
        const bool hiNotDiri = (hiBCType == BCType::Neum || hiBCType == BCType::Periodic);

        if (loNotDiri && hiNotDiri) {
            a_doHorizSolve = true;
            break;
        }
    }

#if CH_MPI
    // If any of the processors reported true, then set all to true.
    int globalDoSolve;
    int localDoSolve = (a_doHorizSolve? 1: 0);

    int ierr = MPI_Allreduce(&localDoSolve, &globalDoSolve, 1, MPI_INT, MPI_SUM, AMRLESMeta::amrComm);
    if (ierr != MPI_SUCCESS) {
        std::ostringstream errmsg;
        errmsg << "MPI_Allreduce failed. Error " << ierr << std::endl;
        MayDay::Error(errmsg.str().c_str());
    }

    a_doHorizSolve = (globalDoSolve > 0);
#endif
}


// -----------------------------------------------------------------------------
// Static utility
// Sets all data referenced by the vector of data indices to zero.
// -----------------------------------------------------------------------------
void LevelLepticSolver::setToZero (BoxLayoutData<FArrayBox>& a_data,
                                   const Vector<DataIndex>&  a_vdi)
{
    for (int idx = 0; idx < a_vdi.size(); ++idx) {
        const DataIndex& di = a_vdi[idx];
        a_data[di].setVal(0.0);
    }
}


// -----------------------------------------------------------------------------
// Static utility
// Brings average[a_phi] to zero.
// -----------------------------------------------------------------------------
void LevelLepticSolver::setZeroAvg (LevelData<FArrayBox>& a_phi)
{
    const DisjointBoxLayout& grids = a_phi.getBoxes();
    const Box& domBox = grids.physDomain().domainBox();
    DataIterator dit = a_phi.dataIterator();

    const int numproc = numProc();
    const int thisproc = procID();
    const int srcproc = uniqueProc(SerialTask::compute);

    Vector<Real> sumVect(numProc());
    Real localSum  = 0.0;
    Real globalSum = 0.0;

    Vector<long> volVect(numProc());
    long localVol  = 0;
    long globalVol = 0;

    Real globalAvg = 0.0;

    // Collect sums on each processor.
    for (dit.reset(); dit.ok(); ++dit) {
        const FArrayBox& phiFAB = a_phi[dit];
        const Box& valid = grids[dit];

        localSum += phiFAB.sum(valid, 0, 1);
        localVol += valid.numPts();
    }

    // Gather / broadcast total average.
    gather(sumVect, localSum, srcproc);
    gather(volVect, localVol, srcproc);
    if(thisproc == srcproc) {
        for(int idx = 0; idx < numproc; ++idx) {
            globalSum += sumVect[idx];
            globalVol += volVect[idx];
        }
        globalAvg = globalSum / Real(globalVol);
    }
    broadcast(globalAvg, srcproc);

    // Remove average from a_phi.
    for (dit.reset(); dit.ok(); ++dit) {
        a_phi[dit] -= globalAvg;
    }
}

