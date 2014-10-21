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
#include "StratUtils.H"
#include "StratUtilsF_F.H"
#include "LevelGeometry.H"
#include "PhysBCUtil.H"
#include "LepticBoxUtils.H"
#include "Subspace.H"
#include "lapack.H"
#include "BoxIterator.H"
#include "Constants.H"
#include "SetValLevel.H"


// -----------------------------------------------------------------------------
// Computes N = sqrt(max|db/dz|) over the entire level.
// -----------------------------------------------------------------------------
Real computeMaxBVFreq (const LevelData<FluxBox>& a_B,
                       const LevelGeometry&      a_levGeo,
                       const Real                a_pow)
{
    CH_TIME("computeMaxBVFreq");

    // Sanity checks
    CH_assert(a_B.nComp() == 1);
    CH_assert(a_B.ghostVect().product() >= 1);

    // Create grid references, etc...
    const RealVect& dx = a_levGeo.getDx();
    const GeoSourceInterface& geoSource = *a_levGeo.getGeoSourcePtr();
    const DisjointBoxLayout& grids = a_B.getBoxes();
    DataIterator dit = grids.dataIterator();

    Real globalMax;
    Real localMax = 0.0;

    // Compute maxN on each grid
    for (dit.reset(); dit.ok(); ++dit) {
        const FluxBox& bFB = a_B[dit];
        const Box& valid = grids[dit];

        // Fill dXi^i/dz
        FArrayBox dXidzFAB(valid, SpaceDim);
        D_TERM(geoSource.fill_dXidx(dXidzFAB, 0, 0, SpaceDim-1, dx);,
               geoSource.fill_dXidx(dXidzFAB, 1, 1, SpaceDim-1, dx);,
               geoSource.fill_dXidx(dXidzFAB, 2, 2, SpaceDim-1, dx);)

        // Compute max(N)^pow
        Real maxN;
        if (SpaceDim == 2) {
            FORT_COMPUTEMAXBVFREQ2D(
                CHF_REAL(maxN),
                CHF_CONST_FRA1(bFB[0],0),
                CHF_CONST_FRA1(bFB[1],0),
                CHF_CONST_FRA(dXidzFAB),
                CHF_CONST_REALVECT(dx),
                CHF_BOX(valid),
                CHF_CONST_REAL(a_pow));
        } else {
            FORT_COMPUTEMAXBVFREQ3D(
                CHF_REAL(maxN),
                CHF_CONST_FRA1(bFB[0],0),
                CHF_CONST_FRA1(bFB[1],0),
                CHF_CONST_FRA1(bFB[2],0),
                CHF_CONST_FRA(dXidzFAB),
                CHF_CONST_REALVECT(dx),
                CHF_BOX(valid),
                CHF_CONST_REAL(a_pow));
        }

        if (localMax < maxN) localMax = maxN;
    }

    // Communicate. Find maxN over all procs.
#ifdef CH_MPI
    int result = MPI_Allreduce(&localMax, &globalMax, 1, MPI_CH_REAL,
                               MPI_MAX, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Communication error in computeDt");
    }
#else
    globalMax = localMax;
#endif

    // Return the max freq.
    return globalMax;
}


// -----------------------------------------------------------------------------
// Computes sqrt(-db/dz)^p over one grid.
// -----------------------------------------------------------------------------
void computeBVFreq (FArrayBox&           a_N,
                    const FluxBox&       a_B,
                    const Box&           a_destBox,
                    const LevelGeometry& a_levGeo,
                    const Real           a_pow)
{
    // Sanity checks
    CH_assert(a_N.nComp() == 1);
    CH_assert(a_B.nComp() == 1);
    CH_assert(a_destBox.type() == a_N.box().type());
    CH_assert(a_N.box().contains(a_destBox));
    CH_assert(a_B.box().contains(a_destBox));

    // Create grid references, etc...
    const RealVect& dx = a_levGeo.getDx();
    const GeoSourceInterface& geoSource = *a_levGeo.getGeoSourcePtr();

    // Fill dXi^i/dz
    FArrayBox dXidzFAB(a_destBox, SpaceDim);
    D_TERM(geoSource.fill_dXidx(dXidzFAB, 0, 0, SpaceDim-1, dx);,
           geoSource.fill_dXidx(dXidzFAB, 1, 1, SpaceDim-1, dx);,
           geoSource.fill_dXidx(dXidzFAB, 2, 2, SpaceDim-1, dx);)

    // Do the calculation
    if (SpaceDim == 2) {
        FORT_COMPUTEBVFREQ2D(
            CHF_FRA1(a_N,0),
            CHF_CONST_FRA1(a_B[0],0),
            CHF_CONST_FRA1(a_B[1],0),
            CHF_CONST_FRA(dXidzFAB),
            CHF_CONST_REALVECT(dx),
            CHF_BOX(a_destBox),
            CHF_CONST_REAL(a_pow));
    } else {
        FORT_COMPUTEBVFREQ3D(
            CHF_FRA1(a_N,0),
            CHF_CONST_FRA1(a_B[0],0),
            CHF_CONST_FRA1(a_B[1],0),
            CHF_CONST_FRA1(a_B[2],0),
            CHF_CONST_FRA(dXidzFAB),
            CHF_CONST_REALVECT(dx),
            CHF_BOX(a_destBox),
            CHF_CONST_REAL(a_pow));
    }
}


// -----------------------------------------------------------------------------
// Computes sqrt(-db/dz)^p over the entire level.
// -----------------------------------------------------------------------------
void computeBVFreq (LevelData<FArrayBox>&     a_N,
                    const LevelData<FluxBox>& a_B,
                    const LevelGeometry&      a_levGeo,
                    const Real                a_pow)
{
    // Sanity checks
    CH_assert(a_N.nComp() == 1);
    CH_assert(a_B.nComp() == 1);
    CH_assert(a_N.getBoxes() == a_B.getBoxes());
    CH_assert(a_N.getBoxes() == a_levGeo.getBoxes());

    // Create grid references, etc...
    const DisjointBoxLayout& grids = a_N.getBoxes();
    DataIterator dit = grids.dataIterator();

    // Do the calculation over each grid
    for (dit.reset(); dit.ok(); ++dit) {
        const Box destBox = grids[dit] & a_N[dit].box();
        computeBVFreq(a_N[dit], a_B[dit], destBox, a_levGeo, a_pow);
    }
}


// -----------------------------------------------------------------------------
// Finds the max speed, c0 (smallest Nsq/c^2) and the corresponding structure
// function, phi0 (corresponding eigenfunction) of phi'' + (Nsq/c^2) phi = 0.
// This level's grids must span the vertical domain.
// -----------------------------------------------------------------------------
void solveVertEigenProblem (LevelData<FArrayBox>&       a_c0,
                            LevelData<FArrayBox>&       a_phi0,
                            const LevelData<FArrayBox>& a_Nsq,
                            const LevelGeometry&        a_levGeo)
{
    CH_TIME("solveVertEigenProblem");

    // Create grid references, etc...
    const ProblemDomain& domain = a_levGeo.getDomain();
    const DisjointBoxLayout& grids = a_phi0.getBoxes();
    CH_assert(grids.physDomain() == domain);

    const IntVect nx = domain.size();
    const int nproc = numProc();

    // Split the domain and create a load-balanced set of grids.
    Vector<Box> vbox;
    LepticBoxUtils::createVerticalSolverGrids(vbox, grids.boxArray(), domain.domainBox());
    DisjointBoxLayout vertGrids;
    vertGrids.defineAndLoadBalance(vbox, NULL, domain);

    // Create data holders
    char dflag[1] = {'S'};
    static const Real abstol = 2.0 * dlamch_(dflag);
    const int N = nx[SpaceDim-1];
    const RealVect& dx = a_levGeo.getDx();
    const Real dz = dx[SpaceDim-1];
    const GeoSourceInterface& geoSource = *a_levGeo.getGeoSourcePtr();

#define ARRAY1D(a)   Box(IntVect(D_DECL((1),(1),(1))), IntVect(D_DECL((a),(1),(1)))), (1)
#define ARRAY2D(a,b) Box(IntVect(D_DECL((1),(1),(1))), IntVect(D_DECL((a),(b),(1)))), (1)
    BaseFab<int>  IFAIL(ARRAY1D(1));    // (1)
    BaseFab<int>  IWORK(ARRAY1D(5*N));  // (5*N)
    BaseFab<Real> AB(ARRAY2D(2,N));     // (2,N)
    BaseFab<Real> BB(ARRAY2D(1,N));     // (1,N)
    BaseFab<Real> Q(ARRAY2D(N,N));      // (N,N)
    BaseFab<Real> W(ARRAY1D(N));        // (N)
    BaseFab<Real> WORK(ARRAY1D(7*N));   // (7*N)
    BaseFab<Real> Z(ARRAY2D(N,N));      // (N,N)
#undef ARRAY1D
#undef ARRAY2D

    LevelData<FArrayBox> vertNsq(vertGrids, 1);
    LevelData<FArrayBox> vertPhi(vertGrids, 1);
    LevelData<FArrayBox> vertC0(vertGrids, 1);
#ifndef NDEBUG
    setValLevel(vertNsq, quietNAN);
    setValLevel(vertPhi, quietNAN);
    setValLevel(vertC0, quietNAN);
#endif

    // Send Nsq to vertGrids.
    a_Nsq.copyTo(vertNsq);

    // Loop over grids and solve the eigenproblem.
    DataIterator vertDit = vertGrids.dataIterator();
    for (vertDit.reset(); vertDit.ok(); ++vertDit) {
        FArrayBox& phiFAB = vertPhi[vertDit];
        FArrayBox& c0FAB = vertC0[vertDit];
        const FArrayBox& NsqFAB = vertNsq[vertDit];
        const Box& valid = vertGrids[vertDit];
        const Box flatValid = flattenBox(valid, SpaceDim-1);
        const Box FCvalid = surroundingNodes(valid, SpaceDim-1);

        // Fill dXi^i/dz
        FArrayBox jacFCFAB(FCvalid, 1);
        geoSource.fill_dXidx(jacFCFAB, 0, SpaceDim-1, SpaceDim-1, dx);

        FArrayBox jacCCFAB(valid, 1);
        geoSource.fill_dXidx(jacCCFAB, 0, SpaceDim-1, SpaceDim-1, dx);

        BoxIterator bit(flatValid);
        for (bit.reset(); bit.ok(); ++bit) {
            Real thisC0;
            IntVect loIV = bit();
            loIV[SpaceDim-1] = valid.smallEnd(SpaceDim-1);

            FORT_SOLVEVERTEIGENPROBLEM(
                CHF_REAL(thisC0),
                CHF_FRA1(phiFAB,0),
                CHF_CONST_FRA1(NsqFAB,0),
                CHF_CONST_FRA1(jacFCFAB,0),
                CHF_CONST_FRA1(jacCCFAB,0),
                CHF_CONST_INTVECT(loIV),
                CHF_CONST_INT(N),
                CHF_CONST_REAL(dz),
                CHF_CONST_REAL(abstol),
                CHF_FIA1(IFAIL,0),
                CHF_FIA1(IWORK,0),
                CHF_FRA1(AB,0),
                CHF_FRA1(BB,0),
                CHF_FRA1(Q,0),
                CHF_FRA1(W,0),
                CHF_FRA1(WORK,0),
                CHF_FRA1(Z,0));

            // Copy thisC0 to all points above loIV.
            Box vertLineBox(loIV,loIV);
            vertLineBox.setBig(SpaceDim-1, valid.bigEnd(SpaceDim-1));
            c0FAB.setVal(thisC0, vertLineBox, 0, 1);
        }
    }

    // Copy data to original grids
    vertPhi.copyTo(a_phi0);
    vertC0.copyTo(a_c0);
}



// These were just copied from old code...

// -----------------------------------------------------------------------------
// Computes N = sqrt(-db/dz) over one grid.
// This is a more general version...N and B can have (almost) any centering.
// TODO: Make a general mapped version that takes the place of all BV funcs.
// -----------------------------------------------------------------------------
void computeBVFreqCartesian (FArrayBox&           a_N,
                             const FArrayBox&     a_B,
                             const Box&           a_destBox,
                             const LevelGeometry& a_levGeo,
                             const Real           a_pow)
{
    // Sanity checks
    CH_assert(a_N.nComp() == 1);
    CH_assert(a_B.nComp() == 1);
    CH_assert(a_destBox.type() == a_N.box().type());
    CH_assert(a_N.contains(a_destBox));
    // CH_assert(strcmp(LevelGeometry::getCoorMapName(), "Cartesian") != 0);

    // Create grid references, etc...
    const IntVect NBoxType = a_N.box().type();
    const IntVect BBoxType = a_B.box().type();
    const Real dz = a_levGeo.getDx()[CH_SPACEDIM-1];

    // Do the calculation
    FORT_COMPUTEBVFREQCARTESIAN(
        CHF_FRA1(a_N,0),
        CHF_CONST_FRA1(a_B,0),
        CHF_CONST_REAL(dz),
        CHF_BOX(a_destBox),
        CHF_CONST_INTVECT(NBoxType),
        CHF_CONST_INTVECT(BBoxType),
        CHF_CONST_REAL(a_pow));
}


// -----------------------------------------------------------------------------
// Finds the max speed, c0 (smallest eigenvalue) and the corresponding structure
// function, phi0 (corresponding eigenfunction) of phi'' + (Nsq/c^2) phi = 0.
// -----------------------------------------------------------------------------
void solveVertEigenProblemCartesian (Real&            a_c0,
                                     FArrayBox&       a_phi0,
                                     const FArrayBox& a_Nsq,
                                     const Box&       a_valid,
                                     const Real       a_dz)
{
    CH_assert(a_valid.size().product() == a_valid.size(SpaceDim-1));
    // CH_assert(strcmp(LevelGeometry::getCoorMapName(), "Cartesian") != 0);

#   define ARRAY1D(a)   Box(IntVect(D_DECL((1),(1),(1))), IntVect(D_DECL((a),(1),(1)))), (1)
#   define ARRAY2D(a,b) Box(IntVect(D_DECL((1),(1),(1))), IntVect(D_DECL((a),(b),(1)))), (1)

    char dflag[1] = {'S'};
    static const Real abstol = 2.0 * dlamch_(dflag);
    const int N = a_valid.size().product();
    BaseFab<int>  IFAIL(ARRAY1D(1));    // (1)
    BaseFab<int>  IWORK(ARRAY1D(5*N));  // (5*N)
    BaseFab<Real> AB(ARRAY2D(2,N));     // (2,N)
    BaseFab<Real> BB(ARRAY2D(1,N));     // (1,N)
    BaseFab<Real> Q(ARRAY2D(N,N));      // (N,N)
    BaseFab<Real> W(ARRAY1D(N));        // (N)
    BaseFab<Real> WORK(ARRAY1D(7*N));   // (7*N)
    BaseFab<Real> Z(ARRAY2D(N,N));      // (N,N)

    FORT_SOLVEVERTEIGENPROBLEMCARTESIAN (
        CHF_REAL(a_c0),
        CHF_FRA1(a_phi0,0),
        CHF_CONST_FRA1(a_Nsq,0),
        CHF_BOX(a_valid),
        CHF_CONST_INT(N),
        CHF_CONST_REAL(a_dz),
        CHF_CONST_REAL(abstol),
        CHF_FIA1(IFAIL,0),
        CHF_FIA1(IWORK,0),
        CHF_FRA1(AB,0),
        CHF_FRA1(BB,0),
        CHF_FRA1(Q,0),
        CHF_FRA1(W,0),
        CHF_FRA1(WORK,0),
        CHF_FRA1(Z,0));

#   undef ARRAY1D
#   undef ARRAY2D
}


// -----------------------------------------------------------------------------
// Fills a_A0 with the horizontal sech^2 structure.
// -----------------------------------------------------------------------------
void fillHorizontalStructure (FArrayBox&           a_A0,
                              StructurePool&       a_pool,
                              const LevelGeometry& a_levGeo,
                              const Real           a_xcenter,
                              const Real           a_amp0)
{
    // Sanity checks
    CH_assert(SpaceDim == 2); // Streamfunction method is different for other dims
    CH_assert(a_A0.nComp() == 1);
    CH_assert(a_levGeo.getDomain() == a_pool.getLevGeo().getDomain());

    // Gather domain data
    const ProblemDomain& domain = a_levGeo.getDomain();
    const Box domBox = domain.domainBox();
    const Real horizPhysDx = a_levGeo.getDx()[0];

    // Get structure data
    const Real c0    = a_pool.c0   (domBox);
    const Real alpha = a_pool.alpha(domBox);
    const Real beta  = a_pool.beta (domBox);
    const Real gamma = a_pool.gamma(domBox);

    // Create calculation regions
    const Box horizCCBox = horizontalDataBox(domain);
    const Box horizFCBox = surroundingNodes(horizCCBox, 0);

    // Fill A0
    if (a_A0.box().type() == IntVect::Zero) {
        // CC data
        CH_assert(a_A0.contains(horizCCBox));

        FORT_FILLHORIZSTRUCTURE (
            CHF_FRA1(a_A0,0),
            CHF_CONST_REAL(alpha),
            CHF_CONST_REAL(beta),
            CHF_CONST_REAL(gamma),
            CHF_CONST_REAL(c0),
            CHF_CONST_REAL(a_xcenter),
            CHF_CONST_REAL(a_amp0),
            CHF_BOX(horizCCBox),
            CHF_CONST_INTVECT(horizCCBox.type()),
            CHF_CONST_REAL(horizPhysDx));

    } else if (a_A0.box().type() == BASISV(0)) {
        // FC data
        CH_assert(a_A0.contains(horizFCBox));

        FORT_FILLHORIZSTRUCTURE (
            CHF_FRA1(a_A0,0),
            CHF_CONST_REAL(alpha),
            CHF_CONST_REAL(beta),
            CHF_CONST_REAL(gamma),
            CHF_CONST_REAL(c0),
            CHF_CONST_REAL(a_xcenter),
            CHF_CONST_REAL(a_amp0),
            CHF_BOX(horizFCBox),
            CHF_CONST_INTVECT(horizFCBox.type()),
            CHF_CONST_REAL(horizPhysDx));
    }

    static bool writtenOnce = false;
    if (!writtenOnce) {
        pout() << "\nStructure:"
               << "\n\tc0    = " << c0
               << "\n\tV     = " << beta * a_amp0 / (3.0 * alpha)
               << "\n\twidth = " << sqrt(6.0 * gamma * c0 / (beta * a_amp0))
               << "\n\talpha = " << alpha
               << "\n\tbeta  = " << beta
               << "\n\tgamma = " << gamma
               << endl;
        writtenOnce = true;
    }
}



// StructurePool...

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
StructurePool::StructurePool ()
: m_levGeoPtr(NULL),
  m_physBCPtr(NULL)
{
    m_deletePtrs.resize(0);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
StructurePool::~StructurePool ()
{
    for (int idx = 0; idx < m_deletePtrs.size(); ++idx) {
        delete m_deletePtrs[idx].CCNsqPtr;
        m_deletePtrs[idx].CCNsqPtr = NULL;

        delete m_deletePtrs[idx].CCphiPtr;
        m_deletePtrs[idx].CCphiPtr = NULL;
    }
    m_deletePtrs.resize(0);
    m_dataMap.erase(m_dataMap.begin(), m_dataMap.end());

    m_levGeoPtr = NULL;
    m_physBCPtr = NULL;
}


// -----------------------------------------------------------------------------
// Sets pointers to objects needed by createNewLevel. If there is any chance
// that createNewLevel will be triggered, be sure to call this function
// first or you'll be sorry.
// -----------------------------------------------------------------------------
void StructurePool::setGeometry (const LevelGeometry* a_levGeoPtr,
                                 const PhysBCUtil*    a_physBCPtr)
{
    m_levGeoPtr = a_levGeoPtr;
    m_physBCPtr = a_physBCPtr;
}


// -----------------------------------------------------------------------------
// The data accessor for internal use. Checks if a call to createNewLevel
// is needed and performs some sanity checks on the pool in debug mode.
// -----------------------------------------------------------------------------
StructurePool::ElementType StructurePool::getStructureData (const Box& a_box)
{
    if (m_dataMap.find(a_box) == m_dataMap.end()) {
        return this->createNewLevel(a_box);
    }

    CH_assert(m_dataMap.count(a_box) == 1);
    return m_dataMap[a_box];
}


// -----------------------------------------------------------------------------
// The data accessor for internal use. Checks if a call to createNewLevel
// is needed and performs some sanity checks on the pool in debug mode.
// (const version)
// -----------------------------------------------------------------------------
const StructurePool::ElementType StructurePool::getStructureData (const Box& a_box) const
{
    CH_assert(m_dataMap.count(a_box) == 1);
    return m_dataMap.find(a_box)->second;
}


// -----------------------------------------------------------------------------
// This is called if a level's data is not yet in the pool. If there is any
// chance that this function will be triggered, be sure to call setGeometry
// first or you'll be sorry. This is the workhorse of the class - it solves
// the eigenprobem using LAPACK.
// -----------------------------------------------------------------------------
StructurePool::ElementType StructurePool::createNewLevel (const Box& a_box) {
    CH_TIME("StructurePool::createNewLevel");

    // Current limitations
    // CH_assert(std::string(m_levGeoPtr->getCoorMapName()) == std::string("Cartesian"));
    CH_assert(SpaceDim == 2);

    // Sanity checks
    CH_assert(m_dataMap.count(a_box) == 0);
    CH_assert(m_levGeoPtr != NULL);
    CH_assert(m_levGeoPtr->getDomain().domainBox() == a_box);

    // Gather domain data
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    const ProblemDomain& domain = m_levGeoPtr->getDomain();
    const Box domBox = domain.domainBox();
    const int Nx = domBox.size(0);
    const int Nz = domBox.size(1);
    const Real physDx = m_levGeoPtr->getDx()[0];
    const Real physDz = m_levGeoPtr->getDx()[1];

    // Create calculation regions
    const Box vertCCBox = verticalDataBox(domain);
    const Box vertFCBox = surroundingNodes(vertCCBox, 1);

    // Get background bouyancy
    FArrayBox CCB(grow(vertCCBox, 2*BASISV(1)), 1);
    m_physBCPtr->setBackgroundScalar(CCB, 0, *m_levGeoPtr, DataIndex(), 0.0);

    FArrayBox FCB(grow(vertFCBox, 2*BASISV(1)), 1);
    m_physBCPtr->setBackgroundScalar(FCB, 0, *m_levGeoPtr, DataIndex(), 0.0);

    // Calculate Brunt-Vaisala frequency
    FArrayBox CCNsq(grow(vertCCBox, 1*BASISV(1)), 1);
    computeBVFreqCartesian(CCNsq, FCB, CCNsq.box(), *m_levGeoPtr, 2.0);

    FArrayBox FCNsq(grow(vertFCBox, 1*BASISV(1)), 1);
    computeBVFreqCartesian(FCNsq, CCB, FCNsq.box(), *m_levGeoPtr, 2.0);

    // Solve vertical structure eigenproblem (fastest mode only)
    Real c0;
    FArrayBox CCphi0(vertCCBox, 1);
    solveVertEigenProblemCartesian(c0, CCphi0, CCNsq, vertCCBox, physDz);

    // Compute KdV coefficients
    Real alpha, beta, gamma;
    FORT_COMPUTEHORIZCOEFFSCARTESIAN (
        CHF_REAL(alpha),
        CHF_REAL(beta),
        CHF_REAL(gamma),
        CHF_CONST_FRA1(CCphi0,0),
        CHF_CONST_FRA1(CCNsq,0),
        CHF_CONST_FRA1(FCNsq,0),
        CHF_CONST_REAL(c0),
        CHF_BOX(vertCCBox),
        CHF_CONST_REAL(physDz));

    // Create a container for the new data.
    ElementType newLevelStructure;
    newLevelStructure.CCNsqPtr = new FArrayBox(vertCCBox, 1);
    newLevelStructure.CCphiPtr = new FArrayBox(grow(vertCCBox, 1*BASISV(1)), 1);

    // Pack data into the new container.
    newLevelStructure.c0 = c0;
    newLevelStructure.alpha = alpha;
    newLevelStructure.beta = beta;
    newLevelStructure.gamma = gamma;
    newLevelStructure.CCNsqPtr->copy(CCNsq);
    newLevelStructure.CCphiPtr->copy(CCphi0);

    // Set homogeneous BCs on phi
    {
        FArrayBox& phiRef = *newLevelStructure.CCphiPtr;
        IntVect bdryIV, ghostIV;

        bdryIV = vertCCBox.smallEnd();
        ghostIV = bdryIV - BASISV(1);
        phiRef(ghostIV) = -phiRef(bdryIV);

        bdryIV = vertCCBox.bigEnd();
        ghostIV = bdryIV + BASISV(1);
        phiRef(ghostIV) = -phiRef(bdryIV);
    }

    // Register the new container full of data with the map.
    m_dataMap[a_box] = newLevelStructure;
    m_deletePtrs.push_back(newLevelStructure);
    return newLevelStructure;
}
