#include "TaylorGreenBCUtil.H"
#include "TaylorGreenBCUtilF_F.H"
#include "ProblemContext.H"
#include "EllipticBCUtils.H"
#include "EllipticBCUtilsF_F.H"
#include "Constants.H"
#include "BoxIterator.H"
#include "LevelGeometryF_F.H"

bool TaylorGreenBCUtil::s_staticParamsSet = false;
int TaylorGreenBCUtil::s_udir;
int TaylorGreenBCUtil::s_vdir;
Real TaylorGreenBCUtil::s_kx;
Real TaylorGreenBCUtil::s_ky;
RealVect TaylorGreenBCUtil::s_L;
Real TaylorGreenBCUtil::s_nu;


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
TaylorGreenBCUtil::TaylorGreenBCUtil ()
{;}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
TaylorGreenBCUtil::~TaylorGreenBCUtil ()
{;}


// -----------------------------------------------------------------------------
// Define constructor
// -----------------------------------------------------------------------------
void TaylorGreenBCUtil::define ()
{
    if (s_staticParamsSet) return;

    const ProblemContext* ctx = ProblemContext::getInstance();

    s_udir = 0;
    s_vdir = 1;
    s_L = ctx->domainLength;
    s_kx = 2.0 * Pi / s_L[s_udir];
    s_ky = 2.0 * Pi / s_L[s_vdir];
    s_nu = ctx->nu;

    PhysBCUtil::define();
    s_staticParamsSet = true;
}


// -----------------------------------------------------------------------------
// Calculates F(t)
// -----------------------------------------------------------------------------
Real TaylorGreenBCUtil::FofT (const Real           a_time,
                              const LevelGeometry& a_levGeo) const
{
    const RealVect& dx = a_levGeo.getDx();
    RealVect L = RealVect(a_levGeo.getDomain().size()) * dx;
    const RealVect kk = 2.0 * Pi / L;
    return exp(-(kk[s_udir]*kk[s_udir] + kk[s_vdir]*kk[s_vdir]) * s_nu * a_time);
}


// -----------------------------------------------------------------------------
// Fills a FAB with the velocity solution
// Results are in a Cartesian basis at mapped locations.
// -----------------------------------------------------------------------------
void TaylorGreenBCUtil::fillVelSoln (FArrayBox&           a_velFAB,
                                     const LevelGeometry& a_levGeo,
                                     const DataIndex&     a_di,
                                     const Real           a_time) const
{
    CH_assert(a_velFAB.nComp() == SpaceDim);

    // Get the Cartesian coordinates of each grid cell.
    const Box& dataBox = a_velFAB.box();
    FArrayBox cartPos(dataBox, SpaceDim);
    a_levGeo.fill_physCoor(cartPos);

    // Fill a FAB with the velocity solution
    const Real F = this->FofT(a_time, a_levGeo);
    for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
        int dirIndicator = -1;
        if (dir == s_udir) dirIndicator = 0;
        if (dir == s_vdir) dirIndicator = 1;

        FORT_FILL_CARTVELSOLN(CHF_FRA1(a_velFAB,dir),
                              CHF_CONST_INT(dirIndicator),
                              CHF_CONST_FRA1(cartPos,s_udir),
                              CHF_CONST_FRA1(cartPos,s_vdir),
                              CHF_CONST_REAL(s_kx),
                              CHF_CONST_REAL(s_ky),
                              CHF_CONST_REAL(F),
                              CHF_CONST_REAL(a_time),
                              CHF_BOX(dataBox));
    }
}


// -----------------------------------------------------------------------------
// Fills a FluxBox with the velocity solution
// Results are in a Cartesian basis at mapped locations.
// -----------------------------------------------------------------------------
void TaylorGreenBCUtil::fillVelSoln (FluxBox&             a_velFB,
                                     const LevelGeometry& a_levGeo,
                                     const DataIndex&     a_di,
                                     const Real           a_time) const
{
    CH_assert(a_velFB.nComp() == 1 || a_velFB.nComp() == SpaceDim);

    for (int FCdir = 0; FCdir < SpaceDim; ++FCdir) {
        // Get the Cartesian coordinates of each grid cell.
        const Box& dataBox = a_velFB[FCdir].box();
        FArrayBox cartPos(dataBox, SpaceDim);
        a_levGeo.fill_physCoor(cartPos);

        // Fill a FAB with the velocity solution
        const Real F = this->FofT(a_time, a_levGeo);
        if (a_velFB.nComp() == 1) {
            int dirIndicator = -1;
            if (FCdir == s_udir) dirIndicator = 0;
            if (FCdir == s_vdir) dirIndicator = 1;

            FORT_FILL_CARTVELSOLN(CHF_FRA1(a_velFB[FCdir],0),
                                  CHF_CONST_INT(dirIndicator),
                                  CHF_CONST_FRA1(cartPos,s_udir),
                                  CHF_CONST_FRA1(cartPos,s_vdir),
                                  CHF_CONST_REAL(s_kx),
                                  CHF_CONST_REAL(s_ky),
                                  CHF_CONST_REAL(F),
                                  CHF_CONST_REAL(a_time),
                                  CHF_BOX(dataBox));

        } else {
            for (int comp = 0; comp < a_velFB.nComp(); ++comp) {
                int dirIndicator = -1;
                if (comp == s_udir) dirIndicator = 0;
                if (comp == s_vdir) dirIndicator = 1;

                FORT_FILL_CARTVELSOLN(CHF_FRA1(a_velFB[FCdir],comp),
                                      CHF_CONST_INT(dirIndicator),
                                      CHF_CONST_FRA1(cartPos,s_udir),
                                      CHF_CONST_FRA1(cartPos,s_vdir),
                                      CHF_CONST_REAL(s_kx),
                                      CHF_CONST_REAL(s_ky),
                                      CHF_CONST_REAL(F),
                                      CHF_CONST_REAL(a_time),
                                      CHF_BOX(dataBox));
            }
        }
    }
}


// -----------------------------------------------------------------------------
// Fills a FAB with the pressure solution
// -----------------------------------------------------------------------------
void TaylorGreenBCUtil::fillPressureSoln (FArrayBox&           a_pressureFAB,
                                          const LevelGeometry& a_levGeo,
                                          const DataIndex&     a_di,
                                          const Real           a_time) const
{
    CH_assert(a_pressureFAB.nComp() == 1);

    // Get the Cartesian coordinates of each grid cell.
    const Box& dataBox = a_pressureFAB.box();
    FArrayBox cartPos(dataBox, SpaceDim);
    a_levGeo.fill_physCoor(cartPos);

    // Fill a FAB with the pressure solution
    const Real F = this->FofT(a_time, a_levGeo);
    FORT_FILL_PRESSURESOLN(CHF_FRA1(a_pressureFAB,0),
                           CHF_CONST_FRA1(cartPos,s_udir),
                           CHF_CONST_FRA1(cartPos,s_vdir),
                           CHF_CONST_REAL(s_kx),
                           CHF_CONST_REAL(s_ky),
                           CHF_CONST_REAL(F),
                           CHF_CONST_REAL(a_time),
                           CHF_BOX(dataBox));
}


// -----------------------------------------------------------------------------
// Computes the velocity error in an FArrayBox.
// -----------------------------------------------------------------------------
void TaylorGreenBCUtil::computeVelError (FArrayBox&           a_errFAB,
                                         const FArrayBox&     a_velFAB,
                                         const LevelGeometry& a_levGeo,
                                         const DataIndex&     a_di,
                                         const Real           a_time,
                                         const bool           a_mapped,
                                         const bool           a_multByJ) const
{
    const Box& valid = a_levGeo.getBoxes()[a_di];
    const int ncomp = SpaceDim;

    // Sanity checks
    CH_assert(a_errFAB.nComp() == ncomp);
    CH_assert(a_velFAB.nComp() == ncomp);
    CH_assert(a_errFAB.box().contains(valid));
    CH_assert(a_velFAB.box().contains(valid));
    CH_assert(a_mapped || !a_multByJ);

    // Get solution data
    const Real scale = -1.0;
    this->fillVelSoln(a_errFAB, a_levGeo, a_di, a_time);

    // Convert solution to match vel.
    if (a_mapped) {
        FArrayBox dXidxFAB(valid, SpaceDim*SpaceDim);
        a_levGeo.fill_dXidx(dXidxFAB);
        FORT_CONTRACTMATRIXVECTORCC(CHF_FRA(a_errFAB),
                                    CHF_CONST_FRA(dXidxFAB),
                                    CHF_BOX(valid));

        if (a_multByJ) a_levGeo.multByJ(a_errFAB, a_di);
    }

    // Compute error
    a_errFAB.minus(a_velFAB, valid, 0, 0, ncomp);
    a_errFAB.negate();
}


// ************************ ICs / background fields ****************************

// -----------------------------------------------------------------------------
// Fills a FAB with the initial velocity
// -----------------------------------------------------------------------------
void TaylorGreenBCUtil::setVelIC (FArrayBox&           a_velFAB,
                                  const int            a_velComp,
                                  const LevelGeometry& a_levGeo,
                                  const DataIndex&     a_di) const
{
    CH_assert(a_velFAB.nComp() == SpaceDim);
    CH_assert(0 <= a_velComp);
    CH_assert(a_velComp < SpaceDim);

    // For now, just set all components.
    this->fillVelSoln(a_velFAB, a_levGeo, a_di, 0.0);
}


// -----------------------------------------------------------------------------
// Fills a FAB with the initial scalars
// -----------------------------------------------------------------------------
void TaylorGreenBCUtil::setScalarIC (FArrayBox&           a_scalarFAB,
                                     const int            a_scalarComp,
                                     const LevelGeometry& a_levGeo,
                                     const DataIndex&     a_di) const
{
    CH_assert(a_scalarFAB.nComp() == 1);

    const Box& dataBox = a_scalarFAB.box();
    BoxIterator bit(dataBox);

    int gradDir = a_scalarComp % SpaceDim;
    Real l = s_L[gradDir];
    Real offset = (1.0 - dataBox.type()[gradDir]) * 0.5;

    for (bit.reset(); bit.ok(); ++bit) {
        const IntVect& cc = bit();
        a_scalarFAB(cc) = (Real(cc[gradDir]) + offset) / l;
    }
}


// -----------------------------------------------------------------------------
// basicVelFuncBC
// Sets physical BCs on velocities.
// -----------------------------------------------------------------------------
BCMethodHolder TaylorGreenBCUtil::basicVelFuncBC (int a_veldir, bool a_isViscous, Real a_scale) const
{
    RefCountedPtr<BCGhostClass> BCPtr;
    BCPtr = RefCountedPtr<BCGhostClass>(
        new TaylorGreenBCGhostClass(a_veldir, s_nu, a_scale)
    );

    BCMethodHolder holder;
    holder.addBCMethod(BCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// TaylorGreenBCGhostClass
// u = sin(k*x) * cos(l*y) * F(t)
// v = -cos(k*x) * sin(l*y) * F(t)
// F(t) = exp(2*nu*t)
// -----------------------------------------------------------------------------
void TaylorGreenBCGhostClass::operator() (FArrayBox&           a_state,
                                          const FArrayBox*     a_extrapPtr,   // Just a dummy
                                          const Box&           a_valid,
                                          const ProblemDomain& a_domain,
                                          const RealVect&      a_dx,
                                          const DataIndex&     a_index,       // Just a dummy
                                          const FluxBox*       a_JgupPtr,     // Just a dummy
                                          bool                 a_homogeneous,
                                          Real                 a_time,
                                          const Interval&      a_interval) const
{
    CH_TIME("TaylorGreenBCGhostClass::operator()");

    // Sanity checks
    const IntVect boxType = a_state.box().type();
    CH_assert(boxType == IntVect::Zero || boxType.sum() == 1);
    CH_assert(boxType == a_valid.type());
    CH_assert(a_JgupPtr != NULL);

    // Alias the state FAB
    const Interval interv = (a_interval == Interval())? a_state.interval(): a_interval;
    const int scomp = interv.begin();
    const int ncomp = interv.size();
    CH_assert(a_state.nComp() > interv.end());
    FArrayBox stateAlias(interv, a_state);

    // Create a domain box with the same centering as the state FAB
    Box thisDom = a_domain.domainBox();
    thisDom.convert(boxType);

    // Calculate k
    const Real twoPi = 2.0 * 3.14159265358979;
    RealVect k(a_domain.size());
    k *= a_dx;
    k = twoPi / k;

    const Real ksq = k[0]*k[0] + k[1]*k[1];

    const Real FofT = exp(-ksq * s_nu * a_time) * m_scale;

    for (int idir = 0; idir < SpaceDim; ++idir) {
        if (a_domain.isPeriodic(idir)) continue;

        const FArrayBox& Jgupi = (*a_JgupPtr)[idir];

        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide iside = sit();
            const int isign = sign(iside);

            // Find the bdry locations
            Box bigDom = thisDom;
            bigDom.grow(4*(IntVect::Unit-BASISV(idir)));    // TODO: This should be done in all BCs!
            Box faceBox = surroundingNodes(a_valid, idir) //bdryBox(a_valid, idir, iside, 1)
                        & bdryBox(bigDom, idir, iside, 1);
            if (faceBox.isEmpty()) continue;

            if (boxType[idir] == 1) {
                // FC data...
                CH_assert(stateAlias.contains(faceBox));
                CH_assert(boxType.sum() == 1);

                if (a_homogeneous) {
                    stateAlias.setVal(0.0, faceBox, 0, stateAlias.nComp());
                } else {
                    FORT_TAYLORGREENFCBC (CHF_FRA(stateAlias),
                                          CHF_BOX(faceBox),
                                          CHF_CONST_INT(idir),
                                          CHF_CONST_INT(isign),
                                          CHF_CONST_REALVECT(a_dx),
                                          CHF_CONST_REALVECT(k),
                                          CHF_CONST_INT(m_velDir),
                                          CHF_CONST_REAL(FofT));
                }
            } else {
                // Find the source data locations
                Box srcBox = faceBox;
                srcBox.shiftHalf(idir, -isign);
                srcBox &= a_state.box();
                CH_assert(!srcBox.isEmpty());

                // Find the ghost locations
                Box ghostBox = faceBox;
                ghostBox.shiftHalf(idir, isign);
                ghostBox &= a_state.box();
                if (ghostBox.isEmpty()) continue;

                CH_assert(stateAlias.contains(srcBox));
                CH_assert(stateAlias.contains(ghostBox));

                // Set the homogeneous ghosts
                a_state.copy(a_state, srcBox, scomp, ghostBox, scomp, ncomp);
                a_state.negate(ghostBox, scomp, ncomp);

                // Increment ghost values to include inhomog contribution
                if (!a_homogeneous) {
                    FORT_TAYLORGREENCCBC (CHF_FRA(stateAlias),
                                          CHF_BOX(ghostBox),
                                          CHF_CONST_INTVECT(boxType),
                                          CHF_CONST_INT(idir),
                                          CHF_CONST_INT(isign),
                                          CHF_CONST_REALVECT(a_dx),
                                          CHF_CONST_REALVECT(k),
                                          CHF_CONST_INT(m_velDir),
                                          CHF_CONST_REAL(FofT));
                }
            }
        } // end loop over sides
    } // end loop over dirs
}
