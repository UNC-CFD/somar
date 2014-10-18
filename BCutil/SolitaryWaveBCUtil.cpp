#include "SolitaryWaveBCUtil.H"
#include "SolitaryWaveBCUtilF_F.H"
#include "StratUtils.H"
#include "ProblemContext.H"
#include "EllipticBCUtils.H"
#include "BoxIterator.H"


// -----------------------------------------------------------------------------
// Static member variable definitions
// -----------------------------------------------------------------------------
RealVect SolitaryWaveBCUtil::s_L;                        // Domain extents
Real SolitaryWaveBCUtil::s_H2;                           // Lower layer depth (undisturbed elevation of the isopycnal)
Real SolitaryWaveBCUtil::s_delta0;                       // Pycnocline thickness
Real SolitaryWaveBCUtil::s_amp0;                         // Amplitude of initial depression
Real SolitaryWaveBCUtil::s_xcenter;                      // Center of solitary wave in physical coordinates

// These are scaled as a reduced gravity. The choice of s_rho0 should be of no consequence.
Real SolitaryWaveBCUtil::s_rho0;                         // Average density
Real SolitaryWaveBCUtil::s_drho;         // Total bouyancy variation

// Vertical/horizontal structure
StructurePool SolitaryWaveBCUtil::s_structure;


// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
SolitaryWaveBCUtil::SolitaryWaveBCUtil ()
{
    // StratifiedFlowUtils::staticDefine();

    static bool paramsRead = false;
    if (!paramsRead) {
        const ProblemContext* ctx = ProblemContext::getInstance();

        s_L = ctx->domainLength;

        s_xcenter = ctx->solitaryWave_xcenter;
        s_H2      = ctx->solitaryWave_H2;
        s_amp0    = ctx->solitaryWave_amp;
        s_rho0    = ctx->solitaryWave_rho0;
        s_drho    = ctx->solitaryWave_drho;
        s_delta0  = ctx->solitaryWave_delta0;

        paramsRead = true;
    }
}


// -----------------------------------------------------------------------------
// Default destructor
// -----------------------------------------------------------------------------
SolitaryWaveBCUtil::~SolitaryWaveBCUtil ()
{;}


// -----------------------------------------------------------------------------
// This object is its own factory
// -----------------------------------------------------------------------------
PhysBCUtil* SolitaryWaveBCUtil::newPhysBCUtil () const
{
    PhysBCUtil* newBCPtr = new SolitaryWaveBCUtil();
    return newBCPtr;
}


// -----------------------------------------------------------------------------
// Fills a FAB with the initial velocity
// Locations are in mapped space, but components are Cartesian.
// -----------------------------------------------------------------------------
void SolitaryWaveBCUtil::setVelIC (FArrayBox&           a_velFAB,
                                   const int            a_velComp,
                                   const LevelGeometry& a_levGeo,
                                   const DataIndex&     a_di) const
{
    // Sanity checks
    CH_assert(SpaceDim == 2); // Streamfunction method is different for other dims
    CH_assert(a_velFAB.nComp() == SpaceDim);
    CH_assert(0 <= a_velComp);
    CH_assert(a_velComp < SpaceDim);

    // Gather domain data
    const ProblemDomain& domain = a_levGeo.getDomain();
    const Box domBox = domain.domainBox();
    const Box valid = a_levGeo.getBoxes()[a_di] & a_velFAB.box();
    const RealVect physDx = a_levGeo.getDx();

    // Get structure data
    s_structure.setGeometry(&a_levGeo, this);
    const FArrayBox& CCphi = s_structure.phi(domBox);

    // Create calculation regions
    const Box horizCCBox = horizontalDataBox(domain);
    const Box horizFCBox = surroundingNodes(horizCCBox, 0);
    const Box vertCCBox = verticalDataBox(domain);

    // Fill A
    FArrayBox CCA(horizCCBox, 1);
    fillHorizontalStructure(CCA, s_structure, a_levGeo, s_xcenter, s_amp0);

    FArrayBox FCA(horizFCBox, 1);
    fillHorizontalStructure(FCA, s_structure, a_levGeo, s_xcenter, s_amp0);

    // Calculate velocity
    FORT_SOLITARYWAVE_SETVELIC (
        CHF_FRA(a_velFAB),
        CHF_CONST_INT(a_velComp),
        CHF_CONST_FRA1(CCA,0),
        CHF_CONST_FRA1(FCA,0),
        CHF_CONST_FRA1(CCphi,0),
        CHF_BOX(valid),
        CHF_BOX(horizCCBox),
        CHF_BOX(vertCCBox),
        CHF_CONST_REALVECT(physDx));
}


// -----------------------------------------------------------------------------
// Fills a FAB with the initial scalars
// -----------------------------------------------------------------------------
void SolitaryWaveBCUtil::setScalarIC (FArrayBox&           a_scalarFAB,
                                      const int            a_scalarComp,
                                      const LevelGeometry& a_levGeo,
                                      const DataIndex&     a_di) const
{
    CH_assert(SpaceDim == 2); // Streamfunction method is different for other dims
    CH_assert(a_scalarFAB.nComp() == 1);

    if (a_scalarComp == 0) {
        // Gather domain data
        const ProblemDomain& domain = a_levGeo.getDomain();
        const Box domBox = domain.domainBox();
        const Box valid = a_levGeo.getBoxes()[a_di] & a_scalarFAB.box();

        // Get structure data
        s_structure.setGeometry(&a_levGeo, this);
        const Real c0 = s_structure.c0(domBox);
        const FArrayBox& CCNsq = s_structure.Nsq(domBox);
        const FArrayBox& CCphi = s_structure.phi(domBox);

        // Create calculation regions
        const Box horizCCBox = horizontalDataBox(domain);
        const Box vertCCBox = verticalDataBox(domain);

        // Fill A
        FArrayBox CCA(horizCCBox, 1);
        fillHorizontalStructure(CCA, s_structure, a_levGeo, s_xcenter, s_amp0);

        // Fill background scalar
        FArrayBox CCBbg(valid, 1);
        this->setBackgroundScalar(CCBbg, 0, a_levGeo, a_di, 0.0);

        // Calculate total initial scalar
        FORT_SOLITARYWAVE_SETSCALARIC (
            CHF_FRA1(a_scalarFAB,0),
            CHF_CONST_FRA1(CCBbg,0),
            CHF_CONST_FRA1(CCA,0),
            CHF_CONST_FRA1(CCphi,0),
            CHF_CONST_FRA1(CCNsq,0),
            CHF_CONST_REAL(c0),
            CHF_BOX(valid),
            CHF_BOX(horizCCBox),
            CHF_BOX(vertCCBox));
    } else {
        MayDay::Error("scalar IC not defined for comp > 0");
    }
}


// -----------------------------------------------------------------------------
// Fills a FAB with the background scalar
// -----------------------------------------------------------------------------
void SolitaryWaveBCUtil::setBackgroundScalar (FArrayBox&           a_scalarFAB,
                                              const int            a_scalarComp,
                                              const LevelGeometry& a_levGeo,
                                              const DataIndex&     a_di,
                                              const Real           a_time) const
{
    CH_assert(a_scalarFAB.nComp() == 1);

    if (s_useBackgroundScalar && a_scalarComp == 0) {
        FArrayBox posFAB(a_scalarFAB.box(), 1);
        const RealVect& dx = a_levGeo.getDx();
        a_levGeo.getGeoSourcePtr()->fill_physCoor(posFAB, 0, SpaceDim-1, dx);

        BoxIterator bit(a_scalarFAB.box());
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& iv = bit();

            Real pos[CH_SPACEDIM];
            pos[CH_SPACEDIM-1] = posFAB(iv);

            int dirDummy;
            Side::LoHiSide sideDummy;
            Real derivScaleDummy = 1.0e300;
            Real value[1];
            this->bscalBCValues(pos, &dirDummy, &sideDummy, value, derivScaleDummy, a_time);
            a_scalarFAB(iv,a_scalarComp) = value[0];
        }

    } else {
        a_scalarFAB.setVal(0.0);
    }
}


// -----------------------------------------------------------------------------
// Sets the (vector scaled) target velocity for the sponge layer. By default,
// this function persuades the velocity field to approach its inflow value.
// -----------------------------------------------------------------------------
void SolitaryWaveBCUtil::fillVelSpongeLayerTarget (FArrayBox&           a_target,
                                                   const int            a_velComp,
                                                   const int            a_spongeDir,
                                                   const Side::LoHiSide a_spongeSide,
                                                   const LevelGeometry& a_levGeo,
                                                   const DataIndex&     a_di,
                                                   const Real           a_time)
{
    // Sanity checks
    CH_assert(useSpongeLayer());
    CH_assert(a_target.nComp() == 1);
    CH_assert(0 <= a_spongeDir);
    CH_assert(a_spongeDir < SpaceDim);

    // This assumes not much is happening in the far-field.
    if (a_spongeDir == 0) {
        a_target.setVal(0.0);
    } else {
        MayDay::Error("BeamGenerationBCUtil::fillVelSpongeLayerTarget "
                      "can only set a sponge target when a_spongeDir = 0");
    }
}


// -----------------------------------------------------------------------------
// basicVelFuncBC
// Sets physical BCs on velocities.
// -----------------------------------------------------------------------------
BCMethodHolder SolitaryWaveBCUtil::basicVelFuncBC (int a_veldir, bool a_isViscous) const
{
    const IntVect hUnit = IntVect::Unit - BASISV(CH_SPACEDIM-1);
    const IntVect vUnit = BASISV(CH_SPACEDIM-1);

    BCMethodHolder holder;

    //             Freeslip
    // u: Neum 0 |==========| Neum 0
    //             Freeslip

    // Low order extrap in horizontal (sponged) directions
    int extrapOrder = 0;
    RefCountedPtr<BCGhostClass> horizBCPtr(
        new EllipticExtrapBCGhostClass(extrapOrder,
                                       hUnit,
                                       hUnit)
    );
    holder.addBCMethod(horizBCPtr);

    RefCountedPtr<BCFluxClass> fluxBCPtr(
        new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                         RealVect::Zero,
                                         BASISV(0),
                                         BASISV(0))
    );
    holder.addBCMethod(fluxBCPtr);

    // Free slip in vertical dir
    RefCountedPtr<BCGhostClass> hiVertBCPtr = RefCountedPtr<BCGhostClass>(
        new BasicVelocityBCGhostClass(0.0,             // inflowVel
                                      -1,              // inflowDir
                                      Side::Lo,        // inflowSide
                                      -1,              // outflowDir
                                      Side::Hi,        // outflowSide
                                      a_veldir,
                                      false,           // isViscous
                                      vUnit,
                                      vUnit)
    );
    holder.addBCMethod(hiVertBCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// basicScalarFuncBC   (Extrapolate BCs)
// Sets physical BCs on a generic passive scalar.
// Chombo uses 1st order extrap
// -----------------------------------------------------------------------------
BCMethodHolder SolitaryWaveBCUtil::basicScalarFuncBC () const
{
    BCMethodHolder holder;

    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticConstDiriBCGhostClass(RealVect::Zero,
                                          RealVect::Zero,
                                          IntVect::Unit,
                                          IntVect::Unit)
    );
    holder.addBCMethod(BCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// Simply sets a_value to the background density at any given point, which
// is all we need for the boundary conditions.
// This function conforms to the EllipticBCValueFunc typedef.
// -----------------------------------------------------------------------------
void SolitaryWaveBCUtil::bscalBCValues (Real*           a_pos,
                                        int*            a_dir,
                                        Side::LoHiSide* a_side,
                                        Real*           a_value,
                                        Real            a_derivScale,
                                        Real            a_time)
{
    Real arg = (a_pos[CH_SPACEDIM-1] - s_H2) / s_delta0;
    a_value[0] = s_rho0 - s_drho*tanh(arg);
}

