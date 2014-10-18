#include "HorizConvBCUtil.H"
#include "EllipticBCUtils.H"
#include "Constants.H"
#include "LevelGeometry.H"


// -----------------------------------------------------------------------------
// Holds Diri BCs for the scalar field at the top of the domain.
// -----------------------------------------------------------------------------
void topBCValueFunc (Real*           a_pos,
                     int*            a_dir,
                     Side::LoHiSide* a_side,
                     Real*           a_value,
                     const RealVect& a_dx,
                     Real            a_time)
{
    // We had better be at the top!
    CH_assert(a_dir[0] == SpaceDim-1);
    CH_assert(a_side[0] == Side::Hi);

    // These params adjust the shape of the hot/cold interface.
    static Real transitionWidth = 1.0e-2;
    static RealVect L = LevelGeometry::getDomainLength();

    static Real xmid = L[0]/2;
    static Real b0 = -0.5;
    static Real deltab = 1.0;

#if CH_SPACEDIM == 2
    Real arg = (a_pos[0] - xmid) / transitionWidth;

#else
    static Real pertAmp = 0.025;
    static Real pertK = 2.0 * (2.0*Pi/L[1]);

    Real xoffset = xmid + pertAmp*cos(pertK * a_pos[1]);
    Real arg = (a_pos[0] - xoffset) / transitionWidth;
#endif

    // Set the BC value
    // a_value[0] = b0 + deltab*tanh(arg);
    a_value[0] = (arg > 0)? 1.0: 0.0;
}


// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
HorizConvBCUtil::HorizConvBCUtil ()
{
    // Nothing is needed here. Setup is done in PhysBCUtil constructor.
}


// -----------------------------------------------------------------------------
// Default destructor
// -----------------------------------------------------------------------------
HorizConvBCUtil::~HorizConvBCUtil ()
{;}


// -----------------------------------------------------------------------------
// This object is its own factory
// -----------------------------------------------------------------------------
PhysBCUtil* HorizConvBCUtil::newPhysBCUtil () const
{

    PhysBCUtil* newBCPtr = new HorizConvBCUtil();
    return newBCPtr;
}


// -----------------------------------------------------------------------------
// Fills a FAB with the initial scalars
// ----------------------------------------------------------------------------- 
void HorizConvBCUtil::setScalarIC (FArrayBox&           a_scalarFAB,
                                   const int            a_scalarComp,
                                   const LevelGeometry& a_levGeo,
                                   const DataIndex&     a_di) const
{
    a_scalarFAB.setVal(0.5);
}


// -----------------------------------------------------------------------------
// diffusiveSolveFuncBC (used in scalar TGA solves)
// This is where the adiabatic/hot/cold BCs are set.
// -----------------------------------------------------------------------------
BCMethodHolder HorizConvBCUtil::diffusiveSolveFuncBC () const
{
    BCMethodHolder holder;

    const IntVect vmask = BASISV(SpaceDim-1);
    const IntVect hmask = IntVect::Unit - vmask;

    // Adiabatic ghost BC
    RefCountedPtr<BCGhostClass> adGhostBCPtr (
        new EllipticConstNeumBCGhostClass(RealVect::Zero,
                                          RealVect::Zero,
                                          IntVect::Unit,
                                          hmask)
    );
    holder.addBCMethod(adGhostBCPtr);

    // Adiabatic flux BC
    RefCountedPtr<BCFluxClass> adFluxBCPtr (
        new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                         RealVect::Zero,
                                         IntVect::Unit,
                                         hmask)
    );
    holder.addBCMethod(adFluxBCPtr);

    // Top BC
    RefCountedPtr<BCGhostClass> topGhostBCPtr (
        new EllipticDiriBCGhostClass(&topBCValueFunc,
                                     IntVect::Zero,
                                     vmask)
    );
    holder.addBCMethod(topGhostBCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// viscousSolveFuncBC (Used in single-component velocity TGA solves)
// -----------------------------------------------------------------------------
BCMethodHolder HorizConvBCUtil::viscousSolveFuncBC (int a_dir) const
{
    BCMethodHolder holder;

    const IntVect normVect = BASISV(a_dir);
    const IntVect transVect = IntVect::Unit - normVect;

    // Transverse, no-slip BCs
    RefCountedPtr<BCGhostClass> transGhostBCPtr (
        new EllipticConstNeumBCGhostClass(RealVect::Zero,
                                          RealVect::Zero,
                                          transVect,
                                          transVect)
    );
    holder.addBCMethod(transGhostBCPtr);

    // Transverse, no-slip BCs (sets fluxes)
    RefCountedPtr<BCFluxClass> transFluxBCPtr (
        new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                         RealVect::Zero,
                                         transVect,
                                         transVect)
    );
    holder.addBCMethod(transFluxBCPtr);

    // Normal, no flux BCs
    RefCountedPtr<BCGhostClass> normGhostBCPtr (
        new EllipticConstDiriBCGhostClass(RealVect::Zero,
                                          RealVect::Zero,
                                          normVect,
                                          normVect)
    );
    holder.addBCMethod(normGhostBCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// diffusiveSourceFuncBC (Used to calculate the diffusive term nu.L[scalar])
// Just use the diffusive solver's BCs.
// -----------------------------------------------------------------------------
BCMethodHolder HorizConvBCUtil::diffusiveSourceFuncBC () const
{
    return HorizConvBCUtil::diffusiveSolveFuncBC();
}


// -----------------------------------------------------------------------------
// basicVelFuncBC
// Sets physical BCs on velocities.
// -----------------------------------------------------------------------------
BCMethodHolder HorizConvBCUtil::basicVelFuncBC (int a_veldir, bool a_isViscous) const
{
    return PhysBCUtil::basicVelFuncBC(a_veldir, false);
}

