#include "CartesianMap.H"
#include "CartesianMapF_F.H"
#include "UsingNamespace.H"



// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
CartesianMap::~CartesianMap ()
{;}


// -----------------------------------------------------------------------------
// 1. Must return the name of the coordinate mapping
// -----------------------------------------------------------------------------
const char* CartesianMap::getCoorMapName () const
{
    return "Cartesian";
}


// -----------------------------------------------------------------------------
// 2. Must return whether or not this metric is diagonal
// -----------------------------------------------------------------------------
bool CartesianMap::isDiagonal () const
{
    return true;
}


// -----------------------------------------------------------------------------
// 3. Must return whether or not this metric is uniform
// -----------------------------------------------------------------------------
bool CartesianMap::isUniform () const
{
    return true;
}


// -----------------------------------------------------------------------------
// 4. Must fill a mapped box with Cartesian locations.
// Typically, only a_dXi[a_mu] will be used, I've included all dXi comps
// just in case we ever need them for some unforseen reason.
// -----------------------------------------------------------------------------
void CartesianMap::fill_physCoor (FArrayBox&      a_dest,
                                  const int       a_destComp,
                                  const int       a_mu,
                                  const RealVect& a_dXi) const
{

    // Sanity checks
    CH_assert(0 <= a_destComp);
    CH_assert(a_destComp < a_dest.nComp());
    CH_assert(0 <= a_mu);
    CH_assert(a_mu < SpaceDim);

    const Real& dXimu = a_dXi[a_mu];
    const Box& destBox = a_dest.box();
    const int destBoxType = destBox.type()[a_mu];

    FORT_CARTESIAN_FILL_PHYSCOOR(
        CHF_FRA1(a_dest,a_destComp),
        CHF_CONST_INT(a_mu),
        CHF_CONST_REAL(dXimu),
        CHF_BOX(destBox),
        CHF_CONST_INT(destBoxType));
}


// -----------------------------------------------------------------------------
// Fills a mapped box with Cartesian locations (a_dest must have SpaceDim comps)
// Mild time saver.
// -----------------------------------------------------------------------------
void CartesianMap::fill_physCoor (FArrayBox&      a_dest,
                                  const RealVect& a_dXi,
                                  const RealVect  a_scale) const
{
    CH_TIME("CartesianMap::fill_physCoor (all comps)");

    const RealVect scaledDXi = a_scale * a_dXi;
    const Box& destBox = a_dest.box();
    const IntVect& destBoxType = destBox.type();

    CH_assert(a_dest.nComp() == SpaceDim);
    CH_assert(destBoxType.product() == 0 || destBoxType.sum() == 1);

    FORT_CARTESIAN_FILL_PHYSCOOR_ALL_COMPS(
        CHF_FRA(a_dest),
        CHF_CONST_REALVECT(scaledDXi),
        CHF_BOX(destBox),
        CHF_CONST_INTVECT(destBoxType));
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the Jacobian matrix elements d[x^mu] / d[xi^nu].
// Major speedup!!!
// -----------------------------------------------------------------------------
void CartesianMap::fill_dxdXi (FArrayBox&      a_dest,
                               const int       a_destComp,
                               const int       a_mu,
                               const int       a_nu,
                               const RealVect& a_dXi,
                               const Real      a_scale) const
{
    CH_TIME("CartesianMap::fill_dxdXi");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));
    CH_assert((0 <= a_nu) && (a_nu < SpaceDim));

    if (a_mu == a_nu) {
        a_dest.setVal(a_scale, a_destComp);
    } else {
        a_dest.setVal(0.0, a_destComp);
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with J = det[Jacobian]
// Major speedup!
// Without the analytic version of this function, machine-epsilon sized
// noise is introduced.
// -----------------------------------------------------------------------------
void CartesianMap::fill_J (FArrayBox&      a_dest,
                           const int       a_destComp,
                           const RealVect& a_dXi,
                           const Real      a_scale) const
{
    CH_TIME("CartesianMap::fill_J");

    // Sanity check
    CH_assert(a_dest.interval().contains(a_destComp));

    a_dest.setVal(a_scale, a_destComp);
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the inverse Jacobian matrix elements d[xi^mu] / d[x^nu].
// -----------------------------------------------------------------------------
void CartesianMap::fill_dXidx (FArrayBox&      a_dest,
                               const int       a_destComp,
                               const int       a_mu,
                               const int       a_nu,
                               const RealVect& a_dXi,
                               const Real      a_scale) const
{
    CH_TIME("CartesianMap::fill_dXidx");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));
    CH_assert((0 <= a_nu) && (a_nu < SpaceDim));

    if (a_mu == a_nu) {
        a_dest.setVal(a_scale, a_destComp);
    } else {
        a_dest.setVal(0.0, a_destComp);
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with 1/J
// -----------------------------------------------------------------------------
void CartesianMap::fill_Jinv (FArrayBox&      a_dest,
                              const int       a_destComp,
                              const RealVect& a_dXi,
                              const Real      a_scale) const
{
    CH_TIME("CartesianMap::fill_Jinv");

    // Sanity check
    CH_assert(a_dest.interval().contains(a_destComp));

    a_dest.setVal(a_scale, a_destComp);
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the covariant metric elements
// gdn_{mu,nu} = Sum over rho [ dx^{rho}/dXi^{mu} * dx^{rho}/dXi^{nu} ]
// -----------------------------------------------------------------------------
void CartesianMap::fill_gdn (FArrayBox&      a_dest,
                             const int       a_destComp,
                             const int       a_mu,
                             const int       a_nu,
                             const RealVect& a_dXi,
                             const Real      a_scale) const
{
    CH_TIME("CartesianMap::fill_gdn");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));
    CH_assert((0 <= a_nu) && (a_nu < SpaceDim));

    if (a_mu == a_nu) {
        a_dest.setVal(a_scale, a_destComp);
    } else {
        a_dest.setVal(0.0, a_destComp);
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the contravariant metric elements
// gup^{mu,nu} = Sum over rho [ dXi^{mu}/dx^{rho} * dXi^{nu}/dx^{rho} ]
// -----------------------------------------------------------------------------
void CartesianMap::fill_gup (FArrayBox&      a_dest,
                             const int       a_destComp,
                             const int       a_mu,
                             const int       a_nu,
                             const RealVect& a_dXi,
                             const Real      a_scale) const
{
    CH_TIME("CartesianMap::fill_gup");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));
    CH_assert((0 <= a_nu) && (a_nu < SpaceDim));

    if (a_mu == a_nu) {
        a_dest.setVal(a_scale, a_destComp);
    } else {
        a_dest.setVal(0.0, a_destComp);
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with detJ * gup
// -----------------------------------------------------------------------------
void CartesianMap::fill_Jgup (FArrayBox&      a_dest,
                              const int       a_destComp,
                              const int       a_mu,
                              const int       a_nu,
                              const RealVect& a_dXi,
                              const Real      a_scale) const
{
    CH_TIME("CartesianMap::fill_Jgup");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_mu) && (a_mu < SpaceDim));
    CH_assert((0 <= a_nu) && (a_nu < SpaceDim));

    if (a_mu == a_nu) {
        a_dest.setVal(a_scale, a_destComp);
    } else {
        a_dest.setVal(0.0, a_destComp);
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the connection elements
// Note that although this function is provided and works, it is never used
// by our NS algorithm.
// -----------------------------------------------------------------------------
void CartesianMap::fill_Gamma (FArrayBox&      a_dest,
                               const int       a_destComp,
                               const int       a_up,
                               const int       a_dn1,
                               const int       a_dn2,
                               const RealVect& a_dXi,
                               const Real      a_scale) const
{
    CH_TIME("CartesianMap::fill_Gamma");

    // Sanity checks
    CH_assert(a_dest.interval().contains(a_destComp));
    CH_assert((0 <= a_up) && (a_up < SpaceDim));
    CH_assert((0 <= a_dn1) && (a_dn1 < SpaceDim));
    CH_assert((0 <= a_dn2) && (a_dn2 < SpaceDim));

    // Start with zero
    a_dest.setVal(0.0, a_destComp);
}
