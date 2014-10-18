/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Copyright (C) 2014 Edward Santilli & Alberto Scotti
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
#include "BeamGenerationBCUtil.H"
#include "EllipticBCUtils.H"
#include "ProblemContext.H"
#include "BoxIterator.H"
#include "Debug.H"


bool BeamGenerationBCUtil::s_tidalParamsRead = false;
Real BeamGenerationBCUtil::s_tidalOmega = 0.0;
Real BeamGenerationBCUtil::s_tidalU0 = 0.0;


// -----------------------------------------------------------------------------
// This not only sets BC values, but also determines the background
// stratification. WARNING: Since this is a static member, do not rely on static
// data set by the constructor such as tidal params.
// -----------------------------------------------------------------------------
void BeamGenerationBCUtil::bscalBCValues (Real*           a_pos,
                                          int*            a_dir,
                                          Side::LoHiSide* a_side,
                                          Real*           a_value,
                                          Real            a_derivScale,
                                          Real            a_time)
{
    // Linear stratification
    const static Real Nsq = 0.00001493;
    a_value[0] = -Nsq * a_pos[CH_SPACEDIM-1];

    // // Quadratic stratification
    // const static Real NsqTop = 7.00;
    // const static Real NsqBottom = 14.93;
    // Real z = a_pos[CH_SPACEDIM-1];
    // Real Nsq = NsqBottom - (NsqBottom - NsqTop) * (z / LevelGeometry::getDomainLength(CH_SPACEDIM-1));
    // a_value[0] = -Nsq * z;
}


// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
BeamGenerationBCUtil::BeamGenerationBCUtil ()
{
    if (s_tidalParamsRead) return;

    const ProblemContext* ctx = ProblemContext::getInstance();

    s_tidalOmega = ctx->tidalOmega;
    s_tidalU0 = ctx->tidalU0;

    s_tidalParamsRead = true;
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
BeamGenerationBCUtil::~BeamGenerationBCUtil ()
{;}


// -----------------------------------------------------------------------------
// Factory
// -----------------------------------------------------------------------------
PhysBCUtil* BeamGenerationBCUtil::newPhysBCUtil () const
{
    PhysBCUtil* newBCPtr = new BeamGenerationBCUtil();
    return newBCPtr;
}


// -----------------------------------------------------------------------------
// This is in case the BC's have an effect on the timestep.
// Pass in currently computed dt, along with the cfl and dx. If the effect
// of the BCs requires a decreased timestep, then the newly reduced timestep
// is returned. In the default case, this just returns a_dt back; however,
// derived classes may actually have an effect.
// -----------------------------------------------------------------------------
void BeamGenerationBCUtil::computeBoundaryDt (Real&                a_dt,
                                              const Real           a_cfl,
                                              const LevelGeometry& a_levGeo) const
{
    // Do nothing.
}


// -----------------------------------------------------------------------------
// Fills a FAB with the initial scalars
// -----------------------------------------------------------------------------
void BeamGenerationBCUtil::setScalarIC (FArrayBox&           a_scalarFAB,
                                        const int            a_scalarComp,
                                        const LevelGeometry& a_levGeo,
                                        const DataIndex&     a_di) const
{
    this->setBackgroundScalar(a_scalarFAB, a_scalarComp, a_levGeo, a_di, 0.0);
}


// -----------------------------------------------------------------------------
// Fills a FAB with the background scalar
// -----------------------------------------------------------------------------
void BeamGenerationBCUtil::setBackgroundScalar (FArrayBox&           a_scalarFAB,
                                                const int            a_scalarComp,
                                                const LevelGeometry& a_levGeo,
                                                const DataIndex&     a_di,
                                                const Real           a_time) const
{
    CH_assert(a_scalarFAB.nComp() == 1);

    if (s_useBackgroundScalar && a_scalarComp == 0) {
        // Get Cartesian coordinates at each point of a_scalarFAB.
        FArrayBox posFAB(a_scalarFAB.box(), 1);
        const RealVect& dx = a_levGeo.getDx();
        LevelGeometry::getGeoSourcePtr()->fill_physCoor(posFAB, 0, SpaceDim-1, dx);

        // Loop over each point of a_scalarFAB and set to background value.
        BoxIterator bit(a_scalarFAB.box());
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& iv = bit();

            Real pos[CH_SPACEDIM];
            pos[CH_SPACEDIM-1] = posFAB(iv);

            int dirDummy;
            Side::LoHiSide sideDummy;
            Real derivScaleDummy = 1.0e300;
            Real value[1];
            bscalBCValues (pos, &dirDummy, &sideDummy, value, derivScaleDummy, a_time);
            a_scalarFAB(iv,a_scalarComp) = value[0];
        }
    } else {
        // We do not have a background scalar.
        a_scalarFAB.setVal(0.0);
    }
}


// -----------------------------------------------------------------------------
// Sets the (vector scaled) target velocity for the sponge layer. By default,
// this function persuades the velocity field to approach its inflow value.
// -----------------------------------------------------------------------------
void BeamGenerationBCUtil::fillVelSpongeLayerTarget (FArrayBox&           a_target,
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

    // Target is the inflow/outflow value.
    if (a_spongeDir == 0) {
        if (a_velComp == 0) {
            a_target.setVal(s_tidalU0 * sin(s_tidalOmega * a_time));
        } else {
            a_target.setVal(0.0);
        }
    } else {
        MayDay::Error("BeamGenerationBCUtil::fillVelSpongeLayerTarget "
                      "can only set a sponge target when a_spongeDir = 0");
    }
}


// -----------------------------------------------------------------------------
// basicVelFuncBC
// Sets physical BCs on velocities.
// -----------------------------------------------------------------------------
BCMethodHolder BeamGenerationBCUtil::basicVelFuncBC (int a_veldir, bool a_isViscous) const
{
    const IntVect hUnit = IntVect::Unit - BASISV(CH_SPACEDIM-1);
    const IntVect vUnit = BASISV(CH_SPACEDIM-1);

    BCMethodHolder holder;

    if (a_veldir == 0) {
        //               Freeslip
        // u: Extrap 0 |==========| Extrap 0
        //                Diri 0

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

        // Free slip on high end in vertical dir
        RefCountedPtr<BCGhostClass> hiVertBCPtr = RefCountedPtr<BCGhostClass>(
            new BasicVelocityBCGhostClass(0.0,             // inflowVel
                                          -1,              // inflowDir
                                          Side::Lo,        // inflowSide
                                          -1,              // outflowDir
                                          Side::Hi,        // outflowSide
                                          a_veldir,
                                          false,           // isViscous?
                                          IntVect(D_DECL(0,0,0)),
                                          vUnit)
        );
        holder.addBCMethod(hiVertBCPtr);


        // No slip on low end in vertical dir
        RefCountedPtr<BCGhostClass> loVertBCPtr = RefCountedPtr<BCGhostClass>(
            new BasicVelocityBCGhostClass(0.0,             // inflowVel
                                          -1,              // inflowDir
                                          Side::Lo,        // inflowSide
                                          -1,              // outflowDir
                                          Side::Hi,        // outflowSide
                                          a_veldir,
                                          a_isViscous,
                                          vUnit,
                                          IntVect(D_DECL(0,0,0)))
        );
        holder.addBCMethod(loVertBCPtr);

    } else {
        //                Diri 0
        // w:   Neum 0 |==========| Neum 0
        //                Diri 0

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

        // Diri 0 in vertical dir
        RefCountedPtr<BCGhostClass> vertBCPtr = RefCountedPtr<BCGhostClass>(
            new BasicVelocityBCGhostClass(0.0,             // inflowVel
                                          -1,              // inflowDir
                                          Side::Lo,        // inflowSide
                                          -1,              // outflowDir
                                          Side::Hi,        // outflowSide
                                          a_veldir,
                                          a_isViscous,
                                          vUnit,
                                          vUnit)
        );
        holder.addBCMethod(vertBCPtr);
    }

    return holder;
}


// -----------------------------------------------------------------------------
// velSlopeBC
// Sets CC bdry slopes (undivided differences) in the tracing scheme.
// NOTE: This needs to set BCs on the velocity in the Cartesian basis.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder,SpaceDim> BeamGenerationBCUtil::velSlopeBC (int a_velComp, bool a_isViscous) const
{
    BCMethodHolder holder;

    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticConstDiriBCGhostClass(RealVect::Zero,
                                          RealVect::Zero,
                                          BASISV(0),
                                          BASISV(0))
    );
    holder.addBCMethod(BCPtr);

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir)
        bcVec[idir] = holder;

    return bcVec;
}


// -----------------------------------------------------------------------------
// diffusiveSourceFuncBC (Used to calculate the diffusive term nu.L[scalar])
// -----------------------------------------------------------------------------
BCMethodHolder BeamGenerationBCUtil::diffusiveSourceFuncBC () const
{
    // return this->basicScalarFuncBC();

    const IntVect loNeumDirs = BASISV(SpaceDim-1);
    const IntVect loDiriDirs = IntVect::Unit - loNeumDirs;
    BCMethodHolder holder;

    // Let buoyancy perturbation approach background value at top and sides.
    RefCountedPtr<BCGhostClass> diriBCPtr(
        new EllipticConstDiriBCGhostClass(RealVect::Zero,
                                          RealVect::Zero,
                                          loDiriDirs,
                                          IntVect::Unit)
    );
    holder.addBCMethod(diriBCPtr);

    // Use zero normal gradient on total buoyancy at lower vertical boundary.
    RefCountedPtr<BCGhostClass> neumBCPtr(
        new BeamGeneratorNeumBCGhostClass(this)
    );
    holder.addBCMethod(neumBCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// diffusiveSolveFuncBC (used in scalar TGA solves)
// -----------------------------------------------------------------------------
BCMethodHolder BeamGenerationBCUtil::diffusiveSolveFuncBC () const
{
    // This works better than 2nd order extrap for at least Diri BCs.
    // extrap exposes the flaws at box boundaries.

    const IntVect loNeumDirs = BASISV(SpaceDim-1);
    const IntVect loDiriDirs = IntVect::Unit - loNeumDirs;
    BCMethodHolder holder;

    // Let buoyancy perturbation approach background value at top and sides.
    RefCountedPtr<BCGhostClass> diriBCPtr(
        new EllipticConstDiriBCGhostClass(RealVect::Zero,
                                          RealVect::Zero,
                                          loDiriDirs,
                                          IntVect::Unit)
    );
    holder.addBCMethod(diriBCPtr);

    // Use zero normal gradient on total buoyancy at lower vertical boundary.
    RefCountedPtr<BCGhostClass> neumBCPtr(
        new BeamGeneratorNeumBCGhostClass(this)
    );
    holder.addBCMethod(neumBCPtr);

    return holder;
}




#include "ExtrapolationUtils.H"
// -----------------------------------------------------------------------------
// Sets nhat.Grad[b'] = -nhat.Grad[bbar].
// Only sets BCs on bottom boundary.
// TODO: Extend this to work at any boundary.
// -----------------------------------------------------------------------------
void BeamGeneratorNeumBCGhostClass::operator() (FArrayBox&           a_state,
                                                const FArrayBox*     a_extrapPtr,
                                                const Box&           a_valid,
                                                const ProblemDomain& a_domain,
                                                const RealVect&      a_dx,
                                                const DataIndex&     a_index,
                                                const FluxBox*       a_JgupPtr,
                                                bool                 a_homogeneous,
                                                Real                 a_time,
                                                const Interval&      a_interval) const
{
    CH_TIME("BeamGeneratorNeumBCGhostClass::operator()");

    // Sanity checks
    CH_assert(a_state.box().type() == IntVect::Zero);
    CH_assert(a_valid.type() == IntVect::Zero);
    CH_assert(a_JgupPtr != NULL);
    CH_assert(m_physBCPtr != NULL);
    CH_assert(a_state.nComp() == 1); // For now, this is all that is needed.

    // Directional info.
    const int idir = SpaceDim-1;
    const Side::LoHiSide iside = Side::Lo;
    const int isign = sign(iside);

    D_TERM(
        ;,
        const int jdir = (idir + 1) % SpaceDim;,
        const int kdir = (idir + 2) % SpaceDim;
    )
    D_TERM(
        ;,
        const IntVect j = BASISV(jdir);,
        const IntVect k = BASISV(kdir);
    )
    CH_assert(idir != jdir);
#if CH_SPACEDIM == 3
    CH_assert(idir != kdir);
    CH_assert(jdir != kdir);
#endif

    // I don't see why the domain would be periodic at the bottom, but we
    // had better check before wasting our time.
    if (a_domain.isPeriodic(idir)) return;

    // Find the boundary locations. If we are not at a boundary, scram.
    const Box faceBox = bdryBox(a_valid, idir, iside, 1)
                      & bdryBox(a_domain.domainBox(), idir, iside, 1);
    if (faceBox.isEmpty()) return;


    // Gather geometric info.
    const FArrayBox& G = (*a_JgupPtr)[idir];

    // Create a LevelGeometry object. This will be used to set the background
    // scalar and compute the bottom topography. As long as this was
    // staticDefined elsewhere, this should not be an expensive operation.
    LevelGeometry levGeo(a_dx);

    // Fill the background scalar at the bottom two rows of cells.
    FArrayBox bkgdState;
    if (!a_homogeneous) {
        Box stencilBox = faceBox;
        stencilBox.shiftHalf(idir, 1);
        stencilBox.grow(1);
        stencilBox.growHi(idir, -1);

        bkgdState.define(stencilBox, 1);
        m_physBCPtr->setBackgroundScalar(bkgdState, 0, levGeo, a_index, a_time);
    }

    // We are done with the geometric stuff.


    // Now, we have to alias and extrapolate the state data.

    // Alias only the components we want to alter
    const Interval interv = (a_interval == Interval())? a_state.interval(): a_interval;
    const int scomp = interv.begin();
    const int ncomp = interv.size();
    CH_assert(a_state.nComp() > interv.end());
    FArrayBox stateAlias(interv, a_state);

    // Allocate space for extrapolated values.
    // Used to compute non-diagonal derivatives at boundaries.
    // This is identical to using one-sided derivative stencils.
    FArrayBox extrapAlias;
    FArrayBox* thisExtrapPtr = NULL;
    if (a_extrapPtr == NULL) {
        thisExtrapPtr = new FArrayBox(a_state.box(), ncomp);
        extrapAlias.define(thisExtrapPtr->interval(), *thisExtrapPtr);
    } else {
        thisExtrapPtr = const_cast<FArrayBox*>(a_extrapPtr);
        extrapAlias.define(interv, *thisExtrapPtr);
    }

    // Extrapolate the boundary values for cross derivatives
    const int extrapOrder = 2;
    ExtrapolateFaceAndCopy(extrapAlias, stateAlias, a_valid, idir, iside, extrapOrder);

    // We are done prepping the state.


    // Find the source data locations
    Box srcBox = faceBox;
    srcBox.shiftHalf(idir, -isign);
    CH_assert(a_state.box().contains(srcBox));

    // Find the ghost locations
    Box ghostBox = srcBox;
    ghostBox.shift(idir, isign);
    CH_assert(a_state.box().contains(ghostBox));


    // Finally, loop over the ghost box and fill it!
    stateAlias.copy(stateAlias, srcBox, 0, ghostBox, 0, ncomp);
    if (!a_homogeneous) {
        for (int icomp = 0; icomp < ncomp; icomp++) {
            stateAlias.plus(bkgdState, ghostBox, ghostBox, -1.0, 0, icomp, 1);
            stateAlias.plus(bkgdState, srcBox  , ghostBox, +1.0, 0, icomp, 1);
        }
    }

    Real bkgdCross = 0.0;
    for (BoxIterator bit(ghostBox); bit.ok(); ++bit) {
        const IntVect& ccGhost = bit();
        const IntVect ccNear = ccGhost - isign*BASISV(idir);
        const IntVect fcBdry = ((iside == Side::Hi)? ccGhost: (ccGhost + BASISV(idir)));

        if (!a_homogeneous) {
            bkgdCross = 0.25 * (D_TERM(
                0.0,

                + (+ bkgdState(ccGhost + j)
                   - bkgdState(ccGhost - j)
                   + bkgdState(ccNear + j)
                   - bkgdState(ccNear - j)) * G(fcBdry,jdir) / a_dx[jdir],

                + (+ bkgdState(ccGhost + k)
                   - bkgdState(ccGhost - k)
                   + bkgdState(ccNear + k)
                   - bkgdState(ccNear - k)) * G(fcBdry,kdir) / a_dx[kdir]
            ));
        }

        for (int icomp = 0; icomp < ncomp; icomp++) {
            Real extrapCross = 0.25 * (D_TERM(
                0.0,

                + (+ extrapAlias(ccGhost + j, icomp)
                   - extrapAlias(ccGhost - j, icomp)
                   + extrapAlias(ccNear + j, icomp)
                   - extrapAlias(ccNear - j, icomp)) * G(fcBdry,jdir) / a_dx[jdir],

                + (+ extrapAlias(ccGhost + k, icomp)
                   - extrapAlias(ccGhost - k, icomp)
                   + extrapAlias(ccNear + k, icomp)
                   - extrapAlias(ccNear - k, icomp)) * G(fcBdry,kdir) / a_dx[kdir]
            ));

            Real cross = (bkgdCross + extrapCross) * a_dx[idir] / G(fcBdry,idir);

            stateAlias(ccGhost, icomp) -= Real(isign) * cross;

        } // end loop over comps
    } // end loop over ghost box


    // Delete extrap space if we created it
    if (a_extrapPtr == NULL) {
        delete thisExtrapPtr;
    }
}
