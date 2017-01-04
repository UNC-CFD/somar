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
#include "PhysBCUtil.H"
#include "PhysBCUtilF_F.H"
#include "ProblemContext.H"
#include "EllipticBCUtils.H"
#include "LevelGeometry.H"
#include "BoxIterator.H"
#include "BGScalarProfiles.H"
#include "ParmParse.H"
#include "Debug.H"


// Declare static members
bool     PhysBCUtil::m_isStaticDefined = false;
RealVect PhysBCUtil::s_domLength = RealVect::Zero;

bool     PhysBCUtil::s_useBackgroundScalar = false;
int      PhysBCUtil::s_bgScalarProfile = ProblemContext::BGScalarProfile::NONE;
RefCountedPtr<EllipticBCValueClass> PhysBCUtil::s_stdProfilePtr;

RealVect PhysBCUtil::s_tidalU0 = RealVect::Zero;
Real     PhysBCUtil::s_tidalOmega = 0.0;
bool     PhysBCUtil::s_doTidalFlow = false;

bool     PhysBCUtil::s_useSpongeLayer = false;
Real     PhysBCUtil::s_spongeWidth[CH_SPACEDIM][2];
Real     PhysBCUtil::s_spongeDtMult[CH_SPACEDIM][2];


// WARNING: Do not add scalars. SOMAR can't handle it yet.
const PhysBCUtil::ScalarMetaData PhysBCUtil::s_scalarMetaData[ScalarIndex::_COUNT] = {
    {"buoyancy", 1, IntVect::Unit}
};


// ************************ Constructors / destructors *************************

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
PhysBCUtil::PhysBCUtil ()
{;}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
PhysBCUtil::~PhysBCUtil ()
{;}


// -----------------------------------------------------------------------------
// Define constructor.
// -----------------------------------------------------------------------------
void PhysBCUtil::define ()
{
    if (!isStaticDefined()) {
        staticDefine();
    }
}


// -----------------------------------------------------------------------------
// Static utility
// Reads BC data from file.
// -----------------------------------------------------------------------------
void PhysBCUtil::staticDefine ()
{
    // Are we already static defined? If so, scram.
    if (isStaticDefined()) return;

    // Set up some basic structures.
    ParmParse ppIBC("ibc");
    const ProblemContext* ctx = ProblemContext::getInstance();
    s_domLength = ctx->domainLength;

    // Background scalar stuff...
    s_useBackgroundScalar = ctx->useBackgroundScalar;
    s_bgScalarProfile = ctx->bgScalarProfile;

    switch (s_bgScalarProfile) {
        case ProblemContext::BGScalarProfile::LINEAR:
        {
            Real Nsq;

            ParmParse ppIBC("ibc");
            ppIBC.get("Nsq", Nsq);

            pout() << "\nStratification profile details:\n";
            pout() << "\tNsq = " << Nsq << endl;

            EllipticBCValueClass* newPtr = new LinearBGScalarProfile(Nsq);
            s_stdProfilePtr = RefCountedPtr<EllipticBCValueClass>(newPtr);
            break;
        }

        case ProblemContext::BGScalarProfile::QUADRATIC:
        {
            Real NsqBottom;
            Real NsqTop;

            ParmParse ppIBC("ibc");
            ppIBC.get("NsqBottom", NsqBottom);
            ppIBC.get("NsqTop", NsqTop);

            pout() << "\nStratification profile details:\n";
            pout() << "\tNsqBottom = " << NsqBottom << endl;
            pout() << "\tNsqTop = " << NsqTop << endl;

            EllipticBCValueClass* newPtr = new QuadraticBGScalarProfile(NsqBottom, NsqTop);
            s_stdProfilePtr = RefCountedPtr<EllipticBCValueClass>(newPtr);
            break;
        }

        case ProblemContext::BGScalarProfile::TANH:
        {
            Real z0;      // Center location
            Real delta;   // Thickness
            Real rho0;    // Median buoyancy value
            Real drho;    // Max buoyancy difference

            ParmParse ppIBC("ibc");
            ppIBC.get("z0", z0);
            ppIBC.get("delta", delta);
            ppIBC.get("rho0", rho0);
            ppIBC.get("drho", drho);

            pout() << "\nStratification profile details:\n";
            pout() << "\tz0 = " << z0 << "\n";
            pout() << "\tdelta = " << delta << "\n";
            pout() << "\trho0 = " << rho0 << "\n";
            pout() << "\tdrho = " << drho << endl;

            EllipticBCValueClass* newPtr = new TanhBGScalarProfile(z0, delta, rho0, drho);
            s_stdProfilePtr = RefCountedPtr<EllipticBCValueClass>(newPtr);
            break;
        }
    }


    // Tidal flow stuff...
    s_tidalU0 = ctx->tidalU0;
    s_tidalOmega = ctx->tidalOmega;
    s_doTidalFlow = (abs(s_tidalU0.sum() * s_tidalOmega) > 1.0e-12);


    // Sponge layer stuff...
    s_useSpongeLayer = ctx->useSpongeLayer;
    for (int dir = 0; dir < SpaceDim; ++dir) {
        // Read in Cartesian sponge widths and dt multipliers.
        s_spongeWidth[dir][0] = ctx->spongeWidth[dir][0];
        s_spongeWidth[dir][1] = ctx->spongeWidth[dir][1];
        s_spongeDtMult[dir][0] = ctx->spongeDtMult[dir][0];
        s_spongeDtMult[dir][1] = ctx->spongeDtMult[dir][1];


        // Next, we need to convert the Cartesian sponge widths
        // to curvilinear sponge widths...

        // Gather coarse level domain info.
        const ProblemDomain& domain = ctx->domain;
        const Box& domBox = domain.domainBox();
        const RealVect dXi = s_domLength / RealVect(domBox.size());
        const RealVect loXi = RealVect(domBox.smallEnd()) * dXi;
        const RealVect hiXi = loXi + s_domLength;

        // Make a strip of cells along this dir.
        Box lineBox;
        {
            IntVect smallEnd = domBox.smallEnd();
            IntVect bigEnd = smallEnd;
            bigEnd[dir] = domBox.bigEnd(dir);
            lineBox.define(smallEnd, bigEnd);
        }

        // Fill that strip with Cartesian locations.
        // Technically, this won't work if the boundaries have irregular shapes.
        FArrayBox fcCartPosFAB(surroundingNodes(lineBox, dir), 1);
        FArrayBox cartPosFAB(lineBox, 1);
        {
            if (!LevelGeometry::isStaticDefined()) {
                LevelGeometry::staticDefine();
            }

            LevelGeometry::getGeoSourcePtr()->fill_physCoor(fcCartPosFAB,
                                                            0,   // dest comp
                                                            dir, // coor index
                                                            dXi);

            LevelGeometry::getGeoSourcePtr()->fill_physCoor(cartPosFAB,
                                                            0,   // dest comp
                                                            dir, // coor index
                                                            dXi);
        }

        // Rescale sponge widths on lo end.
        if (s_spongeWidth[dir][0] > 0.0) {
            const Real xmin = fcCartPosFAB(fcCartPosFAB.box().smallEnd());
            const Real xmax = xmin + s_spongeWidth[dir][0];
            Real xi = loXi[dir];
            Real x, newXi;

            BoxIterator bit(lineBox);
            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& cc = bit();
                x = cartPosFAB(cc);
                if (x <= xmax) {
                    newXi = (Real(cc[dir]) + 0.5) * dXi[dir];
                    xi = max(xi, newXi);
                }
            }
            s_spongeWidth[dir][0] = xi - loXi[dir];
        }

        // Rescale sponge widths on hi end.
        if (s_spongeWidth[dir][1] > 0.0) {
            const Real xmax = fcCartPosFAB(fcCartPosFAB.box().bigEnd());
            const Real xmin = xmax - s_spongeWidth[dir][1];
            Real xi = hiXi[dir];
            Real x, newXi;

            BoxIterator bit(lineBox);
            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& cc = bit();
                x = cartPosFAB(cc);
                if (xmin <= x) {
                    newXi = (Real(cc[dir]) + 0.5) * dXi[dir];
                    xi = min(xi, newXi);
                }
            }
            s_spongeWidth[dir][1] = hiXi[dir] - xi;
        }
    }


    // We are now static defined.
    m_isStaticDefined = true;
}


// ************************ ICs / background fields ****************************

// -----------------------------------------------------------------------------
// Fills a FAB with the initial velocity
// Locations are in mapped space, but components are Cartesian.
// -----------------------------------------------------------------------------
void PhysBCUtil::setVelIC (FArrayBox&           a_velFAB,
                           const int            a_velComp,
                           const LevelGeometry& a_levGeo,
                           const DataIndex&     a_di) const
{
    CH_assert(a_velFAB.nComp() == SpaceDim);
    CH_assert(0 <= a_velComp);
    CH_assert(a_velComp < SpaceDim);

    if (s_doTidalFlow && SpaceDim == 3 && a_velComp == 1) {
        a_velFAB.setVal(s_tidalU0[1], 1);
    } else {
        a_velFAB.setVal(0.0, a_velComp);
    }
}


// -----------------------------------------------------------------------------
// Fills a FAB with the initial scalars
// -----------------------------------------------------------------------------
void PhysBCUtil::setScalarIC (FArrayBox&           a_scalarFAB,
                              const int            a_scalarComp,
                              const LevelGeometry& a_levGeo,
                              const DataIndex&     a_di) const
{
    a_scalarFAB.setVal(0.0);
}


// -----------------------------------------------------------------------------
// Fills a FAB with the background scalar
// -----------------------------------------------------------------------------
void PhysBCUtil::setBackgroundScalar (FArrayBox&           a_scalarFAB,
                                      const int            a_scalarComp,
                                      const LevelGeometry& a_levGeo,
                                      const DataIndex&     a_di,
                                      const Real           a_time) const
{
    CH_assert(a_scalarFAB.nComp() == 1);

    switch (s_bgScalarProfile) {
        case ProblemContext::BGScalarProfile::NONE:
        {
            a_scalarFAB.setVal(0.0);
            break;
        }

        case ProblemContext::BGScalarProfile::USER_DEFINED:
        {
            MayDay::Error("PhysBCUtil::setBackgroundScalar needs to be overridden "
                          "if you plan to use a user-defined background scalar");
            break;
        }

        default:
        {
            CH_assert(!s_stdProfilePtr.isNull());

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
                Real value[1];
                (*s_stdProfilePtr)(pos, &dirDummy, &sideDummy, value, dx, a_time);
                a_scalarFAB(iv,a_scalarComp) = value[0];
            }
            break;
        }
    }
}


// ******************** Miscellaneous utility functions ************************

// -----------------------------------------------------------------------------
// Sets the velocity ICs over an entire level.
// Locations are in mapped space, but components are Cartesian.
// -----------------------------------------------------------------------------
void PhysBCUtil::setVelIC (LevelData<FArrayBox>& a_vel,
                           const int             a_velComp,
                           const LevelGeometry&  a_levGeo) const
{
    DataIterator dit = a_vel.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        setVelIC(a_vel[dit], a_velComp, a_levGeo, dit());
    }
}


// -----------------------------------------------------------------------------
// Sets the scalar ICs over an entire level
// -----------------------------------------------------------------------------
void PhysBCUtil::setScalarIC (LevelData<FArrayBox>& a_scalar,
                              const int             a_scalarComp,
                              const LevelGeometry&  a_levGeo) const
{
    DataIterator dit = a_scalar.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        setScalarIC(a_scalar[dit], a_scalarComp, a_levGeo, dit());
    }
}


// -----------------------------------------------------------------------------
// Sets the background scalars over an entire level
// -----------------------------------------------------------------------------
void PhysBCUtil::setBackgroundScalar (LevelData<FArrayBox>& a_scalar,
                                      const int             a_scalarComp,
                                      const Real            a_time,
                                      const LevelGeometry&  a_levGeo) const
{
    DataIterator dit = a_scalar.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        setBackgroundScalar(a_scalar[dit], a_scalarComp, a_levGeo, dit(), a_time);
    }
}


// -----------------------------------------------------------------------------
// Adds the background scalar to an existing LevelData
// -----------------------------------------------------------------------------
void PhysBCUtil::addBackgroundScalar (LevelData<FArrayBox>& a_scalar,
                                      const int             a_scalarComp,
                                      const Real            a_time,
                                      const LevelGeometry&  a_levGeo) const
{
    CH_assert(0 <= a_scalarComp);

    if (s_useBackgroundScalar) {
        const DisjointBoxLayout& grids = a_scalar.getBoxes();
        DataIterator dit = grids.dataIterator();

        const int numComps = a_scalar.nComp();
        LevelData<FArrayBox> backgroundScal(grids, numComps, a_scalar.ghostVect());
        setBackgroundScalar(backgroundScal, a_scalarComp, a_time, a_levGeo);

        for (dit.reset(); dit.ok(); ++dit) {
            a_scalar[dit].plus(backgroundScal[dit], 0, 0, numComps);
        }
    }
}


// -----------------------------------------------------------------------------
// Subtracts the background scalar from an existing LevelData
// -----------------------------------------------------------------------------
void PhysBCUtil::subtractBackgroundScalar (LevelData<FArrayBox>& a_scalar,
                                           const int             a_scalarComp,
                                           const Real            a_time,
                                           const LevelGeometry&  a_levGeo) const
{
    CH_assert(0 <= a_scalarComp);

    if (s_useBackgroundScalar) {
        const DisjointBoxLayout& grids = a_scalar.getBoxes();
        DataIterator dit = grids.dataIterator();

        const int numComps = a_scalar.nComp();
        LevelData<FArrayBox> backgroundScal(grids, numComps, a_scalar.ghostVect());
        setBackgroundScalar(backgroundScal, a_scalarComp, a_time, a_levGeo);

        for (dit.reset(); dit.ok(); ++dit) {
            a_scalar[dit].minus(backgroundScal[dit], 0, 0, numComps);
        }
    }
}


// -----------------------------------------------------------------------------
// TODO: This should be moved to StratUtils.
// Computes the Brunt–Väisälä frequency on a single grid.
// All ins and outs must be CC.
// a_bFAB is the background buoyancy field.
// a_dXidz needs SpaceDim comps...(dXi/dz, dNu/dz, dZeta/dz).
// This function is a bit raw, but dXidz is expensive to compute and should
// be cached by the user.
// -----------------------------------------------------------------------------
void PhysBCUtil::computeNSq (FArrayBox&       a_NsqFAB,
                             const FArrayBox& a_bFAB,
                             const FArrayBox& a_dXidzFAB,
                             const Box&       a_destBox,
                             const RealVect&  a_dx) const
{
    // Sanity checks
    CH_assert(a_NsqFAB  .nComp() == 1);
    CH_assert(a_bFAB    .nComp() == 1);
    CH_assert(a_dXidzFAB.nComp() == SpaceDim);

    CH_assert(a_NsqFAB  .box().type() == a_destBox.type());
    CH_assert(a_bFAB    .box().type() == a_destBox.type());
    CH_assert(a_dXidzFAB.box().type() == a_destBox.type());
    CH_assert(a_destBox       .type() == IntVect::Zero);

    CH_assert(a_NsqFAB  .box().contains(a_destBox));
    CH_assert(a_bFAB    .box().contains(grow(a_destBox,1)));
    CH_assert(a_dXidzFAB.box().contains(a_destBox));

    // Calculate -dXi^i/dz * dB/dXi^i
    FORT_COMPUTE_CCNSQ (
        CHF_FRA1(a_NsqFAB,0),
        CHF_CONST_FRA1(a_bFAB,0),
        CHF_CONST_FRA(a_dXidzFAB),
        CHF_CONST_REALVECT(a_dx),
        CHF_BOX(a_destBox));
}


// -----------------------------------------------------------------------------
// TODO: This should be moved to StratUtils.
// Computes the Brunt–Väisälä frequency over an entire level.
// a_Nsq must be CC.
// -----------------------------------------------------------------------------
void PhysBCUtil::computeNSq (LevelData<FArrayBox>& a_Nsq,
                             const LevelGeometry&  a_levGeo,
                             const Real            a_time) const
{
    // Sanity checks
    CH_assert(a_levGeo.getBoxes().physDomain() == a_Nsq.getBoxes().physDomain());
    CH_assert(a_levGeo.getBoxes().compatible(a_Nsq.getBoxes()));

    // Declare variables
    const RealVect& dx = a_levGeo.getDx();
    const GeoSourceInterface& geoSource = *(a_levGeo.getGeoSourcePtr());
    DataIterator dit = a_Nsq.dataIterator();

    // Loop over grids and compute Nsq.
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& NsqFAB = a_Nsq[dit];
        const Box& destBox = NsqFAB.box();
        CH_assert(destBox.type() == IntVect::Zero);

        // Fill dXi^i/dz vector.
        FArrayBox dXidzFAB(destBox, SpaceDim);
        for (int dir = 0; dir < SpaceDim; ++dir) {
            geoSource.fill_dXidx(dXidzFAB,
                                 dir,           // fab comp
                                 dir,           // Xi index
                                 SpaceDim-1,    // z index
                                 dx);
        }

        // Fill background bouyancy field.
        FArrayBox BFAB(grow(destBox,1), 1);
        this->setBackgroundScalar(BFAB,
                                  ScalarIndex::BUOYANCY_DEVIATION,
                                  a_levGeo,
                                  dit(),
                                  a_time);

        // Calculate -dXi^i/dz * dB/dXi^i
        this->computeNSq(NsqFAB,
                         BFAB,
                         dXidzFAB,
                         destBox,
                         dx);
    }
}


// -----------------------------------------------------------------------------
// This is in case the BC's have an effect on the timestep.
// Pass in currently computed dt, along with the cfl and dx. If the effect
// of the BCs requires a decreased timestep, then the newly reduced timestep
// is returned. In the default case, this just returns a_dt back; however,
// derived classes may actually have an effect.
// -----------------------------------------------------------------------------
void PhysBCUtil::computeBoundaryDt (Real&                a_dt,
                                    const Real           a_cfl,
                                    const LevelGeometry& a_levGeo) const
{
    // In the default case, nothing is done here.
}


// ************************* Sponge layer functions ****************************


// -----------------------------------------------------------------------------
// Splits the domain into its sponge layers and the interior. The locations of
// the splitting are given as face indices. If the domain is not split, these
// indices will lie outside of the domain.
// -----------------------------------------------------------------------------
void PhysBCUtil::computeSpongeRegions (Tuple<Box, 2>&       a_spongeBox,
                                       Tuple<int, 2>&       a_splitFaceIndex,
                                       Box&                 a_interior,
                                       const LevelGeometry& a_levGeo,
                                       const int            a_dir)
{
    // Sanity checks
    CH_assert(useSpongeLayer());
    CH_assert(0 <= a_dir);
    CH_assert(a_dir < SpaceDim);

    // Gather some needed info
    const Box domBox = a_levGeo.getDomain().domainBox();
    const RealVect& dx = a_levGeo.getDx();
    Box interiorBox[2];

    // Lo side
    int s = 0;
    if (s_spongeWidth[a_dir][s] > 0.0) {
        const int numCells = ceil(s_spongeWidth[a_dir][s] / dx[a_dir]);
        a_splitFaceIndex[s] = domBox.smallEnd(a_dir) + numCells;
        a_spongeBox[s] = domBox;
        interiorBox[s] = a_spongeBox[s].chop(a_dir, a_splitFaceIndex[s]);
    }

    // Hi side
    s = 1;
    if (s_spongeWidth[a_dir][s] > 0.0) {
        const int numCells = ceil(s_spongeWidth[a_dir][s] / dx[a_dir]);
        a_splitFaceIndex[s] = domBox.bigEnd(a_dir) - numCells + 1;
        interiorBox[s] = domBox;
        a_spongeBox[s] = interiorBox[s].chop(a_dir, a_splitFaceIndex[s]);
    }

    // Compute the interior region
    a_interior = interiorBox[0] & interiorBox[1];
}


// -----------------------------------------------------------------------------
// Sets a_srcTerm = (target - state) / (layer time scale).
// To fill velocity sponge layer sources, just use a_comp's default value.
// -----------------------------------------------------------------------------
void PhysBCUtil::fillSpongeLayerSrcTerm (LevelData<FArrayBox>&       a_srcTerm,
                                         const LevelData<FArrayBox>& a_state,
                                         const Real                  a_time,
                                         const Real                  a_dt,
                                         const LevelGeometry&        a_levGeo,
                                         const int                   a_comp)
{
    CH_TIME("PhysBCUtil::fillSpongeLayerSrcTerm");

    // Gather some needed information
    const int ncomp = a_srcTerm.nComp();
    const RealVect& dx = a_levGeo.getDx();
    const DisjointBoxLayout& grids = a_srcTerm.disjointBoxLayout();
    DataIterator dit = grids.dataIterator();

    // Sanity checks
    CH_assert(useSpongeLayer());
    CH_assert(ncomp == a_state.nComp());
    CH_assert(a_comp >=  0 || ncomp == SpaceDim);
    CH_assert(a_comp == -1 || ncomp <= getNumScalars());
    CH_assert(a_state.getBoxes() == a_srcTerm.getBoxes());

    // Start with an inactive sponge layer, then activate the sides we need.
    for (dit.reset(); dit.ok(); ++dit) {
        a_srcTerm[dit].setVal(0.0);
    }

    for (int dir = 0; dir < SpaceDim; ++dir) {
        // Calculate the sponge regions and the interior
        Tuple<Box, 2> spongeBox;
        Tuple<int, 2> splitFaceIndex;
        Box interior;
        this->computeSpongeRegions(spongeBox, splitFaceIndex, interior, a_levGeo, dir);

        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide& iside = sit();
            const int s = ((iside == Side::Lo)? 0: 1);

            // Do we have a sponge on this side?
            if (spongeBox[s].isEmpty()) continue;

            // Calculate the sponge region and it's complement
            const Real splitFaceLoc = Real(splitFaceIndex[s]) * dx[dir];

            // Set sponge layer term coefficient
            const Real coeff = 1.0 / (s_spongeDtMult[dir][s] * a_dt);

            // Loop over the grids.
            for (dit.reset(); dit.ok(); ++dit) {
                FArrayBox& srcTermFAB = a_srcTerm[dit];
                const FArrayBox& stateFAB = a_state[dit];

                const Box valid = grids[dit] & srcTermFAB.box();
                const Box validSpongeBox = valid & spongeBox[s];

                // If the sponge region is empty, just move on.
                if (validSpongeBox.isEmpty()) continue;

                // Compute the target locations
                const Box targetBox = bdryBox(validSpongeBox, dir, iside, 1);
                const int targetFaceIndex = targetBox.smallEnd(dir);
                CH_assert(targetFaceIndex == targetBox.bigEnd(dir));

                FArrayBox targetFAB(targetBox, 1);
                Real ratio, ramp, pos;

                // This is a bit strange. If a_comps is -1, then we use that as a
                // flag to fill the velocity sponge. comp will loop over all
                // SpaceDim elements of the velocity. If a_comps is anything else,
                // then ncomp = 1, comp = 0, and a_comp signals which scalar sponge
                // we need to fill.
                for (int comp = 0; comp < ncomp; ++comp) {
                    // Compute the target field.
                    // TODO: Should we pass dir in?
                    if (a_comp == -1) {
                        this->fillVelSpongeLayerTarget(targetFAB, comp, dir, iside, a_levGeo, dit(), a_time);
                    } else {
                        this->fillScalarSpongeLayerTarget(targetFAB, a_comp, dir, iside, a_levGeo, dit(), a_time);
                    }

                    // Loop over the grid and fill the sponge source term.
                    BoxIterator bit(validSpongeBox);
                    for (bit.reset(); bit.ok(); ++bit) {
                        const IntVect& cc = bit();

                        IntVect targetfc = cc;
                        targetfc[dir] = targetFaceIndex;

                        pos = (Real(cc[dir]) + 0.5) * dx[dir];
                        ratio = Abs(pos - splitFaceLoc) / s_spongeWidth[dir][s];
                        ramp = this->spongeLayerRamp(ratio);

                        srcTermFAB(cc,comp) = coeff * ramp * (targetFAB(targetfc) - stateFAB(cc,comp));

                    } // end loop over sponge box (bit)
                } // end loop over state components (comp)
            } // end loop over grids (dit)
        } // end loop over sides (sit)
    } // end loop over directions (dir)
}


// -----------------------------------------------------------------------------
// Sets the Cartesian-based target velocity for the sponge layer.
// By default, this function throws an error.
// -----------------------------------------------------------------------------
void PhysBCUtil::fillVelSpongeLayerTarget (FArrayBox&           a_target,
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

    if (s_doTidalFlow) {
        // Default velocity is the tidal velocity.
        if (a_velComp == 0) {
            a_target.setVal(s_tidalU0[0] * sin(s_tidalOmega * a_time));
        } else if (SpaceDim == 3 && a_velComp == 1) {
            a_target.setVal(s_tidalU0[1] * cos(s_tidalOmega * a_time));
        } else {
            a_target.setVal(0.0);
        }
    } else {
        // Not using tidal flow. Default value is zero.
        a_target.setVal(0.0);
    }
}


// -----------------------------------------------------------------------------
// Sets the target values for the scalar sponge layer. If we are using a
// background scalar, then this function set the perturbation to zero.
// Otherwise, an error is thrown and this function will need to be overridden.
// -----------------------------------------------------------------------------
void PhysBCUtil::fillScalarSpongeLayerTarget (FArrayBox&           a_target,
                                              const int            a_scalarComp,
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

    if (s_useBackgroundScalar) {
        switch (a_scalarComp) {
            case ScalarIndex::BUOYANCY_DEVIATION:
            {
                // The outgoing scalar should approach its unperturbed value.
                a_target.setVal(0.0, a_scalarComp);
                break;
            }
            default:
            {
                // By default, set to zero and throw a warning.
                MayDay::Warning("PhysBCUtil::fillScalarSpongeLayerTarget has not "
                                "been set up to handle all scalars. Using default "
                                "target value of zero");
                a_target.setVal(0.0, a_scalarComp);
                break;
            }
        }

    } else {
        // I don't know what values the scalar should approach.
        const char* msg = "If you plan to use a sponge layer without a background scalar, "
                          "then you need to override PhysBCUtil::fillScalarSpongeLayerTarget";
        MayDay::Error(msg);
    }
}


// ****************************** Velocity BCs *********************************

// -----------------------------------------------------------------------------
// uStarFuncBC
// Pre-projection velocity BC.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder, SpaceDim> PhysBCUtil::uStarFuncBC (bool a_isViscous) const
{
    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir) {
        bcVec[idir] = this->basicVelFuncBC(idir, a_isViscous);
    }

    return bcVec;
}


// -----------------------------------------------------------------------------
// viscousSourceFuncBC
// Sets ghosts needed to calculate the viscous source term nu.L[vel]
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder, SpaceDim> PhysBCUtil::viscousSourceFuncBC () const
{
    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir) {
        bcVec[idir] = viscousSolveFuncBC(idir);
    }

    return bcVec;
}


// -----------------------------------------------------------------------------
// viscousSolveFuncBC (Used in single-component velocity TGA solves)
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::viscousSolveFuncBC (int a_dir) const
{
    // const bool isViscous = true;
    // return this->basicVelFuncBC(a_dir, isViscous);

    int extrapOrder = 2; // 2 works best in TG test
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticExtrapBCGhostClass(extrapOrder)
    );
    BCMethodHolder protoBC;
    protoBC.addBCMethod(BCPtr);

    return protoBC;
}


// -----------------------------------------------------------------------------
// viscousRefluxBC (Used in sync step by the implicit refluxing op)
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::viscousRefluxBC (int a_dir) const
{
    const bool isViscous = true;
    return this->basicVelFuncBC(a_dir, isViscous);
}


// -----------------------------------------------------------------------------
// viscousVelFuncBC
// Sets BCs on a generic viscous velocity field.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder, SpaceDim> PhysBCUtil::viscousVelFuncBC () const
{
    const bool isViscous = true;

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir) {
        bcVec[idir] = this->basicVelFuncBC(idir, isViscous);
    }

    return bcVec;
}


// -----------------------------------------------------------------------------
// tracingVelFuncBC
// Sets BCs on old-time velocity that will be used to trace characteristics.
// (Used by fill functions)
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder, SpaceDim> PhysBCUtil::tracingVelFuncBC () const
{
    const bool isViscous = false;

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir) {
        bcVec[idir] = this->basicVelFuncBC(idir, isViscous);
    }

    return bcVec;
}


// -----------------------------------------------------------------------------
// uDelUFuncBC
// Sets BCs on the advection term
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder, SpaceDim> PhysBCUtil::uDelUFuncBC (bool a_isViscous) const
{
    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir) {
        bcVec[idir] = this->basicVelFuncBC(idir, a_isViscous);
    }

    return bcVec;
}


// -----------------------------------------------------------------------------
// advectingVelFuncBC
// Sets BCs on the advecting velocity
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder, SpaceDim> PhysBCUtil::advectingVelFuncBC (bool a_isViscous) const
{
    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir) {
        bcVec[idir] = this->basicVelFuncBC(idir, a_isViscous);
    }

    return bcVec;
}


// -----------------------------------------------------------------------------
// vortFuncBC
// Sets BCs on velocity when computing vorticity. We do first-order extrap here.
// The effect is that derivatives at the boundaries become 1-sided differences.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder, SpaceDim> PhysBCUtil::vortFuncBC (bool a_isViscous) const
{
    int extrapOrder = 1;
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticExtrapBCGhostClass(extrapOrder)
    );
    BCMethodHolder protoBC;
    protoBC.addBCMethod(BCPtr);

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir)
        bcVec[idir] = protoBC;

    return bcVec;
}


// -----------------------------------------------------------------------------
// velRiemannBC
// Sets BCs on FC velocity in the Riemann solver.
// NOTE: This needs to set BCs on the velocity in the Cartesian basis.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder,SpaceDim> PhysBCUtil::velRiemannBC (int  a_velComp,
                                                         bool a_isViscous) const
{
    // Really, we don't need to set fluxes because we never request
    // velocity fluxes from the Riemann solver. Also, the BCs for each
    // direction (idir) are the same since predictScalar works on only
    // one velocity component at a time.
    BCMethodHolder protoBC = this->basicVelFuncBC(a_velComp, false);
    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir) {
        bcVec[idir] = protoBC;
    }

    return bcVec;
}


// -----------------------------------------------------------------------------
// velSlopeBC
// Sets CC bdry slopes (undivided differences) in the tracing scheme.
// NOTE: This needs to set BCs on the velocity in the Cartesian basis.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder,SpaceDim> PhysBCUtil::velSlopeBC (int a_velComp, bool a_isViscous) const
{
    // Do nothing BCs
    BCMethodHolder protoBC;

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir)
        bcVec[idir] = protoBC;

    return bcVec;
}


// ******************************* Scalar BCs **********************************

// -----------------------------------------------------------------------------
// diffusiveSourceFuncBC (Used to calculate the diffusive term nu.L[scalar])
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::diffusiveSourceFuncBC () const
{
    return this->diffusiveSolveFuncBC();
}


// -----------------------------------------------------------------------------
// diffusiveSolveFuncBC (used in scalar TGA solves)
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::diffusiveSolveFuncBC () const
{
    // return this->basicScalarFuncBC();

    int extrapOrder = 2;
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticExtrapBCGhostClass(extrapOrder)
    );
    BCMethodHolder protoBC;
    protoBC.addBCMethod(BCPtr);

    return protoBC;
}


// -----------------------------------------------------------------------------
// scalarRefluxSolveBC
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::scalarRefluxSolveBC (int a_scalarType) const
{
    return this->basicScalarFuncBC();
}


// -----------------------------------------------------------------------------
// scalarTraceFuncBC
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::scalarTraceFuncBC (int a_scalarType) const
{
    return this->basicScalarFuncBC();
}


// -----------------------------------------------------------------------------
// scalarRiemannBC
// Sets BCs on FC scalars in the Riemann solver.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder,SpaceDim> PhysBCUtil::scalarRiemannBC (int a_scalarType) const
{
    // Start with standard BC
    BCMethodHolder protoBC = this->basicScalarFuncBC();

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir)
        bcVec[idir] = protoBC;

    return bcVec;
}


// -----------------------------------------------------------------------------
// scalarSlopeBC
// Sets CC bdry slopes (undivided differences) in the tracing scheme.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder,SpaceDim> PhysBCUtil::scalarSlopeBC (int a_scalarType) const
{
    // Do nothing BCs
    BCMethodHolder protoBC;

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir)
        bcVec[idir] = protoBC;

    return bcVec;
}


// ********************** Freestream preservation BCs **************************

// -----------------------------------------------------------------------------
// lambdaFuncBC
// Chombo uses Diri 1.0 BCs at inflow and 1st order extrap elsewhere.
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::lambdaFuncBC () const
{
    BCMethodHolder holder;

    // Solid wall and outflow BCs
    int extrapOrder = 1;
    RefCountedPtr<BCGhostClass> otherBCPtr(
        new EllipticExtrapBCGhostClass(extrapOrder,
                                       IntVect::Unit,
                                       IntVect::Unit)
    );
    holder.addBCMethod(otherBCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// FreestreamCorrFuncBC
// Chombo uses 0 Neum
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::FreestreamCorrFuncBC () const
{
    BCMethodHolder holder;

    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticConstNeumBCGhostClass(IntVect::Zero, IntVect::Zero)
    );
    holder.addBCMethod(BCPtr);

    RefCountedPtr<BCFluxClass> BCFluxPtr(
        new EllipticConstNeumBCFluxClass(IntVect::Zero, IntVect::Zero)
    );
    holder.addBCMethod(BCFluxPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// gradELambdaFuncBC
// Chombo uses basicGradPressureFuncBC (ie, 2nd order extrap)
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::gradELambdaFuncBC () const
{
    return this->basicGradPressureFuncBC();
}


// -----------------------------------------------------------------------------
// lambdaRiemannBC
// Sets BCs on FC lambda in the Riemann solver.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder,SpaceDim> PhysBCUtil::lambdaRiemannBC () const
{
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticConstDiriBCGhostClass(RealVect::Unit,
                                          RealVect::Unit)
    );

    BCMethodHolder protoBC;
    protoBC.addBCMethod(BCPtr);

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir)
        bcVec[idir] = protoBC;

    return bcVec;
}


// -----------------------------------------------------------------------------
// lambdaSlopeBC
// Sets CC bdry slopes (undivided differences) in the tracing scheme.
// -----------------------------------------------------------------------------
Tuple<BCMethodHolder,SpaceDim> PhysBCUtil::lambdaSlopeBC () const
{
    // Do nothing BCs
    BCMethodHolder protoBC;

    Tuple<BCMethodHolder, SpaceDim> bcVec;
    for (int idir = 0; idir < SpaceDim; ++idir)
        bcVec[idir] = protoBC;

    return bcVec;
}


// ****************************** Pressure BCs *********************************

// -----------------------------------------------------------------------------
// MacPressureFuncBC
// Used in levelMacProject solver
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::MacPressureFuncBC () const
{
    bool isHomogeneous = false;
    return this->basicPressureFuncBC(isHomogeneous);
}


// -----------------------------------------------------------------------------
// gradMacPressureFuncBC
// Used to calculate Grad[phi] in MAC projection
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::gradMacPressureFuncBC () const
{
    return this->basicGradPressureFuncBC();
}


// -----------------------------------------------------------------------------
// LevelPressureFuncBC
// Used in LevelProject solver
// Chombo uses homog basicPressureFuncBC
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::LevelPressureFuncBC () const
{
    bool isHomogeneous = false;
    return this->basicPressureFuncBC(isHomogeneous);
}


// -----------------------------------------------------------------------------
// gradPiFuncBC
// Used to calculate CCGrad[Pi] in LevelProject
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::gradPiFuncBC () const
{
    return this->basicGradPressureFuncBC();
}


// -----------------------------------------------------------------------------
// SyncProjFuncBC
// Used in sync projection solver
// Chombo uses inhomog basicPressureFuncBC
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::SyncProjFuncBC () const
{
    bool isHomogeneous = true;
    return this->basicPressureFuncBC(isHomogeneous);
}


// -----------------------------------------------------------------------------
// gradESyncFuncBC
// Used to calculate CCGrad[eSync] in sync projection
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::gradESyncFuncBC () const
{
    return this->basicGradPressureFuncBC();
}


// ************************ Miscellaneous BC functions *************************

// -----------------------------------------------------------------------------
// smoothingSolverBC (used for post-redgrid smoothing)
// Chombo uses 0 Diri at solid walls
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::smoothingSolverBC () const
{
    int extrapOrder = 1;
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticExtrapBCGhostClass(extrapOrder)
    );

    BCMethodHolder holder;
    holder.addBCMethod(BCPtr);

    return holder;
}

// -----------------------------------------------------------------------------
// BCs for the streamfunction solver.
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::streamSolverBC (int comp) const
{
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticConstDiriBCGhostClass(RealVect::Zero,
                                          RealVect::Zero)
    );

    BCMethodHolder holder;
    holder.addBCMethod(BCPtr);

    return holder;
}



// ************************* The basic BC functions ****************************
// ****** These create BCMethodHolders from RefCountedPtr<BC*Class>es. *********

// -----------------------------------------------------------------------------
// basicVelFuncBC
// Sets physical BCs on velocities.
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::basicVelFuncBC (int a_veldir, bool a_isViscous) const
{
    BCMethodHolder holder;

    if (!s_doTidalFlow) {
        RefCountedPtr<BCGhostClass> velBCPtr = RefCountedPtr<BCGhostClass>(
            new BasicVelocityBCGhostClass(0.0,      //s_inflowVel,
                                          -1,       //s_inflowDir,
                                          Side::Lo, //s_inflowSide,
                                          -1,       //s_outflowDir,
                                          Side::Lo, //s_outflowSide,
                                          a_veldir,
                                          a_isViscous)
        );
        holder.addBCMethod(velBCPtr);

    } else {

        if (!s_useSpongeLayer) {
            MayDay::Error("If you are tidally forcing the flow, you should use a sponge. "
                          "basicVelFuncBC only extrapolates in the horizontal. "
                          "It is the sponge that actually enforces the BCs.");
        }

        const IntVect hUnit = IntVect::Unit - BASISV(CH_SPACEDIM-1);
        const IntVect vUnit = BASISV(CH_SPACEDIM-1);

        if (a_veldir < SpaceDim-1) {
            //                 Freeslip
            // u,v: Extrap 0 |==========| Extrap 0
            //                  Diri 0 (or free slip if a_viscous is false)

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
    }

    return holder;
}


// -----------------------------------------------------------------------------
// basicScalarFuncBC   (Extrapolate BCs)
// Sets physical BCs on a generic passive scalar.
// Chombo uses 1st order extrap
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::basicScalarFuncBC () const
{
    BCMethodHolder holder;

    int extrapOrder = 1;
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticExtrapBCGhostClass(extrapOrder)
    );
    holder.addBCMethod(BCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// basicPressureFuncBC   (Zero Neum BCs)
// Sets physical BCs on pressures (used by the Poisson solvers).
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::basicPressureFuncBC (bool a_isHomogeneous) const
{
    BCMethodHolder holder;

    // This sets ghosts so that Grad[CCstate] = Grad[pressure] = 0 at bdry.
    RefCountedPtr<BCGhostClass> ghostBCPtr(
        new EllipticConstNeumBCGhostClass(RealVect::Zero,
                                          RealVect::Zero)
    );
    holder.addBCMethod(ghostBCPtr);

    // This sets face values so that FCstate = Grad[pressure] = 0 at bdry.
    RefCountedPtr<BCFluxClass> fluxBCPtr(
        new EllipticConstNeumBCFluxClass(RealVect::Zero,
                                         RealVect::Zero)
    );
    holder.addBCMethod(fluxBCPtr);

    return holder;
}


// -----------------------------------------------------------------------------
// basicGradPressureFuncBC   (Extrap BCs)
// Sets physical BCs on pressures before taking gradients.
// NOTE: We don't use Neum BCs because many papers say to use Extrap BCs.
// In fact, changing these to extrap cleared up an error in the velocity field!
// -----------------------------------------------------------------------------
BCMethodHolder PhysBCUtil::basicGradPressureFuncBC () const
{
    BCMethodHolder holder;

    int extrapOrder = 2;
    RefCountedPtr<BCGhostClass> BCPtr(
        new EllipticExtrapBCGhostClass(extrapOrder)
    );
    holder.addBCMethod(BCPtr);

    return holder;
}
