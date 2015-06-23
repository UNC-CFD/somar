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
#include "NewBeamGeneratorMap.H"
#include "ProblemContext.H"
#include "CubicSpline.H"
#include "HermiteInterp.H"
#include "BoxIterator.H"
#include "ConvertFABF_F.H"
// #include "MiscUtils.H"
#include "Subspace.H"
#include "LepticMeshRefine.H"
#include "Debug.H"
#include "LoadBalance.H"


// Static variables
bool          NewBeamGeneratorMap::s_staticParamsSet = false;

RealVect      NewBeamGeneratorMap::s_L;
IntVect       NewBeamGeneratorMap::s_lev0Nx;
RealVect      NewBeamGeneratorMap::s_lev0DXi;
ProblemDomain NewBeamGeneratorMap::s_lev0Domain;

RefCountedPtr<BoxLayoutData<NodeFArrayBox> > NewBeamGeneratorMap::s_lev0xPtr;
RefCountedPtr<BoxLayoutData<NodeFArrayBox> > NewBeamGeneratorMap::s_lev0d2xPtr;

Real NewBeamGeneratorMap::s_alpha;

const IntVect NewBeamGeneratorMap::s_ghostVect = 5 * IntVect::Unit;


// -----------------------------------------------------------------------------
// Construtor
// -----------------------------------------------------------------------------
NewBeamGeneratorMap::NewBeamGeneratorMap ()
{
    if (!s_staticParamsSet) {
        // Collect input parameters
        const ProblemContext* ctx = ProblemContext::getInstance();

        s_L = ctx->domainLength;
        s_lev0Nx = ctx->nx;
        s_lev0DXi = s_L / RealVect(s_lev0Nx);
        s_lev0Domain = ctx->domain;

        s_alpha = ctx->beamGenMapAlpha;

        // Cache the cubic spline data.
        BoxLayout lev0Layout = createLev0Grids();

        s_lev0xPtr = RefCountedPtr<BoxLayoutData<NodeFArrayBox> >(new BoxLayoutData<NodeFArrayBox>);
        s_lev0xPtr->define(lev0Layout, SpaceDim);

        s_lev0d2xPtr = RefCountedPtr<BoxLayoutData<NodeFArrayBox> >(new BoxLayoutData<NodeFArrayBox>);
        s_lev0d2xPtr->define(lev0Layout, SpaceDim*SpaceDim);

        // After this point, we only call computeSplineData, which is designed to use
        // only the static members above this point. Go ahead and set the flag.
        s_staticParamsSet = true;

        computeSplineData(this);

#       ifndef NDEBUG
        { // Mow much space does the cache require?
            unsigned long long bytes;

            bytes = s_lev0Domain.domainBox().numPts();
            bytes *= (SpaceDim + 1) * SpaceDim;
            bytes *= sizeof(Real);
            bytes /= (numProc() * 1024);
            pout() << bytes << "Kb should be used by this rank for the NewBeamGeneratorMap "
                   << "cache if load balancing went well." << endl;

            bytes = 0;
            DataIterator dit = lev0Layout.dataIterator();
            for (dit.reset(); dit.ok(); ++dit) {
                bytes += (*s_lev0xPtr)[dit].box().numPts() * s_lev0xPtr->nComp();
                bytes += (*s_lev0d2xPtr)[dit].box().numPts() * s_lev0d2xPtr->nComp();
            }
            bytes *= sizeof(Real);
            bytes /= 1024;
            pout() << bytes << "Kb are being used by this rank." << endl;

#           ifdef CH_MPI
            {
                unsigned long long maxBytes = (unsigned long long)(-1);
                MPI_Allreduce(&bytes, &maxBytes, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, Chombo_MPI::comm);
                bytes = maxBytes;
                pout() << bytes << "Kb is the max across all ranks." << endl;
            }
#           endif
            pout() << endl;
        }
#       endif
    }
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
NewBeamGeneratorMap::~NewBeamGeneratorMap ()
{;}


// -----------------------------------------------------------------------------
// Must return the name of the coordinate mapping
// -----------------------------------------------------------------------------
const char* NewBeamGeneratorMap::getCoorMapName () const
{
    return "NewBeamGeneratorMap";
}


// -----------------------------------------------------------------------------
// Must return whether or not this metric is diagonal
// -----------------------------------------------------------------------------
bool NewBeamGeneratorMap::isDiagonal () const
{
    return false;
}


// -----------------------------------------------------------------------------
// Must return whether or not this metric is uniform
// -----------------------------------------------------------------------------
bool NewBeamGeneratorMap::isUniform () const
{
    return false;
}


// -----------------------------------------------------------------------------
// Fills a mapped box with Cartesian locations.
// This new code tries to generate a grid on AMR level 0, then perform
// a piecewise linear interpolation of the nodal points. This is an attempt to
// reduce boundary face mismatches among different AMR levels.
// -----------------------------------------------------------------------------
void NewBeamGeneratorMap::fill_physCoor (FArrayBox&      a_dest,
                                         const int       a_destComp,
                                         const int       a_mu,
                                         const RealVect& a_dXi) const
{
    CH_TIME("NewBeamGeneratorMap::fill_physCoor");

    // Sanity checks
    CH_assert(0 <= a_destComp && a_destComp < a_dest.nComp());
    CH_assert(0 <= a_mu && a_mu < SpaceDim);
    CH_assert(a_dXi.product() > 0.0);

    if (a_dXi == s_lev0DXi) {
        // Compute the level 0 geometry
        fill_physCoorLev0(a_dest, a_destComp, a_mu);

    } else if (D_TERM(   (a_dXi[0] >= s_lev0DXi[0]),
                      && (a_dXi[1] >= s_lev0DXi[1]),
                      && (a_dXi[2] >= s_lev0DXi[2])   )) {
        // Coarsen the level 0 geometry
        const IntVect refRatio = computeRefToFiner(a_dXi, s_lev0DXi);
        const Box& destBox = a_dest.box();

        // if (destBox.type() == IntVect::Zero) {
        //     Box fineBox = a_dest.box();
        //     fineBox.refine(refRatio);

        //     const Box refBox(IntVect::Zero, refRatio-IntVect::Unit);

        //     FArrayBox fineFAB(fineBox, 1);
        //     fill_physCoorLev0(fineFAB, 0, a_mu);

        //     FArrayBox crseFAB(Interval(a_destComp, a_destComp), a_dest);

        //     FORT_UNMAPPEDAVERAGE(
        //         CHF_FRA(crseFAB),
        //         CHF_CONST_FRA(fineFAB),
        //         CHF_BOX(destBox),
        //         CHF_CONST_INTVECT(refRatio),
        //         CHF_BOX(refBox));

        // } else {
            MayDay::Error("Coarsen the level 0 geometry");
        // }

    } else if (D_TERM(   (a_dXi[0] <= s_lev0DXi[0]),
                      && (a_dXi[1] <= s_lev0DXi[1]),
                      && (a_dXi[2] <= s_lev0DXi[2])   )) {
        // Interpolate the level 0 geometry
        const IntVect refRatio = computeRefToFiner(s_lev0DXi, a_dXi);

        MayDay::Error("Interpolate the level 0 geometry");

    } else {
        MayDay::Error("NewBeamGeneratorMap cannot refine in one direction and coarsen in another");
    }
}


// // -----------------------------------------------------------------------------
// // Fills an FArrayBox with the Jacobian matrix elements d[x^mu] / d[xi^nu].
// // This is a speed bottleneck!!!
// // -----------------------------------------------------------------------------
// void NewBeamGeneratorMap::fill_dxdXi (FArrayBox&      a_dest,
//                                       const int       a_destComp,
//                                       const int       a_mu,
//                                       const int       a_nu,
//                                       const RealVect& a_dXi,
//                                       const Real      a_scale) const
// {
//     CH_TIME("NewBeamGeneratorMap::fill_dxdXi");

//     // Gather spline data
//     const Box& destBox = a_dest.box();
//     const IntVect& destBoxType = destBox.type();
//     const Box destNodeBox = surroundingNodes(destBox);

//     const Tuple<FArrayBox, CH_SPACEDIM> x0FAB;
//     for (int comp = 0; comp < SpaceDim; ++comp) {
//         FArrayBox& x0FABRef = (FArrayBox&)(x0FAB[comp]);
//         BoxLayoutData<NodeFArrayBox>& xCacheRef = (BoxLayoutData<NodeFArrayBox>&)(*s_lev0xPtr);
//         refOrCopy(x0FABRef, destNodeBox, xCacheRef, comp);
//     }

//     const Tuple<Tuple<FArrayBox, CH_SPACEDIM>, CH_SPACEDIM> d2x0FAB;
//     for (int coorComp = 0; coorComp < SpaceDim; ++coorComp) {
//         for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
//             FArrayBox& d2x0FABRef = (FArrayBox&)(d2x0FAB[coorComp][derivDir]);
//             BoxLayoutData<NodeFArrayBox>& d2xCacheRef = (BoxLayoutData<NodeFArrayBox>&)(*s_lev0d2xPtr);
//             refOrCopy(d2x0FABRef, destNodeBox, d2xCacheRef, derivDir+SpaceDim*coorComp);
//         }
//     }

//         // // Use spline data to compute first derivatives
//         // // (If this proved to be slow, we may want to cache this too.)
//         // Tuple<Tuple<FArrayBox, CH_SPACEDIM>, CH_SPACEDIM> dx0FAB;
//         // for (int coorComp = 0; coorComp < SpaceDim; ++coorComp) {
//         //     const FArrayBox& x0FABComp = x0FAB[coorComp];

//         //     for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
//         //         const Real dXi = s_lev0DXi[derivDir];
//         //         const Real invDXi = 1.0 / dXi;
//         //         const Real hoScale = dXi / 6.0;
//         //         const IntVect e = BASISV(derivDir);
//         //         const FArrayBox& d2x0FABComp = d2x0FAB[coorComp][derivDir];

//         //         FArrayBox& dx0FABComp = dx0FAB[coorComp][derivDir];
//         //         dx0FABComp.define(destNodeBox, 1);

//         //         Box hiStencilBox = destNodeBox;
//         //         hiStencilBox.growLo(derivDir, -1);

//         //         Box loStencilBox = destNodeBox;
//         //         loStencilBox.setBig(derivDir, loStencilBox.smallEnd(derivDir));

//         //         CH_assert(destNodeBox.contains(loStencilBox));
//         //         CH_assert(destNodeBox.contains(hiStencilBox));
//         //         CH_assert(destNodeBox.numPts() == loStencilBox.numPts() + hiStencilBox.numPts());

//         //         BoxIterator bit(hiStencilBox);
//         //         for (bit.reset(); bit.ok(); ++bit) {
//         //             const IntVect& hi = bit();
//         //             const IntVect  lo = hi - e;

//         //             Real df = x0FABComp(hi) - x0FABComp(lo);
//         //             dx0FABComp(hi) = df*invDXi
//         //                            + (d2x0FABComp(lo) + 2.0*d2x0FABComp(hi)) * hoScale;
//         //         }

//         //         bit.define(loStencilBox);
//         //         for (bit.reset(); bit.ok(); ++bit) {
//         //             const IntVect& lo = bit();
//         //             const IntVect  hi = lo + e;

//         //             Real df = x0FABComp(hi) - x0FABComp(lo);
//         //             dx0FABComp(lo) = df*invDXi
//         //                            - (2.0*d2x0FABComp(lo) + d2x0FABComp(hi)) * hoScale;
//         //         }
//         //     }
//         // }

//         // Perform interpolation
//         Real fA, fB, fC, fD;
//         Real fxA, fxB, fxC, fxD;
//         Real fyA, fyB, fyC, fyD;
//         Real fAB, fCD, fAC, fBD;

//         const Real u = 0.5;
//         const Real h4u = (u-1.0)*u*u;
//         const Real h3u = ((u-2.0)*u+1.0)*u;
//         const Real h2u = (3.0 - 2.0*u)*u*u;
//         const Real h1u = 1.0 - h2u;

//         const Real v = 0.5;
//         const Real h4v = (v-1.0)*v*v;
//         const Real h3v = ((v-2.0)*v+1.0)*v;
//         const Real h2v = (3.0 - 2.0*v)*v*v;
//         const Real h1v = 1.0 - h2v;

//         BoxIterator bit(destBox);
//         for (bit.reset(); bit.ok(); ++bit) {
//             const IntVect& cc = bit();
//             const IntVect& A = cc;
//             const IntVect  B = cc + BASISV(0);
//             const IntVect  C = cc + BASISV(1);
//             const IntVect  D = B  + BASISV(1);

//             fA = x0FAB[a_mu](A);
//             fB = x0FAB[a_mu](B);
//             fC = x0FAB[a_mu](C);
//             fD = x0FAB[a_mu](D);

//             // Use the Hermite interpolation along edges
//             // with first derivatives from the cubic splines...
//             // fxA = dx0FAB[a_mu][0](A); // These derivs will not be needed if we use
//             // fxB = dx0FAB[a_mu][0](B); // splines for fAB, etc.
//             // fxC = dx0FAB[a_mu][0](C);
//             // fxD = dx0FAB[a_mu][0](D);

//             // fyA = dx0FAB[a_mu][1](A);
//             // fyB = dx0FAB[a_mu][1](B);
//             // fyC = dx0FAB[a_mu][1](C);
//             // fyD = dx0FAB[a_mu][1](D);

//             // // Compute edge interpolations
//             // fAB = fA*h1u + fB*h2u + fxA*h3u + fxB*h4u;
//             // fCD = fC*h1u + fD*h2u + fxC*h3u + fxD*h4u;
//             // fAC = fA*h1v + fC*h2v + fyA*h3v + fyC*h4v;
//             // fBD = fB*h1v + fD*h2v + fyB*h3v + fyD*h4v;


//             // Use splines along edges...
//             // fxA, etc are second derivatives.
//             fxA = d2x0FAB[a_mu][0](A);
//             fxB = d2x0FAB[a_mu][0](B);
//             fxC = d2x0FAB[a_mu][0](C);
//             fxD = d2x0FAB[a_mu][0](D);

//             fyA = d2x0FAB[a_mu][1](A);
//             fyB = d2x0FAB[a_mu][1](B);
//             fyC = d2x0FAB[a_mu][1](C);
//             fyD = d2x0FAB[a_mu][1](D);

//             fAB = 0.5 * (fA + fB) - (fxA + fxB) * s_lev0DXi[0]*s_lev0DXi[0]/16.0;
//             fCD = 0.5 * (fC + fD) - (fxC + fxD) * s_lev0DXi[0]*s_lev0DXi[0]/16.0;
//             fAC = 0.5 * (fA + fC) - (fyA + fyC) * s_lev0DXi[1]*s_lev0DXi[1]/16.0;
//             fBD = 0.5 * (fB + fD) - (fyB + fyD) * s_lev0DXi[1]*s_lev0DXi[1]/16.0;

//             // Compute simplified Hermite interpolant
//             a_dest(cc, a_destComp) = fAB*h1v + fCD*h2v + fAC*h1u + fBD*h2u
//                                    - fA*h1u*h1v - fB*h2u*h1v - fC*h1u*h2v - fD*h2u*h2v;
//         }
// }


// -----------------------------------------------------------------------------
// Fills a mapped box with Cartesian locations.
// Assumes we are on level 0.
// -----------------------------------------------------------------------------
void NewBeamGeneratorMap::fill_physCoorLev0 (FArrayBox& a_dest,
                                             const int  a_destComp,
                                             const int  a_mu) const
{
    CH_TIMERS("NewBeamGeneratorMap::fill_physCoorLev0");
    CH_TIMER("NewBeamGeneratorMap::fill_physCoorLev0...NC", tmrNC);
    CH_TIMER("NewBeamGeneratorMap::fill_physCoorLev0...EC", tmrEC); // EC = FC in 2D.
    CH_TIMER("NewBeamGeneratorMap::fill_physCoorLev0...FC", tmrFC); // Only in 3D.
    CH_TIMER("NewBeamGeneratorMap::fill_physCoorLev0...CC", tmrCC);

    // Sanity checks
    CH_assert(0 <= a_destComp && a_destComp < a_dest.nComp());
    CH_assert(0 <= a_mu && a_mu < SpaceDim);

    // If a_dest is nodal, just fill directly. Otherwise, interpolate carefully.
    const Box& destBox = a_dest.box();
    const IntVect& destBoxType = destBox.type();

    if (destBoxType == IntVect::Unit) {
        // a_dest is nodal.
        CH_START(tmrNC);
        copyLev0NodeCacheData(a_dest, a_destComp, *s_lev0xPtr, a_mu);
        CH_STOP(tmrNC);

    } else if (destBoxType.sum() == 1) {
        // a_dest is an edge-centered box. Just use a cubic spline.
        CH_START(tmrEC);

        // Find the centered (interpolated) direction.
        D_TERM(
        int interpDir = 0;,
        if (destBoxType[1] == 0) interpDir = 1;,
        if (destBoxType[2] == 0) interpDir = 2;)

        // Gather spline data
        const Box destNodeBox = surroundingNodes(destBox);
        const FArrayBox x0FAB, d2x0FAB;
        {
            FArrayBox& x0FABRef = (FArrayBox&)x0FAB;
            BoxLayoutData<NodeFArrayBox>& xCacheRef = (BoxLayoutData<NodeFArrayBox>&)(*s_lev0xPtr);
            refOrCopy(x0FABRef, destNodeBox, xCacheRef, a_mu);

            FArrayBox& d2x0FABRef = (FArrayBox&)d2x0FAB;
            BoxLayoutData<NodeFArrayBox>& d2xCacheRef = (BoxLayoutData<NodeFArrayBox>&)(*s_lev0d2xPtr);
            refOrCopy(d2x0FABRef, destNodeBox, d2xCacheRef, interpDir+SpaceDim*a_mu);
        }

        // Perform the interpolation
        const Real dXi = s_lev0DXi[interpDir];
        const Real halfDXi = 0.5 * dXi;
        const Real invDXi = 1.0 / dXi;
        const Real cubicScale = dXi * dXi / 6.0;
        const IntVect e = BASISV(interpDir);
        Real Xil, Xic, Xir, xpred, A, B;

        BoxIterator bit(destBox);
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& fc = bit();
            const IntVect& ncl = fc;
            const IntVect  ncr = ncl + e;

            Xil = Real(ncl[interpDir]) * dXi;
            Xic = Xil + halfDXi;
            Xir = Xil + dXi;

            // Predict value using linear interpolant.
            A = (Xir - Xic) * invDXi;
            B = (Xic - Xil) * invDXi;
            xpred = A*x0FAB(ncl) + B*x0FAB(ncr);

            // Apply cubic correction.
            A = A*(A*A-1.0);
            B = B*(B*B-1.0);
            a_dest(fc, a_destComp) = xpred
                    + (A*d2x0FAB(ncl) + B*d2x0FAB(ncr)) * cubicScale;
        }

        CH_STOP(tmrEC);

    } else if (destBoxType == IntVect::Zero) {
        // a_dest is CC.
        CH_START(tmrCC);

        // Gather spline data
        // (Hopefully, speeding up the copy function will help here.)
        const Box destNodeBox = surroundingNodes(destBox);

        const Tuple<FArrayBox, CH_SPACEDIM> x0FAB;
        for (int comp = 0; comp < SpaceDim; ++comp) {
            // x0FAB[comp].define(destNodeBox, 1);
            // copyLev0NodeCacheData(x0FAB[comp], 0, *s_lev0xPtr, comp);

            FArrayBox& x0FABRef = (FArrayBox&)(x0FAB[comp]);
            BoxLayoutData<NodeFArrayBox>& xCacheRef = (BoxLayoutData<NodeFArrayBox>&)(*s_lev0xPtr);
            refOrCopy(x0FABRef, destNodeBox, xCacheRef, comp);
        }

        const Tuple<Tuple<FArrayBox, CH_SPACEDIM>, CH_SPACEDIM> d2x0FAB;
        for (int coorComp = 0; coorComp < SpaceDim; ++coorComp) {
            for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
                // d2x0FAB[coorComp][derivDir].define(destNodeBox, 1);

                // const int tensorComp = derivDir + SpaceDim*coorComp;
                // copyLev0NodeCacheData(d2x0FAB[coorComp][derivDir], 0, *s_lev0d2xPtr, tensorComp);

                FArrayBox& d2x0FABRef = (FArrayBox&)(d2x0FAB[coorComp][derivDir]);
                BoxLayoutData<NodeFArrayBox>& d2xCacheRef = (BoxLayoutData<NodeFArrayBox>&)(*s_lev0d2xPtr);
                refOrCopy(d2x0FABRef, destNodeBox, d2xCacheRef, derivDir+SpaceDim*coorComp);
            }
        }

        // // Use spline data to compute first derivatives
        // // (If this proved to be slow, we may want to cache this too.)
        // Tuple<Tuple<FArrayBox, CH_SPACEDIM>, CH_SPACEDIM> dx0FAB;
        // for (int coorComp = 0; coorComp < SpaceDim; ++coorComp) {
        //     const FArrayBox& x0FABComp = x0FAB[coorComp];

        //     for (int derivDir = 0; derivDir < SpaceDim; ++derivDir) {
        //         const Real dXi = s_lev0DXi[derivDir];
        //         const Real invDXi = 1.0 / dXi;
        //         const Real hoScale = dXi / 6.0;
        //         const IntVect e = BASISV(derivDir);
        //         const FArrayBox& d2x0FABComp = d2x0FAB[coorComp][derivDir];

        //         FArrayBox& dx0FABComp = dx0FAB[coorComp][derivDir];
        //         dx0FABComp.define(destNodeBox, 1);

        //         Box hiStencilBox = destNodeBox;
        //         hiStencilBox.growLo(derivDir, -1);

        //         Box loStencilBox = destNodeBox;
        //         loStencilBox.setBig(derivDir, loStencilBox.smallEnd(derivDir));

        //         CH_assert(destNodeBox.contains(loStencilBox));
        //         CH_assert(destNodeBox.contains(hiStencilBox));
        //         CH_assert(destNodeBox.numPts() == loStencilBox.numPts() + hiStencilBox.numPts());

        //         BoxIterator bit(hiStencilBox);
        //         for (bit.reset(); bit.ok(); ++bit) {
        //             const IntVect& hi = bit();
        //             const IntVect  lo = hi - e;

        //             Real df = x0FABComp(hi) - x0FABComp(lo);
        //             dx0FABComp(hi) = df*invDXi
        //                            + (d2x0FABComp(lo) + 2.0*d2x0FABComp(hi)) * hoScale;
        //         }

        //         bit.define(loStencilBox);
        //         for (bit.reset(); bit.ok(); ++bit) {
        //             const IntVect& lo = bit();
        //             const IntVect  hi = lo + e;

        //             Real df = x0FABComp(hi) - x0FABComp(lo);
        //             dx0FABComp(lo) = df*invDXi
        //                            - (2.0*d2x0FABComp(lo) + d2x0FABComp(hi)) * hoScale;
        //         }
        //     }
        // }

        // Perform interpolation
        Real fA, fB, fC, fD;
        Real fxA, fxB, fxC, fxD;
        Real fyA, fyB, fyC, fyD;
        Real fAB, fCD, fAC, fBD;

        const Real u = 0.5;
        const Real h4u = (u-1.0)*u*u;
        const Real h3u = ((u-2.0)*u+1.0)*u;
        const Real h2u = (3.0 - 2.0*u)*u*u;
        const Real h1u = 1.0 - h2u;

        const Real v = 0.5;
        const Real h4v = (v-1.0)*v*v;
        const Real h3v = ((v-2.0)*v+1.0)*v;
        const Real h2v = (3.0 - 2.0*v)*v*v;
        const Real h1v = 1.0 - h2v;

        BoxIterator bit(destBox);
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& cc = bit();
            const IntVect& A = cc;
            const IntVect  B = cc + BASISV(0);
            const IntVect  C = cc + BASISV(1);
            const IntVect  D = B  + BASISV(1);

            fA = x0FAB[a_mu](A);
            fB = x0FAB[a_mu](B);
            fC = x0FAB[a_mu](C);
            fD = x0FAB[a_mu](D);

            // Use the Hermite interpolation along edges
            // with first derivatives from the cubic splines...
            // fxA = dx0FAB[a_mu][0](A); // These derivs will not be needed if we use
            // fxB = dx0FAB[a_mu][0](B); // splines for fAB, etc.
            // fxC = dx0FAB[a_mu][0](C);
            // fxD = dx0FAB[a_mu][0](D);

            // fyA = dx0FAB[a_mu][1](A);
            // fyB = dx0FAB[a_mu][1](B);
            // fyC = dx0FAB[a_mu][1](C);
            // fyD = dx0FAB[a_mu][1](D);

            // // Compute edge interpolations
            // fAB = fA*h1u + fB*h2u + fxA*h3u + fxB*h4u;
            // fCD = fC*h1u + fD*h2u + fxC*h3u + fxD*h4u;
            // fAC = fA*h1v + fC*h2v + fyA*h3v + fyC*h4v;
            // fBD = fB*h1v + fD*h2v + fyB*h3v + fyD*h4v;


            // Use splines along edges...
            // fxA, etc are second derivatives.
            fxA = d2x0FAB[a_mu][0](A);
            fxB = d2x0FAB[a_mu][0](B);
            fxC = d2x0FAB[a_mu][0](C);
            fxD = d2x0FAB[a_mu][0](D);

            fyA = d2x0FAB[a_mu][1](A);
            fyB = d2x0FAB[a_mu][1](B);
            fyC = d2x0FAB[a_mu][1](C);
            fyD = d2x0FAB[a_mu][1](D);

            fAB = 0.5 * (fA + fB) - (fxA + fxB) * s_lev0DXi[0]*s_lev0DXi[0]/16.0;
            fCD = 0.5 * (fC + fD) - (fxC + fxD) * s_lev0DXi[0]*s_lev0DXi[0]/16.0;
            fAC = 0.5 * (fA + fC) - (fyA + fyC) * s_lev0DXi[1]*s_lev0DXi[1]/16.0;
            fBD = 0.5 * (fB + fD) - (fyB + fyD) * s_lev0DXi[1]*s_lev0DXi[1]/16.0;

            // Compute simplified Hermite interpolant
            a_dest(cc, a_destComp) = fAB*h1v + fCD*h2v + fAC*h1u + fBD*h2u
                                   - fA*h1u*h1v - fB*h2u*h1v - fC*h1u*h2v - fD*h2u*h2v;
        }

        CH_STOP(tmrCC);

    } else {
        // a_dest is face-centered in 3D.
        CH_START(tmrFC);

        // TODO
        MayDay::Warning("fill_physCoorLev0 (EC) is using temp code");

        // Collect position data on nodes.
        const Box nodeBox = surroundingNodes(destBox);
        FArrayBox nodeDest(nodeBox, 1);
        copyLev0NodeCacheData(nodeDest, 0, *s_lev0xPtr, a_mu);

        // Collect interpolation data line-by-line so we don't fill memory.
        Box flatBox = flattenBox(nodeBox, a_mu);
        // Compute cubic spline of nodal data.


        // Interpolate to desired centering.
        FORT_CONVERTFAB (
            CHF_FRA1(a_dest, a_destComp),
            CHF_BOX(destBox),
            CHF_CONST_INTVECT(destBoxType),
            CHF_CONST_FRA1(nodeDest, 0),
            CHF_CONST_INTVECT(IntVect::Unit));

        CH_STOP(tmrFC);
    }
}


// -----------------------------------------------------------------------------
// Fills a mapped box with Cartesian locations.
// Assumes we are on level 0 and a_dest is NC.
// Stretching happens here.
// -----------------------------------------------------------------------------
void NewBeamGeneratorMap::fill_physCoorLev0Nodes (FArrayBox& a_dest,
                                                  const int  a_destComp,
                                                  const int  a_mu) const
{
    CH_TIME("NewBeamGeneratorMap::fill_physCoorLev0Nodes");

    // Sanity checks
    CH_assert(s_staticParamsSet);
    CH_assert(a_dest.box().type() == IntVect::Unit);
    CH_assert(0 <= a_destComp && a_destComp < a_dest.nComp());
    CH_assert(0 <= a_mu && a_mu < SpaceDim);

    const Real vertStretch = 1.7;
    const Real horizStretch = 1.7;

    if (a_mu < SpaceDim-1) {
        const Real horizMult = 1.0 / tanh(horizStretch);

        Real x, L, dx;
        L = s_L[a_mu];
        dx = s_lev0DXi[a_mu];

        BoxIterator bit(a_dest.box());
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& nc = bit();

            x = Real(nc[0]) * dx;
            x = 0.5*L * (1.0 - horizMult*tanh(horizStretch*(1.0 - 2.0*abs(x)/L)));
            x = ((nc[0] >= 0)? x: -x);

            a_dest(nc, a_destComp) = x;
        }

    } else {
        const Real vertMult = 1.0 / tanh(vertStretch);

        // Compute x,y at bottom, then elevation of bottom.
        const Box& destBox = a_dest.box();
        const Box bottomNodeBox = flattenBox(destBox, SpaceDim-1);
        const IntVect hmask = IntVect::Unit - BASISV(SpaceDim-1);

        FArrayBox cartPos(bottomNodeBox, SpaceDim);
        D_TERM(,
        fill_physCoorLev0Nodes(cartPos, 0, 0);,
        fill_physCoorLev0Nodes(cartPos, 1, 1);)

        fill_elevationLev0(cartPos, SpaceDim-1, cartPos);

        // Construct vertical grid from bottom data.
        Real z, H, dz, dHat, zetaHat, zHat;
        H = s_L[SpaceDim-1];
        dz = s_lev0DXi[SpaceDim-1];

        BoxIterator bit(destBox);
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& nc = bit();

            dHat = cartPos(nc*hmask, SpaceDim-1) / H;
            zetaHat = Real(nc[SpaceDim-1]) * dz / H;
            zetaHat = (1.0 - vertMult*tanh(vertStretch*(1.0 - zetaHat)));
            zHat = (1.0-zetaHat)*dHat + zetaHat;

            a_dest(nc, a_destComp) = zHat * H;
        }
    }
}


#include "BeamGeneratorMapF_F.H"
// -----------------------------------------------------------------------------
// Returns the topographic elevations at each Cartesian node.
// Assumes we are on level 0 and a_dest is NC.
// -----------------------------------------------------------------------------
void NewBeamGeneratorMap::fill_elevationLev0 (FArrayBox&       a_dest,
                                              const int        a_destComp,
                                              const FArrayBox& a_physCoor) const
{
    // Sanity checks
    CH_assert(a_dest.box().type() == IntVect::Unit);
    CH_assert(a_physCoor.box().type() == IntVect::Unit);
    CH_assert(a_dest.box().size(SpaceDim-1) == 1);
    CH_assert(a_physCoor.box().contains(a_dest.box()));
    CH_assert(a_physCoor.nComp() >= SpaceDim-1);
    CH_assert(0 <= a_destComp);
    CH_assert(a_destComp < a_dest.nComp());

    const Box& destBox = a_dest.box();
    const IntVect& destBoxType = destBox.type();

    FORT_FILL_BEAMGENERATORMAPBATHYMETRY (
        CHF_FRA1(a_dest, a_destComp),
        CHF_BOX(destBox),
        CHF_CONST_INTVECT(destBoxType),
        CHF_CONST_FRA(a_physCoor),
        CHF_CONST_REALVECT(s_lev0DXi),
        CHF_CONST_REALVECT(s_L),
        CHF_CONST_REAL(s_alpha));
}


// -----------------------------------------------------------------------------
// Static utility.
// Computes a refRatio. This will also check if crseDXi >= fineDXi.
// -----------------------------------------------------------------------------
IntVect NewBeamGeneratorMap::computeRefToFiner (const RealVect& a_crseDXi,
                                                const RealVect& a_fineDXi)
{
    IntVect refRatio(D_DECL(
        int(a_crseDXi[0] / a_fineDXi[0]),
        int(a_crseDXi[1] / a_fineDXi[1]),
        int(a_crseDXi[2] / a_fineDXi[2])
    ));

    // Make sure the user didn't switch crse and fine DXi. If they did, refRatio
    // will have been defined with fractions that were truncated to zero.
    D_TERM(
        CH_assert(a_crseDXi[0] == Real(refRatio[0]) * a_fineDXi[0]);,
        CH_assert(a_crseDXi[1] == Real(refRatio[1]) * a_fineDXi[1]);,
        CH_assert(a_crseDXi[2] == Real(refRatio[2]) * a_fineDXi[2]);
    )

    return refRatio;
}


// -----------------------------------------------------------------------------
// Static utility
// Creates lev 0 grids. Does not use a layout suggestion.
// -----------------------------------------------------------------------------
BoxLayout NewBeamGeneratorMap::createLev0Grids ()
{
    // Do not rely on other static member being set.
    const ProblemContext* ctx = ProblemContext::getInstance();
    const ProblemDomain& domain = ctx->domain;
    const IntVect& maxBoxSize = ctx->maxBaseGridSize;
    const int blockFactor = ctx->block_factor;

    Vector<Box> vbox;
    LepticMeshRefine::domainSplit(domain, vbox, maxBoxSize, blockFactor);

    Vector<int> vproc;
    LoadBalance(vproc, vbox);

    BoxLayout layout;
    layout.define(vbox, vproc);
    layout.grow(s_ghostVect);

    return layout;
}


// -----------------------------------------------------------------------------
// This can be used to speed up cache access.
// -----------------------------------------------------------------------------
void NewBeamGeneratorMap::suggestLev0Grids (const DisjointBoxLayout& a_grids)
{
    BoxLayout newLayout;
    newLayout.deepCopy(a_grids);
    newLayout.grow(s_ghostVect);
    newLayout.close();
    DataIterator dit = newLayout.dataIterator();

    // Is this actually a new layout? If now, don't waste time.
    if (newLayout.compatible(s_lev0xPtr->boxLayout())) return;

    // TODO: Use generalCopyTo directly on all comps at once for speed.

    { // Convert s_lev0xPtr
        RefCountedPtr<BoxLayoutData<NodeFArrayBox> > newPtr(new BoxLayoutData<NodeFArrayBox>);
        newPtr->define(newLayout, s_lev0xPtr->nComp());

        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& newFAB = (*newPtr)[dit].getFab();
            for (int comp = 0; comp < newPtr->nComp(); ++comp) {
                copyLev0NodeCacheData(newFAB, comp, *s_lev0xPtr, comp);
            }
        }
        s_lev0xPtr = newPtr;
    }

    { // Convert s_lev0d2xPtr
        RefCountedPtr<BoxLayoutData<NodeFArrayBox> > newPtr(new BoxLayoutData<NodeFArrayBox>);
        newPtr->define(newLayout, s_lev0d2xPtr->nComp());

        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& newFAB = (*newPtr)[dit].getFab();
            for (int comp = 0; comp < newPtr->nComp(); ++comp) {
                copyLev0NodeCacheData(newFAB, comp, *s_lev0d2xPtr, comp);
            }
        }
        s_lev0d2xPtr = newPtr;
    }
}


// -----------------------------------------------------------------------------
// Static utility
// Fills the cubic spline cache.
// -----------------------------------------------------------------------------
void NewBeamGeneratorMap::computeSplineData (const NewBeamGeneratorMap* a_srcPtr)
{
    CH_TIME("NewBeamGeneratorMap::computeSplineData");

    // These need to be well defined prior to call.
    CH_assert(!s_lev0xPtr.isNull());
    CH_assert(!s_lev0d2xPtr.isNull());

    const ProblemContext* ctx = ProblemContext::getInstance();
    const RealVect& dXi = ctx->domainLength / RealVect(ctx->nx);
    const ProblemDomain& domain = ctx->domain;

    DataIterator dit = s_lev0xPtr->dataIterator();

    // This contains all nodes in the domain + ghost nodes.
    // Solves are performed on coordinate lines through these nodes.
    // The ghosts MUST be included in the solves!
    Box allNodeBox = domain.domainBox();
    allNodeBox.surroundingNodes();
    allNodeBox.grow(s_ghostVect);

    // We will compute spline data line-by-line.
    // Loop over each coordinate line.
    for (int interpDir = 0; interpDir < SpaceDim; ++interpDir) {
        const int numNodes = allNodeBox.size(interpDir);

        // These nodes serve as the starting points for each line solve.
        const Box planeNodeBox = bdryLo(allNodeBox, interpDir, 1);
        CH_assert(!planeNodeBox.isEmpty());

        BoxIterator sbit(planeNodeBox);
        for (sbit.reset(); sbit.ok(); ++sbit) {
            const IntVect& startIV = sbit();

            // This is the region of the line solve.
            const Box lineBox(startIV,
                              startIV + (numNodes-1)*BASISV(interpDir),
                              IntVect::Unit);

            // Fill Xi along the line.
            // NOTE: fill_physCoorLev0Nodes will need ALL comps!
            Vector<Vector<Real> > Xi(SpaceDim, Vector<Real>(numNodes));
            IntVect offset = IntVect::Zero;
            for (int idx = 0; idx < numNodes; ++idx) {
                offset[interpDir] = idx;
                D_TERM(Xi[0][idx] = Real(startIV[0] + offset[0]) * dXi[0];,
                       Xi[1][idx] = Real(startIV[1] + offset[1]) * dXi[1];,
                       Xi[2][idx] = Real(startIV[2] + offset[2]) * dXi[2];)
            }

            // Loop over coordinate functions, x^i(\vec{Xi}).
            for (int coorDir = 0; coorDir < SpaceDim; ++coorDir) {
                // Fill x along the line
                Vector<Real> x(numNodes);
                {
                    FArrayBox xFAB(lineBox, 1);
                    a_srcPtr->fill_physCoorLev0Nodes(xFAB, 0, coorDir);
                    IntVect nc = xFAB.smallEnd();
                    for (int i = 0; i < numNodes; ++i) {
                        x[i] = xFAB(nc);
                        ++nc[interpDir];
                    }
                }

                // Solve the spline problem. x^i's dependence on Xi in
                // transverse directions is implicit.
                CubicSpline spline;
                spline.solve(x, Xi[interpDir]);

                // Compute first derivatives
                Vector<Real> d2x(numNodes);
                spline.interpSecondDeriv(d2x, Xi[interpDir]);

                // Copy x^{coorDir} and d^2x^{coorDir}/dXi^{interpDir}^2
                // to the cache.
                // NOTE: x_{coorDir} will be sent to the cache SpaceDim times!
                for (dit.reset(); dit.ok(); ++dit) {
                    FArrayBox& xFAB = (*s_lev0xPtr)[dit].getFab();
                    FArrayBox& d2xFAB = (*s_lev0d2xPtr)[dit].getFab();
                    const int tensorComp = interpDir + SpaceDim*coorDir;

                    Box copyRegion = xFAB.box();
                    copyRegion &= lineBox;

                    BoxIterator cbit(copyRegion);
                    for (cbit.reset(); cbit.ok(); ++cbit) {
                        const IntVect& nc = cbit();
                        int idx = nc[interpDir] - startIV[interpDir];

                        // xFAB(nc, coorDir) = Real(nc[coorDir]) * dXi[coorDir];
                        xFAB(nc, coorDir) = x[idx];
                        d2xFAB(nc, tensorComp) = d2x[idx];
                    } // end loop over line for copy (cbit)
                } // end loop over grids for copy (dit)

            } // end loop over coordinate functions (coorDir)
        } // end loop over line solve starting points (sbit)
    } // end loop over line solve directions (interpDir)
}


// -----------------------------------------------------------------------------
// Static utility
// Copies the Cartesian locations from the cache.
// Assumes we are on level 0 and a_dest is NC.
// -----------------------------------------------------------------------------
void NewBeamGeneratorMap::copyLev0NodeCacheData (FArrayBox&                          a_dest,
                                                 const int                           a_destComp,
                                                 const BoxLayoutData<NodeFArrayBox>& a_srcData,
                                                 const int                           a_srcComp)
{
    CH_TIME("NewBeamGeneratorMap::copy_physCoorLev0Nodes");

    // Sanity checks
    CH_assert(a_dest.box().type() == IntVect::Unit);
    CH_assert(0 <= a_destComp && a_destComp < a_dest.nComp());
    CH_assert(0 <= a_srcComp && a_srcComp < a_srcData.nComp());

    // Gather needed structures
    const ProblemContext* ctx = ProblemContext::getInstance();
    const ProblemDomain& domain = ctx->domain;

    const int np = numProc();
    const int srcProc = uniqueProc(SerialTask::compute);
    const int thisProc = procID();

    const Box& destBox = a_dest.box();
    CH_assert(!destBox.isEmpty());

    // Do we need MPI?
    {
        DataIterator dit = a_srcData.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            const FArrayBox& srcFAB = a_srcData[dit].getFab();

            if (srcFAB.contains(destBox)) {
                // Just copy directly and scram.
                a_dest.copy(srcFAB, a_srcComp, a_destComp, 1);
                return;
            }
        }
    }

    // Construct destination layout for MPI copy.
    Box destCCBox = destBox;
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (destCCBox.size(dir) == 1) {
            destCCBox.grow(dir, 1);
        }
    }
    destCCBox.enclosedCells();
    CH_assert(!destCCBox.isEmpty());

    Vector<Box> vbox(np);
    gather(vbox, destCCBox, srcProc);
    broadcast(vbox, srcProc);

    Vector<int> thisPair(2, -1);
    for (int idx = 0; idx < np; ++idx) {
        if (vbox[idx] != destCCBox) continue;

        thisPair[0] = idx;
        thisPair[1] = thisProc;
        break;
    }
    CH_assert(thisPair[0] >= 0);
    CH_assert(thisPair[1] >= 0);

    Vector<Vector<int> > vpair(np, Vector<int>(2, -1));
    gather(vpair, thisPair, srcProc);

    Vector<int> vproc(np);
    if (thisProc == srcProc) {
        for (int idx = 0; idx < np; ++idx) {
            int sidx;
            for (sidx = 0; sidx < np; ++sidx) {
                if (vpair[sidx][0] == idx) break;
            }
            CH_assert(sidx < np);

            vproc[idx] = vpair[sidx][1];
        }
    }
    broadcast(vproc, srcProc);

    BoxLayout destLayout;
    destLayout.define(vbox, vproc);
    DataIterator dit = destLayout.dataIterator();

    // Copy to temp holder / get MPI copy out of the way.
    LayoutData<Vector<RefCountedPtr<NodeFArrayBox> > > destLDPtrs;
    a_srcData.generalCopyTo(destLayout, destLDPtrs, Interval(a_srcComp, a_srcComp), domain);

    // Send from temp holder to caller's holder.
    for (dit.reset(); dit.ok(); ++dit) {
        const Vector<RefCountedPtr<NodeFArrayBox> >& destNFABPtrs = destLDPtrs[dit];
        for (int idx = 0; idx < destNFABPtrs.size(); ++idx) {
            a_dest.copy(destNFABPtrs[idx]->getFab(), 0, a_destComp, 1);
        }
    }
}


// -----------------------------------------------------------------------------
// Static utility.
// If the cache data is local, we reference it to avoid MPI+copy. Returns
// true if a_dest is a reference, false if we needed to perform the MPI+copy.
// a_dest must enter this function undefined.
// -----------------------------------------------------------------------------
bool NewBeamGeneratorMap::refOrCopy (FArrayBox&                    a_dest,
                                     const Box&                    a_destNodeBox,
                                     BoxLayoutData<NodeFArrayBox>& a_srcData,
                                     const int                     a_srcComp)
{
    CH_assert(a_dest.box().isEmpty());

    // Can we avoid MPI?
    DataIterator dit = a_srcData.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& srcFAB = a_srcData[dit].getFab();
        if (srcFAB.contains(a_destNodeBox)) {
            Interval interv(a_srcComp, a_srcComp);
            a_dest.define(interv, srcFAB);
            return true;
        }
    }

    // Data is not entirely on this rank. Do MPI copy.
    a_dest.define(a_destNodeBox, 1);
    copyLev0NodeCacheData(a_dest, 0, a_srcData, a_srcComp);
    return false;
}
