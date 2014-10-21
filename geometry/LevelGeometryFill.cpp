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
#include "LevelGeometry.H"
#include "LevelGeometryF_F.H"
#include "MappedCoarseAverage.H"


// -----------------------------------------------------------------------------
// Fills a FAB with displacements from Xi to physical locations.
// For use with VisIt's displace operator.
// -----------------------------------------------------------------------------
void LevelGeometry::fill_displacement (FArrayBox& a_dest) const
{
    CH_TIME("LevelGeometry::fill_displacement");

    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_dest.nComp() == SpaceDim);

    // Fill data holder
    this->fill_physCoor(a_dest);

    // Figure out the centering
    const Box destBox(a_dest.box());
    const IntVect boxType(destBox.type());

    // Calculate displacements
    FORT_FILLMAPPINGDISP(
        CHF_FRA(a_dest),
        CHF_CONST_FRA(a_dest),
        CHF_CONST_REALVECT(m_dXi),
        CHF_BOX(destBox),
        CHF_CONST_INTVECT(boxType));
}


// -----------------------------------------------------------------------------
// Fills a mapped box with physical locations
// a_dest must have SpaceDim comps.
// -----------------------------------------------------------------------------
void LevelGeometry::fill_physCoor (FArrayBox& a_dest) const
{
    CH_TIME("LevelGeometry::fill_physCoor");

    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_dest.nComp() == SpaceDim);

    // Fill data holder
    s_metricSourcePtr->fill_physCoor(a_dest, m_dXi);
}


// -----------------------------------------------------------------------------
// Fill a FAB with the Jacobian matrix (mapped to Cartesian)
// J^a_b = dx^a / dXi^b
// a_dest must have SpaceDim^2 comps.
// -----------------------------------------------------------------------------
void LevelGeometry::fill_dxdXi (FArrayBox& a_dest) const
{
    CH_TIME("LevelGeometry::fill_dxdXi");

    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_dest.nComp() == SpaceDim*SpaceDim);

    // Fill data holder
    for (int bdir = 0; bdir < SpaceDim; ++bdir) {
        for (int adir = 0; adir < SpaceDim; ++adir) {
            const int comp = this->tensorCompCC(adir, bdir);
            s_metricSourcePtr->fill_dxdXi(a_dest,
                                          comp,
                                          adir, bdir,
                                          m_dXi);
        }
    }
}


// -----------------------------------------------------------------------------
// Fill a FAB with the inverse Jacobian matrix (Cartesian to mapped)
// Jinv^a_b = dXi^a / dx^b
// a_dest must have SpaceDim^2 comps.
// -----------------------------------------------------------------------------
void LevelGeometry::fill_dXidx (FArrayBox& a_dest) const
{
    CH_TIME("LevelGeometry::fill_dXidx");

    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_dest.nComp() == SpaceDim*SpaceDim);

    // Fill data holder
    for (int bdir = 0; bdir < SpaceDim; ++bdir) {
        for (int adir = 0; adir < SpaceDim; ++adir) {
            const int comp = this->tensorCompCC(adir, bdir);

            s_metricSourcePtr->fill_dXidx(a_dest,
                                          comp,
                                          adir, bdir,
                                          m_dXi);
        }
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with J
// -----------------------------------------------------------------------------
void LevelGeometry::fill_J (FArrayBox& a_dest) const
{
    CH_TIME("LevelGeometry::fill_J");

    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_dest.nComp() == 1);

    // Fill data holder
    s_metricSourcePtr->fill_J(a_dest, 0, m_dXi);
}


// -----------------------------------------------------------------------------
// Fills a FluxBox with J
// -----------------------------------------------------------------------------
void LevelGeometry::fill_J (FluxBox& a_dest) const
{
    CH_TIME("LevelGeometry::fill_J");

    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_dest.nComp() == 1);

    // Fill data holder
    for (int adir = 0; adir < SpaceDim; ++adir) {
        s_metricSourcePtr->fill_J(a_dest[adir], 0, m_dXi);
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with 1/J
// -----------------------------------------------------------------------------
void LevelGeometry::fill_Jinv (FArrayBox& a_dest) const
{
    CH_TIME("LevelGeometry::fill_Jinv");

    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_dest.nComp() == 1);

    // Fill data holder
    s_metricSourcePtr->fill_Jinv(a_dest, 0, m_dXi);
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the covariant metric elements (static version)
// -----------------------------------------------------------------------------
void LevelGeometry::fill_gdn (FArrayBox&      a_dest,
                              const RealVect& a_dXi)
{
    CH_TIME("LevelGeometry::fill_gdn (static version)");

    // Sanity checks
    CH_assert(a_dest.nComp() == (SpaceDim * (SpaceDim + 1)) / 2);

    // Fill data holder
    for (int adir = 0; adir < SpaceDim; ++adir) {
        for (int bdir = adir; bdir < SpaceDim; ++bdir) {
            const int comp = symTensorCompCC(adir, bdir);

            s_metricSourcePtr->fill_gdn(a_dest,
                                        comp,
                                        adir, bdir,
                                        a_dXi);
        }
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the covariant metric elements
// -----------------------------------------------------------------------------
void LevelGeometry::fill_gdn (FArrayBox& a_dest) const
{
    CH_TIME("LevelGeometry::fill_gdn");

    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_dest.nComp() == (SpaceDim * (SpaceDim + 1)) / 2);

    // Fill data holder
    for (int adir = 0; adir < SpaceDim; ++adir) {
        for (int bdir = adir; bdir < SpaceDim; ++bdir) {
            const int comp = this->symTensorCompCC(adir, bdir);

            s_metricSourcePtr->fill_gdn(a_dest,
                                        comp,
                                        adir, bdir,
                                        m_dXi);
        }
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the contravariant metric elements (static version)
// -----------------------------------------------------------------------------
void LevelGeometry::fill_gup (FArrayBox&      a_dest,
                              const RealVect& a_dXi)
{
    CH_TIME("LevelGeometry::fill_gup (static version)");

    // Sanity check
    CH_assert(a_dest.nComp() == (SpaceDim * (SpaceDim + 1)) / 2);

    // Fill data holder
    for (int adir = 0; adir < SpaceDim; ++adir) {
        for (int bdir = adir; bdir < SpaceDim; ++bdir) {
            const int comp = symTensorCompCC(adir, bdir);

            s_metricSourcePtr->fill_gup(a_dest,
                                        comp,
                                        adir, bdir,
                                        a_dXi);
        }
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with the contravariant metric elements
// -----------------------------------------------------------------------------
void LevelGeometry::fill_gup (FArrayBox& a_dest) const
{
    CH_TIME("LevelGeometry::fill_gup");

    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_dest.nComp() == (SpaceDim * (SpaceDim + 1)) / 2);

    // Fill data holder
    for (int adir = 0; adir < SpaceDim; ++adir) {
        for (int bdir = adir; bdir < SpaceDim; ++bdir) {
            const int comp = this->symTensorCompCC(adir, bdir);

            s_metricSourcePtr->fill_gup(a_dest,
                                        comp,
                                        adir, bdir,
                                        m_dXi);
        }
    }
}


// -----------------------------------------------------------------------------
// Fills a FluxBox with the contravariant metric elements
// -----------------------------------------------------------------------------
void LevelGeometry::fill_gup (FluxBox& a_dest) const
{
    CH_TIME("LevelGeometry::fill_gup");

    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_dest.nComp() == SpaceDim);

    // Fill data holder
    for (int adir = 0; adir < SpaceDim; ++adir) {
        for (int bdir = 0; bdir < SpaceDim; ++bdir) {
            s_metricSourcePtr->fill_gup(a_dest[adir], bdir,
                                        adir, bdir,
                                        m_dXi);
        }
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with just one contravariant metric element
// -----------------------------------------------------------------------------
void LevelGeometry::fill_gup (FArrayBox& a_dest, int a_1, int a_2) const
{
    CH_TIME("LevelGeometry::fill_gup");

    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_dest.nComp() == 1);
    CH_assert(0 <= a_1 && a_1 < SpaceDim);
    CH_assert(0 <= a_2 && a_2 < SpaceDim);

    // Fill data holder
    s_metricSourcePtr->fill_gup(a_dest, 0,
                                a_1, a_2,
                                m_dXi);
}


// -----------------------------------------------------------------------------
// Fills a FluxBox with J * the contravariant metric elements
// -----------------------------------------------------------------------------
void LevelGeometry::fill_Jgup (FluxBox& a_dest) const
{
    CH_TIME("LevelGeometry::fill_Jgup");

    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_dest.nComp() == SpaceDim);

    // Fill data holder
    for (int adir = 0; adir < SpaceDim; ++adir) {
        for (int bdir = 0; bdir < SpaceDim; ++bdir) {
            s_metricSourcePtr->fill_Jgup(a_dest[adir], bdir,
                                         adir, bdir,
                                         m_dXi);
        }
    }
}


// -----------------------------------------------------------------------------
// Fills an FArrayBox with Jgup^{a_mu,*}
// This is useful for computing transverse derivatives.
// NOTE: The FArrayBox can have any centering, but must have SpaceDim comps. The
// SpaceDim comps of Jgup^{a_mu,*} will fill the SpaceDim comps of a_dest.
// -----------------------------------------------------------------------------
void LevelGeometry::fill_Jgup (FArrayBox& a_dest, int a_mu) const
{
    CH_TIME("LevelGeometry::fill_Jgup2");

    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_dest.nComp() == SpaceDim);
    CH_assert(0 <= a_mu && a_mu < SpaceDim);

    // Fill data holder
    for (int nu = 0; nu < SpaceDim; ++nu) {
        s_metricSourcePtr->fill_Jgup(a_dest, nu,
                                     a_mu, nu,
                                     m_dXi);
    }
}


// -----------------------------------------------------------------------------
// Fills a single-component FArrayBox with Gamma^{a_1}_{a_2, a_3}.
// This is scaled as h^{a_2} h^{a_3} / h^{a_1} where h = physDxCoarse.
// -----------------------------------------------------------------------------
void LevelGeometry::fill_Gamma (FArrayBox& a_dest, int a_1, int a_2, int a_3) const
{
    CH_TIME("LevelGeometry::fill_Gamma (single-comp)");

    // Sanity checks
    CH_assert(this->isDefined());
    CH_assert(a_dest.nComp() == 1);
    // CH_assert(!isUniform()); // Should be handled as a special case.

    // Fill data holder
    s_metricSourcePtr->fill_Gamma(a_dest,
                                  0,  // a_dest comp
                                  a_1, a_2, a_3,
                                  m_dXi);
}


// -----------------------------------------------------------------------------
// Ensures that metric data is the average of finer data.
// -----------------------------------------------------------------------------
void LevelGeometry::averageMetricsDown ()
{
    CH_TIME("LevelGeometry::averageMetricsDown");

    Vector<LevelGeometry*> amrLevGeos = getAMRLevGeos();
    for (int fineLev = amrLevGeos.size()-1; fineLev > 0; --fineLev) {
        LevelGeometry* const crseLevGeoPtr = amrLevGeos[fineLev-1];

        const LevelGeometry* const fineLevGeoPtr = amrLevGeos[fineLev];
        const DisjointBoxLayout& fineGrids = fineLevGeoPtr->getBoxes();
        const IntVect& refToCrse = fineLevGeoPtr->getCrseRefRatio();

        MappedCoarseAverage crseAvgObj;
        MappedCoarseAverageFace crseAvgFaceObj;
        int ncomps;

        { // CC J
            t_CCJPtr crseMetricPtr = crseLevGeoPtr->m_CCJPtr;
            const t_CCJPtr fineMetricPtr = fineLevGeoPtr->m_CCJPtr;
            ncomps = fineMetricPtr->nComp();
            crseAvgObj.define(fineGrids, ncomps, refToCrse);
            crseAvgObj.averageToCoarse(*crseMetricPtr, *fineMetricPtr);
        }

        { // CC Jinv
            t_CCJinvPtr crseMetricPtr = crseLevGeoPtr->m_CCJinvPtr;
            const t_CCJinvPtr fineMetricPtr = fineLevGeoPtr->m_CCJinvPtr;
            ncomps = fineMetricPtr->nComp();
            crseAvgObj.define(fineGrids, ncomps, refToCrse);
            crseAvgObj.averageToCoarseHarmonic(*crseMetricPtr, *fineMetricPtr);
        }

        { // CC gdn
            t_CCgdnPtr crseMetricPtr = crseLevGeoPtr->m_CCgdnPtr;
            const t_CCgdnPtr fineMetricPtr = fineLevGeoPtr->m_CCgdnPtr;
            ncomps = fineMetricPtr->nComp();
            crseAvgObj.define(fineGrids, ncomps, refToCrse);
            crseAvgObj.averageToCoarse(*crseMetricPtr, *fineMetricPtr);
        }

        { // FC gup
            t_FCgupPtr crseMetricPtr = crseLevGeoPtr->m_FCgupPtr;
            const t_FCgupPtr fineMetricPtr = fineLevGeoPtr->m_FCgupPtr;
            ncomps = fineMetricPtr->nComp();
            crseAvgFaceObj.define(fineGrids, ncomps, refToCrse);
            crseAvgFaceObj.averageToCoarseHarmonic(*crseMetricPtr, *fineMetricPtr);
        }

        { // FC Jgup
            t_FCJgupPtr crseMetricPtr = crseLevGeoPtr->m_FCJgupPtr;
            const t_FCJgupPtr fineMetricPtr = fineLevGeoPtr->m_FCJgupPtr;
            ncomps = fineMetricPtr->nComp();
            crseAvgFaceObj.define(fineGrids, ncomps, refToCrse);
            crseAvgFaceObj.averageToCoarseHarmonic(*crseMetricPtr, *fineMetricPtr);
        }
    }
}
