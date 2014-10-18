#include "RelaxationMethod.H"
#include "HomogeneousCFInterp.H"
#include "ExtrapolationUtils.H"
#include "Debug.H"


// -----------------------------------------------------------------------------
// Constructor -- leaves object unusable.
// -----------------------------------------------------------------------------
RelaxationMethod::RelaxationMethod ()
: m_isDefined(false)
{;}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
RelaxationMethod::~RelaxationMethod ()
{;}


// -----------------------------------------------------------------------------
// Define constructor.
// -----------------------------------------------------------------------------
void RelaxationMethod::define (const Real                                  a_alpha,
                               const Real                                  a_beta,
                               const RealVect&                             a_dx,
                               const RealVect&                             a_crseDx,
                               const RefCountedPtr<LevelData<FluxBox> >&   a_FCJgup,
                               const RefCountedPtr<LevelData<FArrayBox> >& a_CCJinv,
                               const RefCountedPtr<LevelData<FArrayBox> >& a_lapDiag,
                               BCMethodHolder&                             a_bc,
                               const Copier&                               a_exchangeCopier,
                               const CFRegion&                             a_cfRegion,
                               const bool                                  a_isDiagonal,
                               const IntVect&                              a_activeDirs)
{
    // Sanity checks
    CH_assert(!isDefined());

    CH_assert(!a_FCJgup.isNull());
    CH_assert(!a_CCJinv.isNull());
    CH_assert(!a_lapDiag.isNull());

    CH_assert(a_FCJgup->getBoxes().compatible(a_lapDiag->getBoxes()));
    CH_assert(a_CCJinv->getBoxes().compatible(a_lapDiag->getBoxes()));

    CH_assert(a_bc.isDefined());
    CH_assert(a_exchangeCopier.isDefined());

    D_TERM(CH_assert(a_activeDirs[0] == 0 || a_activeDirs[0] == 1);,
           CH_assert(a_activeDirs[1] == 0 || a_activeDirs[1] == 1);,
           CH_assert(a_activeDirs[2] == 0 || a_activeDirs[2] == 1);)

    // Copy inputs
    m_alpha = a_alpha;
    m_beta = a_beta;
    m_dx = a_dx;
    m_crseDx = a_crseDx;
    m_FCJgup = a_FCJgup;
    m_CCJinv = a_CCJinv;
    m_lapDiag = a_lapDiag;
    m_bc = a_bc;
    m_exchangeCopier = a_exchangeCopier;
    m_cfRegion = a_cfRegion;
    m_activeDirs = a_activeDirs;
    m_isDiagonal = a_isDiagonal;

    // Allocate scratch space
    m_extrap.define(a_lapDiag->getBoxes(), 1, a_activeDirs);

    // Create boundary maps. This must be called AFTER members are set!
    this->collectBoundaryData();

    // This object is ready for use.
    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// This collects the boundary box data needed by (simple)boundaryGSRB.
// -----------------------------------------------------------------------------
void RelaxationMethod::collectBoundaryData ()
{
    CH_TIME("RelaxationMethod::collectBoundaryData");

    CH_assert(!m_lapDiag.isNull());

    // Gather grid info
    const DisjointBoxLayout& grids = m_lapDiag->getBoxes();
    DataIterator dit = grids.dataIterator();
    const ProblemDomain& domain = grids.physDomain();
    const Box& domBox = domain.domainBox();
    const BCDescriptor& fluxDesc = m_bc.getFluxDescriptor();

    // Compute domain boundaries and interior
    Box domBdry[CH_SPACEDIM][2];
    for (int dir = 0; dir < SpaceDim; ++dir) {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide iside = sit();
            const int isign = sign(iside);

            domBdry[dir][iside] = bdryBox(domBox, dir, iside, 1);
            domBdry[dir][iside].shiftHalf(dir, -isign);
        }
    }
    const Box domInterior = grow(domBox, -m_activeDirs);

    // Compute domain boundaries first...
    m_boundaryBoxData.resize(0);

    for (dit.reset(); dit.ok(); ++dit) {
        const Box& valid = grids[dit];
        const Box interior = valid & domInterior;

        BoundaryBoxData data;
        data.index = dit();
        data.valid = valid;

#if CH_SPACEDIM == 2
        // Do face sweeps
        Box faceBox, edgeBox;
        for (int fdir = 0; fdir < CH_SPACEDIM; ++fdir) {
            if (m_activeDirs[fdir] == 0) continue;

            for (SideIterator fsit; fsit.ok(); ++fsit) {
                faceBox = adjCellBox(interior, fdir, fsit(), 1) & valid;
                if (!domBdry[fdir][fsit()].intersects(faceBox)) continue;

                data.stencil[0][Side::Lo] = fluxDesc.stencil(faceBox, domain, 0, Side::Lo);
                data.stencil[0][Side::Hi] = fluxDesc.stencil(faceBox, domain, 0, Side::Hi);
                data.stencil[1][Side::Lo] = fluxDesc.stencil(faceBox, domain, 1, Side::Lo);
                data.stencil[1][Side::Hi] = fluxDesc.stencil(faceBox, domain, 1, Side::Hi);
                data.validBdry = faceBox;
                m_boundaryBoxData.push_back(data);

                // Do edge sweeps
                for (int edir = fdir + 1; edir < CH_SPACEDIM; ++edir) {
                    if (m_activeDirs[edir] == 0) continue;

                    for (SideIterator esit; esit.ok(); ++esit) {
                        edgeBox = adjCellBox(faceBox, edir, esit(), 1) & valid;
                        if (!domBdry[edir][esit()].intersects(edgeBox)) continue;

                        data.stencil[0][Side::Lo] = fluxDesc.stencil(edgeBox, domain, 0, Side::Lo);
                        data.stencil[0][Side::Hi] = fluxDesc.stencil(edgeBox, domain, 0, Side::Hi);
                        data.stencil[1][Side::Lo] = fluxDesc.stencil(edgeBox, domain, 1, Side::Lo);
                        data.stencil[1][Side::Hi] = fluxDesc.stencil(edgeBox, domain, 1, Side::Hi);
                        data.validBdry = edgeBox;
                        m_boundaryBoxData.push_back(data);

                    } // end loop over edge sides (esit)
                } // end loop over edge dirs (edir)
            } // end loop over face sides (fsit)
        } // end loop over face dirs (fdir)

#elif CH_SPACEDIM == 3
        // Do face sweeps
        Box faceBox, edgeBox, vertexBox;
        for (int fdir = 0; fdir < CH_SPACEDIM; ++fdir) {
            if (m_activeDirs[fdir] == 0) continue;

            for (SideIterator fsit; fsit.ok(); ++fsit) {
                faceBox = adjCellBox(interior, fdir, fsit(), 1) & valid;
                if (!domBdry[fdir][fsit()].intersects(faceBox)) continue;

                data.stencil[0][Side::Lo] = fluxDesc.stencil(faceBox, domain, 0, Side::Lo);
                data.stencil[0][Side::Hi] = fluxDesc.stencil(faceBox, domain, 0, Side::Hi);
                data.stencil[1][Side::Lo] = fluxDesc.stencil(faceBox, domain, 1, Side::Lo);
                data.stencil[1][Side::Hi] = fluxDesc.stencil(faceBox, domain, 1, Side::Hi);
                data.stencil[2][Side::Lo] = fluxDesc.stencil(faceBox, domain, 2, Side::Lo);
                data.stencil[2][Side::Hi] = fluxDesc.stencil(faceBox, domain, 2, Side::Hi);
                data.validBdry = faceBox;
                m_boundaryBoxData.push_back(data);

                // Do edge sweeps
                for (int edir = fdir + 1; edir < CH_SPACEDIM; ++edir) {
                    if (m_activeDirs[edir] == 0) continue;

                    for (SideIterator esit; esit.ok(); ++esit) {
                        edgeBox = adjCellBox(faceBox, edir, esit(), 1) & valid;
                        if (!domBdry[edir][esit()].intersects(edgeBox)) continue;

                        data.stencil[0][Side::Lo] = fluxDesc.stencil(edgeBox, domain, 0, Side::Lo);
                        data.stencil[0][Side::Hi] = fluxDesc.stencil(edgeBox, domain, 0, Side::Hi);
                        data.stencil[1][Side::Lo] = fluxDesc.stencil(edgeBox, domain, 1, Side::Lo);
                        data.stencil[1][Side::Hi] = fluxDesc.stencil(edgeBox, domain, 1, Side::Hi);
                        data.stencil[2][Side::Lo] = fluxDesc.stencil(edgeBox, domain, 2, Side::Lo);
                        data.stencil[2][Side::Hi] = fluxDesc.stencil(edgeBox, domain, 2, Side::Hi);
                        data.validBdry = edgeBox;
                        m_boundaryBoxData.push_back(data);

                        // Do vertex sweeps
                        const int vdir = SpaceDim - fdir - edir;
                        if (vdir <= edir) continue; // We only want to process each vertex once.
                        if (m_activeDirs[vdir] == 0) continue;

                        for (SideIterator vsit; vsit.ok(); ++vsit) {
                            vertexBox = adjCellBox(edgeBox, vdir, vsit(), 1) & valid;
                            if (!domBdry[vdir][vsit()].intersects(vertexBox)) continue;

                            data.stencil[0][Side::Lo] = fluxDesc.stencil(vertexBox, domain, 0, Side::Lo);
                            data.stencil[0][Side::Hi] = fluxDesc.stencil(vertexBox, domain, 0, Side::Hi);
                            data.stencil[1][Side::Lo] = fluxDesc.stencil(vertexBox, domain, 1, Side::Lo);
                            data.stencil[1][Side::Hi] = fluxDesc.stencil(vertexBox, domain, 1, Side::Hi);
                            data.stencil[2][Side::Lo] = fluxDesc.stencil(vertexBox, domain, 2, Side::Lo);
                            data.stencil[2][Side::Hi] = fluxDesc.stencil(vertexBox, domain, 2, Side::Hi);
                            data.validBdry = vertexBox;
                            m_boundaryBoxData.push_back(data);

                        } // end loop over vertex sides (vsit)
                    } // end loop over edge sides (esit)
                } // end loop over edge dirs (edir)
            } // end loop over face sides (fsit)
        } // end loop over face dirs (fdir)
#else
#   error Bad SpaceDim
#endif
    } // end loop over grids (dit)


    // Next, compute valid boundaries...
    m_simpleBoundaryBoxData.resize(0);

    for (dit.reset(); dit.ok(); ++dit) {
        const FluxBox&   JgupFlub = (*m_FCJgup)[dit];
        const FArrayBox& JinvFAB  = (*m_CCJinv)[dit];
        const Box& valid = grids[dit];
        const Box interior = grow(valid, -m_activeDirs);

        BoundaryBoxData data;
        data.index = dit();
        data.valid = valid;

#if CH_SPACEDIM == 2
        // Do face sweeps
        Box faceBox, edgeBox;
        for (int fdir = 0; fdir < CH_SPACEDIM; ++fdir) {
            if (m_activeDirs[fdir] == 0) continue;

            for (SideIterator fsit; fsit.ok(); ++fsit) {
                faceBox = adjCellBox(interior, fdir, fsit(), 1) & valid;

                data.stencil[0][Side::Lo] = fluxDesc.stencil(faceBox, domain, 0, Side::Lo);
                data.stencil[0][Side::Hi] = fluxDesc.stencil(faceBox, domain, 0, Side::Hi);
                data.stencil[1][Side::Lo] = fluxDesc.stencil(faceBox, domain, 1, Side::Lo);
                data.stencil[1][Side::Hi] = fluxDesc.stencil(faceBox, domain, 1, Side::Hi);
                data.validBdry = faceBox;
                m_simpleBoundaryBoxData.push_back(data);

                // Do edge sweeps
                for (int edir = fdir + 1; edir < CH_SPACEDIM; ++edir) {
                    if (m_activeDirs[edir] == 0) continue;

                    for (SideIterator esit; esit.ok(); ++esit) {
                        edgeBox = adjCellBox(faceBox, edir, esit(), 1) & valid;

                        data.stencil[0][Side::Lo] = fluxDesc.stencil(edgeBox, domain, 0, Side::Lo);
                        data.stencil[0][Side::Hi] = fluxDesc.stencil(edgeBox, domain, 0, Side::Hi);
                        data.stencil[1][Side::Lo] = fluxDesc.stencil(edgeBox, domain, 1, Side::Lo);
                        data.stencil[1][Side::Hi] = fluxDesc.stencil(edgeBox, domain, 1, Side::Hi);
                        data.validBdry = edgeBox;
                        m_simpleBoundaryBoxData.push_back(data);

                    } // end loop over edge sides (esit)
                } // end loop over edge dirs (edir)
            } // end loop over face sides (fsit)
        } // end loop over face dirs (fdir)

#elif CH_SPACEDIM == 3
        // Do face sweeps
        Box faceBox, edgeBox, vertexBox;
        for (int fdir = 0; fdir < CH_SPACEDIM; ++fdir) {
            if (m_activeDirs[fdir] == 0) continue;

            for (SideIterator fsit; fsit.ok(); ++fsit) {
                faceBox = adjCellBox(interior, fdir, fsit(), 1) & valid;

                data.stencil[0][Side::Lo] = fluxDesc.stencil(faceBox, domain, 0, Side::Lo);
                data.stencil[0][Side::Hi] = fluxDesc.stencil(faceBox, domain, 0, Side::Hi);
                data.stencil[1][Side::Lo] = fluxDesc.stencil(faceBox, domain, 1, Side::Lo);
                data.stencil[1][Side::Hi] = fluxDesc.stencil(faceBox, domain, 1, Side::Hi);
                data.stencil[2][Side::Lo] = fluxDesc.stencil(faceBox, domain, 2, Side::Lo);
                data.stencil[2][Side::Hi] = fluxDesc.stencil(faceBox, domain, 2, Side::Hi);
                data.validBdry = faceBox;
                m_simpleBoundaryBoxData.push_back(data);

                // Do edge sweeps
                for (int edir = fdir + 1; edir < CH_SPACEDIM; ++edir) {
                    if (m_activeDirs[edir] == 0) continue;

                    for (SideIterator esit; esit.ok(); ++esit) {
                        edgeBox = adjCellBox(faceBox, edir, esit(), 1) & valid;

                        data.stencil[0][Side::Lo] = fluxDesc.stencil(edgeBox, domain, 0, Side::Lo);
                        data.stencil[0][Side::Hi] = fluxDesc.stencil(edgeBox, domain, 0, Side::Hi);
                        data.stencil[1][Side::Lo] = fluxDesc.stencil(edgeBox, domain, 1, Side::Lo);
                        data.stencil[1][Side::Hi] = fluxDesc.stencil(edgeBox, domain, 1, Side::Hi);
                        data.stencil[2][Side::Lo] = fluxDesc.stencil(edgeBox, domain, 2, Side::Lo);
                        data.stencil[2][Side::Hi] = fluxDesc.stencil(edgeBox, domain, 2, Side::Hi);
                        data.validBdry = edgeBox;
                        m_simpleBoundaryBoxData.push_back(data);

                        // Do vertex sweeps
                        const int vdir = SpaceDim - fdir - edir;
                        if (vdir <= edir) continue; // We only want to process each vertex once.
                        if (m_activeDirs[vdir] == 0) continue;

                        for (SideIterator vsit; vsit.ok(); ++vsit) {
                            vertexBox = adjCellBox(edgeBox, vdir, vsit(), 1) & valid;

                            data.stencil[0][Side::Lo] = fluxDesc.stencil(vertexBox, domain, 0, Side::Lo);
                            data.stencil[0][Side::Hi] = fluxDesc.stencil(vertexBox, domain, 0, Side::Hi);
                            data.stencil[1][Side::Lo] = fluxDesc.stencil(vertexBox, domain, 1, Side::Lo);
                            data.stencil[1][Side::Hi] = fluxDesc.stencil(vertexBox, domain, 1, Side::Hi);
                            data.stencil[2][Side::Lo] = fluxDesc.stencil(vertexBox, domain, 2, Side::Lo);
                            data.stencil[2][Side::Hi] = fluxDesc.stencil(vertexBox, domain, 2, Side::Hi);
                            data.validBdry = vertexBox;
                            m_simpleBoundaryBoxData.push_back(data);

                        } // end loop over vertex sides (vsit)
                    } // end loop over edge sides (esit)
                } // end loop over edge dirs (edir)
            } // end loop over face sides (fsit)
        } // end loop over face dirs (fdir)
#else
#   error Bad SpaceDim
#endif
    } // end loop over grids (dit)


#if 0
    // TEST: The entire splitTest field should be 1.0
    LevelData<FArrayBox> splitTest(grids, 1);
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& testFAB = splitTest[dit];
        testFAB.setVal(0.0);

        Box interior = grow(grids[dit], -1);
        testFAB.plus(1.0, interior);
    }
    for (int idx = 0; idx < m_simpleBoundaryBoxData.size(); ++idx) {
        const BoundaryBoxData data = m_simpleBoundaryBoxData[idx];
        FArrayBox& testFAB = splitTest[data.index];

        testFAB.plus(1.0, data.validBdry);
    }
    writeLevelHDF5(splitTest, 0.0, false);

    LevelData<FArrayBox> splitTest(grids, 1);
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& testFAB = splitTest[dit];
        testFAB.setVal(0.0);

        testFAB.plus(1.0, testFAB.box() & domInterior);
    }
    for (int idx = 0; idx < m_boundaryBoxData.size(); ++idx) {
        const BoundaryBoxData data = m_boundaryBoxData[idx];
        FArrayBox& testFAB = splitTest[data.index];

        testFAB.plus(1.0, data.validBdry);
    }
    writeLevelHDF5(splitTest, 0.0, false);
#endif
}




// -----------------------------------------------------------------------------
// As advertised, this will fill the ghosts of a_phi at the physical
// boundary using the BCs and at the CF interface using a homogeneous
// coarser level. This function will also extrapolate a_phi's valid data to
// help compute mixed second derivatives when the metric is not diagonal.
// This function is not blocking and the exchange is left to the caller.
// -----------------------------------------------------------------------------
void RelaxationMethod::fillGhostsAndExtrapolate (LevelData<FArrayBox>& a_phi,
                                                 const bool            a_doCFInterp,
                                                 const bool            a_doExtrapAndCopy,
                                                 const bool            a_doBCs)
{
    CH_TIME("RelaxationMethod::fillGhostsAndExtrapolate");

    // Sanity checks
    CH_assert(isDefined());
    CH_assert(m_extrap.nComp() == a_phi.nComp()); // If this trips, send comps into define.
    CH_assert(m_extrap.getBoxes().compatible(a_phi.getBoxes()));
    CH_assert(m_FCJgup->getBoxes().compatible(a_phi.getBoxes()));
    CH_assert(m_extrap.getBoxes().physDomain() == a_phi.getBoxes().physDomain());
    CH_assert(m_FCJgup->getBoxes().physDomain() == a_phi.getBoxes().physDomain());

    // Gather needed info
    const int extrapOrder = 1;
    const DisjointBoxLayout& grids = a_phi.disjointBoxLayout();
    DataIterator dit = grids.dataIterator();
    const ProblemDomain& domain = grids.physDomain();

    // Fill ghosts at CF interface
    if (a_doCFInterp) {
        homogeneousCFInterp(a_phi, m_dx, m_crseDx, m_cfRegion, m_activeDirs);
    }

    // Set physical ghosts and extrapolate.
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox&       phiFAB     = a_phi[dit];
        FArrayBox&       extrapFAB  = m_extrap[dit];
        const FluxBox&   JgupFB     = (*m_FCJgup)[dit];
        const Box&       valid      = grids[dit];

        // Extrapolate ghosts for non-diagonal derivatives
        if (!m_isDiagonal && a_doExtrapAndCopy) {
            extrapFAB.copy(phiFAB);

            Box domValid = domain.domainBox();
            domValid &= extrapFAB.box();

            for (int fdir = 0; fdir < CH_SPACEDIM; ++fdir) {
                if (m_activeDirs[fdir] == 0) continue;

                SideIterator fsit;
                for (fsit.reset(); fsit.ok(); ++fsit) {
                    ExtrapolateFaceAndCopy(extrapFAB, extrapFAB, domValid, fdir, fsit(), extrapOrder);
                }
                domValid.grow(fdir, 1);
            }
        }

        // Set physical BCs
        // Note that homogeneous BCs should not need the time.
        if (a_doBCs) {
            const Real bogusTime = 1e-300;
            const bool isHomogeneous = true;
            m_bc.setGhosts(phiFAB, &extrapFAB, valid, domain, m_dx, dit(), &JgupFB, isHomogeneous, bogusTime);
        }
    }
}
