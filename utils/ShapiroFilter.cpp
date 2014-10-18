#include "ShapiroFilter.H"
#include "BoxIterator.H"
#include "EllipticBCUtils.H"
#include "Printing.H"


// -----------------------------------------------------------------------------
// This uses a 3-point stencil on each pass. The parameters and number of passes
// are hard-coded into this funciton, but can be changed easily.
// -----------------------------------------------------------------------------
void ShapiroFilter1D (// Data stuff
                      LevelData<FArrayBox>&       a_data,
                      const Interval              a_interval,
                      // CF BC stuff
                      const LevelData<FArrayBox>* a_crseDataPtr,
                      const MappedQuadCFInterp&   a_cfInterp,
                      // Physical BC stuff
                      VelBCHolder                 a_bc,
                      const bool                  a_isHomogeneous,
                      const LevelGeometry&        a_levGeo,
                      const Real                  a_time,
                      // Filter stuff
                      const int                   a_filterDir)
{
    Vector<Real> S(2);
    S[0] = 0.5;
    S[1] = -0.5;

    const DisjointBoxLayout& grids = a_data.getBoxes();
    DataIterator dit = grids.dataIterator();
    const RealVect& dx = a_levGeo.getDx();
    const LevelData<FluxBox>& Jgup = a_levGeo.getFCJgup();
    const IntVect e = BASISV(a_filterDir);

    for (int pass = 0; pass < S.size(); ++pass) {
        if (a_crseDataPtr != NULL) {
            CH_assert(a_cfInterp.isDefined());
            a_cfInterp.coarseFineInterp(a_data, *a_crseDataPtr);
        }
        a_bc.setGhosts(a_data,
                       NULL,    // extrap ptr
                       dx,
                       &Jgup,
                       a_isHomogeneous,
                       a_time);
        a_data.exchange();

        const Real alpha = 1.0 - S[pass];
        const Real beta = S[pass];

        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox w(a_interval, a_data[dit]);
            BoxIterator bit(grids[dit]);

            for (int comp = 0; comp < w.nComp(); ++comp) {
                for (bit.reset(); bit.ok(); ++bit) {
                    const IntVect& cc = bit();
                    w(cc,comp) = alpha*w(cc,comp) + beta*0.5*(w(cc-e,comp)+w(cc+e,comp));
                }
            }
        }
    }
}


// -----------------------------------------------------------------------------
// This uses a 5-point stencil on each pass. The parameters and number of passes
// are hard-coded into this funciton, but can be changed easily.
// Use this carefully. I found that the 9-point stencil that arises from using
// ShapiroFilter1D twice (it's commutative) produces better results. Shapiro
// talks about this in his paper.
// -----------------------------------------------------------------------------
void ShapiroFilter2D (// Data stuff
                      LevelData<FArrayBox>&       a_data,
                      const Interval              a_interval,
                      // CF BC stuff
                      const LevelData<FArrayBox>* a_crseDataPtr,
                      const MappedQuadCFInterp&   a_cfInterp,
                      // Physical BC stuff
                      VelBCHolder                 a_bc,
                      const bool                  a_isHomogeneous,
                      const LevelGeometry&        a_levGeo,
                      const Real                  a_time,
                      // Filter stuff
                      const int                   a_filterDir0,
                      const int                   a_filterDir1)
{
    Vector<Real> S(2);
    S[0] = 0.5;
    S[1] = -0.5;

    const DisjointBoxLayout& grids = a_data.getBoxes();
    DataIterator dit = grids.dataIterator();
    const RealVect& dx = a_levGeo.getDx();
    const LevelData<FluxBox>& Jgup = a_levGeo.getFCJgup();
    const IntVect e0 = BASISV(a_filterDir0);
    const IntVect e1 = BASISV(a_filterDir1);

    for (int pass = 0; pass < S.size(); ++pass) {
        if (a_crseDataPtr != NULL) {
            CH_assert(a_cfInterp.isDefined());
            a_cfInterp.coarseFineInterp(a_data, *a_crseDataPtr);
        }
        a_bc.setGhosts(a_data,
                       NULL,    // extrap ptr
                       dx,
                       &Jgup,
                       a_isHomogeneous,
                       a_time);
        a_data.exchange();

        const Real alpha = 1.0 - S[pass];
        const Real beta = S[pass];

        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox w(a_interval, a_data[dit]);
            BoxIterator bit(grids[dit]);
            Real avg;

            for (int comp = 0; comp < w.nComp(); ++comp) {
                for (bit.reset(); bit.ok(); ++bit) {
                    const IntVect& cc = bit();
                    avg = 0.25 * (w(cc-e0,comp) + w(cc+e0,comp) + w(cc-e1,comp) + w(cc+e1,comp));
                    w(cc,comp) = alpha * w(cc,comp) + beta * avg;
                }
            }
        }
    }
}

