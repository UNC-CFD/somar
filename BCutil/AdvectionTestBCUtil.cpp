#include "AdvectionTestBCUtil.H"
#include "BoxIterator.H"


// -----------------------------------------------------------------------------
// Default constrictor
// -----------------------------------------------------------------------------
AdvectionTestBCUtil::AdvectionTestBCUtil ()
{;}


// -----------------------------------------------------------------------------
// Default destructor
// -----------------------------------------------------------------------------
AdvectionTestBCUtil::~AdvectionTestBCUtil ()
{;}


// -----------------------------------------------------------------------------
// This object is its own factory
// -----------------------------------------------------------------------------
PhysBCUtil* AdvectionTestBCUtil::newPhysBCUtil () const
{
    PhysBCUtil* newBCPtr = new AdvectionTestBCUtil();
    return newBCPtr;
}


// ************************ ICs / background fields ****************************

// -----------------------------------------------------------------------------
// Fills a FAB with the initial velocity
// Locations are in mapped space, but components are Cartesian.
// -----------------------------------------------------------------------------
void AdvectionTestBCUtil::setVelIC (FArrayBox&           a_velFAB,
                                    const int            a_velComp,
                                    const LevelGeometry& a_levGeo,
                                    const DataIndex&     a_di) const
{
    CH_assert(a_velFAB.nComp() == SpaceDim);
    CH_assert(0 <= a_velComp);
    CH_assert(a_velComp < SpaceDim);

    // a_velFAB.setVal(0.0, a_velComp);

    if (SpaceDim == 2) {
        if (a_velComp == 0) {
            a_velFAB.setVal(1.0, a_velComp);
        } else {
            a_velFAB.setVal(1.0, a_velComp);
        }
    } else if (SpaceDim == 3) {
        if (a_velComp == 0) {
            a_velFAB.setVal(1.0, a_velComp);
        } else if (a_velComp == 1) {
            a_velFAB.setVal(0.0, a_velComp);
        } else {
            a_velFAB.setVal(1.0, a_velComp);
        }
    } else {
        MayDay::Error("AdvectionTestBCUtil::setVelIC: Bad SpaceDim");
    }
}


// -----------------------------------------------------------------------------
// Fills a FAB with the initial scalars
// -----------------------------------------------------------------------------
void AdvectionTestBCUtil::setScalarIC (FArrayBox&           a_scalarFAB,
                                       const int            a_scalarComp,
                                       const LevelGeometry& a_levGeo,
                                       const DataIndex&     a_di) const
{
    CH_assert(a_scalarFAB.nComp() == 1);

    if (a_scalarComp == 0) {
        const Box domBox = a_levGeo.getDomain().domainBox();
        const Box dataBox = a_scalarFAB.box();

        // Get the Cartesian locations
        FArrayBox distFAB(dataBox, SpaceDim);
        a_levGeo.fill_physCoor(distFAB);

        // Convert to distances from the center point.
        for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
            const Real x0 = 0.5;
            distFAB.plus(-x0, dir);
        }

        // Start with zero density
        a_scalarFAB.setVal(0.0, 0);

        // Add the U shape
        BoxIterator bit(dataBox);
        const Real maxRsq = 0.01;
        const Real holeTop = 0.02;
        const Real holeRsq = 0.0004;
        RealVect distSq;

        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& cc = bit();

            D_TERM(distSq[0] = distFAB(cc,0)*distFAB(cc,0);,
                   distSq[1] = distFAB(cc,1)*distFAB(cc,1);,
                   distSq[2] = distFAB(cc,2)*distFAB(cc,2);)
            Real Rsq = D_TERM(distSq[0], +distSq[1], +distSq[2]);

            // Circle
            if (Rsq <= maxRsq) {
                a_scalarFAB(cc,0) = 1.0;
            }

            // Notch
            if (distSq[0] < holeRsq && distFAB(cc,SpaceDim-1) < holeTop) {
                a_scalarFAB(cc,0) = 0.0;
            }
        }

    } else {
        MayDay::Error("AdvectionTestBCUtil::setScalarIC received a_scalarComp > 0");
    }
}

