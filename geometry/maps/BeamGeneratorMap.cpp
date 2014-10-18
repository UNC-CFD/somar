#include "BeamGeneratorMap.H"
#include "BeamGeneratorMapF_F.H"
#include "Subspace.H"
#include "ProblemContext.H"

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
BeamGeneratorMap::BeamGeneratorMap ()
: BathymetricBaseMap()
{
    const ProblemContext* ctx = ProblemContext::getInstance();
    m_alpha = ctx->beamGenMapAlpha;
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
BeamGeneratorMap::~BeamGeneratorMap ()
{;}


// -----------------------------------------------------------------------------
// Must return the name of the coordinate mapping
// -----------------------------------------------------------------------------
const char* BeamGeneratorMap::getCoorMapName () const
{
    return "BeamGeneratorMap";
}


// -----------------------------------------------------------------------------
// Must return whether or not this metric is diagonal
// -----------------------------------------------------------------------------
bool BeamGeneratorMap::isDiagonal () const
{
    return false;
}


// -----------------------------------------------------------------------------
// Fills a NodeFAB with the bathymetric data. a_dest must be flat in the
// vertical. Upon return, each point in the horizontal (Xi,Eta) of a_dest
// will contain the (positive) local elevation from the sea floor.
// NOTE: This vertical distance is measured in a straight line perpendicular
// the the surface. We are measuring this distance along the Cartesian
// vertical coordinate line, not the mapped vertical coordinate line.
// -----------------------------------------------------------------------------
void BeamGeneratorMap::fill_bathymetry (FArrayBox&       a_dest,
                                        const int        a_destComp,
                                        const FArrayBox& a_cartPos,
                                        const RealVect&  a_dXi) const
{
    const Box& destBox = a_dest.box();
    const IntVect destBoxType = destBox.type();

    // The holder needs to be flat and nodal in the vertical.
    CH_assert(destBox == horizontalDataBox(destBox));
    CH_assert(destBoxType[SpaceDim-1] == 1);

    FORT_FILL_BEAMGENERATORMAPBATHYMETRY (
        CHF_FRA1(a_dest, a_destComp),
        CHF_BOX(destBox),
        CHF_CONST_INTVECT(destBoxType),
        CHF_CONST_FRA(a_cartPos),
        CHF_CONST_REALVECT(a_dXi),
        CHF_CONST_REALVECT(m_L),
        CHF_CONST_REAL(m_alpha));
}
