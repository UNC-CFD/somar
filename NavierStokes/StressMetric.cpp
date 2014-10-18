#include "StressMetric.H"
#ifdef USE_STRESSMETRIC

#include "LevelGeometry.H"
#include "ProblemContext.H"
#include "StressMetricF_F.H"


// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
StressMetric::StressMetric ()
: m_geoSourcePtr(NULL)
{
    const ProblemContext* ctx = ProblemContext::getInstance();
    m_dynVisc = ctx->nu;
}


// -----------------------------------------------------------------------------
// Full constructor
// -----------------------------------------------------------------------------
StressMetric::StressMetric (const GeoSourceInterface* a_geoSourcePtr)
{
    this->define(a_geoSourcePtr);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
StressMetric::~StressMetric ()
{
    m_geoSourcePtr = NULL;
}


// -----------------------------------------------------------------------------
// Full constructor
// -----------------------------------------------------------------------------
void StressMetric::define (const GeoSourceInterface* a_geoSourcePtr)
{
    m_geoSourcePtr = a_geoSourcePtr;

    const ProblemContext* ctx = ProblemContext::getInstance();
    m_dynVisc = ctx->nu;
}


// -----------------------------------------------------------------------------
// Set dest = (1+turbVisc/dynVisc)*J*gup^{ij}
// Remember: The solvers will multiply by dynVisc, that's why it's factored out.
// -----------------------------------------------------------------------------
void StressMetric::fill_Jgup (FArrayBox&       a_dest,
                              const int        a_destComp,
                              const int        a_mu,
                              const int        a_nu,
                              const RealVect&  a_dXi,
                              const Real       a_scale) const
{
    // First, fill with the metric tensor
    m_geoSourcePtr->fill_Jgup(a_dest, a_destComp, a_mu, a_nu, a_dXi, a_scale);

    // Then, if this is a viscous fluid, multiply by (1+turbVisc/dynVisc)
    if (m_dynVisc != 0.0) {
        const Box& destBox = a_dest.box();
        const IntVect& destBoxType = destBox.type();

        FArrayBox physCoor(destBox, SpaceDim);
        m_geoSourcePtr->fill_physCoor(physCoor, a_dXi);

        FORT_APPLYVISCOUSFACTOR (
            CHF_FRA1(a_dest, a_destComp),
            CHF_CONST_FRA(physCoor),
            CHF_CONST_REALVECT(a_dXi),
            CHF_CONST_REAL(m_dynVisc),
            CHF_BOX(destBox),
            CHF_CONST_INTVECT(destBoxType));
    }
}


#endif //USE_STRESSMETRIC
