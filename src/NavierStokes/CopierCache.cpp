#include "CopierCache.H"
#include "CornerCopier.H"
#include "CH_Timer.H"


// -----------------------------------------------------------------------------
// Defines a single Copier entry.
// -----------------------------------------------------------------------------
void CopierCache::CopierEntry::define (const DisjointBoxLayout& a_grids,
                                       const IntVect&           a_ghostVect,
                                       const bool               a_isExchange)
{
    this->grids = a_grids;
    this->ghostVect = a_ghostVect;
    this->isExchange = a_isExchange;

    this->cp.define(this->grids,
                    this->grids,
                    this->grids.physDomain(),
                    this->ghostVect,
                    this->isExchange);
}


// -----------------------------------------------------------------------------
// Defines a single CornerCopier entry.
// -----------------------------------------------------------------------------
void CopierCache::CornerCopierEntry::define (const DisjointBoxLayout& a_grids,
                                             const IntVect&           a_ghostVect,
                                             const bool               a_isExchange)
{
    this->grids = a_grids;
    this->ghostVect = a_ghostVect;
    this->isExchange = a_isExchange;

    this->cp.define(this->grids,
                    this->grids,
                    this->grids.physDomain(),
                    this->ghostVect,
                    this->isExchange);
}


// -----------------------------------------------------------------------------
// Constructor.
// -----------------------------------------------------------------------------
CopierCache::CopierCache (const IntVect& a_tracingGhosts)
: m_tracingGhosts(a_tracingGhosts)
{;}


// -----------------------------------------------------------------------------
// Destructor.
// -----------------------------------------------------------------------------
CopierCache::~CopierCache ()
{;}


// -----------------------------------------------------------------------------
// Accessor. Redefines if needed.
// -----------------------------------------------------------------------------
const Copier& CopierCache::getTracingExCopier (const DisjointBoxLayout& a_grids)
{
    CopierEntry& entry = m_tracingExCopierEntry;

    if (!(a_grids == entry.grids)) {
        this->regrid(a_grids);
    }

    return entry.cp;
}


// -----------------------------------------------------------------------------
// Accessor. Redefines if needed.
// -----------------------------------------------------------------------------
const CornerCopier& CopierCache::getTracingExCornerCopier (const DisjointBoxLayout& a_grids)
{
    CornerCopierEntry& entry = m_tracingExCornerCopierEntry;

    if (!(a_grids == entry.grids)) {
        this->regrid(a_grids);
    }

    return entry.cp;
}


// -----------------------------------------------------------------------------
// Accessor. Redefines if needed.
// -----------------------------------------------------------------------------
const Copier& CopierCache::getOneGhostExCopier (const DisjointBoxLayout& a_grids)
{
    CopierEntry& entry = m_oneGhostExCopierEntry;

    if (!(a_grids == entry.grids)) {
        this->regrid(a_grids);
    }

    return entry.cp;
}


// -----------------------------------------------------------------------------
// Defines the copiers with new grids.
// This function re-defines the copiers even if they are already defined
// over the same set of grids.
// -----------------------------------------------------------------------------
void CopierCache::regrid (const DisjointBoxLayout& a_grids)
{
    CH_TIME("CopierCache::regrid");
    // pout() << "CopierCache::regrid" << std::endl;

    m_tracingExCopierEntry.define(a_grids, m_tracingGhosts, true);
    m_tracingExCornerCopierEntry.define(a_grids, m_tracingGhosts, true);
    m_oneGhostExCopierEntry.define(a_grids, IntVect::Unit, true);
}
