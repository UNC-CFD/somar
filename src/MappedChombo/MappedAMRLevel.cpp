#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;

#include "Box.H"
#include "Vector.H"
#include "LayoutIterator.H"
#include "parstream.H"

#include "MappedAMRLevel.H"


int MappedAMRLevel::s_verbosity = 0;

//-----------------------------------------------------------------------
bool MappedAMRLevel::isDefined() const
{
    return m_isDefined;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
MappedAMRLevel::~MappedAMRLevel()
{
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
MappedAMRLevel::MappedAMRLevel()
{
    m_coarser_level_ptr = NULL;
    m_finer_level_ptr = NULL;
    m_isDefined = false;
    m_level = 0;
    m_time = 0;
    m_dt = 0;
    m_initial_dt_multiplier = 0.1;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void MappedAMRLevel::define(MappedAMRLevel*  a_coarser_level_ptr,
                            const Box&       a_problem_domain,
                            int              a_level,
                            const IntVect&   a_ref_ratio)
{
    ProblemDomain physDomain(a_problem_domain);
    define(a_coarser_level_ptr, physDomain, a_level, a_ref_ratio);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void MappedAMRLevel::define(MappedAMRLevel*      a_coarser_level_ptr,
                            const ProblemDomain& a_problem_domain,
                            int                  a_level,
                            const IntVect&       a_ref_ratio)
{
    if (s_verbosity >= 3) {
        pout() << "MappedAMRLevel::define" << endl;
    }

    m_coarser_level_ptr = a_coarser_level_ptr;
    m_problem_domain = a_problem_domain;
    m_level = a_level;
    m_ref_ratio = a_ref_ratio;
    m_finer_level_ptr = NULL;
    m_isDefined = true;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void MappedAMRLevel::finerLevelPtr(MappedAMRLevel* a_finer_level_ptr)
{
    m_finer_level_ptr = a_finer_level_ptr;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void MappedAMRLevel::dt(Real a_dt)
{
    m_dt = a_dt;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real MappedAMRLevel::dt() const
{
    return (m_dt);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
const ProblemDomain& MappedAMRLevel::problemDomain() const
{
    return (m_problem_domain);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Vector<Box> MappedAMRLevel::boxes() const
{
    return (m_level_grids);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool MappedAMRLevel::hasCoarserLevel() const
{
    return ((m_coarser_level_ptr != NULL)
            && (m_coarser_level_ptr->m_isDefined)
            && (m_coarser_level_ptr->m_level_grids.size() > 0));
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
bool MappedAMRLevel::hasFinerLevel() const
{
    return ((m_finer_level_ptr != NULL)
            && (m_finer_level_ptr->m_isDefined)
            && (m_finer_level_ptr->m_level_grids.size() > 0));
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
int MappedAMRLevel::level() const
{
    return m_level;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
const IntVect& MappedAMRLevel::refRatio() const
{
    return (m_ref_ratio);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void MappedAMRLevel::time(Real a_time)
{
    m_time = a_time;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real MappedAMRLevel::time() const
{
    return m_time;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void MappedAMRLevel::initialDtMultiplier(Real a_initial_dt_multiplier)
{
    m_initial_dt_multiplier = a_initial_dt_multiplier;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Real MappedAMRLevel::initialDtMultiplier() const
{
    return m_initial_dt_multiplier;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// static
void MappedAMRLevel::verbosity(int a_verbosity)
{
    s_verbosity = a_verbosity;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
// static
int MappedAMRLevel::verbosity()
{
    return (s_verbosity);
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void MappedAMRLevel::preRegrid(int a_base_level, const Vector<Vector<Box> >& a_new_grids)
{
    if (s_verbosity >= 3) {
        pout() << "MappedAMRLevel::preRegrid" << endl;
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void MappedAMRLevel::postRegrid(int a_base_level)
{
    if (s_verbosity >= 3) {
        pout() << "MappedAMRLevel::postRegrid" << endl;
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void MappedAMRLevel::postInitialGrid(const bool a_restart)
{
    if (s_verbosity >= 3) {
        pout() << "MappedAMRLevel::postInitialGrid" << endl;
    }
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
Vector<MappedAMRLevel*> MappedAMRLevel::getAMRLevelHierarchy()
{
    Vector<MappedAMRLevel*> retval;
    // First go to level 0
    MappedAMRLevel* levelPtr = this;
    while (levelPtr->hasCoarserLevel()) {
        levelPtr = levelPtr->m_coarser_level_ptr;
    }

    // Now can accumulate the pointers by chasing finer level
    retval.push_back(levelPtr);
    while (levelPtr->hasFinerLevel()) {
        levelPtr = levelPtr->m_finer_level_ptr;
        retval.push_back(levelPtr);
    }

    return retval;
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MappedAMRLevel::writeCustomPlotFile(const std::string& a_prefix,
                                    int a_step) const
{
    // By default, this does nothing.
}
//-----------------------------------------------------------------------

//-----------------------------------------------------------------------
void
MappedAMRLevel::conclude(int a_step) const
{
    // By default, this does nothing.
}
//-----------------------------------------------------------------------

