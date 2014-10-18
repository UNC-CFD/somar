/*******************************************************************************
 *    FILE: FCDataFactory.H
 *    DESC: Create one comp of a LevelData<FluxBox> (useful for saving memory).
 *    DATE: Thu 27 Jun 2013 01:50:19 PM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __FCDataFactory_H__INCLUDED__
#define __FCDataFactory_H__INCLUDED__


// -----------------------------------------------------------------------------
class FCDataFactory: public DataFactory<FArrayBox>
{
public:
    // Constructor
    FCDataFactory (const int a_FCdir)
    : m_FCdir(a_FCdir)
    {
        CH_assert(0 <= m_FCdir);
        CH_assert(m_FCdir < SpaceDim);
    }

    // Factory function.
    // Creates a new FArray object and returns a pointer to it. Responsiblitly
    // for calling delete on this pointer is passed to the user.
    virtual FArrayBox* create (const Box&       a_box,
                               int              a_ncomps,
                               const DataIndex& a_datInd) const
    {
        CH_assert(0 < a_ncomps);

        Box FCregion = a_box;
        FCregion.surroundingNodes(m_FCdir);

        return new FArrayBox(FCregion, a_ncomps);
    }

protected:
    int m_FCdir;
};


// -----------------------------------------------------------------------------
class NodeDataFactory: public DataFactory<FArrayBox>
{
public:
    // Constructor
    NodeDataFactory ()
    {;}

    // Factory function.
    // Creates a new FArray object and returns a pointer to it. Responsiblitly
    // for calling delete on this pointer is passed to the user.
    virtual FArrayBox* create (const Box&       a_box,
                               int              a_ncomps,
                               const DataIndex& a_datInd) const
    {
        CH_assert(0 < a_ncomps);
        return new FArrayBox(surroundingNodes(a_box), a_ncomps);
    }
};


// -----------------------------------------------------------------------------
class FlatNodeDataFactory: public DataFactory<FArrayBox>
{
public:
    // Constructor
    FlatNodeDataFactory (const IntVect& a_flatDirs)
    : m_flatDirs(a_flatDirs)
    {
        D_TERM(CH_assert(m_flatDirs[0] == 0 || m_flatDirs[0] == 1);,
               CH_assert(m_flatDirs[1] == 0 || m_flatDirs[1] == 1);,
               CH_assert(m_flatDirs[2] == 0 || m_flatDirs[2] == 1);)
    }

    // Factory function.
    // Creates a new FArray object and returns a pointer to it. Responsiblitly
    // for calling delete on this pointer is passed to the user.
    virtual FArrayBox* create (const Box&       a_box,
                               int              a_ncomps,
                               const DataIndex& a_datInd) const
    {
        CH_assert(0 < a_ncomps);

        Box region = a_box;
        region.surroundingNodes();
        for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
            if (m_flatDirs[dir] != 0) {
                region.setBig(dir, region.smallEnd(dir));
                region.shift(dir, -region.smallEnd(dir));
            }
        }

        return new FArrayBox(region, a_ncomps);
    }

protected:
    IntVect m_flatDirs;
};


// -----------------------------------------------------------------------------
inline void aliasLevelData (LevelData<FArrayBox>& a_alias,
                            LevelData<FluxBox>*   a_original,
                            const Interval&       a_comps,
                            const int             a_FCdir)
{
    FABAliasFlBxDataFactory factory(a_original, a_comps, a_FCdir);
    a_alias.define(a_original->disjointBoxLayout(),
                   a_comps.size(),
                   a_original->ghostVect(),
                   factory);
}


#endif //!__FCDataFactory_H__INCLUDED__
