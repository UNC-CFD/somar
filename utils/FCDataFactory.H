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
#ifndef __FCDataFactory_H__INCLUDED__
#define __FCDataFactory_H__INCLUDED__


// -----------------------------------------------------------------------------
// Creates a face-centered data holder. This avoids the creation of a
// LevelData<FluxBox> when only one component is needed.
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
// Creates a node-centered data holder.
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
// Creates a node-centered data holder that has been projected (flattened) in
// the directions designated by a_flatDirs.
//
// If a_flatDirs[d] = 0, then the grids will not be flattened in direction d.
// If a_flatDirs[d] = 1, then the grids will be flattened in direction d.
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
// Allows us to alias just one face-centered component of a LevelData<FluxBox>.
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
