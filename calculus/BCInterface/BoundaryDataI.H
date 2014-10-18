/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Copyright (C) 2014 Edward Santilli & Alberto Scotti
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


// **************** The BoundaryData template implementation *******************

// -----------------------------------------------------------------------------
// Static member definitions
// -----------------------------------------------------------------------------
template <class T>
BaseFab<T> BoundaryData<T>::s_emptyFab;


// -----------------------------------------------------------------------------
// Default constructor (leaves unusable)
// -----------------------------------------------------------------------------
template <class T>
BoundaryData<T>::BoundaryData ()
{
    this->clear();
}


// -----------------------------------------------------------------------------
// Full constructor (calls define)
// -----------------------------------------------------------------------------
template <class T>
BoundaryData<T>::BoundaryData (const DisjointBoxLayout& a_grids,
                               const ProblemDomain&     a_domain,
                               const int                a_ncomp)
{
    this->define(a_grids, a_domain, a_ncomp);
}


// -----------------------------------------------------------------------------
// Copy constructor
// -----------------------------------------------------------------------------
template <class T>
BoundaryData<T>::BoundaryData (const BoundaryData<T>& a_src)
{
    this->deepCopy(a_src);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
template <class T>
BoundaryData<T>::~BoundaryData ()
{
    this->clear();
}


// -----------------------------------------------------------------------------
// Full define
// -----------------------------------------------------------------------------
template <class T>
void BoundaryData<T>::define (const DisjointBoxLayout& a_grids,
                              const ProblemDomain&     a_domain,
                              const int                a_ncomp)
{
    CH_assert(a_ncomp >= 1);

    // Remember the source layout so we can call check() later.
    m_grids = a_grids;

    // This object will not be flat until we call a function that flattens it.
    m_isFlat = false;

    const Box& domBox = a_domain.domainBox();
    DataIterator dit = a_grids.dataIterator();

    for (int dir = 0; dir < SpaceDim; ++dir) {
        // Only process non-periodic boundaries
        if (a_domain.isPeriodic()) continue;

        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide& iside = sit();
            const int s = int(iside);
            const Box domBdryFaces = bdryBox(domBox, dir, iside, 1);
            BoundaryDataMapType& thisMap = m_bdryFaceMap[dir][s];

            for (dit.reset(); dit.ok(); ++dit) {
                // Compute valid faces at the physical boundary.
                const Box& valid = a_grids[dit];
                Box validBdryFaces = bdryBox(valid, dir, iside, 1);
                validBdryFaces &= domBdryFaces;

                // If there are no cells at the boundary, just move on.
                if (validBdryFaces.isEmpty()) continue;

                // If we get here, then we have faces to store in the map.
                BaseFab<T>* newDataPtr = new BaseFab<T>(validBdryFaces, a_ncomp);
                thisMap[dit()] = DataPtrType(newDataPtr);

            } // end loop over grids (dit)
        } // end loop over sides (sit)
    } // end loop over directions (dir)
}


// -----------------------------------------------------------------------------
// Copy operator
// -----------------------------------------------------------------------------
template <class T>
BoundaryData<T>& BoundaryData<T>::operator= (const BoundaryData<T>& a_src)
{
    this->deepCopy(a_src);
    return *this;
}


// -----------------------------------------------------------------------------
// Copies a_src entirely without sharing pointers.
// NOTE: m_grids will still be a shallow copy.
// -----------------------------------------------------------------------------
template <class T>
void BoundaryData<T>::deepCopy (const BoundaryData<T>& a_src)
{
    this->clear();

    m_grids = a_src.m_grids;
    m_isFlat = a_src.m_isFlat;

    const ProblemDomain& domain = m_grids.physDomain();
    const Box& domBox = domain.domainBox();
    DataIterator dit = m_grids.dataIterator();

    for (int dir = 0; dir < SpaceDim; ++dir) {
        // Only process non-periodic boundaries
        if (domain.isPeriodic()) continue;

        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide& iside = sit();
            const int s = int(iside);

            const_iterator cit = a_src.begin(dir, iside);
            const_iterator citEnd = a_src.end(dir, iside);
            for (; cit != citEnd; ++cit) {
                const DataIndex srcDi = cit->first;
                const RefCountedPtr<BaseFab<T> >& srcFABPtr = cit->second;
                const Box& srcBox = srcFABPtr->box();
                const int srcComps = srcFABPtr->nComp();

                RefCountedPtr<BaseFab<T> > newFABPtr(new BaseFab<T>(srcBox, srcComps));
                m_bdryFaceMap[dir][s][srcDi] = newFABPtr;
            }
        } // end loop over sides (sit)
    } // end loop over directions (dir)
}


// -----------------------------------------------------------------------------
// Clears the std::maps and leaves this object unusable.
// -----------------------------------------------------------------------------
template <class T>
void BoundaryData<T>::clear ()
{
    // We no longer need the grids for index checking.
    m_grids = DisjointBoxLayout();

    m_isFlat = false;

    // Loop over directions and sides, freeing memory.
    for (int dir = 0; dir < SpaceDim; ++dir) {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const int s = int(sit());
            m_bdryFaceMap[dir][s].clear();
        }
    }
}


// -----------------------------------------------------------------------------
// Checks if a DataIndex is compatible with the grids that defined this object.
// -----------------------------------------------------------------------------
template <class T>
bool BoundaryData<T>::check (const DataIndex& a_di) const
{
    if (!m_grids.isClosed()) return false;
    return m_grids.check(a_di);
}


// -----------------------------------------------------------------------------
// Checks if the grids are compatible and data lies on the exact same boxes.
// -----------------------------------------------------------------------------
template <class T>
bool BoundaryData<T>::sameLayout (const BoundaryData<T>& a_src) const
{
    // Both grids need to be well-defined.
    if (!m_grids.isClosed()) return false;
    if (!a_src.m_grids.isClosed()) return false;

    // Are the grids compatible?
    if (!m_grids.compatible(a_src.m_grids)) return false;

    // Are they both flat / not flat?
    if (isFlat() != a_src.isFlat()) return false;

    // Loop over directions and side. Check if maps are identical.
    for (int dir = 0; dir < SpaceDim; ++dir) {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const int s = int(sit());
            const BoundaryDataMapType& thisMap = m_bdryFaceMap[dir][s];
            const BoundaryDataMapType& srcMap = a_src.m_bdryFaceMap[dir][s];

            if (thisMap.size() != srcMap.size()) return false;

            const_iterator thisIt = thisMap.begin();
            const_iterator srcIt = srcMap.begin();
            const const_iterator endIt = thisMap.end();

            while (thisIt != endIt) {
                const DataIndex& thisDit = thisIt->first;
                const DataPtrType& thisFABPtr = thisIt->second;

                const DataIndex& srcDit = srcIt->first;
                const DataPtrType& srcFABPtr = srcIt->second;

                if (thisDit != srcDit) return false;
                if (thisFABPtr->box() != srcFABPtr->box()) return false;

                ++thisIt;
                ++srcIt;
            }
        } // end loop over sides (sit)
    } // end loop over directions (dir)

    // If we get here, the layouts are the same.
    return true;
}


// -----------------------------------------------------------------------------
// Checks if a_src is a flattened version of this object.
// -----------------------------------------------------------------------------
template <class T>
bool BoundaryData<T>::flattenedLayout (const BoundaryData<T>& a_flatSrc) const
{
    // If this object is flat or the src is not flat, why would we be here?!
    CH_assert(!isFlat());
    CH_assert(a_flatSrc.isFlat());

    // Both grids need to be well-defined.
    if (!m_grids.isClosed()) return false;
    if (!a_flatSrc.m_grids.isClosed()) return false;

    // Are the grids compatible?
    if (!m_grids.compatible(a_flatSrc.m_grids)) return false;

    // Loop over directions and side. Check if maps are identical.
    // NOTE: This does not check the vertical boundaries.
    for (int dir = 0; dir < SpaceDim-1; ++dir) {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const int s = int(sit());
            const BoundaryDataMapType& thisMap = m_bdryFaceMap[dir][s];
            const BoundaryDataMapType& srcMap = a_flatSrc.m_bdryFaceMap[dir][s];

            if (thisMap.size() != srcMap.size()) return false;

            const_iterator thisIt = thisMap.begin();
            const_iterator srcIt = srcMap.begin();
            const const_iterator endIt = thisMap.end();

            while (thisIt != endIt) {
                const DataIndex& thisDit = thisIt->first;
                const DataIndex& srcDit = srcIt->first;
                if (thisDit != srcDit) return false;

                const Box& thisFlatBox = flattenBox(thisIt->second->box(), SpaceDim-1);
                const Box& srcBox = srcIt->second->box();
                if (thisFlatBox != srcBox) return false;

                ++thisIt;
                ++srcIt;
            }
        } // end loop over sides (sit)
    } // end loop over directions (dir)

    // If we get here, the layouts are the same.
    return true;
}


// -----------------------------------------------------------------------------
// Checks if this object can be flattened
// -----------------------------------------------------------------------------
template <class T>
bool BoundaryData<T>::flattenable () const
{
    CH_assert(m_grids.isClosed());
    CH_assert(!isFlat());

    const Box& domBox = m_grids.physDomain().domainBox();
    const int vertSize = domBox.size(SpaceDim-1);
    DataIterator dit = m_grids.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        if (m_grids[dit].size(SpaceDim-1) != vertSize) return false;
    }
    return true;
}


// -----------------------------------------------------------------------------
// Checks if this object is vertically flat.
// -----------------------------------------------------------------------------
template <class T>
bool BoundaryData<T>::isFlat () const
{
    return m_isFlat;
}


// -----------------------------------------------------------------------------
// Returns an iterator at the start of the appropriate map. (const version)
// -----------------------------------------------------------------------------
template <class T>
typename BoundaryData<T>::const_iterator
BoundaryData<T>::begin (const int             a_dir,
                        const Side::LoHiSide& a_side) const
{
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    CH_assert(!isFlat() || a_dir != SpaceDim-1);

    const int s = int(a_side);
    const_iterator it = m_bdryFaceMap[a_dir][s].begin();
    return it;
}


// -----------------------------------------------------------------------------
// Returns an iterator at the start of the appropriate map.
// -----------------------------------------------------------------------------
template <class T>
typename BoundaryData<T>::iterator
BoundaryData<T>::begin (const int             a_dir,
                        const Side::LoHiSide& a_side)
{
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    CH_assert(!isFlat() || a_dir != SpaceDim-1);

    const int s = int(a_side);
    iterator it = m_bdryFaceMap[a_dir][s].begin();
    return it;
}


// -----------------------------------------------------------------------------
// Returns an iterator at the end of the appropriate map. (const version)
// -----------------------------------------------------------------------------
template <class T>
typename BoundaryData<T>::const_iterator
BoundaryData<T>::end (const int             a_dir,
                      const Side::LoHiSide& a_side) const
{
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    CH_assert(!isFlat() || a_dir != SpaceDim-1);

    const int s = int(a_side);
    const_iterator it = m_bdryFaceMap[a_dir][s].end();
    return it;
}


// -----------------------------------------------------------------------------
// Returns an iterator at the end of the appropriate map.
// -----------------------------------------------------------------------------
template <class T>
typename BoundaryData<T>::iterator
BoundaryData<T>::end (const int             a_dir,
                      const Side::LoHiSide& a_side)
{
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    CH_assert(!isFlat() || a_dir != SpaceDim-1);

    const int s = int(a_side);
    iterator it = m_bdryFaceMap[a_dir][s].end();
    return it;
}


// -----------------------------------------------------------------------------
// Returns the appropriate FC boundary data holder. (const version)
// -----------------------------------------------------------------------------
template <class T>
const BaseFab<T>& BoundaryData<T>::getData (const DataIndex&      a_di,
                                            const int             a_dir,
                                            const Side::LoHiSide& a_side) const
{
    // Sanity checks
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    CH_assert(!isFlat() || a_dir != SpaceDim-1);
    CH_assert(this->check(a_di));

    // Grab a reference to the appropriate map.
    const int s = int(a_side);
    const BoundaryDataMapType& thisMap = m_bdryFaceMap[a_dir][s];

    // Find the DataIndex in the map.
    const_iterator it;
    it = thisMap.find(a_di);

    // If we didn't find anything, return the default empty holder.
    if (it == thisMap.end()) return s_emptyFab;

    // We found a holder at the boundary. Return it.
    const DataPtrType thisDataPtr = it->second;
    CH_assert(!thisDataPtr.isNull());
    return *thisDataPtr;
}


// -----------------------------------------------------------------------------
// Returns the appropriate FC boundary data holder.
// -----------------------------------------------------------------------------
template <class T>
BaseFab<T>& BoundaryData<T>::getData (const DataIndex&      a_di,
                                      const int             a_dir,
                                      const Side::LoHiSide& a_side)
{
    // Sanity checks
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    CH_assert(!isFlat() || a_dir != SpaceDim-1);
    CH_assert(this->check(a_di));

    // Grab a reference to the appropriate map.
    const int s = int(a_side);
    BoundaryDataMapType& thisMap = m_bdryFaceMap[a_dir][s];

    // Find the DataIndex in the map.
    iterator it;
    it = thisMap.find(a_di);

    // If we didn't find anything, return the default empty holder.
    if (it == thisMap.end()) return s_emptyFab;

    // We found a holder at the boundary. Return it.
    DataPtrType thisDataPtr = it->second;
    CH_assert(!thisDataPtr.isNull());
    return *thisDataPtr;
}


// -----------------------------------------------------------------------------
// Sets all data to a uniform value.
// -----------------------------------------------------------------------------
template <class T>
void BoundaryData<T>::setVal (T a_val)
{
    int numDirs = SpaceDim;
    if (isFlat()) --numDirs;

    for (int dir = 0; dir < numDirs; ++dir) {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const int s = int(sit());
            BoundaryDataMapType& thisMap = m_bdryFaceMap[dir][s];

            iterator it = thisMap.begin();
            const iterator endIt = thisMap.end();

            while (it != endIt) {
                it->second->setVal(a_val);
                ++it;
            }
        } // end loop over sides (sit)
    } // end loop over directions (dir)
}

