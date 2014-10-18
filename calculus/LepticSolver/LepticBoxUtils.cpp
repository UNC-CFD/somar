#include "LepticBoxUtils.H"
#include "AnisotropicRefinementTools.H"
#include "Box.H"
#include "CH_Timer.H"
#include "MergeBoxesOnLines.H"
#include "Subspace.H"


// -----------------------------------------------------------------------------
// Flattens all grids in a BoxLayout. This preserves proc assignments, but not
// disjointness.
//
// a_pos is the desired vertical index of the flattened grids.
// -----------------------------------------------------------------------------
Box LepticBoxUtils::FlattenTransform::operator() (const Box& a_inputBox)
{
    const int vdir = SpaceDim-1;

    Box retBox = a_inputBox;
    retBox.setBig(vdir, a_inputBox.smallEnd(vdir));
    retBox.shift(vdir, m_pos - a_inputBox.smallEnd(vdir));

    CH_assert(!retBox.isEmpty());
    return retBox;
}


// -----------------------------------------------------------------------------
// Returns the largest dimensions exhibited in a vector of boxes.
// This throws an error if the vector is empty or if any of the boxes are empty.
// -----------------------------------------------------------------------------
IntVect LepticBoxUtils::getMaxBoxSize (const Vector<Box>& a_boxes)
{
    const int numBoxes = a_boxes.size();
    CH_assert(numBoxes > 0);

    CH_assert(!a_boxes[0].isEmpty());
    IntVect maxBoxSize = a_boxes[0].size();

    for (int idx = 1; idx < numBoxes; ++idx) {
        CH_assert(!a_boxes[idx].isEmpty());
        const IntVect& boxSize = a_boxes[idx].size();

        D_TERM(
        maxBoxSize[0] = Max(maxBoxSize[0], boxSize[0]);,
        maxBoxSize[1] = Max(maxBoxSize[1], boxSize[1]);,
        maxBoxSize[2] = Max(maxBoxSize[2], boxSize[2]);)
    }

    return maxBoxSize;
}


// -----------------------------------------------------------------------------
// Create an array of grids suitable for the vertical solver.
// -----------------------------------------------------------------------------
void LepticBoxUtils::createVerticalSolverGrids (Vector<Box>&       a_vertBoxes,
                                                const Vector<Box>& a_origBoxes,
                                                const Box&         a_domBox)
{
    // If there are no boxes, we have nothing to do.
    if (a_origBoxes.size() == 0) return;

    // Create grids that are unsplit in the vertical.
    a_vertBoxes = a_origBoxes;
    MergeBoxesOnLines().mergeBoxes(a_vertBoxes, SpaceDim-1);

    // TODO: we should reorganize the covered area so that
    // load balancing will be easier.
}


// -----------------------------------------------------------------------------
// Create an array of grids suitable for the horizontal solver.
// -----------------------------------------------------------------------------
void LepticBoxUtils::createHorizontalSolverGrids (Vector<Box>&       a_horizBoxes,
                                                  const Vector<Box>& a_vertBoxes,
                                                  const Box&         a_vertDomBox,
                                                  const int          a_blockFactor)
{
    // For now, just do the naive thing and use flatBoxes
    // with the unused boxes removed.
    a_horizBoxes.clear();

    for (int idx = 0; idx < a_vertBoxes.size(); ++idx) {
        const Box& b = a_vertBoxes[idx];

        if (LepticBoxUtils::vertSpanCheck(b, a_vertDomBox)) {
            CH_assert(coarsenable(b, IntVect(D_DECL(a_blockFactor, a_blockFactor, a_blockFactor))));
            a_horizBoxes.push_back(flattenBox(b, SpaceDim-1));
        }
    }
}


// -----------------------------------------------------------------------------
// Vertically flattens a set of grids.
// -----------------------------------------------------------------------------
void LepticBoxUtils::createLevelGrids (Vector<Box>&       a_flatBoxes,         // Out: For the horizontal solver
                                       Vector<Box>&       a_inflatedBoxes,     // Out: For the vertical solver
                                       Vector<Box>&       a_preInflatedBoxes,  // Out: For the excess function.
                                       const Vector<Box>& a_origBoxes,         // In: The initial grids that aren't suitable for leptic solves.
                                       const Box&         a_domainBox,
                                       const int          a_blockFactor)
{
    CH_TIME("LepticBoxUtils::flattenLevelGrids");
    using std::list;
    using std::copy;

    // If there are no boxes, we have nothing to do.
    int numBoxes = a_origBoxes.size();
    if (numBoxes == 0) return;

    // We will work with coarsened boxes
    IntVect vBlockFactor = a_blockFactor * IntVect::Unit;

    Box coarsenedDomBox = a_domainBox;
    CH_assert(coarsenable(coarsenedDomBox, vBlockFactor));
    coarsenedDomBox.coarsen(vBlockFactor);

    a_inflatedBoxes.resize(numBoxes);
    for (int idx = 0; idx < numBoxes; ++idx) {
        Box& b = a_inflatedBoxes[idx];
        b.define(a_origBoxes[idx]);
        CH_assert(coarsenable(b, vBlockFactor));
        b.coarsen(vBlockFactor);
    }

    // Create grids that are unsplit in the vertical.
    MergeBoxesOnLines merger;
    merger.mergeBoxes(a_inflatedBoxes, SpaceDim-1);
    numBoxes = a_inflatedBoxes.size();

    // Flatten the merged boxes that span the domain vertically.
    // TODO: Should also check for Diri BCs.
    a_preInflatedBoxes.reserve(numBoxes);
    for (int idx = 0; idx < numBoxes; ++idx) {
        const int Nz = coarsenedDomBox.size(SpaceDim-1);
        if (a_inflatedBoxes[idx].size(SpaceDim-1) != Nz) continue;

        Box b = a_inflatedBoxes[idx];
        b.setBig(SpaceDim-1, b.smallEnd(SpaceDim-1));
        b.shift(SpaceDim-1, -b.smallEnd(SpaceDim-1));
        a_preInflatedBoxes.push_back(b);
    }

    // Attempt to reduce the number of boxes.
    IntVect joinDirs(D_DECL(1,1,1));
    joinDirs[SpaceDim-1] = 0;

    IntVect flatBlockFactor = vBlockFactor;
    flatBlockFactor[SpaceDim-1] = 1;

    list<Box> disjointBoxList;
    createDisjointBoxList(disjointBoxList, a_preInflatedBoxes);
    joinLevelGrids(disjointBoxList, joinDirs);

    // Refine boxes
    for (int idx = 0; idx < numBoxes; ++idx) {
        a_inflatedBoxes[idx].refine(vBlockFactor);
    }
    for (int idx = 0; idx < a_preInflatedBoxes.size(); ++idx) {
        a_preInflatedBoxes[idx].refine(vBlockFactor);
    }

    a_flatBoxes.resize(disjointBoxList.size());
    list<Box>::const_iterator dbit = disjointBoxList.begin();
    for (int idx = 0; dbit != disjointBoxList.end(); ++dbit, ++idx) {
        a_flatBoxes[idx] = refine(*dbit, flatBlockFactor);
    }
}


// -----------------------------------------------------------------------------
// Attempts to reduce the number of level grids.
// a_fndVector is a flat, non-disjoint vector of pre-joined boxes that serve
// as a set of suggestions for the joiner. This vector can be empty.
// The elements of a_joinDirs can be 0 (skip dir) or 1 (join dir).
//
// TODO: respect maxBoxSize.
// -----------------------------------------------------------------------------
void LepticBoxUtils::joinLevelGrids (std::list<Box>&    a_joinList,
                                     const IntVect&     a_joinDirs,
                                     const Vector<Box>& a_fndVector)
{
    CH_TIME("LepticBoxUtils::joinLevelGrids");

    typedef std::list<Box>::iterator listIterator;
    using std::back_inserter;
    using std::copy;
    using std::list;
    listIterator it1;

    // Do we have anything to do?
    if (a_joinList.size() == 0) return;

    // ***** Pass #1: Collect boxes from a_fndVector *****
    if (a_fndVector.size() > 0) {
        // Copy the flat non-disjoint boxes to a locally modifiable list.
        list<Box> fndList;
        const std::vector<Box>& fndStdVector = a_fndVector.constStdVector();
        copy(fndStdVector.begin(), fndStdVector.end(), back_inserter(fndList));

        for (it1 = a_joinList.begin(); it1 != a_joinList.end(); ++it1) {
            // Fill matchIts with iterators to elements of fndList containing *it1.
            std::vector<listIterator> matchIts;
            listIterator it2;
            for (it2 = fndList.begin(); it2 != fndList.end(); ++it2) {
                if (it2->contains(*it1)) {
                    matchIts.push_back(it2);
                }
            }

            // If there are no matches, just move on. We may have already
            // removed anything that would have matched.
            const int numMatches = matchIts.size();
            if (numMatches == 0) continue;

            // Choose one of the matches.
            Box& curWinnerRef = *(matchIts[0]);
            for (int idx = 1; idx < numMatches; ++idx) {
                Box& testRef = *(matchIts[idx]);

                // Which box covers more area?
                if (testRef.numPts() > curWinnerRef.numPts()) {
                    curWinnerRef = testRef;
                }
            }
            const Box winner = curWinnerRef;

            // Remove all entries of fndList that are in matchIts.
            for (int idx = 0; idx < numMatches; ++idx) {
                fndList.erase(matchIts[idx]);
            }
            matchIts.resize(0);

            // Remove all elements of fndList that intersect winner.
            for (it2 = fndList.begin(); it2 != fndList.end(); ++it2) {
                if (it2->intersects(winner)) {
                    it2 = fndList.erase(it2);
                    --it2;
                }
            }

            // Remove all boxes in a_joinList that are contained in winner.
            bool isFirst = true;
            for (it2 = a_joinList.begin(); it2 != a_joinList.end(); ++it2) {
                CH_assert(winner.contains(*it2) || !winner.intersects(*it2));
                if (winner.contains(*it2)) {
                    it2 = a_joinList.erase(it2);
                    --it2;

                    // Reset the iterator since we clobbered a_joinList.
                    if (isFirst) {
                        isFirst = false;
                        it1 = it2;
                    }
                }
            }

            // Add the lucky winner to the end of a_joinList.
            a_joinList.push_back(winner);
        }
    }

    // ***** Pass #2: Try to join neighboring boxes *****
    int changeCount;
    do {
        changeCount = 0;
        for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
            if (a_joinDirs[dir] == 0) continue;

            for (it1 = a_joinList.begin(); it1 != a_joinList.end(); ++it1) {
                listIterator it2 = it1;
                ++it2;

                for (; it2 != a_joinList.end(); ++it2) {
                    if (fastJoin(*it1, *it2, dir)) {
                        // it1 is still valid and points to the joined box.
                        // Just remove *it2 and move on the next *it1.
                        a_joinList.erase(it2);
                        ++changeCount;
                        break;
                    }
                } // end loop over a_joinList (it2 > it1)
            } // end loop over a_joinList (it1)
        } // end loop over join directions (dir)
    } while (changeCount > 0);
}


// -----------------------------------------------------------------------------
// Vertically unflattens a set of grids.
//
// a_inflatedBoxes (OUT): spans a_originalBoxes, but its boxes are not split in
//      the vertical.
// a_flatBoxList (IN): The disjoint, flattened version of a_originalBoxes. All
//      of these boxes should be square with a side length of the block factor.
// a_originalBoxes (IN): The full set of disjoint grids. These grids are not
//      suitable for leptic solves because they can be broken in the vertical.
// -----------------------------------------------------------------------------
void LepticBoxUtils::inflateLevelGrids (std::list<Box>&       a_inflatedList,
                                        const std::list<Box>& a_flatBoxList,
                                        const Vector<Box>&    a_originalBoxes)
{
    using std::list;
    const int numOriginalBoxes = a_originalBoxes.size();

    // Iterate over the flat boxes, inflating as we go.
    list<Box>::const_iterator flatIt = a_flatBoxList.begin();
    for (; flatIt != a_flatBoxList.end(); ++flatIt) {
        const IntVect smallEnd = flatIt->smallEnd();
        const IntVect bigEnd = flatIt->bigEnd();

        // Create a list of all originalBoxes that lie over *flatIt.
        list<Box> hoverList;
        for (int idx = 0; idx < numOriginalBoxes; ++idx) {
            const Box& ob = a_originalBoxes[idx];

            // Does ob contain *flatIt? If so, add ob to the list.
            if (D_TERM(,
                smallEnd[0] >= ob.smallEnd(0) && bigEnd[0] <= ob.bigEnd(0), &&
                smallEnd[1] >= ob.smallEnd(1) && bigEnd[1] <= ob.bigEnd(1))) {

                hoverList.push_back(ob);
            }
        }

        // Create boxes that span *flatIt horizontally and spans each box in
        // hoverList vertically. Put these boxes in inflatedList.
        list<Box>::const_iterator hovIt = hoverList.begin();
        for (; hovIt != hoverList.end(); ++hovIt) {
            Box newBox(smallEnd, bigEnd);
            newBox.setBig(CH_SPACEDIM-1, hovIt->bigEnd(CH_SPACEDIM-1));
            newBox.setSmall(CH_SPACEDIM-1, hovIt->smallEnd(CH_SPACEDIM-1));

            a_inflatedList.push_back(newBox);
        }
    }

    // Join as many boxes as possible in the vertical.
    const IntVect joinDirs = BASISV(CH_SPACEDIM-1);
    joinLevelGrids(a_inflatedList, joinDirs);
}


// -----------------------------------------------------------------------------
// For each element of a_destBoxArray, we assign the same processor that is
// assigned to the corresponding box in a_srcAssignments. Correspondence is
// assumed if the src box contains the dest box.
// -----------------------------------------------------------------------------
void LepticBoxUtils::copyProcArray (Vector<int>&                          a_destProcArray,
                                    const Vector<Box>&                    a_destBoxArray,
                                    const std::list<std::pair<int,Box> >& a_srcAssignments)
{
    const int numDestBoxes = a_destBoxArray.size();
    a_destProcArray.resize(numDestBoxes, -1);

    list<pair<int,Box> >::const_iterator it = a_srcAssignments.begin();
    for (; it != a_srcAssignments.end(); ++it) {
        const Box& assignedBox = it->second;

        for (int idx = 0; idx < numDestBoxes; ++idx) {
            const Box& unassignedBox = a_destBoxArray[idx];
            if (!assignedBox.contains(unassignedBox)) continue;

            // Assign the proc
            a_destProcArray[idx] = it->first;
        } // end loop over assigned boxes (idx)
    } // end loop over a_srcAssignments (it)

#ifndef NDEBUG
    for (int idx = 0; idx < numDestBoxes; ++idx) {
        if (a_destProcArray[idx] < 0) {
            pout() << "LepticBoxUtils::copyProcArray: a_destProcArray = " << a_destProcArray << endl;
            MayDay::Error("a_destProcArray not properly assigned in LepticBoxUtils::copyProcArray. See pout");
        }
    }
#endif
}


// -----------------------------------------------------------------------------
// Given: (&disjointBoxList, boxVector, depth)
// If depth = 0:
//   Empty disjointBoxList.
//   Compute minBox of boxVector.
//   Stretch each side of minBox so it is a power of 2 and square.
//   Initialize disjointBoxList = {minBox}.
// Set crseBoxVector = boxVector coarsened by 2.
// Loop over each disjointBox in disjointBoxList:
//   If disjointBox is coarsenable by 4:
//     Set crseDBList = {disjointBox coarsened by 2}.
//     Recurse with (&crseDBList, crseBoxVector, depth+1).     When function returns, crseDBList will be the resolved set of coarsened boxes.
//     Refine crseDBList by 2.
//     Perform quadrantCheck on each crseDB in crseDBList -> resolvedDBList
//   else:
//     Perform quadrantCheck on disjointBox -> resolvedDBList
//   Replace disjointBox in list with resolvedDBList and set iterator to last new element.
//
// Creates a disjoint list of boxes using a binary (quadtree) search in the
// horizontal domain. NOTE: a_boxArray must be flat!
// -----------------------------------------------------------------------------
void LepticBoxUtils::createDisjointBoxList (std::list<Box>&    a_disjointBoxList,
                                            const Vector<Box>& a_boxArray,
                                            const int          a_depth)
{
    // I don't think a timer will work in a recursive function.

    typedef std::list<Box>::iterator listIterator;
    using std::list;

    if (a_depth == 0) {
        // Compute square minBox of boxVector. This will be flat.
        Box minSqBox = minBox(a_boxArray);

        // Initialize disjointBoxList = {minBox}.
        a_disjointBoxList.clear();
        a_disjointBoxList.push_back(minSqBox);
    }

    // Set crseBoxVector = boxVector coarsened by 2.
    int numBoxes = a_boxArray.size();
    Vector<Box> crseBoxVector = a_boxArray;
    for (int idx = 0; idx < numBoxes; ++idx) {
        crseBoxVector[idx].coarsen(2);
    }

    // Loop over each disjointBox in disjointBoxList and resolve.
    listIterator dbit;
    for (dbit = a_disjointBoxList.begin(); dbit != a_disjointBoxList.end(); ++dbit) {
        Box& disjointBox = *dbit;

        list<Box> resolvedDBList;

        IntVect ref(D_DECL(4,4,4));
        ref[CH_SPACEDIM-1] = 1;

        bool needsReplacement = false;

        // Create the resolved list of boxes over disjointBox.
        if (coarsenable(disjointBox, ref)) {
            // Set crseDBList = {disjointBox coarsened by 2}.
            ref = IntVect(D_DECL(2,2,2));
            ref[CH_SPACEDIM-1] = 1;

            list<Box> crseDBList;
            crseDBList.push_back(disjointBox);
            crseDBList.back().coarsen(ref);

            // Recurse with (&crseDBList, crseBoxVector, depth+1). When function returns, crseDBList will be the resolved set of coarsened boxes.
            createDisjointBoxList(crseDBList, crseBoxVector, a_depth + 1);

            // Refine crseDBList by 2.
            listIterator cdbit = crseDBList.begin();
            for (; cdbit != crseDBList.end(); ++cdbit) {
                cdbit->refine(ref);
            }

            // Perform quadrantCheck on each crseDB in crseDBList and append to resolvedDBList
            for (cdbit = crseDBList.begin(); cdbit != crseDBList.end(); ++cdbit) {
                needsReplacement = resolveQuadrants(resolvedDBList, *cdbit, a_boxArray);
            }
        } else {
            // Perform quadrantCheck on each crseDB in crseDBList and append to resolvedDBList
            needsReplacement = resolveQuadrants(resolvedDBList, disjointBox, a_boxArray);
        }

        // Replace disjointBox in list with resolvedDBList and set iterator to last new element.
        dbit = replace(a_disjointBoxList, dbit, resolvedDBList);

    } // end loop over disjointBoxList (dbit)
}


// -----------------------------------------------------------------------------
// Splits a_testBox into quadrants and throws away quadrants that are not
// used in a_boxArray. Appends the result to a_quadrantList. If all
// quadrants are used this function returns false.
// -----------------------------------------------------------------------------
bool LepticBoxUtils::resolveQuadrants (std::list<Box>&   a_quadrantList,
                                       const Box&        a_testBox,
                                       const Vector<Box> a_boxArray)
{
    CH_assert(SpaceDim == 2 || SpaceDim == 3);
    CH_assert(!a_testBox.isEmpty());
    CH_assert(a_boxArray.size() > 0);

    // Split the test box into quadrants.
    const int numBoxes = D_TERM(,2,*2);
    Vector<Box> quadrants(numBoxes, a_testBox);

    int choppt = a_testBox.smallEnd(0) + a_testBox.size(0) / 2;
    quadrants[1] = quadrants[0].chop(0, choppt);

#if CH_SPACEDIM == 3
    quadrants[2] = quadrants[0];
    quadrants[3] = quadrants[1];

    choppt = a_testBox.smallEnd(1) + a_testBox.size(1) / 2;
    quadrants[2] = quadrants[0].chop(1, choppt);
    quadrants[3] = quadrants[1].chop(1, choppt);
#endif

    // Loop over quadrants
    int numQuadrantsNeeded = 0;
    for (int idx = 0; idx < numBoxes; ++idx) {
        const Box& quadBox = quadrants[idx];
        CH_assert(!quadBox.isEmpty());

        // Check if quadBox covers cells used by a_boxArray
        bool overlap = checkOverlap(quadBox, a_boxArray);
        if (overlap) {
            a_quadrantList.push_back(quadBox);
            ++numQuadrantsNeeded;
        }
    }

    // If all quadrants were added to the list, return false.
    return (numQuadrantsNeeded != numBoxes);
}


// -----------------------------------------------------------------------------
// Returns true if any points in a_testBox lie in a_boxArray.
// -----------------------------------------------------------------------------
bool LepticBoxUtils::checkOverlap (const Box&        a_testBox,
                                   const Vector<Box> a_boxArray)
{
    const int numBoxes = a_boxArray.size();
    for (int idx = 0; idx < numBoxes; ++idx) {
        if (a_testBox.intersects(a_boxArray[idx])) return true;
    }
    return false;
}


// -----------------------------------------------------------------------------
// Computes the smallest box that contains all boxes in a_boxArray.
// -----------------------------------------------------------------------------
Box LepticBoxUtils::minBox (const Vector<Box>& a_boxArray,
                            const bool         a_makeSquare)
{
    // If there are no boxes to enclose, return an empty box.
    const int numBoxes = a_boxArray.size();
    if (numBoxes == 0) return Box();

    // Begin with the first box in the array. Then, loop over the array,
    // stretching retBox as we go.
    Box retBox = a_boxArray[0];
    for (int idx = 1; idx < numBoxes; ++idx) {
        retBox.minBox(a_boxArray[idx]);
    }

    // This code is for the leptic horizontal solver.
    CH_assert(retBox.size(CH_SPACEDIM-1) == 1);

    if (a_makeSquare) {
        // Stretch each side of minBox so it is a power of 2 and square.
        int maxLenPow2 = 0;
        for (int dir = 0; dir < CH_SPACEDIM-1; ++dir) {
            int len = retBox.size(dir);
            int lenPow2 = 1;
            while (lenPow2 < len) {
                lenPow2 *= 2;
            }
            maxLenPow2 = Max(maxLenPow2, lenPow2);
        }
        for (int dir = 0; dir < CH_SPACEDIM-1; ++dir) {
            retBox.setBig(dir, retBox.smallEnd(dir) + maxLenPow2 - 1);
        }
    }

    return retBox;
}


// -----------------------------------------------------------------------------
// Replaces an element in a list with the elements of another list.
//
// a_list is the list being modified.
// a_pos points to the element being removed.
// a_replacementList's elements will be added to a_list.
//
// Upon exit, a_pos will point to the first element inserted and the return
// iterator will point to the last element inserted. If the replacement list
// is empty, a_pos will be invalidated and the return iterator will point to
// the element before the removal point.
// -----------------------------------------------------------------------------
template<class T>
typename std::list<T>::iterator
LepticBoxUtils::replace (std::list<T>&                    a_list,
                         typename std::list<T>::iterator& a_pos,
                         const std::list<T>&              a_replacementList)
{
    using std::list;

    // Sanity checks
    CH_assert(a_list.size() > 0);
    CH_assert(a_pos != a_list.end());

    list<Box>::iterator listIt = a_pos;

    // Is there anything to insert?
    if (a_replacementList.size() == 0) {
        listIt = a_list.erase(a_pos);
        return --listIt;
    }

    // Instead of removing the old element then making N insertions,
    // redefine the old element to be the first replacement, then make
    // N-1 insertions.

    // Perform the redefinition.
    list<Box>::const_iterator replIt = a_replacementList.begin();
    *listIt = *replIt;

    // Increment iterators.
    // The replacement iterator should always point to the element being
    // inserted. The list iterator should always point to the element AFTER
    // the insertion.
    ++replIt;
    ++listIt;

    // Perform the insertions.
    for (; replIt != a_replacementList.end(); ++replIt) {
        a_list.insert(listIt, *replIt);
    }

    // Set the return iterators.
    --listIt;
    CH_assert(*a_pos == a_replacementList.front());
    CH_assert(*listIt == a_replacementList.back());

    // Done!
    return listIt;
}


// -----------------------------------------------------------------------------
// If the area covered by box1 and box2 can be exactly covered by a single
// box, this single box is returned. Otherwise, the empty box is returned.
// -----------------------------------------------------------------------------
Box LepticBoxUtils::join (const Box& a_box1,
                          const Box& a_box2)
{
    Box retBox = a_box1;
    retBox.minBox(a_box2);

    int sumPts = a_box1.numPts() + a_box2.numPts();
    if (a_box1.intersects(a_box2)) {
        Box overlap = a_box1;
        overlap &= a_box2;
        sumPts -= overlap.numPts();
    }

    if (retBox.numPts() == sumPts) {
        return retBox;
    }

    return Box();
}


// -----------------------------------------------------------------------------
// This is the same as join, except box1 and box2 are assumed to be disjoint
// and non-empty. If the boxes can be joined, the result is placed in box1
// and true is returned. If the boxes cannot be joined, box1 will be
// minBox(box1, box2) and false is returned.
// -----------------------------------------------------------------------------
bool LepticBoxUtils::fastJoin (Box&       a_box1,
                               const Box& a_box2)
{
    const int sumPts = a_box1.numPts() + a_box2.numPts();
    a_box1.minBox(a_box2);
    return (a_box1.numPts() == sumPts);
}


// -----------------------------------------------------------------------------
// If extents of box1 in the dir direction can be adjusted so that the new
// box1 covers the exact same area as the original box1 and box2, then this
// new box1 is returned along with true. Otherwise, false will be returned
// and box1 will remain untouched.
//
// The "fast" in the function name implies that you should not expect error
// checking in this function. Both input boxes must be non-empty and disjoint.
// -----------------------------------------------------------------------------
bool LepticBoxUtils::fastJoin (Box&       a_box1,
                               const Box& a_box2,
                               const int  a_dir)
{
    int transDir = (a_dir + 1) % CH_SPACEDIM;
    if (   a_box1.smallEnd(transDir) != a_box2.smallEnd(transDir)
        || a_box1  .bigEnd(transDir) != a_box2  .bigEnd(transDir)   )
        return false;

#if CH_SPACEDIM == 3
    transDir = (a_dir + 2) % CH_SPACEDIM;
    if (   a_box1.smallEnd(transDir) != a_box2.smallEnd(transDir)
        || a_box1  .bigEnd(transDir) != a_box2  .bigEnd(transDir)   )
        return false;
#endif

    if (a_box1.bigEnd(a_dir) + 1 == a_box2.smallEnd(a_dir)) {
        a_box1.setBig(a_dir, a_box2.bigEnd(a_dir));
        return true;
    } else if (a_box2.bigEnd(a_dir) + 1 == a_box1.smallEnd(a_dir)) {
        a_box1.setSmall(a_dir, a_box2.smallEnd(a_dir));
        return true;
    }

    return false;
}
