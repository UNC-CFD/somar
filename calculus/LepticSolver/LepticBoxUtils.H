/*******************************************************************************
 *    FILE: LepticBoxUtils.H
 *    DESC: Utilities used by the leptic solver for creating grids.
 *    DATE: Tue 02 Jul 2013 03:52:00 PM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __LepticBoxUtils_H__INCLUDED__
#define __LepticBoxUtils_H__INCLUDED__

#include "Vector.H"
#include "Box.H"
#include "BoxLayout.H"


class LepticBoxUtils {
public:

    // Returns the largest dimensions exhibited in a vector of boxes.
    // This throws an error if the vector is empty or if any of the boxes are empty.
    static IntVect getMaxBoxSize (const Vector<Box>& a_boxes);

    // Create an array of grids suitable for the vertical solver.
    static void createVerticalSolverGrids (Vector<Box>&       a_vertBoxes,
                                           const Vector<Box>& a_origBoxes,
                                           const Box&         a_domBox);

    // Flattens grids that span the domain vertically and erases (makes empty)
    // those grids that do not. This preserves proc assignments, etc.
    class FlattenTransform;

    // Create an array of grids suitable for the horizontal solver.
    static void createHorizontalSolverGrids (Vector<Box>&       a_horizBoxes,
                                             const Vector<Box>& a_origBoxes,
                                             const Box&         a_origDomBox,
                                             const int          a_blockFactor);


    // Vertically flattens a set of grids. This assumes all boxes in the array
    // have the same type.
    static void createLevelGrids (Vector<Box>&       a_flatBoxes,       // Out: For the horizontal solver
                                  Vector<Box>&       a_inflatedBoxes,   // Out: For the vertical solver
                                  Vector<Box>&       a_preInflatedBoxes,// Out: For the excess function
                                  const Vector<Box>& a_origBoxes,       // In: The initial grids that aren't suitable for leptic solves.
                                  const Box&         a_domainBox,
                                  const int          a_blockFactor);

    // Attempts to reduce the number of level grids.
    // a_fndVector is a flat, non-disjoint vector of pre-joined boxes that serve
    // as a set of suggestions for the joiner. This vector can be empty.
    // The elements of a_joinDirs can be 0 (skip dir) or 1 (join dir).
    static void joinLevelGrids (std::list<Box>&    a_joinList,
                                const IntVect&     a_joinDirs,
                                const Vector<Box>& a_fndVector = Vector<Box>());

    // Vertically unflattens a set of grids.
    //
    // a_inflatedBoxes (OUT): spans a_originalBoxes, but its boxes are not split in
    //      the vertical.
    // a_flatBoxList (IN): The disjoint, flattened version of a_originalBoxes. All
    //      of these boxes should be square with a side length of the block factor.
    // a_originalBoxes (IN): The full set of disjoint grids. These grids are not
    //      suitable for leptic solves because they can be broken in the vertical.
    static void inflateLevelGrids (std::list<Box>&       a_inflatedList,
                                   const std::list<Box>& a_flatBoxList,
                                   const Vector<Box>&    a_originalBoxes);

    // For each element of a_destBoxArray, we assign the same processor that is
    // assigned to the corresponding box in a_srcAssignments. Correspondence is
    // assumed if the src box contains the dest box.
    static void copyProcArray (Vector<int>&                          a_destProcArray,
                               const Vector<Box>&                    a_destBoxArray,
                               const std::list<std::pair<int,Box> >& a_srcAssignments);

    // Checks if the given box vertically spans the domain.
    inline static bool vertSpanCheck (const Box&           a_box,
                                      const ProblemDomain& a_domain);

    // Checks if the given box vertically spans the domain.
    inline static bool vertSpanCheck (const Box& a_box,
                                      const Box& a_domBox);

    // Returns true if container contains element in the horizontal. Vertical
    // containment is ignored.
    inline static bool horizontallyContains (const Box& a_container,
                                             const Box& a_element);

    // I could make everything public since these are just stateless utilities,
    // but I want to distinguish these functions from what is meant for public
    // consumption.
private:

    // Creates a disjoint list of boxes using a binary (quadtree) search in the
    // horizontal domain. NOTE: a_boxArray must be flat!
    static void createDisjointBoxList (std::list<Box>&    a_disjointBoxList,
                                       const Vector<Box>& a_boxArray,
                                       const int          a_depth = 0);

    // Splits a_testBox into quadrants and throws away quadrants that are not
    // used in a_boxArray. Appends the result to a_quadrantList. If all
    // quadrants are used this function returns false.
    static bool resolveQuadrants (std::list<Box>&   a_quadrantList,
                                  const Box&        a_testBox,
                                  const Vector<Box> a_boxArray);

    // Returns true if any points in a_testBox lie in a_boxArray.
    static bool checkOverlap (const Box&        a_testBox,
                              const Vector<Box> a_boxArray);

    // Computes the smallest box that contains all boxes in a_boxArray.
    static Box minBox (const Vector<Box>& a_boxArray,
                       const bool         a_makeSquare = true);

    // Replaces an element in a list with the elements of another list.
    // a_list is the list being modified.
    // a_pos points to the element being removed.
    // a_replacementList's elements will be added to a_list.
    // Upon exit, a_pos will point to the first element inserted
    // and the return iterator will point to the last element inserted.
    template<class T>
    static typename std::list<T>::iterator
    replace (std::list<T>&                    a_list,
             typename std::list<T>::iterator& a_pos,
             const std::list<T>&              a_replacementList);

    // If the area covered by box1 and box2 can be exactly covered by a single
    // box, this single box is returned. Otherwise, the empty box is returned.
    static Box join (const Box& a_box1,
                     const Box& a_box2);

    // This is the same as join, except box1 and box2 are assumed to be disjoint
    // and non-empty. If the boxes can be joined, the result is placed in box1
    // and true is returned. If the boxes cannot be joined, box1 will be
    // minBox(box1, box2) and false is returned.
    static bool fastJoin (Box&       a_box1,
                          const Box& a_box2);

    // If extents of box1 in the dir direction can be adjusted so that the new
    // box1 covers the exact same area as the original box1 and box2, then this
    // new box1 is returned along with true. Otherwise, false will be returned
    // and box1 will remain untouched.
    //
    // The "fast" in the function name implies that you should not expect error
    // checking in this function. Both input boxes must be non-empty and disjoint.
    static bool fastJoin (Box&       a_box1,
                          const Box& a_box2,
                          const int  a_dir);
};


// -----------------------------------------------------------------------------
// Flattens grids that span the domain vertically and erases (makes empty)
// those grids that do not. This preserves proc assignments, etc.
// -----------------------------------------------------------------------------
class LepticBoxUtils::FlattenTransform: public BaseTransform
{
public:
    // a_pos is the desired vertical index of the flattened grids.
    FlattenTransform (const int a_pos)
    : m_pos(a_pos)
    {;}

    virtual ~FlattenTransform ()
    {;}

    virtual Box operator() (const Box& a_inputBox);

protected:
    int m_pos;
};



// -----------------------------------------------------------------------------
// Checks if the given box vertically spans the domain.
// -----------------------------------------------------------------------------
bool LepticBoxUtils::vertSpanCheck (const Box&           a_box,
                                    const ProblemDomain& a_domain)
{
    return LepticBoxUtils::vertSpanCheck(a_box, a_domain.domainBox());
}


// -----------------------------------------------------------------------------
// Checks if the given box vertically spans the domain.
// -----------------------------------------------------------------------------
bool LepticBoxUtils::vertSpanCheck (const Box& a_box,
                                    const Box& a_domBox)
{
    const int vdir = SpaceDim-1;
    if (a_box.smallEnd(vdir) != a_domBox.smallEnd(vdir)) return false;
    if (a_box.  bigEnd(vdir) != a_domBox.  bigEnd(vdir)) return false;
    return true;
}


// -----------------------------------------------------------------------------
// Returns true if container contains element in the horizontal. Vertical
// containment is ignored.
// -----------------------------------------------------------------------------
bool LepticBoxUtils::horizontallyContains (const Box& a_container,
                                           const Box& a_element)
{
    CH_assert(a_container.sameType(a_element));
    CH_assert(!a_container.isEmpty());
    CH_assert(!a_element.isEmpty());
    return D_TERM(,
                  (a_container.smallEnd(0) <= a_element.smallEnd(0) && a_element.bigEnd(0) <= a_container.bigEnd(0)), &&
                  (a_container.smallEnd(1) <= a_element.smallEnd(1) && a_element.bigEnd(1) <= a_container.bigEnd(1)));
}



#endif //!__LepticBoxUtils_H__INCLUDED__
