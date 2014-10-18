#include "MiscUtils.H"
#include "MergeBoxesOnLines.H"
#include "LepticMeshRefine.H"


// -----------------------------------------------------------------------------
// Takes a DBL and returns a new one whose boxes are merged in a_dir.
// -----------------------------------------------------------------------------
void mergeLayout (DisjointBoxLayout&       a_newLayout,
                  const DisjointBoxLayout& a_origLayout,
                  const int                a_dir)
{
    // Sanity checks
    CH_assert(0 <= a_dir);
    CH_assert(a_dir < SpaceDim);

    // Merge the boxes
    Vector<Box> vbox = a_origLayout.boxArray();
    MergeBoxesOnLines().mergeBoxes(vbox, a_dir);

    // Create the merged layout
    a_newLayout.defineAndLoadBalance(vbox, NULL, a_origLayout.physDomain());
}


// -----------------------------------------------------------------------------
// Splits a box into a load balanced set of boxes that are not split in a_dir.
// -----------------------------------------------------------------------------
void lineLayout (DisjointBoxLayout&   a_newLayout,
                 const ProblemDomain& a_domain,
                 const int            a_dir)
{
    // Sanity checks
    CH_assert(0 <= a_dir);
    CH_assert(a_dir < SpaceDim);

    const int np = numProc();
    IntVect maxBoxSize = a_domain.size() / np;
    maxBoxSize[a_dir] = 0;

    Vector<Box> vbox;
    LepticMeshRefine::domainSplit(a_domain.domainBox(), vbox, maxBoxSize, 1);

    a_newLayout.defineAndLoadBalance(vbox, NULL, a_domain);
}


// -----------------------------------------------------------------------------
// Define a dbl with just one box.
// This function must be called on ALL procs, but a_box only needs to be
// defined on a_srcProc.
// -----------------------------------------------------------------------------
void defineOneProcGrids (DisjointBoxLayout&   a_grids,
                         const ProblemDomain& a_domain,
                         Box                  a_box,
                         const int            a_srcProc)
{
    broadcast(a_box, a_srcProc);
    Vector<Box> boxArray(1, a_box);
    Vector<int> procArray(1, a_srcProc);
    a_grids.define(boxArray, procArray, a_domain);
}


// -----------------------------------------------------------------------------
// This will define a copier that does not care about valid vs invalid data -
// it will just copy everything.
// -----------------------------------------------------------------------------
void defineImpartialCopier (Copier&                  a_copier,
                            const DisjointBoxLayout& a_srcGrids,
                            const DisjointBoxLayout& a_destGrids,
                            const IntVect&           a_ghostVect,
                            const IntVect&           a_shift)
{
    BoxLayout srcLayout;
    srcLayout.deepCopy(a_srcGrids);
    srcLayout.grow(a_ghostVect);
    srcLayout.close();

    a_copier.define(srcLayout,
                    a_destGrids,
                    a_srcGrids.physDomain(),
                    a_ghostVect,
                    false,    // exchange copier?
                    a_shift);
}
