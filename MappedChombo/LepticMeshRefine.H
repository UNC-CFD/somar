/*******************************************************************************
 *    FILE: LepticMeshRefine.H
 *    DESC: Override of Chombo's BRMeshRefine for leptic grids.
 *    NOTE: Much of this code was taken from Chombo's MeshRefine/BRMeshRefine.
 *    DATE: Mon 27 May 2013 02:03:01 PM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __LepticMeshRefine_H__INCLUDED__
#define __LepticMeshRefine_H__INCLUDED__

#include "BRMeshRefine.H"


// Override of Chombo's BRMeshRefine for leptic grids.
class LepticMeshRefine: public BRMeshRefine
{
public:
    // Default constructor -- leaves object in an unusable state
    LepticMeshRefine ();

    // Full constructor -- leaves object in usable state
    LepticMeshRefine (const Box&             a_baseDomain,  // Level 0 domain
                      const Vector<IntVect>& a_refRatios,   // Refinement ratios -- refRatio[0] is btwn levels 0 and 1
                      const Real             a_fillRatio,   // Measure of how efficiently tagged cells will be covered
                      const int              a_blockFactor, // Amount by which grids are guaranteed to be coarsenable
                      const int              a_bufferSize,  // Proper nesting buffer amount
                      const IntVect&         a_maxSize,     // Maximum grid length in any direction -- 0 means no limit.
                      const IntVect&         a_spanDirs = IntVect::Zero);   // Set to 1 for new boxes to span the dim.
    // Full constructor -- leaves object in usable state
    LepticMeshRefine (const ProblemDomain&   a_baseDomain,  // Level 0 domain
                      const Vector<IntVect>& a_refRatios,   // Refinement ratios -- refRatio[0] is btwn levels 0 and 1
                      const Real             a_fillRatio,   // Measure of how efficiently tagged cells will be covered
                      const int              a_blockFactor, // Amount by which grids are guaranteed to be coarsenable
                      const int              a_bufferSize,  // Proper nesting buffer amount
                      const IntVect&         a_maxSize,     // Maximum grid length in any direction -- 0 means no limit.
                      const IntVect&         a_spanDirs = IntVect::Zero);   // Set to 1 for new boxes to span the dim.

    // Destructor
    virtual ~LepticMeshRefine ();

    // Define function -- size of RefRatios will define maximum number of levels
    void define (const Box&             a_baseDomain,       // Level 0 domain
                 const Vector<IntVect>& a_refRatios,        // Refinement ratios -- refRatio[0] is btwn levels 0 and 1
                 const Real             a_fillRatio,        // Measure of how efficiently tagged cells will be covered
                 const int              a_blockFactor,      // Amount by which grids are guaranteed to be coarsenable
                 const int              a_bufferSize,       // Proper nesting buffer amount
                 const IntVect&         a_maxSize,          // Maximum grid length in any direction -- 0 means no limit.
                 const IntVect&         a_spanDirs = IntVect::Zero);   // Set to 1 for new boxes to span the dim.

    // Define function -- size of RefRatios will define maximum number of levels
    // This is the function that all other full constructors call.
    void define (const ProblemDomain&   a_baseDomain,       // Level 0 domain
                 const Vector<IntVect>& a_refRatios,        // Refinement ratios -- refRatio[0] is btwn levels 0 and 1
                 const Real             a_fillRatio,        // Measure of how efficiently tagged cells will be covered
                 const int              a_blockFactor,      // Amount by which grids are guaranteed to be coarsenable
                 const int              a_bufferSize,       // Proper nesting buffer amount
                 const IntVect&         a_maxSize,          // Maximum grid length in any direction -- 0 means no limit.
                 const IntVect&         a_spanDirs = IntVect::Zero);   // Set to 1 for new boxes to span the dim.

    // Create hierarchy of grids from a single level of tags.
    //  This function creates a hierarchy of grids from a single level of
    //  tags on BaseLevel. If tags exist, then all levels will have grids.
    //  Returns the new finest level of grids.
    virtual int regrid (Vector<Vector<Box> >&       a_newMeshes,  // new set of grids at every level
                        const IntVectSet&           a_tags,       // tagged cells on baseLevel
                        const int                   a_baseLevel,  // index of base mesh level (finest unchanged level)
                        const int                   a_topLevel,   // top level to refine (one less than finest possible level)
                        const Vector<Vector<Box> >& a_oldMeshes); // existing grids (if no previous grids, set to domains)

    // Create hierarchy of grids from tags at all levels.
    //  This function creates a hierarchy of grids from tags at all
    //  refinement levels.  It is possible that not all levels will
    //  return with grids, since there may not be tags at all levels.
    //  Returns the new finest level of grids.
    virtual int regrid (Vector<Vector<Box> >&       a_newMeshes,  // new set of grids at every level
                        Vector<IntVectSet>&         a_tags,       // tagged cells on each existing level
                        const int                   a_baseLevel,  // index of base mesh level (finest unchanged level)
                        const int                   a_topLevel,   // top level to refine (one less than finest possible level)
                        const Vector<Vector<Box> >& a_oldMeshes); // existing grids (if no previous grids, set to domains)

    // Gets/sets maximum box size in any dimension
    inline virtual const IntVect& getMaxSize () const;
    inline virtual int getMaxSize (const int a_dir) const;
    inline virtual void setMaxSize (const IntVect& a_maxSize);
    inline virtual void setMaxSize (const int a_dir,
                                    const int a_maxSize);

    // Gets/sets the refinement ratios
    inline virtual const Vector<IntVect>& getRefRatios () const;
    inline virtual void setRefRatios (const Vector<IntVect>& a_refRatios);

    // Splits domain into vector of disjoint boxes with max size maxsize.
    // This version does not split the domain in planes perpendicular to the
    // vertical. This means the resulting grids will be suitable for leptic solves.
    // If a_maxBoxSize[i] == 0, the domain will not be split in the ith direction.
    static void domainSplit (const ProblemDomain& a_domain,
                             Vector<Box>&         a_vbox,
                             const IntVect&       a_maxBoxSize,
                             int                  a_blockFactor);

    // Splits domain into vector of disjoint boxes with max size maxsize.
    // This version does not split the domain in planes perpendicular to the
    // vertical. This means the resulting grids will be suitable for leptic solves.
    // If a_maxBoxSize[i] == 0, the domain will not be split in the ith direction.
    static void domainSplit (const Box&     a_domain,
                             Vector<Box>&   a_vbox,
                             const IntVect& a_maxBoxSize,
                             int            a_blockFactor);

    // Constructs a set of boxes which covers a set of tagged cells.
    virtual void makeBoxes (Vector<Box>&         a_mesh,        // Output: refined boxes at each new level
                            const IntVectSet&    a_tags,        // Input: set of tagged cells to cover
                            const IntVectSet&    a_pnd,         // Input: proper nesting domain in which mesh boxes must live
                            const ProblemDomain& a_domain,      // Input: physical domain
                            const IntVect&       a_maxSize,     // Input: largest number of cells in any dimension for any box
                            const IntVect&       a_totalBufferSize) const;

    // Free the buffer needed by receiveBoxesParallel. Call this before the program
    // ends to solve the BRMeshRefine leak problem. Calling this in serial does nothing.
    static void deleteBuffer ();

protected:
    // Constructs a set of boxes which covers a set of tagged cells.
    virtual void makeBoxes (std::list<Box>&      a_mesh,        // Output: refined boxes at each new level
                            IntVectSet&          a_tags,        // Input: set of tagged cells to cover
                            const IntVectSet&    a_pnd,         // Input: proper nesting domain in which mesh boxes must live
                            const ProblemDomain& a_domain,      // Input: physical domain
                            const IntVect&       a_maxSize,     // Input: largest number of cells in any dimension for any box
                            const int            a_depth,       // Input: depth of this recursion in the algorithm
                            const IntVect&       a_totalBufferSize) const;

    // Does the same thing as makeBoxes, but across multiple processors.
    virtual void makeBoxesParallel (std::list<Box>&      a_mesh,
                                    IntVectSet&          a_tags,
                                    const IntVectSet&    a_pnd,
                                    const ProblemDomain& a_domain,
                                    const IntVect&       a_maxBoxSize,
                                    const int            a_depth,
                                    const IntVect&       a_totalBufferSize,
                                    const int            a_minSize,
                                    const Interval&      a_procInterval) const;

    // Simply checks if a_pnd contains points that are properly nested in a_box.
    virtual bool properlyNested (const Box&           a_box,
                                 const ProblemDomain& a_domain,
                                 const IntVectSet&    a_pnd,
                                 const IntVect&       a_totalBuffer) const;

    // Compute the traces (signatures) of the minbox in each direction, and
    // find a hole in the trace (zero value) and an inflection point (zero
    // Laplacian) in each direction; keep the best of each.
    virtual void splitTagsInBestDimension (IntVectSet&    a_tags_inout_lo,
                                           IntVectSet&    a_tags_hi,
                                           const IntVect& a_maxBoxSize) const;

    // Takes an element of the a_boxes list and splits it as needed.
    virtual void splitBox (std::list<Box>&                 a_boxes,
                           const std::list<Box>::iterator& a_boxindex,
                           const IntVect&                  a_maxboxsize) const;

    // Computes local blockFactors used internally to enforce the BlockFactor.
    // This function computes values for m_local_blockfactors array, which is
    // the amount that tags on a level are coarsened in order to guarantee that
    // the grids on the next finer level are coarsenable by the BlockFactor.
    virtual void computeLocalBlockFactors ();

    // Computes proper nesting domains.
    // This should only be called by refine. It assumes that everything has
    // already been coarsened by the local blocking factor.
    virtual void makePNDs (Vector<IntVectSet>&          a_pnds,
                           Vector<IntVect>&             a_totalBufferSize,
                           const int                    a_baseLevel,
                           const int                    a_topLevel,
                           const Vector<ProblemDomain>& a_domains,
                           const IntVectSet&            a_baseMesh,
                           const Vector<IntVect>&       a_bufferSize) const;

    // Computes proper nesting domains.
    // This should only be called by refine. It assumes that everything has
    // already been coarsened by the local blocking factor.
    virtual void makePNDs (Vector<IntVectSet>&          a_pnds,
                           Vector<IntVect>&             a_totalBufferSize,
                           const int                    a_baseLevel,
                           const int                    a_topLevel,
                           const Vector<ProblemDomain>& a_domains,
                           const Vector<Box>&           a_baseMesh,
                           const Vector<IntVect>&       a_bufferSize) const;

    // Recursive function to enforce max size of boxes in a given direction.
    static void breakBoxes (Vector<Box>& a_vboxin,
                            const int&   a_maxSize,
                            const int&   a_idir);

    // This streches boxes to conform to m_spanDir and fixes any overlaps
    // that arise. This function returns true if boxes were removed or fixed to
    // avoid an overlap.
    static bool spanBoxes (Vector<Box>&   a_boxes,
                           const Box&     a_domBox,
                           const IntVect& a_spanDirs);

    // Checks if a_i is a power of 2.
    inline static bool isPower2 (const int a_i);

    // If a comp is 1, that dir will have its boxes span the dimension.
    IntVect         m_spanDirs;

    // Member variable overrides
    IntVect         m_maxSize;
    Vector<IntVect> m_nRefVect;
    Vector<IntVect> m_level_blockfactors;

    // BRMeshRefine does not delete its version of this, so we will keep
    // our own version of it and the user can call deleteBuffer when the
    // program ends.
    static int* s_recBuffer;

    // The same as the BRMeshRefine version, except this used s_recBuffer.
    void receiveBoxesParallel (const Interval& a_from,
                               const Interval& a_to,
                               std::list<Box>& a_mesh,
                               int             tag) const;

private:
    // Hidden functions
    void define (const Box&         a_baseDomain,
                 const Vector<int>& a_refRatios,
                 const Real         a_fillRatio,
                 const int          a_blockFactor,
                 const int          a_bufferSize,
                 const int          a_maxSize);

    void define (const ProblemDomain& a_baseDomain,
                 const Vector<int>&   a_refRatios,
                 const Real           a_fillRatio,
                 const int            a_blockFactor,
                 const int            a_bufferSize,
                 const int            a_maxSize);

    int maxSize () const;
    void maxSize (const int);

    const Vector<int>& refRatios () const;
    void refRatios (const Vector<int>&);

    virtual void makeBoxes (Vector<Box>&         a_mesh,
                            const IntVectSet&    a_tags,
                            const IntVectSet&    a_pnd,
                            const ProblemDomain& a_domain,
                            const int            a_maxSize,
                            const int            a_totalBufferSize) const
    {
        MayDay::Error("LepticMeshRefine::makeBoxes: This version should not be called");
    }

    void makeBoxes (std::list<Box>&      a_mesh,
                    IntVectSet&          a_tags,
                    const IntVectSet&    a_pnd,
                    const ProblemDomain& a_domain,
                    const int            a_maxSize,
                    const int            a_depth,
                    const int            a_totalBufferSize) const
    {
        MayDay::Error("LepticMeshRefine::makeBoxes: This version should not be called");
    }

    void makeBoxesParallel (std::list<Box>&      a_mesh,
                            IntVectSet&          a_tags,
                            const IntVectSet&    a_pnd,
                            const ProblemDomain& a_domain,
                            const int            a_maxBoxSize,
                            const int            a_depth,
                            const int            a_totalBufferSize,
                            const int            a_minSize,
                            const Interval&      a_procInterval) const;

    bool properlyNested (const Box&           a_box,
                         const ProblemDomain& a_domain,
                         const IntVectSet&    a_pnd,
                         int                  a_totalBuffer) const;

    void splitTagsInBestDimension (IntVectSet&    a_tags_inout_lo,
                                   IntVectSet&    a_tags_hi,
                                   const int a_maxBoxSize) const;

    void splitBox (std::list<Box>&                 a_boxes,
                   const std::list<Box>::iterator& a_boxindex,
                   const int                       a_maxboxsize) const;
};



// ********************** Inline function definitions **************************

// -----------------------------------------------------------------------------
// Returns maximum box size in any dimension
// -----------------------------------------------------------------------------
const IntVect& LepticMeshRefine::getMaxSize () const
{
    return m_maxSize;
}


// -----------------------------------------------------------------------------
// Returns maximum box size in any dimension
// -----------------------------------------------------------------------------
int LepticMeshRefine::getMaxSize (const int a_dir) const
{
    return m_maxSize[a_dir];
}


// -----------------------------------------------------------------------------
// Sets maximum box size in any dimension
// -----------------------------------------------------------------------------
void LepticMeshRefine::setMaxSize (const IntVect& a_maxSize)
{
    D_TERM(CH_assert(a_maxSize[0] >= 0);,
           CH_assert(a_maxSize[1] >= 0);,
           CH_assert(a_maxSize[2] >= 0);)

    m_maxSize = a_maxSize;
}


// -----------------------------------------------------------------------------
// Sets maximum box size in any dimension
// -----------------------------------------------------------------------------
void LepticMeshRefine::setMaxSize (const int a_dir,
                                   const int a_maxSize)
{
    CH_assert(0 <= a_dir);
    CH_assert(a_dir < SpaceDim);
    CH_assert(a_maxSize >= 0);

    m_maxSize[a_dir] = a_maxSize;
}


// -----------------------------------------------------------------------------
// Returns the refinement ratios
// -----------------------------------------------------------------------------
const Vector<IntVect>& LepticMeshRefine::getRefRatios () const
{
    return m_nRefVect;
}


// -----------------------------------------------------------------------------
// Sets the refinement ratios
// -----------------------------------------------------------------------------
void LepticMeshRefine::setRefRatios (const Vector<IntVect>& a_refRatios)
{
    m_nRefVect = a_refRatios;
    this->computeLocalBlockFactors();
}


// -----------------------------------------------------------------------------
// Static utility.
// Checks if a_i is a power of 2.
// -----------------------------------------------------------------------------
bool LepticMeshRefine::isPower2 (const int a_i)
{
    // Catch <= 0
    if (a_i <= 0) return false;

    // Find first bit
    unsigned i = a_i;
    while (!(i & 1)) {
        i >>= 1;
    }

    // Can't have any other bits
    return (i >> 1) ? false : true;
}


#endif //!__LepticMeshRefine_H__INCLUDED__
