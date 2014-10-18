#include "Subspace.H"
#include "UsingNamespace.H"


// -----------------------------------------------------------------------------
// Takes an N-dimensional Box and projects it to a subspace defined by a_mask.
// a_mask should be 0 in flattened directions and 1 in unmodified directions.
// In other words, a_mask flags the tangential directions of the subspace.
// -----------------------------------------------------------------------------
Box flattenBox (const Box&     a_box,
                const IntVect& a_mask)
{
    D_TERM(CH_assert(a_mask[0] == 0 || a_mask[0] == 1);,
           CH_assert(a_mask[1] == 0 || a_mask[1] == 1);,
           CH_assert(a_mask[2] == 0 || a_mask[2] == 1);)

    const IntVect boxType = a_box.type();
    const IntVect smallIV = a_mask * a_box.smallEnd();
    const IntVect bigIV   = a_mask * a_box.bigEnd();

    Box retBox(smallIV, bigIV, boxType);
    return retBox;
}


// -----------------------------------------------------------------------------
// Takes an N-dimensional Box and projects it to the (N_1)-dimensional surface
// whose normal direction is a_normDir.
// -----------------------------------------------------------------------------
Box flattenBox (const Box& a_box,
                const int  a_normDir)
{
    CH_assert(0 <= a_normDir);
    CH_assert(a_normDir < SpaceDim);

    IntVect mask = IntVect::Unit - BASISV(a_normDir);
    return flattenBox(a_box, mask);
}


// -----------------------------------------------------------------------------
// Returns a box that can define holders for functions that only depend on the
// vertical coordinate.
// -----------------------------------------------------------------------------
Box verticalDataBox (const ProblemDomain& a_domain)
{
    return flattenBox(a_domain.domainBox(), BASISV(SpaceDim-1));
}


// -----------------------------------------------------------------------------
// Returns a box that can define holders for functions that only depend on the
// vertical coordinate.
// -----------------------------------------------------------------------------
Box verticalDataBox (const Box& a_box)
{
    return flattenBox(a_box, BASISV(SpaceDim-1));
}


// -----------------------------------------------------------------------------
// Returns a box that can define holders for functions that only depend on the
// horizontal coordinates.
// -----------------------------------------------------------------------------
Box horizontalDataBox (const ProblemDomain& a_domain)
{
    return flattenBox(a_domain.domainBox(), SpaceDim-1);
}


// -----------------------------------------------------------------------------
// Returns a box that can define holders for functions that only depend on the
// horizontal coordinates.
// -----------------------------------------------------------------------------
Box horizontalDataBox (const Box& a_box)
{
    return flattenBox(a_box, SpaceDim-1);
}

