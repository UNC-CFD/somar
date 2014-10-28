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
#include "Subspace.H"


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

