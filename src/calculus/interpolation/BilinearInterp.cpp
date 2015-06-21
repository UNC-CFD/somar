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
#include "BilinearInterp.H"
#include "BilinearInterpF_F.H"


// -----------------------------------------------------------------------------
// The nodes will be labeled as A, B, C, and D in this order:   C---D
//                                                              |   |
//                                                              A---B
// Simple bilinear interp where
// f_{u,v}=fA*(1-u)*(1-v)+FB*u*(1-v)+FC*(1-v)*u+FD*u*v
// where (u,v) are the coordinate in the middle of the box
// -----------------------------------------------------------------------------
void BilinearInterp2D (FArrayBox&          a_fInterp,
                       const FArrayBox&    a_xInterp,
                       const FArrayBox&    a_yInterp,
                       const Box&          a_interpBox,
                       const int           a_xdir,
                       const int           a_ydir,
                       const Vector<Real>& a_x,
                       const Vector<Real>& a_y,
                       const FArrayBox&    a_f)
{
#ifndef NDEBUG
    {
        // Check centerings.
        CH_assert(a_fInterp.box().type() == a_interpBox.type());
        CH_assert(a_xInterp.box().type() == a_interpBox.type());
        CH_assert(a_yInterp.box().type() == a_interpBox.type());

        // Check FAB regions.
        CH_assert(a_fInterp.box().contains(a_interpBox));
        CH_assert(a_xInterp.box().contains(a_interpBox));
        CH_assert(a_yInterp.box().contains(a_interpBox));
        CH_assert(a_f.box().size(a_xdir) == a_x.size());
        CH_assert(a_f.box().size(a_ydir) == a_y.size());

        // Check number of comps.
        CH_assert(a_xInterp.nComp() == 1);
        CH_assert(a_yInterp.nComp() == 1);
        CH_assert(a_f      .nComp() == a_fInterp.nComp());

        // Check dirs
        CH_assert(0 <= a_xdir);
        CH_assert(a_xdir < SpaceDim);
        CH_assert(0 <= a_ydir);
        CH_assert(a_ydir < SpaceDim);
        CH_assert(a_xdir != a_ydir);
    }
#endif

    // The fortran function can only handle right-handed
    // permutations of the x, y, and z directions.
    if ((a_xdir == 0 && a_ydir == 1) ||
        (a_xdir == 1 && a_ydir == 2) ||
        (a_xdir == 2 && a_ydir == 0)) {

        // The fortran function requires the vector indices to
        // coincide with the constraint indices.
        IntVect shift = IntVect::Zero;
        shift[a_xdir] = a_f.box().smallEnd(a_xdir);
        shift[a_ydir] = a_f.box().smallEnd(a_ydir);

        // Interpolate
        FORT_BILINEARINTERP2DF (
            CHF_FRA(a_fInterp),
            CHF_CONST_FRA1(a_xInterp,0),
            CHF_CONST_FRA1(a_yInterp,0),
            CHF_BOX(a_interpBox),
            CHF_CONST_INT(a_xdir),
            CHF_CONST_INT(a_ydir),
            CHF_CONST_VR(a_x),
            CHF_CONST_VR(a_y),
            CHF_CONST_FRA_SHIFT(a_f,shift));
    } else {

        // The fortran function requires the vector indices to
        // coincide with the constraint indices.
        IntVect shift = IntVect::Zero;
        shift[a_ydir] = a_f.box().smallEnd(a_xdir);
        shift[a_xdir] = a_f.box().smallEnd(a_ydir);

        // Interpolate

        FORT_BILINEARINTERP2DF (
            CHF_FRA(a_fInterp),
            CHF_CONST_FRA1(a_yInterp,0),
            CHF_CONST_FRA1(a_xInterp,0),
            CHF_BOX(a_interpBox),
            CHF_CONST_INT(a_ydir),
            CHF_CONST_INT(a_xdir),
            CHF_CONST_VR(a_y),
            CHF_CONST_VR(a_x),
            CHF_CONST_FRA_SHIFT(a_f,shift));

    }
}
