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
#ifndef __PeriodicLoHiCenter_H__INCLUDED__
#define __PeriodicLoHiCenter_H__INCLUDED__

#include "Box.H"
#include "ProblemDomain.H"


// This is the same as Chombo's LoHiCenter, but considers ghosts in periodic
// directions as valid.
//
// NOTE: These functions assume _FOUR_ ghost cells! It would be better to
// create a more general version of these functions that assumes nothing.


///
/**
   This function is used when in direction a_dir a 3-point stencil of
   cell-centered data is being used to compute something at the
   center of the central cell of the stencil.

   Inputs:
   * a_dir is the direction of the 3-point stencil;
   * a_inBox is the cell-centered box on which we are to compute a
     stencil (2-point or 3-point), expanded by 1 in a_dir direction,
     and not intersected with the domain;
   * a_domain is the problem domain.

   Output boxes are all cell-centered subboxes of a_inBox and
   are also contained in a_domain:
   * a_loBox is where a 2-point stencil must be used on the low side;
   * a_centerBox is where the full 3-point stencil can be used;
   * a_hiBox is where a 2-point stencil must be used on the high side;
   * a_entireBox is the union of a_loBox and a_centerBox and a_hiBox.

   The boxes a_loBox and a_hiBox will be at most 1 cell wide;
   here the 2-point stencil consists of a cell from this box and
   the neighboring cell from a_centerBox in direction a_dir.

   a_entireBox will be a_inBox with one layer of cells removed
   from both sides in dimension a_dir.

   Output flags:
   * a_hasLo:  1 or 0, according to whether a_loBox is defined or not;
   * a_hasHi:  1 or 0, according to whether a_hiBox is defined or not.
 */
void PeriodicLoHiCenter (Box&                 a_loBox,
                         int&                 a_hasLo,
                         Box&                 a_hiBox,
                         int&                 a_hasHi,
                         Box&                 a_centerBox,
                         Box&                 a_entireBox,
                         const Box&           a_inBox,
                         const ProblemDomain& a_domain,
                         const int&           a_dir);


///
/**
   This function is used when in direction a_dir a 2-point stencil of
   cell-centered data is being used to compute something at the
   face between cells of the stencil.

   Inputs:
   * a_dir is the direction of the 2-point stencil;
   * a_inBox is the cell-centered box (on whose faces we are to
     compute a stencil), expanded by 1 in a_dir direction, and not
     intersected with the domain.
   * a_domain is the problem domain.

   Output boxes are all face-centered subboxes of the faces of a_inBox
   in direction a_dir, and are also contained in the faces of a_domain:
   * a_loBox is where a 1-point stencil must be used on the low side;
   * a_centerBox is where the full 2-point stencil can be used;
   * a_hiBox is where a 1-point stencil must be used on the high side;
   * a_entireBox is the union of a_loBox and a_centerBox and a_hiBox.

   The boxes a_loBox and a_hiBox will be at most 1 face wide;
   here the 1-point stencil consists of only 1 of the adjacent cells
   in direction a_dir.

   Output flags:
   * a_hasLo:  1 or 0, according to whether a_loBox is defined or not;
   * a_hasHi:  1 or 0, according to whether a_hiBox is defined or not.
 */
void PeriodicLoHiCenterFace (Box&                 a_loBox,
                             int&                 a_hasLo,
                             Box&                 a_hiBox,
                             int&                 a_hasHi,
                             Box&                 a_centerBox,
                             Box&                 a_entireBox,
                             const Box&           a_inBox,
                             const ProblemDomain& a_domain,
                             const int&           a_dir);


///
/**
   This function is used when in direction a_dir a 4-point stencil of
   cell-centered data is being used to compute something at the
   face between cells of the stencil.

   Inputs:
   * a_dir is the direction of the 4-point stencil;
   * a_inBox is the cell-centered box on which we have data;
   * a_domain is the problem domain.

   Output boxes are all face-centered subboxes of the faces of a_inBox
   in direction a_dir, and are also contained in the faces of a_domain:
   * a_loBox is where a stencil must be used on the low side
     with all four points higher;
   * a_nextLoBox is where a stencil must be used on the low side
     with one point lower and three points higher;
   * a_hiBox is where a stencil must be used on the high side
     with all four points lower;
   * a_nextHiBox is where a stencil must be used on the high side
     with one point higher and three points lower;
   * a_centerBox is union of a_innerCenterBox, a_nextLoBox, and a_nextHiBox,
     where the regular 2-point stencil can be used;
   * a_innerCenterBox is where the regular 4-point stencil can be used;
   * a_entireBox is the union of a_(lo|nextLo|center|hi|nextHi)Box.

   Each of the boxes a_loBox, a_nextLoBox, a_hiBox, a_nextHiBox
   will be at most 1 face wide.

   a_loBox and a_nextLoBox will both be defined or both be undefined.
   a_hiBox and a_nextHiBox will both be defined or both be undefined.

   Output flags:
   * a_hasLo:  1 or 0, according to whether a_loBox and a_nextLoBox
     are defined or not;
   * a_hasHi:  1 or 0, according to whether a_hiBox and a_nextHiBox
     are defined or not.
**/
void PeriodicLoHiCenterFace4 (Box&                 a_loBox,
                              Box&                 a_nextLoBox,
                              int&                 a_hasLo,
                              Box&                 a_hiBox,
                              Box&                 a_nextHiBox,
                              int&                 a_hasHi,
                              Box&                 a_centerBox,
                              Box&                 a_innerCenterBox,
                              Box&                 a_entireBox,
                              const Box&           a_inBox,
                              const ProblemDomain& a_domain,
                              const int&           a_dir);


///
/**
   This function is used when in direction a_dir a 5-point stencil of
   cell-centered data is being used to compute something at the
   center of the central cell of the stencil.

   Inputs:
   * a_dir is the direction of the 5-point stencil;
   * a_inBox is the cell-centered box on which we are to compute a
     stencil (2-point or 3-point or 5-point), expanded by 1 in a_dir direction,
     and not intersected with the domain;
   * a_domain is the problem domain.

   Output boxes are all cell-centered subboxes of a_inBox and
   are also contained in a_domain:
   * a_loBox is where a 2-point stencil must be used on the low side;
   * a_nextLoBox is where a 3-point stencil must be used on the low side
     with one point lower and three points higher;
   * a_hiBox is where a 2-point stencil must be used on the high side;
   * a_nextHiBox is where a 3-point stencil must be used on the high side
     with one point higher and three points lower;
   * a_centerBox is union of a_innerCenterBox, a_nextLoBox, and a_nextHiBox,
     where a 3-point stencil may be used;
   * a_innerCenterBox is where the regular 5-point stencil can be used;
     where a 5-point stencil may be used;
   * a_entireBox is the union of a_(lo|nextLo|center|hi|nextHi)Box.

   Each of the boxes a_loBox, a_nextLoBox, a_hiBox, a_nextHiBox
   will be at most 1 cell wide.

   a_loBox and a_nextLoBox will both be defined or both be undefined.
   a_hiBox and a_nextHiBox will both be defined or both be undefined.

   Output flags:
   * a_hasLo:  1 or 0, according to whether a_loBox and a_nextLoBox
     are defined or not;
   * a_hasHi:  1 or 0, according to whether a_hiBox and a_nextHiBox
     are defined or not.
**/
void PeriodicLoHiCenter5 (Box&                 a_loBox,
                          Box&                 a_nextLoBox,
                          int&                 a_hasLo,
                          Box&                 a_hiBox,
                          Box&                 a_nextHiBox,
                          int&                 a_hasHi,
                          Box&                 a_centerBox,
                          Box&                 a_innerCenterBox,
                          Box&                 a_entireBox,
                          const Box&           a_inBox,
                          const ProblemDomain& a_domain,
                          const int&           a_dir);


#endif //!__PeriodicLoHiCenter_H__INCLUDED__
