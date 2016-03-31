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
                            const BoxLayout&         a_destGrids,
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
