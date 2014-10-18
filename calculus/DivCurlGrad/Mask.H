#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#ifndef __MASK_H__INCLUDED__
#define __MASK_H__INCLUDED__

#include "LevelData.H"
#include "MASKVAL.H"


// Class to determine coarse-fine validity info
class Mask
{
public:
    enum cellType
    {
        maskCovered  = MASKCOVERED, // Invalid cells. These are covered by valid finer level data.
        maskPhysical = MASKPHYS,    // Coarse-level invalid ghosts outside the domain.
        maskCopy     = MASKCOPY,    // Valid cells and exchange ghosts NOT covered by valid finer level data.
        maskCoarse   = MASKCOARSE   // Coarse-level invalid ghosts inside the domain.
    };

    static void buildMask (BaseFab<int>&        a_mask,
                           const ProblemDomain& a_dProblem,
                           const BoxLayout&     a_grids,
                           const BoxLayout*     a_fineGridsPtr = NULL,
                           const IntVect&       a_nRefFine = IntVect(D_DECL(-1,-1,-1)));

    static void buildMasks (LevelData<BaseFab<int> >& a_masks,
                            const ProblemDomain&      a_dProblem,
                            const BoxLayout&          a_grids,
                            const BoxLayout*          a_fineGridsPtr = NULL,
                            const IntVect&            a_nRefFine = IntVect(D_DECL(-1,-1,-1)));
};


#endif //!__MASK_H__INCLUDED__
