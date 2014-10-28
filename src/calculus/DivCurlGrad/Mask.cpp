#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "Mask.H"


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Mask::buildMask (BaseFab<int>&        a_mask,
                      const ProblemDomain& a_dProblem,
                      const BoxLayout&     a_grids,
                      const BoxLayout*     a_fineGridsPtr,
                      const IntVect&       a_nRefFine)
{
    // Start with entire box set to "physical BC"
    a_mask.setVal(maskPhysical);

    // Set all of domain interior to "coarse"
    Box domainInterior(a_mask.box());
    domainInterior &= a_dProblem;
    a_mask.setVal(maskCoarse, domainInterior, 0);

    // Loop over this level's boxes and set them to "copy"
    LayoutIterator lit = a_grids.layoutIterator();
    for (lit.reset(); lit.ok(); ++lit) {
        Box intersectBox = a_grids.get(lit());
        intersectBox &= a_mask.box();

        if (!intersectBox.isEmpty()) {
            a_mask.setVal(maskCopy, intersectBox, 0);
        }
    }

    // If finer grids exist, set them to "covered"
    if (a_fineGridsPtr != NULL) {
        D_TERM(CH_assert(a_nRefFine[0] >= 1);,
               CH_assert(a_nRefFine[1] >= 1);,
               CH_assert(a_nRefFine[2] >= 1);)

        LayoutIterator litFine = a_fineGridsPtr->layoutIterator();
        for (litFine.reset(); litFine.ok(); ++litFine) {
            Box coarsenedBox(a_fineGridsPtr->get(litFine()));
            coarsenedBox.coarsen(a_nRefFine);
            coarsenedBox &= a_mask.box();

            if (!coarsenedBox.isEmpty()) {
                a_mask.setVal(maskCovered, coarsenedBox, 0);
            }
        }
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Mask::buildMasks (LevelData<BaseFab<int> >& a_masks,
                       const ProblemDomain&      a_dProblem,
                       const BoxLayout&          a_grids,
                       const BoxLayout*          a_fineGridsPtr,
                       const IntVect&            a_nRefFine)
{
    DataIterator dit = a_masks.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        buildMask(a_masks[dit()],
                  a_dProblem,
                  a_grids,
                  a_fineGridsPtr,
                  a_nRefFine);
    }
}

