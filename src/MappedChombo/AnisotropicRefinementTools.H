/*******************************************************************************
 *    FILE: AnisotropicRefinementTools.H
 *    DESC: Tools to deal with anisotropic refinement and coarsening.
 *    DATE: Mon 03 Jun 2013 12:59:01 AM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __AnisotropicRefinementTools_H__INCLUDED__
#define __AnisotropicRefinementTools_H__INCLUDED__

class BoxLayout;
class DisjointBoxLayout;
class IntVect;
class IntVectSet;
class ProblemDomain;
class Copier;
class CFRegion;
class CFIVS;


// -----------------------------------------------------------------------------
// Box tools...
// -----------------------------------------------------------------------------

bool coarsenable (const Box&     a_box,
                  const IntVect& a_refRatio);

IntVect calculateRefinementRatio (const Box& a_src,
                                  const Box& a_dest);



// -----------------------------------------------------------------------------
// ProblemDomain tools...
// -----------------------------------------------------------------------------

void coarsen (ProblemDomain&       a_crseDomain,
              const ProblemDomain& a_fineDomain,
              const IntVect&       a_refRatio);

void refine (ProblemDomain&       a_fineDomain,
             const ProblemDomain& a_crseDomain,
             const IntVect&       a_refRatio);


// -----------------------------------------------------------------------------
// BoxLayout tools...
// -----------------------------------------------------------------------------

bool coarsenable (const BoxLayout& a_layout,
                  const IntVect&   a_refRatio);

void coarsen (BoxLayout&       a_output,
              const BoxLayout& a_input,
              const IntVect&   a_refRatio);

void refine (BoxLayout&       a_output,
             const BoxLayout& a_input,
             const IntVect&   a_refinement);

void coarsen (DisjointBoxLayout&       a_output,
              const DisjointBoxLayout& a_input,
              const IntVect&           a_refRatio);

void refine (DisjointBoxLayout&       a_output,
             const DisjointBoxLayout& a_input,
             const IntVect&           a_refinement);


// -----------------------------------------------------------------------------
// Copier tools...
// -----------------------------------------------------------------------------
void coarsen (Copier&        a_copier,
              const IntVect& a_refRatio);


// -----------------------------------------------------------------------------
// IntVectSet tools...
// -----------------------------------------------------------------------------

void coarsen (IntVectSet&    a_ivs,
              const IntVect& a_refRatio);

void refine (IntVectSet&    a_ivs,
             const IntVect& a_refRatio);


// -----------------------------------------------------------------------------
// CFRegion tools...
// -----------------------------------------------------------------------------
void coarsen (CFRegion&      a_cfRegion,
              const IntVect& a_refRatio);


// -----------------------------------------------------------------------------
// CFIVS tools...
// -----------------------------------------------------------------------------
void coarsen (CFIVS&         a_cfivs,
              const IntVect& a_refRatio);



#endif //!__AnisotropicRefinementTools_H__INCLUDED__


