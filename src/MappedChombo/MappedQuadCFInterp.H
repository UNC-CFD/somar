/*******************************************************************************
 *    FILE: MappedQuadCFInterp.H
 *    DESC: Fills ghost cells on coarse-fine interface via quadratic interpolation.
 *    NOTE: Much of this code was taken from Chombo's QuadCFInterp.
 *    DATE: Wed 05 Jun 2013 11:05:57 PM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __MappedQuadCFInterp_H__INCLUDED__
#define __MappedQuadCFInterp_H__INCLUDED__

#include "LevelData.H"
#include "RealVect.H"
class CFIVS;
class IntVectSet;
class TensorFineStencilSet;
class MappedQuadCFStencil;


// Quadratic coarse-fine interpolation utility
class MappedQuadCFInterp {
public:
    // Static member
    static bool newCFInterMode;

    // Default constructor leaves QCFI undefined.
    MappedQuadCFInterp ();

    // Full Constructor. Makes all coarse-fine information and sets internal
    // variables calls full define.
    MappedQuadCFInterp (const DisjointBoxLayout& a_fineBoxes,
                        const DisjointBoxLayout* a_coarBoxes,
                        const RealVect&          a_dxFine,
                        const IntVect&           a_refRatio,
                        const int                a_nComp,
                        const Box&               a_domf);

    // Full Constructor. Makes all coarse-fine information and sets internal
    // variables calls full define.
    MappedQuadCFInterp (const DisjointBoxLayout& a_fineBoxes,
                        const DisjointBoxLayout* a_coarBoxes,
                        const RealVect&          a_dxFine,
                        const IntVect&           a_refRatio,
                        const int                a_nComp,
                        const ProblemDomain&     a_domf);

    // Destructor
    virtual ~MappedQuadCFInterp ();

    // Full define function. Makes all coarse-fine information and sets
    // internal variables.
    virtual void define (const DisjointBoxLayout& a_fineBoxes,
                         const DisjointBoxLayout* a_coarBoxes,
                         const RealVect&          a_dxFine,
                         const IntVect&           a_refRatio,
                         const int                a_nComp,
                         const ProblemDomain&     a_domf);


    // Has full define function been called?
    virtual bool isDefined () const;

    // Return QCFI to undefined state
    virtual void clear ();

    // Coarse / Fine interpolation operator.
    virtual void coarseFineInterp (LevelData<FArrayBox>&       a_phif,
                                   const LevelData<FArrayBox>& a_phic) const;

protected:
    // Get extended phi (lives next to interpivs)
    virtual void getPhiStar (BaseFab<Real>& a_phistar,
                             const BaseFab<Real>& a_phic,
                             const MappedQuadCFStencil& a_qcfs,
                             const Side::LoHiSide a_hiorlo,
                             const int a_idir,
                             const Interval& a_variables) const;

    // Interpolate over correct intvectset
    virtual void interpOnIVS (BaseFab<Real>& a_phif,
                              const BaseFab<Real>& a_phiStar,
                              const MappedQuadCFStencil& a_qcfs,
                              const Side::LoHiSide a_hiorlo,
                              const int a_idir,
                              const Interval& a_variables) const;

    // Coarse / Fine interpolation operator. Using one sided differencing
    // only at a particular face.
    virtual void coarseFineInterp (BaseFab<Real>& a_phif,
                                   const BaseFab<Real>& a_phic,
                                   const MappedQuadCFStencil& a_qcfs,
                                   const Side::LoHiSide a_hiorlo,
                                   const int a_idir,
                                   const Interval& a_variables) const;

    // Boxes for coarse buffers
    BoxLayout m_coarBoxes;

    // Fake level identifier.
    int m_level;

    // Number of components
    int m_nComp;

    // Refinement ratio between levels
    IntVect m_refRatio;

    // Mesh spacing at fine level
    RealVect m_dxFine;

    // cfstencils on low-side faces
    // Lives on fine input grid layout.
    LayoutData<MappedQuadCFStencil> m_loQCFS[SpaceDim];

    // cfstencils on high-side faces
    // Lives on fine input grid layout.
    LayoutData<MappedQuadCFStencil> m_hiQCFS[SpaceDim];

    mutable BoxLayoutData<FArrayBox> m_coarBuffer;
    Copier m_copier;

    DisjointBoxLayout m_inputFineLayout;
    DisjointBoxLayout m_inputCoarLayout;

    // Has full define function been called?
    bool m_isDefined;

    bool m_fineCoversCoarse;
    ProblemDomain m_domainFine;

    bool m_isFlat;
};


#endif //!__MappedQuadCFInterp_H__INCLUDED__
