/*******************************************************************************
 *    FILE: MappedCFStencil.H
 *    DESC: C/F derivative utilities for anisotropic refinement.
 *    NOTE: Much of this code was taken from Chombo's CFStencil/QuadCFStencil.
 *    DATE: Thu 06 Jun 2013 11:02:44 AM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __MappedCFStencil_H__INCLUDED__
#define __MappedCFStencil_H__INCLUDED__

#include "BaseFab.H"
#include "DerivStencil.H"
#include "DisjointBoxLayout.H"
#include "LoHiSide.H"
#include "IntVectSet.H"
class RealVect;


/// Class to encapsulate coarse-fine information for AMR operations.
/**
    MappedCFStencil is a wrapper around a bunch of things
    that allow one to do coarse-fine interpolation.
    This class should be considered internal to QuadCFInterp
    and should not be considered part of the Chombo API.
  */
class MappedCFStencil {

private:

protected:
    //for internal use
    void  setDefaultValues();

    //the following variables are mainly for self-identification

    // direction normal to interface
    int m_direction;

    // interface is on high or low side of grid
    Side::LoHiSide m_hiorlo;

    //
    DataIndex m_dataIndex;

    //the following variables are to be used by others

    /* fine intvects which need to be interpolated */
    IntVectSet m_fineIVS;

    /* coar intvects where slopes can be computed
        == coarsened version of fiinterp ivs*/
    IntVectSet m_coarIVS;

    bool m_isDefined;

    bool m_isPacked;

    Box  m_packedBox;

public:
    /// {\bf Access functions}

    ///
    /** returns true if MappedCFStencil is defined */
    bool isDefined() const;

    ///
    /**
       are there any interpolation points?  \\
       returns false if so.*/
    bool isEmpty() const;

    ///
    bool isPacked() const {
        return  m_isPacked;
    }
    ///
    const Box& packedBox() const

    {
        return m_packedBox;
    }

    ///
    /** get fine intvects which need to be interpolated  \\
        This will be empty if isEmpty() returns true*/
    const IntVectSet& getFineIVS() const;

    ///
    /** get coarse intvects that underly fiinterpivs.
        This will be empty if isEmpty() returns true
        this is the same as getCoarInterpIVS
    */
    const IntVectSet& getCoarIVS() const ;

    ///
    MappedCFStencil& operator= (const MappedCFStencil& cfs_in);

    ///
    MappedCFStencil(const MappedCFStencil& cfs_in);

    /*dummy constructor so i can use this with
       LocalBoxLayoutData
     */
    MappedCFStencil(const Box& a_b, int a_nComp);

    void define(const Box& a_b, int a_nComp);

    ///
    MappedCFStencil();

    ///
    ~MappedCFStencil();

    ///
    MappedCFStencil(
        const Box& FineDomain,
        const Box& grid,
        const DisjointBoxLayout& Levboxes,
        const DisjointBoxLayout& LevCBoxes,
        const IntVect& RefRatio,
        int Direction,
        Side::LoHiSide hiorlo);

    ///
    MappedCFStencil(
        const ProblemDomain& FineDomain,
        const Box& grid,
        const DisjointBoxLayout& Levboxes,
        const DisjointBoxLayout& LevCBoxes,
        const IntVect& RefRatio,
        int Direction,
        Side::LoHiSide hiorlo);

    ///
    void   define(
        const Box& FineDomain,
        const Box& grid,
        const DisjointBoxLayout& fine_boxes,
        const DisjointBoxLayout& coar_boxes,
        const IntVect& Refratio,
        int Direction,
        Side::LoHiSide hiorlo);

    ///
    void   define(
        const ProblemDomain& FineDomain,
        const Box& grid,
        const DisjointBoxLayout& fine_boxes,
        const DisjointBoxLayout& coar_boxes,
        const IntVect& Refratio,
        int Direction,
        Side::LoHiSide hiorlo);

    ///  faster define function that uses a sorted periodic vector
    /** user calls the helper function "buildPeriodicVector", then
        uses that returned Vector<Box> in the subsequent calls to this define.
    */
    void   define(
        const ProblemDomain& a_fineDomain,
        const Box& a_grid,
        const Vector<Box>& a_periodicVector,
        const IntVect& a_refRatio,
        int a_direction,
        Side::LoHiSide a_hiorlo);
    ///
    void clear();


    static void buildPeriodicVector(Vector<Box>& a_periodicVector,
                                    const ProblemDomain& a_fineDomain,
                                    const DisjointBoxLayout& a_fineBoxes);

};


/// class to encapsulate CF info for quadratic interpolation
/**
    MappedQuadCFStencil is a wrapper around a bunch of things
    that allow LevelOp to do coarse-fine interpolation
    on a particular face of a grid.
  */
class MappedQuadCFStencil {

public:
    /// {\bf Derivative functions}

    ///
    /**
       compute second derivative in devdir_a direction
       at coarse point a_ivin
       Uses centered finite diff approximation if
       ivs_standard_m.contains(a_ivin ).
       Otherwise, it uses the stencil from second_m
       Asserts that a_derivdir != direction of face
    */
    Real computeSecondDerivative(
        const BaseFab<Real>& a_phic,
        int a_derivdir,
        int a_ivar,
        const IntVect& a_ivin,
        Real a_dxDir) const;

    ///
    /**
       compute first derivative in devdir_a direction
       at coarse point a_ivin
       Uses centered finite diff approximation if
       ivs_standard_m.contains(a_ivin ).
       Otherwise, it uses the stencil from firstd_m
       Asserts that a_derivdir != direction of face
    */
    Real computeFirstDerivative(
        const BaseFab<Real>& a_phic,
        int a_derivdir,
        int a_ivar,
        const IntVect& a_ivin,
        Real a_dxDir) const;

    ///
    /**
       compute mixed derivative (direction is unambiguous.
       x and y are the two directions tangential to face)
       at coarse point a_ivin. \\
       In two dimensions, always returns 0. \\
       Uses centered finite diff approximation if
       ivs_standard_m.contains(a_ivin ).
       It uses the stencil from mixed_sten_m.\\
       Returns 0 if SpaceDim != 3
    */
    Real computeMixedDerivative(
        const BaseFab<Real>& a_phic,
        int a_ivar,
        const IntVect& a_ivin,
        const RealVect& a_dx) const;

    ///
    MappedQuadCFStencil();

    ///
    ~MappedQuadCFStencil();

    ///
    MappedQuadCFStencil(
        const Box& a_fine_domain,
        const Box& a_grid,
        const DisjointBoxLayout& a_fineBoxes,
        const DisjointBoxLayout& a_coarBoxes,
        const IntVect& a_refRatio,
        int a_direction,
        Side::LoHiSide a_hiorlo,
        bool a_isFlat = false);

    ///
    MappedQuadCFStencil(
        const ProblemDomain& a_fine_domain,
        const Box& a_grid,
        const DisjointBoxLayout& a_fineBoxes,
        const DisjointBoxLayout& a_coarBoxes,
        const IntVect& a_refRatio,
        int a_direction,
        Side::LoHiSide a_hiorlo,
        bool a_isFlat = false);

    ///
    void   define(
        const Box& a_fine_domain,
        const Box& a_grid,
        const DisjointBoxLayout& a_fineBoxes,
        const DisjointBoxLayout& a_coarBoxes,
        const IntVect& a_refRatio,
        int a_direction,
        Side::LoHiSide a_hiorlo,
        bool a_isFlat = false);

    ///
    void   define(
        const ProblemDomain& a_fine_domain,
        const Box& a_grid,
        const DisjointBoxLayout& a_fineBoxes,
        const DisjointBoxLayout& a_coarBoxes,
        const IntVect& a_refRatio,
        int a_direction,
        Side::LoHiSide a_hiorlo,
        bool a_isFlat = false);


    ///  faster constructor that requires the use of CFStencil::buildPeriodicVector
    void   define(
        const ProblemDomain& a_fine_domain,
        const Box& a_grid,
        const Vector<Box>& a_periodicBoxes,
        const Vector<Box>& a_coarsenedPeriodicBoxes,
        const DisjointBoxLayout& a_coarBoxes,
        const IntVect& a_refRatio,
        int a_direction,
        Side::LoHiSide a_hiorlo,
        bool a_isFlat);

    /// {\bf Access functions}

    ///
    /**
       has full define function been called?  return true if so
    */
    bool isDefined() const;

    ///
    /**
       are there any interpolation points?  \\
       returns false if so.*/
    bool isEmpty() const {
        CH_assert(m_baseCFS.isDefined());
        return m_baseCFS.isEmpty();
    }

    ///
    /** get fine intvects which need to be interpolated  \\
        This will be empty if isEmpty() returns true*/
    const IntVectSet& getFineIVS() const {
        CH_assert(m_baseCFS.isDefined());
        return m_baseCFS.getFineIVS();
    }

    ///
    /** get coarse intvects that underly fiinterpivs.
        This will be empty if isEmpty() returns true*/
    const IntVectSet& getCoarIVS() const {
        CH_assert(m_baseCFS.isDefined());
        return m_baseCFS.getCoarIVS();
    }

    bool finePacked() const {
        return m_baseCFS.isPacked();
    }

    const Box& packedBox() const {
        return m_baseCFS.packedBox();
    }


protected:

    //has define fcn been called?
    bool m_isDefined;

    //direction of face
    int m_direction;

    ///
    /**
       set of all points which have ALL their
       derivatives taken normally.  This means
       that no stencil either leaves the domain or
       intersects the next finer level.
    */
    IntVectSet m_ivsStandard;

    ///
    /**
       set of all points which have at least
       one derivative shifted
    */
    IntVectSet m_ivsQuadd;

    /**
       Derivative stencils.  Only defined at intvects where there
       is at least one non-standard derivative to take
       (ie. at all points of ivs_Quadd)
    */

    ///
    /**
       stencil information for mixed derivatives
       undefined if SpaceDim != 3
    */
    BaseFab<DerivStencil> m_mixedSten;

    ///
    /**
       stencil information for unmixed derivatives
       (second derivatives)
    */
    BaseFab<DerivStencil> m_second[SpaceDim];

    ///
    /**
       stencil information for unmixed derivatives
       (first derivatives)
    */
    BaseFab<DerivStencil> m_firstD[SpaceDim];

    ///
    /**
       whether to drop order at any
       one-sided point
    */
    BaseFab<bool> m_dropOrd;

    ///
    MappedCFStencil m_baseCFS;

    bool m_isFlat;


private:
    // internally useful but not for public consumption

    //
    void  setDefaultValues();

    /*
      add fab to the stencil.  For internal use only.
    */
    void addFabToSten(const BaseFab<Real>& fabin_a,
                      DerivStencil&  sten_a);

    // there is no assignment operator for this class
    void operator= (const MappedQuadCFStencil& cfs_in) {
    }

    // there is no copy constructor for this class
    MappedQuadCFStencil(const MappedQuadCFStencil& cfs_in) {
    }

    void buildStencils(const IntVectSet& ivsAllGood);

};




#endif //!__MappedCFStencil_H__INCLUDED__
