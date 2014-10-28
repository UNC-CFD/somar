/*******************************************************************************
 *    FILE: MappedLevelFluxRegister.H
 *    DESC: A remake of Chombo's LevelFluxRegister to handle mapped coordinates
 *          and anisotropic refinement.
 *    DATE: Sun 09 Jun 2013 11:20:14 AM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __MappedLevelFluxRegister_H__INCLUDED__
#define __MappedLevelFluxRegister_H__INCLUDED__

#include "LevelData.H"
class LevelGeometry;


///MappedLevelFluxRegister-A class to encapsulate a levels worth of flux registers.
/**
   A MappedLevelFluxRegister  handles all the data choreography
   necessary to create a levels worth of flux registers.
*/
class MappedLevelFluxRegister {
protected:

    /// Define status
    enum FluxRegDefineStatus {
        FluxRegUndefined      =        0,
        FluxRegDefined        = (1 << 0),      // Indicates define() was called
        FluxRegFineDefined    = (1 << 1),      // m_fineFlux was defined
        FluxRegCoarseDefined  = (1 << 2),      // m_coarFlux was defined
        FluxRegAllDefined     = (1 << 3) - 1   // Everything
    };

public:
    static bool s_verbose;
    ///
    /**
       Default constructor.  Creates an uninitialized levelfluxregister.
    */
    MappedLevelFluxRegister();

    ///
    /**
       Full constructor.  Calls the define function which creates
       a levels worth of flux registers.  a_dProblem is fine-level
       domain.

       If a_scaleFineFluxes is true (default then all fine fluxes are
       scaled by 1/a_nRefine^(SpaceDim-1).  Otherwise, no scaling is done.
    */
    MappedLevelFluxRegister(const DisjointBoxLayout& a_dbl,
                            const DisjointBoxLayout& a_dblCoarse,
                            const ProblemDomain&     a_dProblem,
                            const IntVect&           a_nRefine,
                            int                      a_nComp,
                            bool                     a_scaleFineFluxes = true);

    /// This is an old style constructor - do not use if at all possible.
    /**
       Full constructor.  Calls the define function which creates
       a levels worth of flux registers.  a_dProblem is fine-level
       domain.

       If a_scaleFineFluxes is true (default then all fine fluxes are
       scaled by 1/a_nRefine^(SpaceDim-1).  Otherwise, no scaling is done.
    */
    MappedLevelFluxRegister(const DisjointBoxLayout& a_dbl,
                            const DisjointBoxLayout& a_dblCoarse,
                            const Box&               a_dProblem,
                            const IntVect&           a_nRefine,
                            int                      a_nComp,
                            bool                     a_scaleFineFluxes = true);

    ///
    virtual ~MappedLevelFluxRegister();

    ///
    /**
       Full define function.  Creates a levels worth of flux registers.
       The values in the flux registers are still undefined, however.
       To zero the fluxregisters, you must call setToZero().  a_dProblem
       is fine-level domain. This method may be overridden by a subclass.

       If a_scaleFineFluxes is true then all fine fluxes are
       scaled by 1/a_nRefine^(SpaceDim-1).  Otherwise, no scaling is done.
    */
    virtual void define(const DisjointBoxLayout& a_dbl,
                        const DisjointBoxLayout& a_dblCoarse,
                        const ProblemDomain&     a_dProblem,
                        const IntVect&           a_nRefine,
                        int                      a_nComp,
                        bool                     a_scaleFineFluxes);

    ///
    /**
       This define function calls the above define(), setting a_scaleFineFluxes to true.
    */
    void define(const DisjointBoxLayout& a_dbl,
                const DisjointBoxLayout& a_dblCoarse,
                const ProblemDomain&     a_dProblem,
                const IntVect&           a_nRefine,
                int                      a_nComp);

    /// This is an old style define - do not use if at all possible.
    /**
       Full define function.  Creates a levels worth of flux registers.
       The values in the flux registers are still undefined, however.
       To zero the fluxregisters, you must call setToZero().  a_dProblem
       is fine-level domain.

       If a_scaleFineFluxes is true (default then all fine fluxes are
       scaled by 1/a_nRefine^(SpaceDim-1).  Otherwise, no scaling is done.
    */
    void define(const DisjointBoxLayout& a_dbl,
                const DisjointBoxLayout& a_dblCoarse,
                const Box&               a_dProblem,
                const IntVect&           a_nRefine,
                int                      a_nComp,
                bool                     a_scaleFineFluxes = true);

    ///
    /**
       Modifies this MappedLevelFluxRegister so that it is returned to the
       uninitialized state.  User must now call the full define() before
       using it.

    */
    void
    undefine();

    ///
    /**
        Initialize values of registers to zero.
    */
    void setToZero();

    ///
    /**
       increments the register with data from coarseFlux, multiplied by scale.
       coarseFlux must contain the coarse fluxes in the dir direction
       for the grid   m_coarseLayout[coarseDataIndex].  This operation is local.
       Increment coarse does register += -sign(side)*incrScale*coarseFlux
       for both Lo and Hi "side".
    */
    void incrementCoarse(
        const FArrayBox& a_coarseFlux,
        Real a_scale,
        const DataIndex& a_coarseDataIndex,
        const Interval& a_srcInterval,
        const Interval& a_dstInterval,
        int a_dir);

    ///
    /**
       This method implements the register's coarse increment on a single side.
       It may be overridden by a subclass.
    */
    virtual void incrementCoarse(
        const FArrayBox& a_coarseFlux,
        Real a_scale,
        const DataIndex& a_coarseDataIndex,
        const Interval& a_srcInterval,
        const Interval& a_dstInterval,
        int a_dir,
        Side::LoHiSide a_sd);

    ///
    /**
       increments the register with data from fineFlux, multiplied by scale.
       fineFlux must contain the coarse fluxes in the dir direction for the
       grid m_fineLayout[fineDataIndex].  This operation is local.
       Increment fine does register +=  sign(side)*incrScale*Ave(fineFlux)
       for both Lo and Hi "side".
    */
    void incrementFine(
        const FArrayBox& a_fineFlux,
        Real a_scale,
        const DataIndex& a_fineDataIndex,
        const Interval& a_srcInterval,
        const Interval& a_dstInterval,
        int a_dir);

    ///
    /**
       This method implements the register's fine increment on one side.
       It may be overridden by a subclass.
    */
    virtual void incrementFine(
        const FArrayBox& a_fineFlux,
        Real a_scale,
        const DataIndex& a_fineDataIndex,
        const Interval& a_srcInterval,
        const Interval& a_dstInterval,
        int a_dir,
        Side::LoHiSide a_sd);

    // Reflux does cellValue -= refluxScale * registerContents / J
    virtual void reflux (LevelData<FArrayBox>& a_levelDiv,
                         const LevelGeometry&  a_crseLevGeo,
                         const Real            a_scale = 1.0);

    // Reflux does cellValue -= refluxScale * registerContents / J
    virtual void reflux (LevelData<FArrayBox>&       a_levelDiv,
                         const LevelData<FArrayBox>& a_crseCCJinv,
                         const Real                  a_scale = 1.0);

    // // Increments uCoarse with the reflux divergence of the
    // // contents of the flux register.
    // // Note that there is no srccomp etc crap here.
    // // this is done for all components so uCoarse has
    // // to have the same number of components as input nComp.
    // // This operation is global and blocking.
    // //
    // // Reflux does cellValue -= refluxScale*registerContents
    // virtual void reflux (LevelData<FArrayBox>& a_uCoarse,
    //                      Real                  a_scale);

    // // Interface for multiblock algorithms where, if an interval
    // // of the components represents a vector, it must be specially
    // // handled.  If called from a singleblock code, the vector
    // // interval is just ignored.
    // virtual void reflux (LevelData<FArrayBox>& a_uCoarse,
    //                      const Interval&       a_coarseVectorIntv,
    //                      Real                  a_scale);


    bool hasCF(const DataIndex& a_fineDataIndex, Side::LoHiSide) const;
    bool hasCF(const DataIndex& a_coarseIndex) const;

    ///
    /**
       has full define function been called?  return true if so.
    */
    bool isDefined() const;

    bool isAllDefined() const;

    ///
    void poutCoarseRegisters() const;

    ///
    void poutFineRegisters() const;

    LayoutData< Vector<Box> >& getCoarseLocations(int a_idir, Side::LoHiSide a_sd) {
        CH_assert(isDefined());
        int index = a_idir + a_sd * CH_SPACEDIM;
        return m_coarseLocations[index];
    }

    Copier& getReverseCopier() {
        return m_reverseCopier;
    }

protected:
    // Increments uCoarse with the reflux divergence of the
    // contents of the flux register. This can be overridden by
    // a subclass. This operation is global and blocking.
    // Reflux does cellValue -= refluxScale*registerContents
    virtual void reflux (LevelData<FArrayBox>& a_uCoarse,
                         const Interval&       a_coarse_interval,
                         const Interval&       a_flux_interval,
                         Real                  a_scale);

    // Same as above with a variable scale multiplied in
    void reflux (LevelData<FArrayBox>&       a_uCoarse,
                 Real                        a_scale,
                 const Interval&             a_coarseInterv,
                 const Interval&             a_fluxInterval,
                 const LevelData<FArrayBox>& a_beta);

    // both of these are defined in the coarse index space,
    // but, m_coarFlux uses the coarse grid layout, and
    // m_fineFlux uses the fine grid layout.
    LevelData<FArrayBox> m_coarFlux;
    LevelData<FArrayBox> m_fineFlux;

    //LayoutData<IntVectSet> m_coarseLocations[CH_SPACEDIM*2];
    LayoutData< Vector<Box> > m_coarseLocations[CH_SPACEDIM * 2];

    ProblemDomain m_domain;

    int m_isDefined;

    IntVect m_nRefine;

    bool m_scaleFineFluxes;

    Copier m_reverseCopier;

private:

    ///there is no operator= for this class
    void operator= (const MappedLevelFluxRegister&);
    ///there is no copy constructor for this class
    MappedLevelFluxRegister(const MappedLevelFluxRegister&);
};



#endif //!__MappedLevelFluxRegister_H__INCLUDED__
