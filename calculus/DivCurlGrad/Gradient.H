#ifndef __GRADIENT_H__INCLUDED__
#define __GRADIENT_H__INCLUDED__

#include "LevelData.H"
class BCMethodHolder;
class IntVectSet;
class LevelGeometry;
class MappedQuadCFInterp;
class RealVect;


// Class to encapsulate Gradient functions (both CC and face-centered)
class Gradient
{
public:
    // The gradient functions will use setFluxes() to override the calculated
    // values at the physical boundaries. This is especially useful in mapped
    // coordinates when homogeneous Neumann BCs must be enforced.
    // Also, if a setGhosts() method is attached to the holder, this will be
    // called to fill ghosts on a_phi before calculating the gradient. This
    // saves a call to ExtrapolateFaceAndCopy().
    //
    // NOTE: The inhomogeneous version of all BC methods will be called.
    // NOTE: No exchanges are made in the gradient functions!
    typedef BCMethodHolder BC_type;


    // Single-level Gradient functions -----------------------------------------

    // Computes FC gradient of CC phi.
    // If phiCrsePtr != NULL, does quadratic interpolation to compute the
    // coarse-fine BCs for phi.
    static void levelGradientMAC (LevelData<FluxBox>&           a_edgeGrad,
                                  LevelData<FArrayBox>&         a_phi,
                                  const LevelData<FArrayBox>*   a_phiCrsePtr,
                                  const LevelGeometry&          a_levGeo,
                                  const Real                    a_time,
                                  const BC_type*                a_fluxBC);

    /// computes edge-centered level-operator gradient of cell-centered phi
    /**
       if phiCrsePtr != NULL, does quadratic interpolation
       to compute coarse-fine boundary conditions for phi.  A predefined
       MappedQuadCFInterp object is passed in for efficiency */
    static void levelGradientMAC (LevelData<FluxBox>&           a_edgeGrad,
                                  LevelData<FArrayBox>&         a_phi,
                                  const LevelData<FArrayBox>*   a_phiCrsePtr,
                                  const MappedQuadCFInterp&     a_cfInterpCrse,
                                  const LevelGeometry&          a_levGeo,
                                  const Real                    a_time,
                                  const BC_type*                a_fluxBC);


    /// computes edge-centered level-operator gradient of cell-centered phi
    /**
        assumes _ALL_ ghost cell values have been preset (so phi can be const)
    */
    static void levelGradientMAC (LevelData<FluxBox>&           a_edgeGrad,
                                  const LevelData<FArrayBox>&   a_phi,
                                  const LevelGeometry&          a_levGeo,
                                  const Real                    a_time,
                                  const BC_type*                a_fluxBC);


    /// computes cell-centered level-operator gradient of cell-centered phi
    /** if phiCrse != NULL, does coarse-fine boundary
        conditions for phi (quadratic interpolation) */
    static void levelGradientCC (LevelData<FArrayBox>&          a_grad,
                                 LevelData<FArrayBox>&          a_phi,
                                 const LevelData<FArrayBox>*    a_phiCrsePtr,
                                 const LevelGeometry&           a_levGeo,
                                 const Real                     a_time,
                                 const BC_type*                 a_fluxBC);


    /// computes cell-centered level-operator gradient of cell-centered phi
    /** if phiCrse != NULL, does coarse-fine boundary
        conditions for phi (quadratic C/F interpolation); predefined
        MappedQuadCFInterp object is passed in for efficiency */
    static void levelGradientCC (LevelData<FArrayBox>&          a_grad,
                                 LevelData<FArrayBox>&          a_phi,
                                 const LevelData<FArrayBox>*    a_phiCrsePtr,
                                 MappedQuadCFInterp&            a_cfInterp,
                                 const LevelGeometry&           a_levGeo,
                                 const Real                     a_time,
                                 const BC_type*                 a_fluxBC);

    /// computes cell-centered, level-operator gradient of cell-centered phi
    /** in this case, assume that all relevant
        BC's (coarse-fine and physical) have already been set, so
        phi can be a const variable */
    static void levelGradientCC (LevelData<FArrayBox>&          a_grad,
                                 const LevelData<FArrayBox>&    a_phi,
                                 const LevelGeometry&           a_levGeo,
                                 const Real                     a_time,
                                 const BC_type*                 a_fluxBC);


    // Composite Gradient functions --------------------------------------------

    /// computes edge-centered composite gradient of cell-centered phi
    /**  if phiCrsePtr != NULL, does quadratic interpolation to compute
         coarse-fine boundary conditions for phi; since the edge between
         this level and finer levels is not considered to be a part of
         this level, there is no fine-level coarse-fine BC.  because of
         this, this function produces exactly the same results as
         LevelGradientMAC -- it's just included for completeness... */
    static void compGradientMAC (LevelData<FluxBox>&            a_edgeGrad,
                                 LevelData<FArrayBox>&          a_phi,
                                 const LevelData<FArrayBox>*    a_phiCrse,
                                 const LevelData<FArrayBox>*    a_phiFine,
                                 const LevelGeometry&           a_levGeo,
                                 const Real                     a_time,
                                 const BC_type*                 a_fluxBC);

    /// computes edge-centered composite gradient of cell-centered phi
    /** if phiCrsePtr != NULL, does quadratic interpolation to compute
        coarse-fine boundary conditions for phi; since the edge between
        this level and finer levels is not considered to be a part of
        this level, there is no fine-level coarse-fine BC.  because of
        this, this function produces exactly the same results as
        LevelGradientMAC -- it's just included for completeness... In
        this one, a predefined MappedQuadCFInterp object is passed in.  This
        (deprecated) interface uses a Box instead of a ProblemDomain object */
    static void compGradientMAC (LevelData<FluxBox>&            a_edgeGrad,
                                 LevelData<FArrayBox>&          a_phi,
                                 const LevelData<FArrayBox>*    a_phiCrse,
                                 const LevelData<FArrayBox>*    a_phiFine,
                                 MappedQuadCFInterp&            a_cfInterpCrse,
                                 const LevelGeometry&           a_levGeo,
                                 const Real                     a_time,
                                 const BC_type*                 a_fluxBC);


    /// computes cell-centered composite gradient of cell centered phi
    /*
        Uses same coarse-level C/F BC's as LevelGradientCC, if phiFinePtr !=
        NULL, then also uses one-sided differencing to compute gradient
        on coarse side of coarse-fine interface. */
    static void compGradientCC (LevelData<FArrayBox>&       a_grad,
                                LevelData<FArrayBox>&       a_phi,
                                const LevelData<FArrayBox>* a_phiCrsePtr,
                                const LevelData<FArrayBox>* a_phiFinePtr,
                                const LevelGeometry&        a_levGeo,
                                const Real                  a_time,
                                const BC_type*              a_fluxBC);


    /// computes cell-centered composite gradient of cell-centered phi
    /**
        uses same coarse-level C/F BC's as LevelGradientCC, if phiFinePtr !=
        NULL, then also uses one-sided differencing to compute gradient
        on coarse side of coarse-fine interface. A predefined MappedQuadCFInterp
        is also passed in for efficiency.  Note that no fine MappedQuadCFInterp
        is necessary because we use one-sided differencing for coarse
        cells adjacent to finer-level coarse-fine interfaces. Note also
        that gradient is only really defined in valid regions of grids. */
    static void compGradientCC (LevelData<FArrayBox>&           a_grad,
                                LevelData<FArrayBox>&           a_phi,
                                const LevelData<FArrayBox>*     a_phiCrsePtr,
                                const LevelData<FArrayBox>*     a_phiFinePtr,
                                MappedQuadCFInterp&             a_cfInterpCrse,
                                const LevelGeometry&            a_levGeo,
                                const Real                      a_time,
                                const BC_type*                  a_fluxBC);


    /// computes cell-centered composite gradient of cell-centered phi
    /**
        this one assumes that all ghost-cell values have already been set;
        if phiFinePtr != NULL, then also uses one-sided differencing to compute
        gradient on coarse side of corarse-fine interface.  note that gradient
        is only really defined in valid regions of grids. */
    static void compGradientCC (LevelData<FArrayBox>&       a_grad,
                                const LevelData<FArrayBox>& a_phi,
                                const LevelData<FArrayBox>* a_phiFinePtr,
                                const LevelGeometry&        a_levGeo,
                                const Real                  a_time,
                                const BC_type*              a_fluxBC);


    // Additional utilities ----------------------------------------------------

    /// Use this to compute (u.Del)U...that is, CC.grad(FC)
    static void levelCCDotGradFC (LevelData<FArrayBox>&       a_CCDotGradFC,
                                  const LevelData<FArrayBox>& a_CC,
                                  const LevelData<FluxBox>&   a_FC,
                                  const LevelGeometry&        a_levGeo);

    /// utility function for internal use
    static void singleBoxMacGrad (FArrayBox&             a_gradFab,
                                  const FArrayBox&       a_phiFab,
                                  int                    a_gradComp,
                                  int                    a_phiComp,
                                  int                    a_numComp,
                                  const Box&             a_edgeBox,
                                  const Box&             a_validPhi,
                                  int                    a_dir,
                                  int                    a_edgeDir,
                                  const DataIndex&       a_di,
                                  const LevelGeometry&   a_levGeo,
                                  const Real             a_time,
                                  const BC_type*         a_fluxBC);
};


#endif //!__GRADIENT_H__INCLUDED__
