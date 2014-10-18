/*******************************************************************************
 *  FILE: VelocityAMRPoissonOp.H
 *  DESC: Operator for solving (alpha I + beta*Laplacian)(phi) = rho over an
 *        AMR hierarchy in mapped coordiantes.
 *  DATE: Some time around Thu 24 May 2012, I think.
 *  MAIL: santilli@physics.unc.edu
 *  NOTE: This is a rework of VelocityAMRPoissonOp in Chombo's AMRINS example.
 ******************************************************************************/
#ifndef __VELOCITYAMRPOISSONOP_H__INCLUDED__
#define __VELOCITYAMRPOISSONOP_H__INCLUDED__

#include "VelBCHolder.H"
#include "MappedAMRPoissonOp.H"


// The reason for this class -- want to be able to apply operator
// with different BCs for each component in a_phi.
class VelocityAMRPoissonOp: public MappedAMRPoissonOp
{
public:
    // Default constructor
    VelocityAMRPoissonOp()
    {;}

    // Destructor
    virtual ~VelocityAMRPoissonOp()
    {;}

    // Defines this op to simply calculate Laplacians and residuals.
    void define(const LevelGeometry&                        a_levGeo,
                const Real                                  a_alpha,
                const Real                                  a_beta,
                const RefCountedPtr<LevelData<FArrayBox> >& a_lapDiagsPtr = RefCountedPtr<LevelData<FArrayBox> >(NULL));

    //
    virtual void applyOpI(LevelData<FArrayBox>&        a_lhs,
                          const LevelData<FArrayBox>&  a_phi,
                          bool                         a_homogeneous,
                          VelBCHolder*                 a_velBCPtr);

    ///
    virtual void AMROperatorNF(LevelData<FArrayBox>&        a_LofPhi,
                               const LevelData<FArrayBox>&  a_phi,
                               const LevelData<FArrayBox>&  a_phiCoarse,
                               bool                         a_homogeneousPhysBC,
                               VelBCHolder*                 a_velBCPtr);

    ///
    virtual void applyOp(LevelData<FArrayBox>&        a_lhs,
                         const LevelData<FArrayBox>&  a_phi,
                         const LevelData<FArrayBox>*  a_phiCoarsePtr,
                         bool                         a_homogeneous,
                         VelBCHolder*                 a_velBCPtr);
};

#endif //!__VELOCITYAMRPOISSONOP_H__INCLUDED__

