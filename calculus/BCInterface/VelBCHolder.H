/*******************************************************************************
 *  FILE: VelBCHolder.H
 *  DESC: BC interface for the VelocityAMRPoissonOp class.
 *  DATE: Some time around Thu 24 May 2012, I think.
 *  MAIL: santilli@physics.unc.edu
 *  NOTE: This is a rework of VelBCHolder in Chombo's AMRINS example.
 *        Make sure this class DOES NOT own level-dependent data.
 ******************************************************************************/
#ifndef __VELBCHOLDER_H__INCLUDED__
#define __VELBCHOLDER_H__INCLUDED__

#include "EllipticBCInterface.H"
#include "LevelData.H"
#include "Tuple.H"


// This is a physical BC class designed to handle all velocity components.
// It is used by VelocityAMRPoissonOp.
//
// TODO: I don't like this. Either BCMethodHolder should be altered to work
// properly on multi-component datasets, or this should specialize (be a child
// of) BCMethodHolder. When SOMAR gets revamped, this class will not survive.
class VelBCHolder
{
public:

	///
	VelBCHolder();

	///
	VelBCHolder(const Tuple<BCMethodHolder, SpaceDim>& a_componentBC);

	///
	virtual ~VelBCHolder();

	///
	virtual void setGhosts(FArrayBox&           a_state,
	                       const FArrayBox*     a_extrapPtr,
                           const Box&           a_valid,
                           const ProblemDomain& a_domain,
                           const RealVect&      a_dx,
                           const DataIndex&     a_index,
                           const FluxBox*       a_JgupPtr,
                           bool                 a_homogeneous,
                           Real                 a_time = BOGUS_TIME);

	///
	virtual void setGhosts(LevelData<FArrayBox>&       a_state,
	                       const LevelData<FArrayBox>* a_extrapPtr,
                           const RealVect&             a_dx,
                           const LevelData<FluxBox>*   a_JgupPtr,
                           bool                        a_homogeneous,
                           Real                        a_time = BOGUS_TIME);

	///
	virtual void setFluxes(FluxBox&             a_state,
	                       const FluxBox*       a_extrapPtr,
                           const Box&           a_validCC,
                           const ProblemDomain& a_domain,
                           const RealVect&      a_dx,
                           const DataIndex&     a_index,
                           const FluxBox*       a_JgupPtr,
                           bool                 a_homogeneous,
                           Real                 a_time = BOGUS_TIME);

protected:
	///
	Tuple<BCMethodHolder, SpaceDim> m_componentBC;
};


// This is a physical BC class designed to handle velocities on edges
// This class handles velocities, which are by their nature multicomponent
class EdgeVelBCHolder
{
public:

	///
	EdgeVelBCHolder();

	///
	EdgeVelBCHolder(const Tuple<BCMethodHolder, SpaceDim>& a_componentBC);

	///
	virtual ~EdgeVelBCHolder();

	///
	virtual void setGhosts(LevelData<FluxBox>&       a_state,
	                       const LevelData<FluxBox>* a_extrapPtr,
                           const RealVect&           a_dx,
	                       const LevelData<FluxBox>* a_JgupPtr,
	                       bool                      a_homogeneous,
                           Real                      a_time = BOGUS_TIME);

protected:
	///
	Tuple<BCMethodHolder, SpaceDim> m_componentBC;
};


#endif //!__VELBCHOLDER_H__INCLUDED__

