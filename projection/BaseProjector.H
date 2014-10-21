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
#ifndef BaseProjector_H__INCLUDED__
#define BaseProjector_H__INCLUDED__

#include "AMRPressureSolver.H"

extern const Real BOGUS_TIME;


// -----------------------------------------------------------------------------
// An interface for all projectors.
// -----------------------------------------------------------------------------
template <class FluxType>
class BaseProjector
{
public:
    // Default constructor
    BaseProjector ();

    // Destructor
    virtual ~BaseProjector ();

    // Sets pointers to the AMR pressure fields.
    virtual void setAMRPressure (Vector<LevelData<FArrayBox>*>& a_amrPhi);

    // Sets pointers to the pressure fields for single-level solves.
    virtual void setLevelPressure (LevelData<FArrayBox>*       a_phiPtr,
                                   const LevelData<FArrayBox>* a_crsePhiPtr = NULL);

    // Do the projection!
    // a_lmin is the coarsest projected level -- not the base as in Chombo!
    // Set a_velIsFlux if a_amrVel has been multiplied by J.
    // a_amrVel MUST be indexed to match m_pressure.
    virtual void project (Vector<LevelData<FluxType>*> a_amrVel,
                          const LevelGeometry&         a_levGeo,
                          const int                    a_lmin,
                          const int                    a_lmax,
                          const Real                   a_newTime,
                          const Real                   a_dt,
                          const bool                   a_velIsFlux,
                          const bool                   a_zeroPressure,
                          const bool                   a_forceHomogeneous);

    // Do the projection!
    // a_lmin is the coarsest projected level -- not the base as in Chombo!
    // Set a_velIsFlux if a_amrVel has been multiplied by J.
    // a_amrVel and a_amrLevGeo MUST be indexed to match m_pressure.
    virtual void project (Vector<LevelData<FluxType>*> a_amrVel,
                          Vector<const LevelGeometry*> a_amrLevGeo,
                          const int                    a_lmin,
                          const int                    a_lmax,
                          const Real                   a_newTime,
                          const Real                   a_dt,
                          const bool                   a_velIsFlux,
                          const bool                   a_zeroPressure,
                          const bool                   a_forceHomogeneous);

    // Do the projection!
    // Set a_velIsFlux if a_amrVel has been multiplied by J.
    virtual void levelProject (Vector<LevelData<FluxType>*> a_amrVel,
                               const LevelGeometry*         a_levGeoPtr,
                               const Real                   a_newTime,
                               const Real                   a_dt,
                               const bool                   a_velIsFlux,
                               const bool                   a_zeroPressure,
                               const bool                   a_forceHomogeneous);

    // Computes Jg^{i,j}\partial_j(phi) on a single level. This simply calls
    // the more general AMR version, which needs to be overridden. If your fluid
    // is compressible (for testing purposes), set a_scale to zero.
    // NOTE: m_pressure.size() must be 1 or 2!
    virtual void computeLevelGradPressure (LevelData<FluxType>& a_flux,
                                           const Real           a_scale = 1.0) const;

    // AMRNavierStokes::computeInitDt may simply call computeDt. If that
    // happens, the projectors will not have been defined yet and
    // grad[pressure] is not available. This function simply looks to see
    // if m_pressure contains any pointers. It does not test for the
    // validity of those pointers or the data they point to.
    virtual bool isPressureAvail () const;

protected:
    // Sets the time used to evaluate BCs.
    virtual inline void setTime (const Real a_time);

    // Gets the time used to evaluate BCs.
    virtual inline Real getTime () const;

    // Computes J^{-1}\partial_i(J u^i) over an AMR hierarchy.
    // This must be overriden or an error will be thrown.
    virtual void computeDiv (Vector<LevelData<FArrayBox>*>&      a_div,
                             const Vector<LevelData<FluxType>*>& a_flux,
                             const int                           a_lmin,
                             const int                           a_lmax) const = 0;

    // Computes Jg^{i,j}\partial_j(phi) over an AMR hierarchy.
    // This must be overriden or an error will be thrown.
    virtual void computeGrad (Vector<LevelData<FluxType>*>&        a_flux,
                              const Vector<LevelData<FArrayBox>*>& a_phi,
                              const int                            a_lmin,
                              const int                            a_lmax) const = 0;

    // Applies vel = vel - dt*corr.
    virtual void applyCorrection (Vector<LevelData<FluxType>*>&       a_amrVel,
                                  const Vector<LevelData<FluxType>*>& a_amrCorr,
                                  const Real                          a_dt,
                                  const int                           a_lmin,
                                  const int                           a_lmax) const = 0;

    // The time used to evaluate BCs.
    Real m_time;

    // The pressure field
    Vector<LevelData<FArrayBox>*> m_pressure;

    // This solves the Poisson eq.
    AMRPressureSolver m_solver;

private:
    // Copy and assignment not allowed
    BaseProjector (const BaseProjector&);
    void operator= (const BaseProjector&);
};


#include "BaseProjectorI.H"

#endif //!BaseProjector_H__INCLUDED__
