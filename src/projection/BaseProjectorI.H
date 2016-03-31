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
#include "Debug.H"


// -----------------------------------------------------------------------------
// Default constructor -- sets the default solver parameters.
// -----------------------------------------------------------------------------
template <class FluxType>
BaseProjector<FluxType>::
BaseProjector ()
: m_time(BOGUS_TIME),
  m_pressure(Vector<LevelData<FArrayBox>*>(0))
{;}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
template <class FluxType>
BaseProjector<FluxType>::
~BaseProjector ()
{
    m_pressure.clear();
    m_solver.undefine();
}


// -----------------------------------------------------------------------------
// Sets pointers to the AMR pressure fields.
// -----------------------------------------------------------------------------
template <class FluxType>
void
BaseProjector<FluxType>::
setAMRPressure (Vector<LevelData<FArrayBox>*>& a_amrPhi)
{
    m_pressure = a_amrPhi;
}


// -----------------------------------------------------------------------------
// Sets pointers to the pressure fields for single-level solves.
// -----------------------------------------------------------------------------
template <class FluxType>
void
BaseProjector<FluxType>::
setLevelPressure (LevelData<FArrayBox>*       a_phiPtr,
                  const LevelData<FArrayBox>* a_crsePhiPtr)
{
    m_pressure.clear();
    if (a_crsePhiPtr != NULL) {
        m_pressure.push_back((LevelData<FArrayBox>*)a_crsePhiPtr);
    }
    m_pressure.push_back(a_phiPtr);
}


// -----------------------------------------------------------------------------
// Computes Jg^{i,j}\partial_j(phi) on a single level. This simply calls
// the more general AMR version, which needs to be overridden. If your fluid
// is compressible (for testing purposes), set a_scale to zero.
// NOTE: m_pressure.size() must be 1 or 2!
// -----------------------------------------------------------------------------
template <class FluxType>
void
BaseProjector<FluxType>::
computeLevelGradPressure (LevelData<FluxType>& a_gradPressure,
                          const Real           a_scale) const
{
    // Sanity checks and initialization for debugging.
    CH_assert(m_pressure.size() == 1 || m_pressure.size() == 2);
    CH_assert(m_pressure[0] != NULL);
    CH_assert(m_pressure.back() != NULL);
    debugInitLevel(a_gradPressure);

    // For compressible flows, m_pressure can be undefined.
    DataIterator dit = a_gradPressure.dataIterator();
    if (a_scale == 0.0) {
        for (dit.reset(); dit.ok(); ++dit) {
            a_gradPressure[dit].setVal(0.0);
        }
        return;
    }

    // Compute grad
    const int lmax = m_pressure.size() - 1;
    Vector<LevelData<FluxType>*> vGrad(lmax, NULL);
    vGrad.push_back(&a_gradPressure);

    computeGrad(vGrad, m_pressure, lmax, lmax);

    // Scale
    for (dit.reset(); dit.ok(); ++dit) {
        a_gradPressure[dit] *= a_scale;
    }
}


// -----------------------------------------------------------------------------
// AMRNavierStokes::computeInitDt may simply call computeDt. If that
// happens, the projectors will not have been defined yet and
// grad[pressure] is not available. This function simply looks to see
// if m_pressure contains any pointers. It does not test for the
// validity of those pointers or the data they point to.
// -----------------------------------------------------------------------------
template <class FluxType>
bool
BaseProjector<FluxType>::
isPressureAvail () const
{
    return (m_pressure.size() > 0);
}



// -----------------------------------------------------------------------------
// Do the projection!
// a_lmin is the coarsest projected level -- not the base as in Chombo!
// Set a_velIsFlux if a_amrVel has been multiplied by J.
// a_amrVel MUST be indexed to match m_pressure.
// -----------------------------------------------------------------------------
template <class FluxType>
void
BaseProjector<FluxType>::
project (Vector<LevelData<FluxType>*> a_amrVel,
         const LevelGeometry&         a_levGeo,
         const int                    a_lmin,
         const int                    a_lmax,
         const Real                   a_newTime,
         const Real                   a_dt,
         const bool                   a_velIsFlux,
         const bool                   a_zeroPressure,
         const bool                   a_forceHomogeneous)
{
    project(a_amrVel,
            a_levGeo.getAMRLevGeos(),
            a_lmin,
            a_lmax,
            a_newTime,
            a_dt,
            a_velIsFlux,
            a_zeroPressure,
            a_forceHomogeneous);
}


// -----------------------------------------------------------------------------
// Do the projection!
// a_lmin is the coarsest projected level -- not the base as in Chombo!
// Set a_velIsFlux if a_amrVel has been multiplied by J.
// -----------------------------------------------------------------------------
template <class FluxType>
void
BaseProjector<FluxType>::
project (Vector<LevelData<FluxType>*> a_amrVel,
         Vector<const LevelGeometry*> a_amrLevGeo,
         const int                    a_lmin,
         const int                    a_lmax,
         const Real                   a_newTime,
         const Real                   a_dt,
         const bool                   a_velIsFlux,
         const bool                   a_zeroPressure,
         const bool                   a_forceHomogeneous)
{
    // Compute the base level. This is the level needed for CF-BCs.
    int lbase = a_lmin;
    if (lbase > 0) --lbase;

#ifndef NDEBUG
    // Do we have all of the data we need to project?
    CH_assert(lbase >= 0);
    CH_assert(a_lmin <= a_lmax);
    CH_assert(a_lmax < a_amrVel.size());
    CH_assert(a_lmax < a_amrLevGeo.size());
    CH_assert(a_lmax < m_pressure.size());
    CH_assert(m_solver.isDefined());
    CH_assert(m_solver.getNumLevels() > a_lmax-lbase);

    for (int lev = a_lmin; lev <= a_lmax; ++lev) {
        // Everything must be well-defined
        CH_assert(a_amrVel[lev] != NULL);
        CH_assert(a_amrVel[lev]->isDefined());

        CH_assert(a_amrLevGeo[lev] != NULL);
        CH_assert(a_amrLevGeo[lev]->isDefined());
        CH_assert(a_amrLevGeo[lev]->getBoxes() == a_amrVel[lev]->getBoxes());

        CH_assert(m_pressure[lev] != NULL);
        CH_assert(m_pressure[lev]->isDefined());
        CH_assert(m_pressure[lev]->getBoxes() == a_amrVel[lev]->getBoxes());
    }

    if (lbase < a_lmin) {
        // We only need the pressure here.
        CH_assert(m_pressure[lbase] != NULL);
        CH_assert(m_pressure[lbase]->isDefined());

        // If the velocity is defined here, then check everything.
        // NOTE: The only reason we would ever need the vel at lbase is to
        // compute CFBCs on a CC velocity before taking the divergence. And
        // to do that, we may need the levGeo to mult/div by J.
        if (a_amrVel[lbase] != NULL) {
            CH_assert(a_amrLevGeo[lbase] != NULL);
            CH_assert(a_amrLevGeo[lbase]->isDefined());
            CH_assert(m_pressure[lbase]->getBoxes() == a_amrLevGeo[lbase]->getBoxes());

            CH_assert(a_amrVel[lbase]->isDefined());
            CH_assert(m_pressure[lbase]->getBoxes() == a_amrVel[lbase]->getBoxes());
        }
    }
#endif

    // Make sure a_amrVel is scaled by J
    if (!a_velIsFlux) {
        if ((lbase < a_lmin) && (a_amrVel[lbase] != NULL)) {
            a_amrLevGeo[lbase]->multByJ(*a_amrVel[lbase]);
        }
        for (int lev = a_lmin; lev <= a_lmax; ++lev) {
            a_amrLevGeo[lev]->multByJ(*a_amrVel[lev]);
        }
    }

    // Compute: amrRhs = Div[amrVel/amrDt]
    Vector<LevelData<FArrayBox>*> amrRhs(a_amrVel.size(), NULL);
    for (int lev = a_lmin; lev <= a_lmax; ++lev) {
        amrRhs[lev] = new LevelData<FArrayBox>(a_amrVel[lev]->getBoxes(), 1);
    }
    computeDiv(amrRhs, a_amrVel, a_lmin, a_lmax);
    if (a_dt != 0.0) {
        for (int lev = a_lmin; lev <= a_lmax; ++lev) {
            DataIterator dit = amrRhs[lev]->dataIterator();
            for (dit.reset(); dit.ok(); ++dit) {
                (*amrRhs[lev])[dit] /= a_dt;
            }
        }
    }

    // Solve: L[pressure] = amrRhs
    m_solver.solve(m_pressure,
                   amrRhs,
                   a_lmin,
                   a_lmax,
                   a_zeroPressure,
                   a_forceHomogeneous);

    // Free rhs
    for (int lev = a_lmin; lev <= a_lmax; ++lev) {
        delete amrRhs[lev];
        amrRhs[lev] = NULL;
    }
    amrRhs.clear();

    // Correct: amrVel = amrVel - amrDt*Grad[pressure]
    Vector<LevelData<FluxType>*> amrCorr(a_amrVel.size(), NULL);
    for (int lev = a_lmin; lev <= a_lmax; ++lev) {
        amrCorr[lev] = new LevelData<FluxType>(a_amrVel[lev]->getBoxes(), a_amrVel[lev]->nComp());
    }

    computeGrad(amrCorr, m_pressure, a_lmin, a_lmax);
    applyCorrection(a_amrVel, amrCorr, a_dt, a_lmin, a_lmax);

    // Free grad[pressure]
    for (int lev = a_lmin; lev <= a_lmax; ++lev) {
        delete amrCorr[lev];
        amrCorr[lev] = NULL;
    }
    amrCorr.clear();

    // Undo the scaling, if necessary.
    if (!a_velIsFlux) {
        if ((lbase < a_lmin) && (a_amrVel[lbase] != NULL)) {
            a_amrLevGeo[lbase]->divByJ(*a_amrVel[lbase]);
        }
        for (int lev = a_lmin; lev <= a_lmax; ++lev) {
            a_amrLevGeo[lev]->divByJ(*a_amrVel[lev]);
        }
    }
}


// -----------------------------------------------------------------------------
// Do the projection!
// Set a_velIsFlux if a_amrVel has been multiplied by J.
// -----------------------------------------------------------------------------
template <class FluxType>
void
BaseProjector<FluxType>::
levelProject (Vector<LevelData<FluxType>*> a_amrVel,
              const LevelGeometry*         a_levGeoPtr,
              const Real                   a_newTime,
              const Real                   a_dt,
              const bool                   a_velIsFlux,
              const bool                   a_zeroPressure,
              const bool                   a_forceHomogeneous)
{
    CH_assert(a_amrVel.size() == 1 || a_amrVel.size() == 2);
    CH_assert(a_amrVel.back() != NULL);
    CH_assert(a_amrVel.back()->isDefined());

    CH_assert(a_levGeoPtr != NULL);
    CH_assert(a_levGeoPtr->isDefined());
    CH_assert(a_levGeoPtr->getBoxes() == a_amrVel.back()->getBoxes());

    Vector<const LevelGeometry*> amrLevGeo(0);
    int lmin, lmax;

    if (m_pressure.size() == 1) {
        CH_assert(a_amrVel.size() == 1);
        amrLevGeo.push_back(a_levGeoPtr);
        lmin = 0;
        lmax = 0;

    } else {
        CH_assert(m_pressure.size() == 2);

        amrLevGeo.push_back(a_levGeoPtr->getCoarserPtr());
        amrLevGeo.push_back(a_levGeoPtr);

        lmin = 1;
        lmax = 1;

        if (a_amrVel.size() == 1) {
#ifndef NDEBUG
            // TODO: Why did I write this and what does it mean?!
            // I'll allow this only for FC velocity holders since we won't need
            // to interpolate ghosts.
            FluxType fluxDummy;
            CH_assert(dynamic_cast<FluxBox*>(&fluxDummy) != NULL);
#endif

            a_amrVel.push_back(a_amrVel[0]);
            a_amrVel[0] = NULL;
        }
    }

    project(a_amrVel,
            amrLevGeo,
            lmin,
            lmax,
            a_newTime,
            a_dt,
            a_velIsFlux,
            a_zeroPressure,
            a_forceHomogeneous);
}


// -----------------------------------------------------------------------------
// Sets the time used to evaluate BCs.
// -----------------------------------------------------------------------------
template <class FluxType>
void
BaseProjector<FluxType>::
setTime (const Real a_time)
{
    m_time = a_time;
}


// -----------------------------------------------------------------------------
// Gets the time used to evaluate BCs.
// -----------------------------------------------------------------------------
template <class FluxType>
Real
BaseProjector<FluxType>::
getTime () const
{
    return m_time;
}
