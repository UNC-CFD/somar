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
#include "AMRNavierStokes.H"
#include "AMRNSF_F.H"
#include "SetValLevel.H"
#include "Constants.H"


// -----------------------------------------------------------------------------
// Encapsulates the finer-level-pointer casting
// -----------------------------------------------------------------------------
AMRNavierStokes* AMRNavierStokes::fineNSPtr () const
{

    AMRNavierStokes* fine_ns_ptr;
    if (m_finer_level_ptr != NULL) {
        CH_assert(m_finer_level_ptr != NULL);

        fine_ns_ptr = dynamic_cast <AMRNavierStokes*> (m_finer_level_ptr);
        if (fine_ns_ptr == NULL) {
            MayDay::Error ("AMRNavierStokes::fineNSPtr: fineptr not castable to AMRNavierStokes*");
        }
    } else {
        fine_ns_ptr = NULL;
    }

    return fine_ns_ptr;
}


// -----------------------------------------------------------------------------
// Encapsulates the coarser-level pointer casting
// -----------------------------------------------------------------------------
AMRNavierStokes* AMRNavierStokes::crseNSPtr () const
{
    AMRNavierStokes* crse_ns_ptr;
    if (m_coarser_level_ptr != NULL) {
        CH_assert(m_coarser_level_ptr != NULL);
        crse_ns_ptr = dynamic_cast <AMRNavierStokes*> (m_coarser_level_ptr);
        if (crse_ns_ptr == NULL) {
            MayDay::Error ("AMRNavierStokes::crseNSPtr: crseptr not castable to AMRNavierStokes*");
        }
    } else {
        crse_ns_ptr = NULL;
    }

    return crse_ns_ptr;
}


// -----------------------------------------------------------------------------
// Returns a pointer to the finest level
// -----------------------------------------------------------------------------
AMRNavierStokes* AMRNavierStokes::finestNSPtr ()
{
    AMRNavierStokes* levelNSPtr = this;
    while (!levelNSPtr->finestLevel()) {
        levelNSPtr = levelNSPtr->fineNSPtr();
    }

    CH_assert(levelNSPtr->finestLevel());
    CH_assert(!(levelNSPtr->m_is_empty));

    return levelNSPtr;
}


// -----------------------------------------------------------------------------
// Returns a pointer to the finest level (const version)
// -----------------------------------------------------------------------------
const AMRNavierStokes* AMRNavierStokes::finestNSPtr () const
{
    const AMRNavierStokes* levelNSPtr = this;
    while (!levelNSPtr->finestLevel()) {
        levelNSPtr = levelNSPtr->fineNSPtr();
    }

    CH_assert(levelNSPtr->finestLevel());
    CH_assert(!(levelNSPtr->m_is_empty));

    return levelNSPtr;
}


// -----------------------------------------------------------------------------
// Returns a pointer to the coarsest level (const version)
// -----------------------------------------------------------------------------
const AMRNavierStokes* AMRNavierStokes::coarsestNSPtr () const
{
    const AMRNavierStokes* levelNSPtr = this;
    while (levelNSPtr->m_level != 0) {
        levelNSPtr = levelNSPtr->crseNSPtr();
    }
    return levelNSPtr;
}


// -----------------------------------------------------------------------------
// Returns true if this level is not in use.
// -----------------------------------------------------------------------------
bool AMRNavierStokes::isEmpty () const
{
    return m_is_empty;
}


// -----------------------------------------------------------------------------
// Returns true if this is the finest level.
// -----------------------------------------------------------------------------
bool AMRNavierStokes::finestLevel () const
{
    return m_finest_level;
}


// -----------------------------------------------------------------------------
// Sets whether this is the finest level or not.
// -----------------------------------------------------------------------------
void AMRNavierStokes::finestLevel (bool a_finest_level)
{

    m_finest_level = a_finest_level;
    // m_projection.isFinestLevel(a_finest_level);
}


// -----------------------------------------------------------------------------
// Fill this level's data holders with bogus data
// -----------------------------------------------------------------------------
void AMRNavierStokes::setAllBogus()
{
    if (s_set_bogus_values) {
        setValLevel(*m_vel_old_ptr, s_bogus_value);
        setValLevel(*m_vel_new_ptr, s_bogus_value);

        setValLevel(*m_lambda_old_ptr, s_bogus_value);
        setValLevel(*m_lambda_new_ptr, s_bogus_value);

        setValLevels(m_scal_old, 0, s_num_scal_comps-1, s_bogus_value);
        setValLevels(m_scal_new, 0, s_num_scal_comps-1, s_bogus_value);

        setValLevel(m_macPressure, s_bogus_value);
        setValLevel(m_ccPressure, s_bogus_value);
        setValLevel(m_syncPressure, s_bogus_value);

        setValLevel(m_eLambda, s_bogus_value);
        setValLevel(m_gradELambda, s_bogus_value);

        m_ccPressureState = CCPressureState::BOGUS;
        m_syncPressureState = SyncPressureState::BOGUS;
        m_eLambdaState = ELambdaState::BOGUS;
        m_gradELambdaState = GradELambdaState::BOGUS;
    }
}


// -----------------------------------------------------------------------------
// newVel
// -----------------------------------------------------------------------------
LevelData<FArrayBox>& AMRNavierStokes::newVel ()
{
    CH_assert(m_vel_new_ptr != NULL);

    return *(m_vel_new_ptr);
}


// -----------------------------------------------------------------------------
// const newVel
// -----------------------------------------------------------------------------
const LevelData<FArrayBox>& AMRNavierStokes::newVel () const
{
    CH_assert(m_vel_new_ptr != NULL);

    return *(m_vel_new_ptr);
}


// -----------------------------------------------------------------------------
// oldVel
// -----------------------------------------------------------------------------
LevelData<FArrayBox>& AMRNavierStokes::oldVel ()
{
    CH_assert(m_vel_old_ptr != NULL);

    return *(m_vel_old_ptr);
}


// -----------------------------------------------------------------------------
// const oldVel
// -----------------------------------------------------------------------------
const LevelData<FArrayBox>& AMRNavierStokes::oldVel () const
{
    CH_assert(m_vel_old_ptr != NULL);

    return *(m_vel_old_ptr);
}


// -----------------------------------------------------------------------------
// newLambda
// -----------------------------------------------------------------------------
LevelData<FArrayBox>& AMRNavierStokes::newLambda ()
{
    CH_assert(m_lambda_new_ptr != NULL);

    return *(m_lambda_new_ptr);
}


// -----------------------------------------------------------------------------
// const newLambda
// -----------------------------------------------------------------------------
const LevelData<FArrayBox>& AMRNavierStokes::newLambda () const
{
    CH_assert(m_lambda_new_ptr != NULL);

    return *(m_lambda_new_ptr);
}


// -----------------------------------------------------------------------------
// oldLambda
// -----------------------------------------------------------------------------
LevelData<FArrayBox>& AMRNavierStokes::oldLambda ()
{
    CH_assert(m_lambda_old_ptr != NULL);

    return *(m_lambda_old_ptr);
}


// -----------------------------------------------------------------------------
// const oldLambda
// -----------------------------------------------------------------------------
const LevelData<FArrayBox>& AMRNavierStokes::oldLambda () const
{
    CH_assert(m_lambda_old_ptr != NULL);

    return *(m_lambda_old_ptr);
}


// -----------------------------------------------------------------------------
// newScal
// -----------------------------------------------------------------------------
LevelData<FArrayBox>& AMRNavierStokes::newScal (const int a_comp)
{
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < m_scal_new.size());
    CH_assert(m_scal_new[a_comp] != NULL);

    return *(m_scal_new[a_comp]);
}


// -----------------------------------------------------------------------------
// const newScal
// -----------------------------------------------------------------------------
const LevelData<FArrayBox>& AMRNavierStokes::newScal (const int a_comp) const
{
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < m_scal_new.size());
    CH_assert(m_scal_new[a_comp] != NULL);

    return *(m_scal_new[a_comp]);
}

// -----------------------------------------------------------------------------
// oldScal
// -----------------------------------------------------------------------------
LevelData<FArrayBox>& AMRNavierStokes::oldScal (const int a_comp)
{
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < m_scal_old.size());
    CH_assert(m_scal_old[a_comp] != NULL);

    return *(m_scal_old[a_comp]);
}


// -----------------------------------------------------------------------------
// const oldScal
// -----------------------------------------------------------------------------
const LevelData<FArrayBox>& AMRNavierStokes::oldScal (const int a_comp) const
{
    CH_assert(0 <= a_comp);
    CH_assert(a_comp < m_scal_old.size());
    CH_assert(m_scal_old[a_comp] != NULL);

    return *(m_scal_old[a_comp]);
}


// -----------------------------------------------------------------------------
// Collects the new velocities from each level at and above a_startNSPtr.
// Lower levels are set to NULL. If a_startPtr is NULL, all data is retrieved.
// WARNING: It's up to you to make sure the base vels are at the correct time!
// -----------------------------------------------------------------------------
Vector<LevelData<FArrayBox>*> AMRNavierStokes::gatherNewVelStartingWith (AMRNavierStokes* a_startNSPtr)
{
    // Translate NULL (which is 0) to coarsest level.
    if (a_startNSPtr == NULL) {
        a_startNSPtr = this;
        while (a_startNSPtr->m_coarser_level_ptr != NULL) {
            a_startNSPtr = a_startNSPtr->crseNSPtr();
        }
    }

    // Begin with NULL pointers below starting point.
    Vector<LevelData<FArrayBox>*> retVect(a_startNSPtr->m_level, NULL);

    // Work up the levels, gathering pointers.
    while (1) {
        retVect.push_back(&(a_startNSPtr->newVel()));
        a_startNSPtr = a_startNSPtr->fineNSPtr();

        if (a_startNSPtr == NULL) break;
        if (a_startNSPtr->isEmpty()) break;
    };
    return retVect;
}


// -----------------------------------------------------------------------------
// Collects the sync pressures from each level at and above this level-1.
// Lower levels are set to NULL.
// -----------------------------------------------------------------------------
Vector<LevelData<FArrayBox>*> AMRNavierStokes::gatherSyncPressure ()
{
    Vector<LevelData<FArrayBox>*> retVect(0);

    // Collect a base level pointer.
    AMRNavierStokes* levPtr = this;
    if (m_level > 0) {
        levPtr = crseNSPtr();
    }
    CH_assert(levPtr != NULL);

    // Begin with NULL pointers below levPtr.
    int lev = 0;
    for (; lev < levPtr->m_level; ++lev) {
        retVect.push_back(NULL);
    }

    // Work up the levels, gathering pointers.
    while (1) {
        // The field will be set to zero before the solve, so as long as it's
        // well-defined, we are happy.
        CH_assert(levPtr->m_syncPressureState != SyncPressureState::UNDEFINED);

        retVect.push_back(&(levPtr->m_syncPressure));

        levPtr = levPtr->fineNSPtr();
        if (levPtr == NULL) break;
        if (levPtr->isEmpty()) break;
    };

    return retVect;
}


// -----------------------------------------------------------------------------
// Sets m_syncPressureState on this and all higher levels.
// -----------------------------------------------------------------------------
void AMRNavierStokes::setSyncPressureStates (int a_state)
{
    AMRNavierStokes* levPtr = this;
    while (1) {
        levPtr->m_syncPressureState = a_state;

        levPtr = levPtr->fineNSPtr();
        if (levPtr == NULL) break;
        if (levPtr->isEmpty()) break;
    };
}


// -----------------------------------------------------------------------------
// Collects eLambda from each level at and above this level-1.
// Lower levels are set to NULL.
// -----------------------------------------------------------------------------
Vector<LevelData<FArrayBox>*> AMRNavierStokes::gatherELambda ()
{
    Vector<LevelData<FArrayBox>*> retVect(0);

    // Collect a base level pointer.
    AMRNavierStokes* levPtr = this;
    if (m_level > 0) {
        levPtr = crseNSPtr();
    }
    CH_assert(levPtr != NULL);

    // Begin with NULL pointers below levPtr.
    int lev = 0;
    for (; lev < levPtr->m_level; ++lev) {
        retVect.push_back(NULL);
    }

    // Work up the levels, gathering pointers.
    while (1) {
        // The field will be set to zero before the solve, so as long as it's
        // well-defined, we are happy.
        CH_assert(levPtr->m_eLambdaState != ELambdaState::UNDEFINED);

        retVect.push_back(&(levPtr->m_eLambda));

        levPtr = levPtr->fineNSPtr();
        if (levPtr == NULL) break;
        if (levPtr->isEmpty()) break;
    };

    return retVect;
}


// -----------------------------------------------------------------------------
// Sets m_eLambdaState on this and all higher levels.
// -----------------------------------------------------------------------------
void AMRNavierStokes::setELambdaStates (int a_state)
{
    AMRNavierStokes* levPtr = this;
    while (1) {
        levPtr->m_eLambdaState = a_state;

        levPtr = levPtr->fineNSPtr();
        if (levPtr == NULL) break;
        if (levPtr->isEmpty()) break;
    };
}


// -----------------------------------------------------------------------------
// Sets m_gradELambdaState on this and all higher levels.
// -----------------------------------------------------------------------------
void AMRNavierStokes::setGradELambdaStates (int a_state)
{
    AMRNavierStokes* levPtr = this;
    while (1) {
        levPtr->m_gradELambdaState = a_state;

        levPtr = levPtr->fineNSPtr();
        if (levPtr == NULL) break;
        if (levPtr->isEmpty()) break;
    };
}


// -----------------------------------------------------------------------------
// Move new-time state into old-time state -- also advances time by dt
// -----------------------------------------------------------------------------
void AMRNavierStokes::swapOldAndNewStates ()
{

    LevelData<FArrayBox>* tempPtr;

    // swap both velocity and scalar data
    tempPtr = m_vel_new_ptr;
    m_vel_new_ptr = m_vel_old_ptr;
    m_vel_old_ptr = tempPtr;

    tempPtr = m_lambda_new_ptr;
    m_lambda_new_ptr = m_lambda_old_ptr;
    m_lambda_old_ptr = tempPtr;

    // loop over scalar components
    CH_assert(m_scal_new.size() == m_scal_old.size());
    CH_assert(m_scal_new.size() == s_num_scal_comps);

    for (int comp = 0; comp < m_scal_new.size(); ++comp) {
        CH_assert(m_scal_new[comp] != NULL);
        CH_assert(m_scal_old[comp] != NULL);

        tempPtr = m_scal_new[comp];
        m_scal_new[comp] = m_scal_old[comp];
        m_scal_old[comp] = tempPtr;
    }

    time(m_time + m_dt);
}


// -----------------------------------------------------------------------------
// Swap old and new states, resets time to a_time
// -----------------------------------------------------------------------------
void AMRNavierStokes::resetStates (const Real a_time)
{

    LevelData<FArrayBox>* tempPtr;
    tempPtr = m_vel_new_ptr;
    m_vel_new_ptr = m_vel_old_ptr;
    m_vel_old_ptr = tempPtr;

    tempPtr = m_lambda_new_ptr;
    m_lambda_new_ptr = m_lambda_old_ptr;
    m_lambda_old_ptr = tempPtr;

    // loop over scalar components
    CH_assert (m_scal_new.size() == m_scal_old.size());
    CH_assert(m_scal_new.size() == s_num_scal_comps);

    for (int comp = 0; comp < m_scal_new.size(); ++comp) {
        CH_assert(m_scal_new[comp] != NULL);
        CH_assert(m_scal_old[comp] != NULL);

        tempPtr = m_scal_new[comp];
        m_scal_new[comp] = m_scal_old[comp];
        m_scal_old[comp] = tempPtr;
    }

    time(a_time);
}


// -----------------------------------------------------------------------------
// computeInitialDt
// -----------------------------------------------------------------------------
Real AMRNavierStokes::computeInitialDt ()
{
    if (s_verbosity >= 5) {
        pout() << "AMRNavierStokes::computeInitialDt " << m_level << endl;
    }

    return s_init_shrink * computeDt();
}


// -----------------------------------------------------------------------------
// computeDt
// -----------------------------------------------------------------------------
Real AMRNavierStokes::computeDt ()
{
    if (s_verbosity >= 5) {
        pout() << "AMRNavierStokes::computeDt " << m_level << endl;
    }

    ostringstream limFactor;
    const int verbThresh = 3;
    pout() << std::fixed;

    Real dt;
    const Real old_dt = m_dt;

    // if prescribedDt > 0 then use that
    // (don't worry about not being on level 0 -- Amr class
    // will take care of that
    if (s_prescribedDt > 0.0) {
        dt = s_prescribedDt;

    } else {
        // Compute min(dx), max(speed), etc.
        const RealVect dx = m_levGeoPtr->getDx();
        RealVect maxVel;
        int maxVelDir = 0;
        Real maxSpeed = 0.0;
        Real minDx = dx[0];

        for (int dir = 0; dir < SpaceDim; ++dir) {
            maxVel[dir] = norm(newVel(), Interval(dir,dir), 0);

            if (maxSpeed < maxVel[dir]) {
                maxSpeed = maxVel[dir];
                maxVelDir = dir;
            }

            if (dx[dir] < minDx) {
                minDx = dx[dir];
            }
        }

        // Initialize to max allowed dt.
        limFactor << "max_dt";
        dt = s_max_dt;
        CH_assert(dt > 0.0);


        // Advective stability limit...
        for (int dir = 0; dir < SpaceDim; ++dir) {
            if (maxVel[dir] != 0.0) {
                Real advDt = m_cfl * dx[dir] / maxVel[dir];
                if (s_verbosity > verbThresh) {
                    pout() << "dt limit by advection in dir " << dir << " = " << advDt << endl;
                }
                CH_assert(advDt >= 0.0);

                if (dt > advDt) {
                    limFactor.clear();
                    limFactor.str("");
                    limFactor << "advection in dir " << dir;
                    dt = advDt;
                }
            }
        }

        // Viscous stability limit...
        if (s_limitDtViaViscosity) {
            if (s_nu > 0.0) {
                Real viscDt = m_cfl * minDx * minDx / s_nu;
                if (s_verbosity > verbThresh) {
                    pout() << "dt limit by viscosity = " << viscDt << endl;
                }
                CH_assert(viscDt >= 0.0);

                if (dt > viscDt) {
                    limFactor.clear();
                    limFactor.str("");
                    limFactor << "viscosity";
                    dt = viscDt;
                }
            }
        }

        // Diffusive stability limit...
        if (s_limitDtViaDiffusion) {
            for (int comp = 0; comp < s_num_scal_comps; ++comp) {
                if (s_scal_coeffs[comp] > 0.0) {
                    // Real diffDt = m_cfl * minDx * minDx / s_scal_coeffs[comp];
                    Real diffDt = m_cfl * minDx * minDx / s_scal_coeffs[comp];
                    if (s_verbosity > verbThresh) {
                        pout() << "dt limit by diffusion on scalar " << comp << " = " << diffDt << endl;
                    }
                    CH_assert(diffDt >= 0.0);

                    if (dt > diffDt) {
                        limFactor.clear();
                        limFactor.str("");
                        limFactor << "diffusion";
                        dt = diffDt;
                    }
                }
            }
        }

        // Limit by acceleration. This normally doesn't kick in, but if the
        // velocity is very small (initially zero), this will prevent a large
        // timestep when the velocity forcing is large. If the pressure is not
        // available, then we just skip this entire step.
        //
        // TODO: Should we add the sync pressure?
        if (s_limitDtViaPressureGradient    &&
            m_ccProjector.isPressureAvail() &&
            m_ccPressureState == CCPressureState::VALID) {

            const DisjointBoxLayout& grids = newVel().getBoxes();
            DataIterator dit = grids.dataIterator();

            // Initialize with -Grad[Pressure]
            LevelData<FArrayBox> denom(grids, SpaceDim);
            gradCCPressure(denom, -1.0);

            // -b
            LevelData<FArrayBox> gravSource(grids, SpaceDim);
            this->fillGravSource(gravSource, m_time);                       // TODO: Use all vel sources.
            m_levGeoPtr->sendToMappedBasis(gravSource);

            for (dit.reset(); dit.ok(); ++dit) {
                denom[dit] += gravSource[dit];
            }

            // Now, limit dt.
            for (int dir = 0; dir < SpaceDim; ++dir) {
                Real newMinDt = 1.0e100;
                const Real p = 0;

                Real denomNorm = norm(denom, Interval(dir, dir), p);
                if (denomNorm > 0.0) {
                    newMinDt = m_cfl * sqrt(2.0 * dx[0] / denomNorm);
                    if (s_verbosity > verbThresh) {
                        pout() << "dt limit by pressure gradient in dir "
                               << dir << " = " << newMinDt << endl;
                    }
                    CH_assert(newMinDt >= 0.0);
                } else {
                    if (s_verbosity > verbThresh) {
                        pout() << "dt limit by pressure gradient in dir "
                               << dir << " = infinity" << endl;
                    }
                }

                if (newMinDt < dt) {
                    limFactor.clear();
                    limFactor.str("");
                    limFactor << "pressure gradient";
                    dt = newMinDt;
                }

            }
        } // end if s_limitDtViaPressureGradient

        // Limit by internal wave speed
        if (s_limitDtViaInternalWaveSpeed) {
            // Create grid references, etc...
            const LevelData<FArrayBox>& vel = this->newVel();
            const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
            DataIterator dit = grids.dataIterator();

            // Get the internal wave speed
            LevelData<FArrayBox> c0i(grids, SpaceDim);
            this->fillInternalWaveSpeed(c0i);

            // Loop over grids and compute minDt.
            Real newMinDt = dt;
            for (dit.reset(); dit.ok(); ++dit) {
                const FArrayBox& c0iFAB = c0i[dit];
                const FArrayBox& velFAB = vel[dit];
                const Box& valid = grids[dit];

                FORT_COMPUTEMINBVDT(
                    CHF_REAL(newMinDt),
                    CHF_CONST_FRA(c0iFAB),
                    CHF_CONST_FRA(velFAB),
                    CHF_CONST_REALVECT(dx),
                    CHF_BOX(valid));

                CH_assert(newMinDt >= 0.0);
            }
            newMinDt = m_cfl * newMinDt;

            if (s_verbosity > verbThresh) {
                pout() << "dt limit by internal wave speed = " << newMinDt << endl;
            }

            if (newMinDt < dt) {
                limFactor.clear();
                limFactor.str("");
                limFactor << "internal wave speed";
                dt = newMinDt;
            }
        } // end if s_limitDtViaInternalWaveSpeed

        // check to see if boundary conditions restrict dt
        m_physBCPtr->computeBoundaryDt(dt, m_cfl, *m_levGeoPtr);
        CH_assert(dt >= 0.0);
    }

    // Limit by old_dt * s_max_dt_grow
    if (old_dt > 0.0) {
        if (s_verbosity > verbThresh) {
            pout() << "dt limit by old_dt * max_dt_grow = " << old_dt * s_max_dt_grow << endl;
        }

        if (dt > old_dt * s_max_dt_grow) {
            limFactor.clear();
            limFactor.str("");
            limFactor << "max_dt_grow";
            dt = old_dt * s_max_dt_grow;
        }
    }

    // If this fails, check your dt limiting code!
    CH_assert(dt <= s_max_dt);

#ifdef CH_MPI
    // Communicate: compute the minimum over all procs.
    Real recv;
    int result = MPI_Allreduce(&dt, &recv, 1, MPI_CH_REAL,
                               MPI_MIN, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeDt");
    }

    dt = recv;
#endif

    if (s_verbosity >= verbThresh) {
        pout() << "-- dt limiting factor: " << limFactor.str() << endl;
    }

    CH_assert(dt > 0.0);
    return dt;
}
