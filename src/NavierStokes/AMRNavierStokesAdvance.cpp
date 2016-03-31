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
#include <iomanip>
#include <fstream>

#include "AMRNavierStokes.H"
#include "memusage.H"
#include "ProblemContext.H"

#ifndef CH_NTIMER
#   include "OldTimer.H"
#endif



// -----------------------------------------------------------------------------
// Advance solution by one timestep. Return max safe timestep.
// -----------------------------------------------------------------------------
Real AMRNavierStokes::advance ()
{
    if (s_verbosity >= 5) {
        pout() << "AMRNavierStokes::advance" << endl;
    }

    {// Has the halt flag been thrown? If so, abort!
        ifstream ifile("halt");
        if (ifile) {
            remove("halt");
            MayDay::Abort("Halt flag was set. Terminating run...");
        }
    }

#   ifndef NDEBUG
        // We are about to update the dynamical variables. At this point, we
        // better have valid pressure data.
        if (s_initial_pressure_iters > 0) {
            if (m_ccPressureState != CCPressureState::VALID) {
                ostringstream msg;
                msg << "Pressure on level " << m_level << " is in state " << m_ccPressureState;
                MayDay::Error(msg.str().c_str());
            }
            if (m_level > 0) {
                const AMRNavierStokes* crsePtr = crseNSPtr();
                if (crsePtr->m_ccPressureState != CCPressureState::VALID) {
                    ostringstream msg;
                    msg << "Pressure on level " << crsePtr->m_level << " is in state " << m_ccPressureState;
                    MayDay::Error(msg.str().c_str());
                }
            }
        }

        if (s_applyFreestreamCorrection) {
            if (m_gradELambdaState != GradELambdaState::VALID) {
                ostringstream msg;
                msg << "GradELambda on level " << m_level << " is in state " << m_gradELambdaState;
                MayDay::Error(msg.str().c_str());
            }
        }
#   endif

#   ifndef CH_NTIMER
        OldTimer TimeAdvanceStep;
        TimeAdvanceStep.start();
#   endif

    if (s_verbosity >= 2) {
        pout() << "AMRNavierStokes::advance " << m_level
               << ", starting time = "
               << setiosflags(ios::fixed) << setprecision(6)
               << setw(12) << m_time
               << ", dt = "
               << setiosflags(ios::fixed) << setprecision(6) << setw(12) << dt()
               << endl;
    }

    if (s_verbosity >= 4) {
        const DisjointBoxLayout& grids = newVel().getBoxes();
        DataIterator dit = grids.dataIterator();
        Vector<Box> vbox = grids.boxArray();
        Vector<int> vproc = grids.procIDs();

        pout () << "processor map: " << endl;
        for (int idx = 0; idx < vbox.size(); ++idx) {
            pout() << idx << ": " << vproc[idx] << "\t" << vbox[idx] << endl;
        }
        pout() << endl;

        pout() << "This proc = " << procID() << ":" << endl;
        for (dit.reset(); dit.ok(); ++dit) {
            pout() << grids[dit] << endl;
        }
        pout() << endl;
    }

    // The "new" state is now the old state.
    const Real old_time = m_time;
    this->swapOldAndNewStates();
    // const Real new_time = m_time; // unused, so commented.

    // Initialize all flux registers
    if (!finestLevel()) {
        m_vel_flux_reg.setToZero();
        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            m_scal_fluxreg_ptrs[comp]->setToZero();
        }
        m_lambda_flux_reg.setToZero();
    }

    if (s_updateScheme == ProblemContext::UpdateScheme::FiniteVolume) {
        if (s_gravityMethod == ProblemContext::GravityMethod::IMPLICIT) {
            this->PPMIGTimeStep(old_time,
                                m_dt,
                                true,      // updatePassiveScalars
                                true);     // doLevelProj
        } else {
            this->PPMTimeStep(old_time,
                              m_dt,
                              true,      // updatePassiveScalars
                              true);     // doLevelProj
        }

    } else if (s_updateScheme == ProblemContext::UpdateScheme::RK3) {
        this->RK3TimeStep(old_time,
                          m_dt,
                          true,      // updatePassiveScalars
                          true);     // doLevelProj
    } else {
        MayDay::Error("Unrecognized update scheme");
    }


    // compute maximum safe timestep for next iteration
    Real newDt = this->computeDt();

#   ifndef CH_NTIMER
        TimeAdvanceStep.stop();
#   endif

    if (s_verbosity >= 2) {
        pout() << "AMRNavierStokes::advance " << m_level
               << ",      end time = " << setiosflags(ios::fixed) << setprecision(6) << setw(12) << m_time
               << ", dt = " << setiosflags(ios::fixed) << setprecision(6) << setw(12) << m_dt
#              ifndef CH_NTIMER
                   << "  wc = " << setiosflags(ios::fixed) << setprecision(4) << setw(12) << TimeAdvanceStep.wc_time()
#              endif
               << setiosflags(ios::scientific) << setprecision(8) << endl;
    }

    if (s_verbosity >= 2) {
        const Real memory = get_memory_usage_from_OS();
        pout() << "Current memory usage = " << memory << "MB\n" << endl;
    }

    return newDt;
}


// -----------------------------------------------------------------------------
// Increments the coarse and fine registers
// -----------------------------------------------------------------------------
void AMRNavierStokes::updateLambdaFluxRegister (LevelData<FluxBox>& a_flux,
                                                const Real          a_scale)
{
    CH_TIME("AMRNavierStokes::updateLambdaFluxRegister");

    // Sanity checks and basic references
    CH_assert(s_advective_lambda_reflux);
    CH_assert(m_levGeoPtr != NULL);

    DataIterator dit = a_flux.dataIterator();
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    CH_assert(a_flux.getBoxes().compatible(grids));

    // Reference the flux registers
    MappedLevelFluxRegister* fluxRegPtr = &m_lambda_flux_reg;
    MappedLevelFluxRegister* crseFluxRegPtr = NULL;
    if (m_level > 0) {
        crseFluxRegPtr = &(crseNSPtr()->m_lambda_flux_reg);
    }

    for (dit.reset(); dit.ok(); ++dit) {
        FluxBox& flux = a_flux[dit];
        const Interval interv(0,0);

        // Increment the coarse side of the fine register
        if (!finestLevel()) {
            CH_assert(fluxRegPtr != NULL);
            CH_assert(fluxRegPtr->isDefined());
            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real dx = m_levGeoPtr->getDx()[dir];
                fluxRegPtr->incrementCoarse(flux[dir],
                                            -a_scale / dx,
                                            dit(),
                                            interv,
                                            interv,
                                            dir);
            }
        }

        // Increment the fine side of the coarse register
        if (m_level > 0) {
            CH_assert(crseFluxRegPtr != NULL);
            CH_assert(crseFluxRegPtr->isDefined());
            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real dx = crseNSPtr()->m_levGeoPtr->getDx()[dir];
                crseFluxRegPtr->incrementFine(flux[dir],
                                              -a_scale / dx,
                                              dit(),
                                              interv,
                                              interv,
                                              dir);
            }
        }
    } // end loop over grids
}


// -----------------------------------------------------------------------------
// Increments the coarse and fine registers
// -----------------------------------------------------------------------------
void AMRNavierStokes::updateScalarFluxRegister (LevelData<FluxBox>& a_flux,
                                                const Real          a_dtScale,
                                                const int           a_scalComp)
{
    CH_TIME("AMRNavierStokes::updateScalarFluxRegister");

    // Sanity checks and basic references
    CH_assert(s_advective_scalar_reflux);
    CH_assert(m_levGeoPtr != NULL);
    CH_assert(0 <= a_scalComp);
    CH_assert(a_scalComp < s_num_scal_comps);

    DataIterator dit = a_flux.dataIterator();
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    CH_assert(a_flux.getBoxes().compatible(grids));

    // Reference the flux registers
    MappedLevelFluxRegister* fluxRegPtr = m_scal_fluxreg_ptrs[a_scalComp];
    MappedLevelFluxRegister* crseFluxRegPtr = NULL;
    if (m_level > 0) {
        crseFluxRegPtr = crseNSPtr()->m_scal_fluxreg_ptrs[a_scalComp];
    }

    for (dit.reset(); dit.ok(); ++dit) {
        FluxBox& flux = a_flux[dit];
        const Interval interv(0,0);

        // Increment the coarse side of the fine register
        if (!finestLevel()) {
            CH_assert(fluxRegPtr != NULL);
            CH_assert(fluxRegPtr->isDefined());
            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real dx = m_levGeoPtr->getDx()[dir];
                fluxRegPtr->incrementCoarse(flux[dir],
                                            -a_dtScale / dx,
                                            dit(),
                                            interv,
                                            interv,
                                            dir);
            }
        }

        // Increment the fine side of the coarse register
        if (m_level > 0) {
            CH_assert(crseFluxRegPtr != NULL);
            CH_assert(crseFluxRegPtr->isDefined());
            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real dx = crseNSPtr()->m_levGeoPtr->getDx()[dir];
                crseFluxRegPtr->incrementFine(flux[dir],
                                              -a_dtScale / dx,
                                              dit(),
                                              interv,
                                              interv,
                                              dir);
            }
        }
    } // end loop over grids
}


// -----------------------------------------------------------------------------
// Increments the coarse and fine registers
// -----------------------------------------------------------------------------
void AMRNavierStokes::updateVelFluxRegister (LevelData<FluxBox>& a_flux,
                                             const Real          a_scale)
{
    CH_TIME("AMRNavierStokes::updateVelFluxRegister");

    // Sanity checks and basic references
    CH_assert(s_advective_momentum_reflux);
    CH_assert(m_levGeoPtr != NULL);

    DataIterator dit = a_flux.dataIterator();
    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();
    CH_assert(a_flux.getBoxes().compatible(grids));

    // Reference the flux registers
    MappedLevelFluxRegister* fluxRegPtr = &(m_vel_flux_reg);
    MappedLevelFluxRegister* crseFluxRegPtr = NULL;
    if (m_level > 0) {
        crseFluxRegPtr = &(crseNSPtr()->m_vel_flux_reg);
    }

    for (dit.reset(); dit.ok(); ++dit) {
        FluxBox& flux = a_flux[dit];
        const Interval interv(0,SpaceDim-1);
        const bool reflux_normal_momentum = true;

        // Increment the coarse side of the fine register
        if (!finestLevel()) {
            CH_assert(fluxRegPtr != NULL);
            CH_assert(fluxRegPtr->isDefined());

            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real dx = m_levGeoPtr->getDx()[dir];

                if (reflux_normal_momentum) {
                    fluxRegPtr->incrementCoarse(flux[dir],
                                                -a_scale / dx,
                                                dit(),
                                                interv,
                                                interv,
                                                dir);
                } else {
                    // if we're not refluxing normal momentum component,
                    // then do this component by component
                    for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                        if (velComp == dir) continue;

                        Interval velInt(velComp, velComp);
                        fluxRegPtr->incrementCoarse(flux[dir],
                                                    -a_scale / dx,
                                                    dit(),
                                                    velInt,
                                                    velInt,
                                                    dir);
                    } // end loop over velocity components
                } // end if we're not refluxing normal momentum
            } // end loop over directions
        }

        // Increment the fine side of the coarse register
        if (m_level > 0) {
            CH_assert(crseFluxRegPtr != NULL);
            CH_assert(crseFluxRegPtr->isDefined());

            for (int dir = 0; dir < SpaceDim; ++dir) {
                const Real dx = crseNSPtr()->m_levGeoPtr->getDx()[dir];

                if (reflux_normal_momentum) {
                    crseFluxRegPtr->incrementFine(flux[dir],
                                                  -a_scale / dx,
                                                  dit(),
                                                  interv,
                                                  interv,
                                                  dir);
                } else {
                    // for case where we only want to reflux tangential
                    // momentum components, do this component by component
                    // so that we can isolate the normal component.
                    for (int velComp = 0; velComp < SpaceDim; ++velComp) {
                        if (velComp == dir) continue;

                        Interval velInt(velComp, velComp);
                        crseFluxRegPtr->incrementFine(flux[dir],
                                                      -a_scale / dx,
                                                      dit(),
                                                      velInt,
                                                      velInt,
                                                      dir);
                    } // end loop over velocity components
                } // end if we're not refluxing normal momentum
            } // end loop over directions
        } // end if not base level
    } // end loop over grids
}

