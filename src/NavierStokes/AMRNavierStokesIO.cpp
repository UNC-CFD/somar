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
#include "Divergence.H"
#include "ExtrapolationUtils.H"
#include "AMRNSF_F.H"
#include "BoxIterator.H"
#include "EdgeToCell.H"
#include "computeMappedNorm.H"
#include "SetValLevel.H"
#include "TaylorGreenBCUtil.H"

#ifdef CH_USE_HDF5
#include "CH_HDF5.H"


// -----------------------------------------------------------------------------
// Writes the checkpoint metadata.
// The metadata written by this function will only be used to perform a sanity
// check on the checkpoint file.
// -----------------------------------------------------------------------------
void AMRNavierStokes::writeCheckpointHeader(HDF5Handle& a_handle) const
{
    if (s_verbosity >= 5) {
        pout() << "AMRNavierStokes::writeCheckpointHeader" << endl;
    }

    // Create the header. This will only store metadata
    // about the number of fields and thier names.
    HDF5HeaderData header;

    // Scalar metadata...
    header.m_int["num_components"] = s_num_scal_comps;
    char comp_str[30];
    for (int comp = 0; comp < s_num_scal_comps; ++comp) {
        sprintf (comp_str, "component_%d", comp);
        header.m_string[comp_str] = s_scal_names[comp];
    }

    // Lambda metadata...
    header.m_string["lambda_component"] = "lambda";

    // Velocity metadata...
    for (int comp = 0; comp < CH_SPACEDIM; ++comp) {
        sprintf (comp_str, "vel_component_%d", comp);
        header.m_string[comp_str] = s_vel_names[comp];
    }

    // Pressure
    header.m_string["ccPressure_component"] = "ccPressure";

    // eLambda
    header.m_string["eLambda_component"] = "eLambda";

    // Write the metadata to HDF5 and pout.*
    header.writeToFile(a_handle);
    if (s_verbosity >= 3) {
          pout () << header << endl;
    }
}


// -----------------------------------------------------------------------------
// Writes the checkpoint data.
// This function writes the field data to the checkpoint file. It also calls
// writeCheckpointLevel() of all state objects owned by this level.
//
// NOTE: Do not write static data that comes from ParmParse.
// -----------------------------------------------------------------------------
void AMRNavierStokes::writeCheckpointLevel(HDF5Handle& a_handle) const
{
    if (s_verbosity >= 5) {
        pout () << "AMRNavierStokes::writeCheckpointLevel" << endl;
    }

    // Set group for this level.
    char level_str[20];
    sprintf (level_str, "%d", m_level);
    const std::string label = std::string("level_") + level_str;
    a_handle.setGroup(label);

    // Create the header for this group.
    // This will store the level's metadata.
    HDF5HeaderData header;

    // Collect all metadata that will be needed at restart.
    header.m_int     ["step_number"] = s_step_number;
    header.m_intvect ["ref_ratio"]   = m_ref_ratio;
    header.m_realvect["vec_dx"]      = m_levGeoPtr->getDx();
    header.m_real    ["dt"]          = m_dt;
    header.m_real    ["time"]        = m_time;
    header.m_real    ["cfl"]         = m_cfl;

    header.m_int["finest_level"] = m_finest_level;
    header.m_int["is_empty"]     = m_is_empty;

    header.m_box["prob_domain"]  = m_problem_domain.domainBox();

    D_TERM(
        header.m_int["is_periodic_0"] = (m_problem_domain.isPeriodic(0)? 1: 0);,
        header.m_int["is_periodic_1"] = (m_problem_domain.isPeriodic(1)? 1: 0);,
        header.m_int["is_periodic_2"] = (m_problem_domain.isPeriodic(2)? 1: 0);
    )

    // Write the metadata to file and pout.*
    header.writeToFile(a_handle);
    if (s_verbosity >= 3) {
        pout() << header << endl;
    }

    // If this level has valid data, we need to write it to HDF5.
    if (!isEmpty()) {
        // First, write the grids
        write (a_handle, m_vel_new_ptr->boxLayout());

        // Then, the velocity and lambda
        write (a_handle, *m_vel_new_ptr, "new_velocity");
        write (a_handle, *m_vel_old_ptr, "old_velocity");

        write (a_handle, *m_lambda_new_ptr, "new_lambda");
        write (a_handle, *m_lambda_old_ptr, "old_lambda");

        // Then, all of the scalars
        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            ostringstream old_scal_str;
            old_scal_str << "new_scalar_component_" << comp;
            write (a_handle, *m_scal_new[comp], old_scal_str.str().c_str());

            ostringstream new_scal_str;
            new_scal_str << "old_scalar_component_" << comp;
            write (a_handle, *m_scal_old[comp], new_scal_str.str().c_str());
        }

        // Pressure and VD correction stuff
        write (a_handle, m_ccPressure, "ccPressure");
        write (a_handle, m_eLambda, "eLambda");
    }
}


// -----------------------------------------------------------------------------
// Reads the checkpoint metadata.
// This function just performs a sanity check on the checkpoint file.
// -----------------------------------------------------------------------------
void AMRNavierStokes::readCheckpointHeader(HDF5Handle& a_handle)
{
    if (s_verbosity >= 5) {
        pout () << "AMRNavierStokes::readCheckpointHeader" << endl;
    }

    // Get the checkpoint's metadata
    HDF5HeaderData header;
    header.readFromFile(a_handle);

    // Write the metadata to pout.*
    if (s_verbosity >= 3) {
        pout() << "hdf5 header data:" << endl;
        pout() << header << endl;
    }

    // Read and check number of scalars
    if (header.m_int.find("num_components") == header.m_int.end()) {
        MayDay::Error("AMRNavierStokes::readCheckpointHeader: checkfile does not have num_components");
    }
    int num_comps = header.m_int["num_components"];
    if (num_comps != s_num_scal_comps) {
        MayDay::Error("AMRNavierStokes::readCheckpointHeader: num_components in input file does not match solver");
    }

    // Read and check scalar names
    std::string state_name;
    char comp_str[60];
    for (int comp = 0; comp < s_num_scal_comps; ++comp) {
        sprintf (comp_str, "component_%d", comp);
        if (header.m_string.find(comp_str) == header.m_string.end()) {
            MayDay::Error("AMRNavierStokes::readCheckpointHeader: checkfile does not have enough component names");
        }
        state_name = header.m_string[comp_str];
        if (state_name != s_scal_names[comp]) {
            MayDay::Error("AMRNavierStokes::readCheckpointHeader: state name in checkpoint does not match solver");
        }
    }

    // Read and check lambda name
    sprintf(comp_str, "lambda_component");
    if (header.m_string.find(comp_str) == header.m_string.end()) {
        MayDay::Error("AMRNavierStokes::readCheckpointHeader: checkfile does not have lambda name");
    }
    state_name = header.m_string[comp_str];
    if (state_name != "lambda") {
        MayDay::Error("AMRNavierStokes::readCheckpointHeader: lambda name in checkfile does not match that in solver");
    }

    // Read and check number and names of velocity components
    for (int comp = 0; comp < SpaceDim; ++comp) {
        sprintf(comp_str, "vel_component_%d", comp);
        if (header.m_string.find(comp_str) == header.m_string.end()) {
            MayDay::Error("AMRNavierStokes::readCheckpointHeader: checkfile does not have enough velocity names");
        }
        state_name = header.m_string[comp_str];
        if (state_name != s_vel_names[comp]) {
            MayDay::Error("AMRNavierStokes::readCheckpointHeader: vel name in checkfile does not match solver");
        }
    }

    // Check pressure
    sprintf(comp_str, "ccPressure_component");
    if (header.m_string.find(comp_str) == header.m_string.end()) {
        MayDay::Error("AMRNavierStokes::readCheckpointHeader: checkfile does not have CC pressure name");
    }
    state_name = header.m_string[comp_str];
    if (state_name != "ccPressure") {
        MayDay::Error("AMRNavierStokes::readCheckpointHeader: ccPressure name in checkfile does not match that in solver");
    }

    // Check eLambda
    sprintf(comp_str, "eLambda_component");
    if (header.m_string.find(comp_str) == header.m_string.end()) {
        MayDay::Error("AMRNavierStokes::readCheckpointHeader: checkfile does not have eLambda name");
    }
    state_name = header.m_string[comp_str];
    if (state_name != "eLambda") {
        MayDay::Error("AMRNavierStokes::readCheckpointHeader: eLambda name in checkfile does not match that in solver");
    }
}


// -----------------------------------------------------------------------------
// Reads the checkpoint data.
// -----------------------------------------------------------------------------
void AMRNavierStokes::readCheckpointLevel (HDF5Handle& a_handle)
{
    if (s_verbosity >= 5) {
        pout () << "AMRNavierStokes::readCheckpointLevel " << m_level << endl;
    }

    // Open this level's group.
    char level_str[20];
    sprintf(level_str, "%d", m_level);
    const std::string label = std::string("level_") + level_str;
    a_handle.setGroup(label);

    // Read this group's metadata
    HDF5HeaderData header;
    header.readFromFile(a_handle);

    // Write the metadata to pout.*
    if (s_verbosity >= 3) {
        pout() << "hdf5 header data:" << endl;
        pout() << header << endl;
    }

    // Read current step number
    if (header.m_int.find("step_number") == header.m_int.end()) {
        MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain step_number");
    }
    int tmp_step_number = header.m_int["step_number"];
    if (m_level == 0) {
        s_step_number = tmp_step_number;
    } else {
        if (s_step_number != tmp_step_number) {
            pout() << "Checkpoint step number = " << tmp_step_number
                   << " while s_step_number = " << s_step_number
                   << endl;
            MayDay::Error("Checkpoint is at wrong step number");
        }
    }

    // Read refinement ratio
    if (header.m_intvect.find("ref_ratio") == header.m_intvect.end()) {
        MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain ref_ratio");
    }
    m_ref_ratio = header.m_intvect["ref_ratio"];

    // Read dt
    if (header.m_real.find("dt") == header.m_real.end()) {
        MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain dt");
    }
    m_dt = header.m_real["dt"];

    // Read time
    if (header.m_real.find("time") == header.m_real.end()) {
        MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain time");
    }
    m_time = header.m_real["time"];

    // Read cfl
    if (header.m_real.find("cfl") == header.m_real.end()) {
        MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain cfl");
    }
    m_cfl = header.m_real["cfl"];

    // Read finest_level
    if (header.m_int.find("finest_level") == header.m_int.end()) {
        MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain finest_level");
    }
    m_finest_level = header.m_int["finest_level"];

    // Read is_empty
    if (header.m_int.find("is_empty") == header.m_int.end()) {
        MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain is_empty");
    }
    m_is_empty = header.m_int["is_empty"];

    // Read domain box
    if (header.m_box.find("prob_domain") == header.m_box.end()) {
        MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain prob_domain");
    }
    Box domainBox = header.m_box["prob_domain"];

    // Read in periodicity info
    bool is_periodic[SpaceDim];
    D_TERM(
        if (!(header.m_int.find("is_periodic_0") == header.m_int.end())) {
            is_periodic[0] = (header.m_int["is_periodic_0"] == 1);
        } else {
            MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain is_periodic_0");
        },
        if (!(header.m_int.find("is_periodic_1") == header.m_int.end())) {
            is_periodic[1] = (header.m_int["is_periodic_1"] == 1);
        } else {
            MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain is_periodic_1");
        },
        if (!(header.m_int.find("is_periodic_2") == header.m_int.end())) {
            is_periodic[2] = (header.m_int["is_periodic_2"] == 1);
        } else {
            MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain is_periodic_2");
        }
    );

    // Create level domain
    m_problem_domain = ProblemDomain(domainBox, is_periodic);

    // Read grids
    Vector<Box> boxArrayFromFile;
    const int grid_status = read(a_handle, boxArrayFromFile);
    if (grid_status != 0) {
        MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain a Vector<Box>");
    }

    // Create level grids
    const DisjointBoxLayout grids = loadBalance(boxArrayFromFile);

    // I am not using boxArrayFromFile because loadBalance may have changed things.
    m_level_grids = grids.boxArray();

    // Write level grids to pout.*
    if (s_verbosity >= 4) {
        pout() << "read level domain: " << endl;
        LayoutIterator lit = grids.layoutIterator();
        for (lit.begin(); lit.ok(); ++lit) {
            const Box& b = grids[lit];
            pout() << lit().intCode() << ": " << b << endl;
        }
        pout() << endl;
    }

    // If this level has valid data, read it now.
    if (!isEmpty()) {
        // Allocate and define fields.
        // This includes allocating the scalar Vector<LevelFluxRegister*>,
        // but do not define it. That will be done in this->levelSetup().

        const IntVect ghostVect = IntVect::Unit;

        m_vel_new_ptr = new LevelData<FArrayBox>(grids, CH_SPACEDIM, ghostVect);
        m_vel_old_ptr = new LevelData<FArrayBox>(grids, CH_SPACEDIM, ghostVect);

        m_lambda_new_ptr = new LevelData<FArrayBox>(grids, 1, ghostVect);
        m_lambda_old_ptr = new LevelData<FArrayBox>(grids, 1, ghostVect);

        m_scal_new.resize(s_num_scal_comps, NULL);
        m_scal_old.resize(s_num_scal_comps, NULL);
        m_scal_fluxreg_ptrs.resize(s_num_scal_comps,NULL);

        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            m_scal_new[comp] = new LevelData<FArrayBox>(grids, 1, ghostVect);
            m_scal_old[comp] = new LevelData<FArrayBox>(grids, 1, ghostVect);
            m_scal_fluxreg_ptrs[comp] = new MappedLevelFluxRegister;
        }

        m_macPressure.define(grids, 1, IntVect::Unit);
        m_ccPressure.define(grids, 1, IntVect::Unit);
        m_syncPressure.define(grids, 1, IntVect::Unit);

        // Fill all fields with bogus data
        this->setAllBogus();

        // Read velocity data
        {
            LevelData<FArrayBox>& new_vel = *m_vel_new_ptr;
            const int velData_status = read<FArrayBox>(a_handle,
                                                       new_vel,
                                                       "new_velocity",
                                                       grids);
            if (velData_status != 0) {
                MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain new_velocity data");
            }
        }
        {
            LevelData<FArrayBox>& old_vel = *m_vel_old_ptr;
            const int velData_status = read<FArrayBox>(a_handle,
                                                       old_vel,
                                                       "old_velocity",
                                                       grids);
            if (velData_status != 0) {
                MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain old_velocity data");
            }
        }

        // Read lambda data
        {
            LevelData<FArrayBox>& new_lambda = *m_lambda_new_ptr;
            const int lambdaData_status = read<FArrayBox>(a_handle,
                                                          new_lambda,
                                                          "new_lambda",
                                                          grids);
            if (lambdaData_status != 0) {
                MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain new_lambda data");
            }
        }
        {
            LevelData<FArrayBox>& old_lambda = *m_lambda_old_ptr;
            const int lambdaData_status = read<FArrayBox>(a_handle,
                                                          old_lambda,
                                                          "old_lambda",
                                                          grids);
            if (lambdaData_status != 0) {
                MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain old_lambda data");
            }
        }

        // Read scalar data
        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            {
                char scal_str[20];
                sprintf(scal_str, "new_scalar_component_%d", comp);

                LevelData<FArrayBox>& new_scal = *m_scal_new[comp];
                const int scalData_status = read<FArrayBox>(a_handle,
                                                            new_scal,
                                                            scal_str,
                                                            grids);
                if (scalData_status != 0) {
                    MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain new_scalar data");
                }
            }
            {
                char scal_str[20];
                sprintf(scal_str, "old_scalar_component_%d", comp);

                LevelData<FArrayBox>& old_scal = *m_scal_old[comp];
                const int scalData_status = read<FArrayBox>(a_handle,
                                                            old_scal,
                                                            scal_str,
                                                            grids);
                if (scalData_status != 0) {
                    MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain old_scalar data");
                }
            }
        }

    } else {
        // This level is empty.
        m_vel_new_ptr = NULL;
        m_vel_old_ptr = NULL;
        m_lambda_new_ptr = NULL;
        m_lambda_old_ptr = NULL;

        // The pressure fields.
        m_macPressure.clear();
        m_ccPressure.clear();
        m_syncPressure.clear();
    }

    // Now that this level has a domain, grids, and filled data holders,
    // we can set up all other objects on this level as if we just finished
    // regridding.
    this->levelSetup(grids);

    // We read the pressure data *after* calling levelSetup so that it
    // doesn't get clobbered.
    if (!isEmpty()) {
        const int ccPressureData_status = read<FArrayBox>(a_handle,
                                                          m_ccPressure,
                                                          "ccPressure",
                                                          grids);
        if (ccPressureData_status != 0) {
            MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain ccPressure data");
        }
        m_ccPressureState = CCPressureState::VALID;

        // Set CFBCs. The grad functions will do the rest.
        if (m_level > 0) {
            const AMRNavierStokes* crsePtr = crseNSPtr();
            if (crsePtr->m_ccPressureState != CCPressureState::VALID) {
                MayDay::Error("AMRNavierStokes::readCheckpointLevel: cannot set CF BCs");
            }

            // CFBCs
            const LevelData<FArrayBox>& crsePressure = crsePtr->m_ccPressure;
            MappedQuadCFInterp interpObj(grids,
                                         &(crsePressure.getBoxes()),
                                         m_levGeoPtr->getDx(),
                                         m_levGeoPtr->getCrseRefRatio(),
                                         1, // ncomp
                                         m_problem_domain);
            interpObj.coarseFineInterp(m_ccPressure, crsePressure);

            // Finish the job by extrapolating the edges and vertices.
            CFRegion cfregion(grids, m_problem_domain);
            ExtrapolateCFEV(m_ccPressure, cfregion, 2);
        }

        // eLambda
        if (s_etaLambda > 0.0) {
            // Read eLambda from file.
            const int eLambdaData_status = read<FArrayBox>(a_handle,
                                                           m_eLambda,
                                                           "eLambda",
                                                           grids);
            if (ccPressureData_status != 0) {
                MayDay::Error("AMRNavierStokes::readCheckpointLevel: file does not contain eLambda data");
            }

            // Set physical BCs
            BCMethodHolder bcHolder = m_physBCPtr->gradELambdaFuncBC();

            const DisjointBoxLayout& levelGrids = m_eLambda.getBoxes();
            const ProblemDomain& levelDomain = levelGrids.physDomain();
            const RealVect& levelDx = m_levGeoPtr->getDx();

            DataIterator ditLev = m_eLambda.dataIterator(); // TODO: May not be needed.
            for (ditLev.reset(); ditLev.ok(); ++ditLev) {
                bcHolder.setGhosts(m_eLambda[ditLev],   // stateFAB
                                   NULL,                // &extrapFAB
                                   levelGrids[ditLev],  // valid box
                                   levelDomain,         // ProblemDomain
                                   levelDx,             // dx
                                   ditLev(),            // DataIndex
                                   NULL,                // &JgupFAB
                                   false,               // isHomogeneous
                                   m_time);             // time
            }

            // CF BCs
            const AMRNavierStokes* thisCrseNSPtr = crseNSPtr();
            if (thisCrseNSPtr != NULL) {
                const LevelData<FArrayBox>& crseELambda = thisCrseNSPtr->m_eLambda;

                MappedQuadCFInterp interpObj(levelGrids,
                                             &(crseELambda.getBoxes()),
                                             levelDx,
                                             m_levGeoPtr->getCrseRefRatio(),
                                             1, // ncomp
                                             levelDomain);
                interpObj.coarseFineInterp(m_eLambda, crseELambda);

                // Don't miss the corner cells excluded by the interpolator!
                CFRegion cfregion(levelGrids, levelDomain);
                ExtrapolateCFEV(m_eLambda, cfregion, 2, IntVect::Unit);
            }

            // Do exchanges.
            Copier exc(levelGrids, levelGrids, levelDomain, IntVect::Unit, true);
            m_eLambda.exchange(exc);

            CornerCopier excc(levelGrids, levelGrids, levelDomain, IntVect::Unit, true);
            m_eLambda.exchange(excc);

        } else {
            // Just set to zero.
            m_eLambda.define(grids, 1, IntVect::Unit);
            setValLevel(m_eLambda, 0.0);
        }
        m_eLambdaState = ELambdaState::VALID;

        // Grad[eLambda]
        // This must be computed from the finest extant level, down.
        // I handle this a bit redundantly, but we only need to do this once.
        m_gradELambda.define(grids, 1, IntVect::Unit);
        if (s_applyFreestreamCorrection) {
            AMRNavierStokes* thisNSPtr = this;
            AMRNavierStokes* thisFineNSPtr = NULL;

            while (thisNSPtr != NULL) {
                const AMRNavierStokes* thisCrseNSPtr = thisNSPtr->crseNSPtr();

                // Sanity checks.
                if (thisCrseNSPtr != NULL) {
                    if (thisNSPtr->m_eLambdaState != ELambdaState::VALID) {
                        MayDay::Error("AMRNavierStokes::readCheckpointLevel: eLambda is not valid");
                    }
                    if (thisCrseNSPtr->m_eLambdaState != ELambdaState::VALID) {
                        MayDay::Error("AMRNavierStokes::readCheckpointLevel: Coarse level eLambda is not valid");
                    }
                    if (thisCrseNSPtr->m_gradELambdaState != GradELambdaState::VALID) {
                        MayDay::Error("AMRNavierStokes::readCheckpointLevel: Coarse level gradELambda is not valid");
                    }
                }

                // Compute gradient.
                // Recall that composite MAC gradient is the same as the level
                // gradient, since finer level is not considered to be part
                // of this level (take care of covered regions with avgDown)
                BCMethodHolder fluxBC = thisNSPtr->m_physBCPtr->gradELambdaFuncBC();
                const LevelData<FArrayBox>* crseELambdaPtr = (thisCrseNSPtr != NULL)?
                                                             &(thisCrseNSPtr->m_eLambda):
                                                             NULL;
                Gradient::levelGradientMAC(thisNSPtr->m_gradELambda,
                                           thisNSPtr->m_eLambda,
                                           crseELambdaPtr,
                                           *(thisNSPtr->m_levGeoPtr),
                                           m_time,
                                           &fluxBC);

                // Average gradient down from finer level.
                if (thisFineNSPtr != NULL) {
                    const LevelGeometry* fineLevGeoPtr = thisFineNSPtr->m_levGeoPtr;

                    const DisjointBoxLayout& fineGrids = fineLevGeoPtr->getBoxes();
                    const IntVect& nRefFine = fineLevGeoPtr->getCrseRefRatio();
                    const int nComp = 1;
                    const LevelData<FluxBox>& fineEdgeGrad = thisFineNSPtr->m_gradELambda;

                    MappedCoarseAverageFace avgDownObj(fineGrids, nComp, nRefFine);
                    avgDownObj.averageToCoarse(thisNSPtr->m_gradELambda, fineEdgeGrad);
                }

                // Exchanges.
                {
                    const DisjointBoxLayout& levelGrids = thisNSPtr->m_gradELambda.getBoxes();
                    const ProblemDomain& levelDomain = levelGrids.physDomain();

                    Copier exc(levelGrids, levelGrids, levelDomain, IntVect::Unit, true);
                    thisNSPtr->m_gradELambda.exchange(exc);

                    CornerCopier excc(levelGrids, levelGrids, levelDomain, IntVect::Unit, true);
                    thisNSPtr->m_gradELambda.exchange(excc);
                }

                // Move on to the next coarser level.
                thisFineNSPtr = thisNSPtr;
                thisNSPtr = thisNSPtr->crseNSPtr();
            } // end while loop down levels

        } else {
            // Just set to zero.
            setValLevel(m_gradELambda, 0.0);

        } // end if/if not applying VD correction.

        m_gradELambdaState = GradELambdaState::VALID;

    } else {
        m_macPressure.clear();
        m_ccPressure.clear();
        m_syncPressure.clear();
    }
}


// -----------------------------------------------------------------------------
// writePlotHeader
// -----------------------------------------------------------------------------
void AMRNavierStokes::writePlotHeader (HDF5Handle& a_handle) const
{
    if (s_verbosity >= 5) {
        pout() << "AMRNavierStokes::writePlotHeader " << m_level << endl;
    }

    HDF5HeaderData header;

    int numcomp = numPlotComps();
    header.m_int ["num_components"] = numcomp;

    char comp_str[30];
    int comp = 0;

    sprintf(comp_str, "component_%d", comp);
    header.m_string[comp_str] = "x_Vel";
    comp++;

    if (SpaceDim > 1) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "y_Vel";
        comp++;

        if (SpaceDim > 2) {
            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "z_Vel";
            comp++;
        }

        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "mag_vel";
        comp++;
    }

    if (s_write_divergence) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "divergence";
        comp++;
    }

    if (s_write_lambda) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "lambda-1";
        comp++;
    }

    if (s_write_grad_eLambda) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "x_Grad_eLambda";
        comp++;

        if (SpaceDim >1) {
            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "y_Grad_eLambda";
            comp++;

            if (SpaceDim >2) {
                sprintf(comp_str, "component_%d", comp);
                header.m_string[comp_str] = "z_Grad_eLambda";
                comp++;
            }
        }
    } // end if writing grad_eLambda

    if (s_write_pressure) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "pressure";
        comp++;
    }

    if (s_write_vorticity) {
        if (SpaceDim == 2) {
            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "vorticity";
            comp++;
        } else if (SpaceDim == 3) {
            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "x_vort";
            comp++;

            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "y_vort";
            comp++;

            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "z_vort";
            comp++;

            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "mag_vort";
            comp++;
        } // end if 3d
    } // end if writing vorticity

    if (s_write_streamfunction) {
        if  (SpaceDim == 2) {
            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "streamfunction";
            comp++;
        } else if (SpaceDim == 3) {
            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "x_streamfunction";
            comp++;

            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "y_streamfunction";
            comp++;

            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "z_streamfunction";
            comp++;
        }
    }

    if (s_write_scalars) {
        for (int scal = 0; scal < s_num_scal_comps; ++scal) {
            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = s_scal_names[scal];
            comp++;
        }
    }

    if (s_write_scalarsMinusBackground) {
        for (int scal = 0; scal < s_num_scal_comps; ++scal) {
            sprintf(comp_str, "component_%d", comp);
            std::ostringstream thisName;
            thisName << s_scal_names[scal] << "_pert";
            header.m_string[comp_str] = thisName.str().c_str();
            comp++;
        }
    }

    // write procID's
    if (s_write_proc_ids) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "procIDs";
        comp++;
    }

    // write level ID's
    if (s_write_level_ids) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "levelIDs";
        comp++;
    }

    // write displacement
    if (s_write_displacement) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "x_Displacement";
        comp++;

        if (SpaceDim >1) {
            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "y_Displacement";
            comp++;

            if (SpaceDim >2) {
                sprintf(comp_str, "component_%d", comp);
                header.m_string[comp_str] = "z_Displacement";
                comp++;
            }
        }
    }

    // Write geometry
    if (s_write_geometry) {
        // physCoor
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "Geo_x_physCoor";
        comp++;

        if (SpaceDim >1) {
            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "Geo_y_physCoor";
            comp++;

            if (SpaceDim >2) {
                sprintf(comp_str, "component_%d", comp);
                header.m_string[comp_str] = "Geo_z_physCoor";
                comp++;
            }
        }

        // dxdXi
        for (int adir = 0; adir < SpaceDim; ++adir) {
            for (int bdir = 0; bdir < SpaceDim; ++bdir) {
                sprintf(comp_str, "component_%d", comp);
                if (adir == 0) {
                    if (bdir == 0) header.m_string[comp_str] = "Geo_xx_dxdXi";
                    if (bdir == 1) header.m_string[comp_str] = "Geo_xy_dxdXi";
                    if (bdir == 2) header.m_string[comp_str] = "Geo_xz_dxdXi";
                } else if (adir == 1) {
                    if (bdir == 0) header.m_string[comp_str] = "Geo_yx_dxdXi";
                    if (bdir == 1) header.m_string[comp_str] = "Geo_yy_dxdXi";
                    if (bdir == 2) header.m_string[comp_str] = "Geo_yz_dxdXi";
                } else {
                    if (bdir == 0) header.m_string[comp_str] = "Geo_zx_dxdXi";
                    if (bdir == 1) header.m_string[comp_str] = "Geo_zy_dxdXi";
                    if (bdir == 2) header.m_string[comp_str] = "Geo_zz_dxdXi";
                }
                comp++;
            }
        }

        // dXidx
        for (int adir = 0; adir < SpaceDim; ++adir) {
            for (int bdir = 0; bdir < SpaceDim; ++bdir) {
                sprintf(comp_str, "component_%d", comp);
                if (adir == 0) {
                    if (bdir == 0) header.m_string[comp_str] = "Geo_xx_dXidx";
                    if (bdir == 1) header.m_string[comp_str] = "Geo_xy_dXidx";
                    if (bdir == 2) header.m_string[comp_str] = "Geo_xz_dXidx";
                } else if (adir == 1) {
                    if (bdir == 0) header.m_string[comp_str] = "Geo_yx_dXidx";
                    if (bdir == 1) header.m_string[comp_str] = "Geo_yy_dXidx";
                    if (bdir == 2) header.m_string[comp_str] = "Geo_yz_dXidx";
                } else {
                    if (bdir == 0) header.m_string[comp_str] = "Geo_zx_dXidx";
                    if (bdir == 1) header.m_string[comp_str] = "Geo_zy_dXidx";
                    if (bdir == 2) header.m_string[comp_str] = "Geo_zz_dXidx";
                }
                comp++;
            }
        }

        // J
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "Geo_J";
        comp++;

        // Jinv
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "Geo_Jinv";
        comp++;

        // gdn
        for (int adir = 0; adir < SpaceDim; ++adir) {
            for (int bdir = 0; bdir <= adir; ++bdir) {
                sprintf(comp_str, "component_%d", comp);
                if (adir == 0) {
                    if (bdir == 0) header.m_string[comp_str] = "Geo_xx_gdn";
                    if (bdir == 1) header.m_string[comp_str] = "Geo_xy_gdn";
                    if (bdir == 2) header.m_string[comp_str] = "Geo_xz_gdn";
                } else if (adir == 1) {
                    if (bdir == 0) header.m_string[comp_str] = "Geo_yx_gdn";
                    if (bdir == 1) header.m_string[comp_str] = "Geo_yy_gdn";
                    if (bdir == 2) header.m_string[comp_str] = "Geo_yz_gdn";
                } else {
                    if (bdir == 0) header.m_string[comp_str] = "Geo_zx_gdn";
                    if (bdir == 1) header.m_string[comp_str] = "Geo_zy_gdn";
                    if (bdir == 2) header.m_string[comp_str] = "Geo_zz_gdn";
                }
                comp++;
            }
        }

        // gup
        for (int adir = 0; adir < SpaceDim; ++adir) {
            for (int bdir = 0; bdir <= adir; ++bdir) {
                sprintf(comp_str, "component_%d", comp);
                if (adir == 0) {
                    if (bdir == 0) header.m_string[comp_str] = "Geo_xx_gup";
                    if (bdir == 1) header.m_string[comp_str] = "Geo_xy_gup";
                    if (bdir == 2) header.m_string[comp_str] = "Geo_xz_gup";
                } else if (adir == 1) {
                    if (bdir == 0) header.m_string[comp_str] = "Geo_yx_gup";
                    if (bdir == 1) header.m_string[comp_str] = "Geo_yy_gup";
                    if (bdir == 2) header.m_string[comp_str] = "Geo_yz_gup";
                } else {
                    if (bdir == 0) header.m_string[comp_str] = "Geo_zx_gup";
                    if (bdir == 1) header.m_string[comp_str] = "Geo_zy_gup";
                    if (bdir == 2) header.m_string[comp_str] = "Geo_zz_gup";
                }
                comp++;
            }
        }
    }

    // Stress tensor
#   ifdef USE_STRESSMETRIC
    if (m_stressMetric.isDefined()) {
        for (int adir = 0; adir < SpaceDim; ++adir) {
            for (int bdir = 0; bdir <= adir; ++bdir) {
                sprintf(comp_str, "component_%d", comp);
                if (adir == 0) {
                    if (bdir == 0) header.m_string[comp_str] = "xx_StressTensor";
                    if (bdir == 1) header.m_string[comp_str] = "xy_StressTensor";
                    if (bdir == 2) header.m_string[comp_str] = "xz_StressTensor";
                } else if (adir == 1) {
                    if (bdir == 0) header.m_string[comp_str] = "yx_StressTensor";
                    if (bdir == 1) header.m_string[comp_str] = "yy_StressTensor";
                    if (bdir == 2) header.m_string[comp_str] = "yz_StressTensor";
                } else {
                    if (bdir == 0) header.m_string[comp_str] = "zx_StressTensor";
                    if (bdir == 1) header.m_string[comp_str] = "zy_StressTensor";
                    if (bdir == 2) header.m_string[comp_str] = "zz_StressTensor";
                }
                comp++;
            }
        }
    }
#   endif //USE_STRESSMETRIC

    // Write Taylor-Green vortex solution
    if (dynamic_cast<TaylorGreenBCUtil*>(m_physBCPtr)) {
        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "x_Vel_Sol";
        comp++;

        if (SpaceDim >1) {
            sprintf(comp_str, "component_%d", comp);
            header.m_string[comp_str] = "y_Vel_Sol";
            comp++;

            if (SpaceDim >2) {
                sprintf(comp_str, "component_%d", comp);
                header.m_string[comp_str] = "z_Vel_Sol";
                comp++;
            }
        }

        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "Pressure_Sol";
        comp++;

        sprintf(comp_str, "component_%d", comp);
        header.m_string[comp_str] = "FofT";
        comp++;
    }

    header.writeToFile(a_handle);

    if (s_verbosity >= 5) {
        pout () << header << endl;
    }
}


// -----------------------------------------------------------------------------
// writePlotLevel
// -----------------------------------------------------------------------------
void AMRNavierStokes::writePlotLevel(HDF5Handle& a_handle) const
{
    if (s_verbosity >= 5) {
        pout() << "AMRNavierStokes::writePlotLevel " << m_level << endl;
    }

    char level_str[20];
    sprintf (level_str, "%d", m_level);
    const std::string label = std::string("level_") + level_str;

    a_handle.setGroup(label);

    HDF5HeaderData header;

    header.m_intvect ["ref_ratio"]   = m_ref_ratio;
    header.m_realvect["vec_dx"]      = m_levGeoPtr->getDx();
    header.m_real    ["dt"]          = m_dt;
    header.m_real    ["time"]        = m_time;
    header.m_box     ["prob_domain"] = m_problem_domain.domainBox();
    header.writeToFile(a_handle);

    if (s_verbosity >= 3) {
        pout () << header << endl;
    }

    const DisjointBoxLayout& levelGrids = m_vel_new_ptr->getBoxes();
    int numcomp = numPlotComps();

    // We will need one ghost in order for the displacement field to
    // do its job properly.
    LevelData<FArrayBox> plotData(levelGrids, numcomp, IntVect::Unit);
    getPlotData(plotData);

    // This will fill all ghosts of plotData, including edges and vertices so
    // that VisIt will be able to generate contour plots cleanly.
#   define EXTRAPOLATE_GHOSTS

#   ifdef EXTRAPOLATE_GHOSTS
    {
        extrapAllGhosts(plotData, 2);

        // New code
        ProblemDomain nonPeriodicDomain(m_problem_domain.domainBox());
        Copier nonPeriodicExCopier(levelGrids, levelGrids, nonPeriodicDomain, plotData.ghostVect(), true);
        plotData.exchange(nonPeriodicExCopier);
    }
#   endif

    write (a_handle, levelGrids);
    write (a_handle, plotData, "data", plotData.ghostVect());
}

#endif //CH_USE_HDF5


// -----------------------------------------------------------------------------
// numPlotComps
// -----------------------------------------------------------------------------
int AMRNavierStokes::numPlotComps () const
{

    // velocity
    int numcomp = SpaceDim;
    // also include mag(vel)
    ++numcomp;

    // divergence
    if (s_write_divergence) {
        ++numcomp;
    }

    // lambda
    if (s_write_lambda) {
        ++numcomp;
    }

    // add in grad_eLambda
    if (s_write_grad_eLambda) {
        numcomp += SpaceDim;
    }

    // pressure
    if (s_write_pressure) {
        ++numcomp;
    }

    // vorticity
    if (s_write_vorticity) {
        if (SpaceDim == 2) {
            ++numcomp;
        } else if (SpaceDim == 3) {
            numcomp += SpaceDim;
            // also include mag(vort)
            ++numcomp;
        }
    }

    // streamfunction
    if (s_write_streamfunction) {
        if (SpaceDim == 2) {
            ++numcomp;
        } else if (SpaceDim == 3) {
            numcomp += SpaceDim;
        }
    }

    // scalars
    if (s_write_scalars) {
        numcomp += s_num_scal_comps;
    }

    if (s_write_scalarsMinusBackground) {
        numcomp += s_num_scal_comps;
    }

    // procID's
    if (s_write_proc_ids) {
        ++numcomp;
    }

    // level ID's
    if (s_write_level_ids) {
        ++numcomp;
    }

    // displacement
    if (s_write_displacement) {
        numcomp += SpaceDim;
    }

    // geometry
    if (s_write_geometry) {
        numcomp += SpaceDim;                    // physCoor
        numcomp += SpaceDim*SpaceDim;           // dxdXi
        numcomp += SpaceDim*SpaceDim;           // dXidx
        numcomp += 1;                           // J
        numcomp += 1;                           // Jinv
        numcomp += (SpaceDim*(SpaceDim+1))/2;   // gdn
        numcomp += (SpaceDim*(SpaceDim+1))/2;   // gup
    }

#   ifdef USE_STRESSMETRIC
    if (m_stressMetric.isDefined()) {
        numcomp += (SpaceDim*(SpaceDim+1))/2;
    }
#   endif

    // Taylor-Green vortex solution
    if (dynamic_cast<TaylorGreenBCUtil*>(m_physBCPtr)) {
        numcomp += SpaceDim;    // velocity
        ++numcomp;              // pressure
        ++numcomp;              // FofT
    }

    return numcomp;
}


// -----------------------------------------------------------------------------
// getPlotData
// -----------------------------------------------------------------------------
void AMRNavierStokes::getPlotData (LevelData<FArrayBox>& a_plot_data) const
{

    // Sanity checks
    CH_assert(a_plot_data.nComp() == numPlotComps());

    // Gather data
    const DisjointBoxLayout& grids = a_plot_data.getBoxes();
    DataIterator dit = grids.dataIterator();
    const ProblemDomain& domain = grids.physDomain();

    const Interval srcComp(0, 0);
    const Interval velComps(0, SpaceDim-1);
    const IntVect ghostVect = IntVect::Unit;

    LevelData<FArrayBox>* fineVelPtr = NULL;
    LevelData<FArrayBox>* crseVelPtr = NULL;
    IntVect nRefCrse(D_DECL(-1,-1,-1));

    // Grab coarse/fine level data
    if (m_level > 0) {
        crseVelPtr = crseNSPtr()->m_vel_new_ptr;
        nRefCrse = m_coarser_level_ptr->refRatio();
    }

    if (!finestLevel()) {
        fineVelPtr = fineNSPtr()->m_vel_new_ptr;
    }

    // Initialize plot counter
    int plot_data_counter = 0;

    // Create a copier
    Copier copier(grids, grids, domain, ghostVect, false);

    // Velocity
    const Interval destVelComps(plot_data_counter, plot_data_counter+SpaceDim-1);
    {
        // Copy
        m_vel_new_ptr->copyTo(velComps, a_plot_data, destVelComps, copier);

        // Convert to Cartesian basis
        LevelData<FArrayBox> velAlias;
        aliasLevelData(velAlias, &a_plot_data, destVelComps);
        m_levGeoPtr->sendToCartesianBasis(velAlias, true);

        // Increment plot count
        plot_data_counter += SpaceDim;
    }

    // Norm of velocity
    {
        // Calculate norm
        LevelData<FArrayBox> magVel(grids, 1);
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& magVelFAB = magVel[dit];
            const FArrayBox& velFAB = (*m_vel_new_ptr)[dit];
            const Box& valid = grids[dit];
            const Real half = 0.5;

            m_levGeoPtr->contractVectors(magVelFAB, velFAB, velFAB, dit());
            FORT_POWFAB(CHF_FRA(magVelFAB),
                        CHF_BOX(valid),
                        CHF_CONST_REAL(half));
        }

        // Copy
        Interval magVelComps(plot_data_counter, plot_data_counter);
        magVel.copyTo(srcComp, a_plot_data, magVelComps);

        // Increment plot count
        ++plot_data_counter;
    }

    // Divergence
    if (s_write_divergence) {
        // Set up Ju on this level
        LevelData<FArrayBox> Ju(grids, SpaceDim, m_vel_new_ptr->ghostVect());
        m_vel_new_ptr->copyTo(velComps, Ju, velComps, copier);
        m_levGeoPtr->multByJ(Ju);
        Ju.exchange(velComps, m_oneGhostExCopier);

        // Set up coarser Ju
        LevelData<FArrayBox>* crseJuPtr = NULL;
        if (m_level > 0) {
            crseJuPtr = new LevelData<FArrayBox>(crseVelPtr->getBoxes(),
                                                 SpaceDim,
                                                 ghostVect);
            crseVelPtr->copyTo(velComps, *crseJuPtr, velComps);
            crseNSPtr()->m_levGeoPtr->multByJ(*crseJuPtr);
        }

        // Set up finer Ju
        LevelData<FArrayBox>* fineJuPtr = NULL;
        if (!finestLevel()) {
            fineJuPtr = new LevelData<FArrayBox>(fineVelPtr->getBoxes(),
                                                 SpaceDim,
                                                 fineVelPtr->ghostVect());
            fineVelPtr->copyTo(velComps, *fineJuPtr, velComps);
            fineNSPtr()->m_levGeoPtr->multByJ(*fineJuPtr);
        }

        // Compute divergence
        bool isViscousDummy = true;
        Tuple<BCMethodHolder, SpaceDim> uStarBCHolder = m_physBCPtr->uStarFuncBC(isViscousDummy);

        LevelData<FArrayBox> divU(grids, 1);
        Divergence::compDivergenceCC(divU,
                                     Ju,
                                     crseJuPtr,
                                     fineJuPtr,
                                     true,
                                     *m_levGeoPtr,
                                     m_time,
                                     &uStarBCHolder);

        // Free memory
        delete crseJuPtr;
        delete fineJuPtr;

        // Copy
        Interval divDestComp(plot_data_counter, plot_data_counter);
        divU.copyTo(srcComp, a_plot_data, divDestComp);

        // Increment plot counter
        plot_data_counter += 1;
    } // end divergence

    // Lambda - 1.0
    if (s_write_lambda) {
        // Copy
        Interval lambdaDestComp(plot_data_counter, plot_data_counter);
        m_lambda_new_ptr->copyTo(srcComp, a_plot_data, lambdaDestComp, copier);

        for (dit.reset(); dit.ok(); ++dit) {
            a_plot_data[dit].plus(-1.0, plot_data_counter, 1);
        }

        // Increment plot counter
        plot_data_counter += 1;
    }

    // VD correction
    if (s_write_grad_eLambda) {
        // Write the EdgeToCell average right to the plot holder
        for (dit.begin(); dit.ok(); ++dit) {
            FArrayBox& thisPlotData = a_plot_data[dit];
            const FluxBox& thisGradE = m_gradELambda[dit];

            for (int dir = 0; dir < SpaceDim; ++dir) {
                int edgeComp = 0;
                int plotComp = plot_data_counter + dir;

                EdgeToCell(thisGradE,
                           edgeComp,
                           thisPlotData,
                           plotComp,
                           dir);

                m_levGeoPtr->divByJ(thisPlotData, dit(), plotComp);
            }
        }

        // Increment plot counter
        plot_data_counter += SpaceDim;
    }

    // Pressure
    Interval piDestComp; // This will be needed for the TG vortex test output.
    if (s_write_pressure) {
        // Copy
        piDestComp.define(plot_data_counter, plot_data_counter);
        m_ccPressure.copyTo(srcComp, a_plot_data, piDestComp, copier);

        // Add eSync to Pi. (Not available during startup from a checkpoint.)
        if (m_syncPressureState == SyncPressureState::SYNC) {
            m_syncPressure.addTo(srcComp, a_plot_data, piDestComp, domain);
        }

        // Increment plot counter
        plot_data_counter += 1;
    }

    // Vorticity and norm(vorticity)
    if (s_write_vorticity) {
        const int numVortComps = D_TERM(0,+1,+2);

        // Compute vorticity
        LevelData<FArrayBox> vort(grids, numVortComps);
        this->computeVorticity(vort);

        // Copy vorticity
        Interval srcVortComps(0, numVortComps-1);
        Interval destVortComps(plot_data_counter, plot_data_counter+numVortComps-1);
        vort.copyTo(srcVortComps, a_plot_data, destVortComps);

        // Increment plot counter
        plot_data_counter += numVortComps;

        if (SpaceDim == 3) {
            // Compute norm(vorticity)
            LevelData<FArrayBox> normVort(grids, 1);
            for (dit.reset(); dit.ok(); ++dit) {
                m_levGeoPtr->contractVectors(normVort[dit],
                                             vort[dit],
                                             vort[dit],
                                             dit());

                BoxIterator bit(normVort[dit].box());
                for (bit.reset(); bit.ok(); ++bit) {
                    normVort[dit](bit(),0) = sqrt(normVort[dit](bit(),0));
                }
            }

            // Copy norm(vorticity)
            Interval magVortComps(plot_data_counter, plot_data_counter);
            normVort.copyTo(srcComp, a_plot_data, magVortComps);

            // Increment plot counter
            ++plot_data_counter;
        }
    }

    // Streamfunction
    if (s_write_streamfunction) {
        const int numStreamComps = D_TERM(0,+1,+2);

        // Compute streamfunction
        LevelData<FArrayBox> streamFunction(grids, numStreamComps, IntVect::Unit);
        this->computeStreamFunction(streamFunction);

        // Copy streamfunction
        Interval srcStreamComps(0, numStreamComps-1);
        Interval destStreamComps(plot_data_counter, plot_data_counter+numStreamComps-1);
        streamFunction.copyTo(srcStreamComps, a_plot_data, destStreamComps);

        // Increment plot counter
        plot_data_counter += numStreamComps;
    }

    // Scalars
    if (s_write_scalars) {
        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            const LevelData<FArrayBox>& src = newScal(comp);
            const int srcComps = src.nComp();
            const Interval srcInterval = src.interval();

            const Interval destInterval(plot_data_counter,
                                        plot_data_counter + srcComps - 1);

            // Copy
            src.copyTo(srcInterval, a_plot_data, destInterval, copier);

            // Add background scalar to a_plot_data.
            LevelData<FArrayBox> dest;
            aliasLevelData(dest, &a_plot_data, destInterval);
            m_physBCPtr->addBackgroundScalar(dest, comp, m_time, *m_levGeoPtr);

            // Increment plot counter
            plot_data_counter += srcComps;
        }
    }

    if (s_write_scalarsMinusBackground) {
        for (int comp = 0; comp < s_num_scal_comps; ++comp) {
            const LevelData<FArrayBox>& src = newScal(comp);
            const int srcComps = src.nComp();
            const Interval srcInterval = src.interval();

            const Interval destInterval(plot_data_counter,
                                        plot_data_counter + srcComps - 1);

            // Copy
            src.copyTo(srcInterval, a_plot_data, destInterval, copier);

            // Increment plot counter
            plot_data_counter += srcComps;
        }
    }

    // procID's
    if (s_write_proc_ids) {
        // Write proc id directly to plot holder
        const Real procVal = (Real) procID();
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& thisPlotData = a_plot_data[dit];
            thisPlotData.setVal(procVal, plot_data_counter);
        }

        // Increment plot counter
        ++plot_data_counter;
    }

    // level ID's
    if (s_write_level_ids) {
        // Write level id directly to plot holder
        const Real levVal = (Real) m_level;
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& thisPlotData = a_plot_data[dit];
            thisPlotData.setVal(levVal, plot_data_counter);
        }

        // Increment plot counter
        ++plot_data_counter;
    }

    // Displacement (used to visualize geometries in VisIt w/ displace operator)
    if (s_write_displacement) {
        // Write displacements directly to plot holder
        Interval interv(plot_data_counter, plot_data_counter + SpaceDim - 1);
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox dispFAB(interv, a_plot_data[dit]);
            m_levGeoPtr->fill_displacement(dispFAB);
        }

        // Increment plot counter
        plot_data_counter += SpaceDim;
    }

    // Geometry
    if (s_write_geometry) {
        { // physCoor
            Interval interv(plot_data_counter, plot_data_counter + SpaceDim - 1);
            for (dit.reset(); dit.ok(); ++dit) {
                FArrayBox physCoorFAB(interv, a_plot_data[dit]);
                m_levGeoPtr->fill_physCoor(physCoorFAB);
            }
            plot_data_counter += SpaceDim;
        }

        { // dxdXi
            Interval interv(plot_data_counter, plot_data_counter + SpaceDim*SpaceDim - 1);
            for (dit.reset(); dit.ok(); ++dit) {
                FArrayBox dxdXiFAB(interv, a_plot_data[dit]);
                m_levGeoPtr->fill_dxdXi(dxdXiFAB);
            }
            plot_data_counter += SpaceDim*SpaceDim;
        }

        { // dXidx
            Interval interv(plot_data_counter, plot_data_counter + SpaceDim*SpaceDim - 1);
            for (dit.reset(); dit.ok(); ++dit) {
                FArrayBox dXidxFAB(interv, a_plot_data[dit]);
                m_levGeoPtr->fill_dXidx(dXidxFAB);
            }
            plot_data_counter += SpaceDim*SpaceDim;
        }

        { // J
            Interval interv(plot_data_counter, plot_data_counter);
            for (dit.reset(); dit.ok(); ++dit) {
                FArrayBox JFAB(interv, a_plot_data[dit]);
                m_levGeoPtr->fill_J(JFAB);
            }
            plot_data_counter += 1;
        }

        { // Jinv
            Interval interv(plot_data_counter, plot_data_counter);
            for (dit.reset(); dit.ok(); ++dit) {
                FArrayBox JinvFAB(interv, a_plot_data[dit]);
                m_levGeoPtr->fill_Jinv(JinvFAB);
            }
            plot_data_counter += 1;
        }

        { // gdn
            const int ncomps = (SpaceDim * (SpaceDim + 1)) / 2;
            Interval interv(plot_data_counter, plot_data_counter + ncomps - 1);
            for (dit.reset(); dit.ok(); ++dit) {
                FArrayBox gdnFAB(interv, a_plot_data[dit]);
                m_levGeoPtr->fill_gdn(gdnFAB);
            }
            plot_data_counter += ncomps;
        }

        { // gup
            const int ncomps = (SpaceDim * (SpaceDim + 1)) / 2;
            Interval interv(plot_data_counter, plot_data_counter + ncomps - 1);
            for (dit.reset(); dit.ok(); ++dit) {
                FArrayBox gupFAB(interv, a_plot_data[dit]);
                m_levGeoPtr->fill_gup(gupFAB);
            }
            plot_data_counter += ncomps;
        }
    }

    // Stress tensor
#   ifdef USE_STRESSMETRIC
    if (m_stressMetric.isDefined()) {
        const int ncomps = (SpaceDim * (SpaceDim + 1)) / 2;
        Interval interv(plot_data_counter, plot_data_counter + ncomps - 1);
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox stressFAB(interv, a_plot_data[dit]);
            const RealVect& dx = m_levGeoPtr->getDx();

            for (int adir = 0; adir < SpaceDim; ++adir) {
                for (int bdir = adir; bdir < SpaceDim; ++bdir) {
                    const int comp = LevelGeometry::symTensorCompCC(adir, bdir);

                    m_stressMetric.fill_Jgup(stressFAB,
                                             comp,
                                             adir, bdir,
                                             dx);
                }
            }
        }
        plot_data_counter += ncomps;
    }
#   endif //USE_STRESSMETRIC

    // Taylor-Green vortex solution
    if (dynamic_cast<TaylorGreenBCUtil*>(m_physBCPtr)) {
        // Write solution directly to plot holder
        Interval velInt(plot_data_counter, plot_data_counter + SpaceDim - 1);
        Interval pressureInt(plot_data_counter + SpaceDim, plot_data_counter + SpaceDim);
        Interval FofTInt(plot_data_counter + SpaceDim + 1, plot_data_counter + SpaceDim + 1);

        TaylorGreenBCUtil* tgBC = dynamic_cast<TaylorGreenBCUtil*>(m_physBCPtr);
        const Real FofTval = tgBC->FofT(m_time, *m_levGeoPtr);

        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox solFAB(velInt, a_plot_data[dit]);
            FArrayBox pressureFAB(pressureInt, a_plot_data[dit]);
            FArrayBox FofTFAB(FofTInt, a_plot_data[dit]);

            tgBC->fillVelSoln(solFAB, *m_levGeoPtr, dit(), m_time);
            tgBC->fillPressureSoln(pressureFAB, *m_levGeoPtr, dit(), m_time-0.5*m_dt);
            FofTFAB.setVal(FofTval);
        }

        LevelData<FArrayBox> pCalcRef;
        CH_assert(!(piDestComp == Interval()));
        aliasLevelData(pCalcRef, &a_plot_data, piDestComp);

        LevelData<FArrayBox> pSolRef;
        aliasLevelData(pSolRef, &a_plot_data, pressureInt);

        LevelData<FArrayBox> velSolRef;
        aliasLevelData(velSolRef, &a_plot_data, velInt);

        LevelData<FArrayBox> velCalc(grids, SpaceDim);
        newVel().copyTo(velCalc);
        m_levGeoPtr->sendToCartesianBasis(velCalc);

        LevelData<FArrayBox> pError(grids, 1);
        LevelData<FArrayBox> velError(grids, 1);
        LevelData<FArrayBox> velSolMag(grids, 1);
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& pErrorFAB = pError[dit];
            const FArrayBox& pCalcFAB = pCalcRef[dit];
            const FArrayBox& pSolFAB = pSolRef[dit];

            FArrayBox& velErrorFAB = velError[dit];
            FArrayBox& velSolMagFAB = velSolMag[dit];
            const FArrayBox& velCalcFAB = velCalc[dit];
            const FArrayBox& velSolFAB = velSolRef[dit];

            BoxIterator bit(grids[dit]);
            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& cc = bit();

                pErrorFAB(cc) = pCalcFAB(cc) - pSolFAB(cc);

                RealVect ev(D_DECL(
                    velCalcFAB(cc,0)-velSolFAB(cc,0),
                    velCalcFAB(cc,1)-velSolFAB(cc,1),
                    velCalcFAB(cc,2)-velSolFAB(cc,2)));
                velErrorFAB(cc) = ev.vectorLength();

                RealVect sm(D_DECL(velSolFAB(cc,0),
                                   velSolFAB(cc,1),
                                   velSolFAB(cc,2)));
                velSolMagFAB(cc) = sm.vectorLength();
            }
        }

        Real uSolNorm0 = computeMappedNorm(velSolMag, NULL, *m_levGeoPtr, 0);
        Real uSolNorm1 = computeMappedNorm(velSolMag, NULL, *m_levGeoPtr, 1);
        Real uSolNorm2 = computeMappedNorm(velSolMag, NULL, *m_levGeoPtr, 2);

        Real unorm0 = computeMappedNorm(velError, NULL, *m_levGeoPtr, 0) / uSolNorm0;
        Real unorm1 = computeMappedNorm(velError, NULL, *m_levGeoPtr, 1) / uSolNorm1;
        Real unorm2 = computeMappedNorm(velError, NULL, *m_levGeoPtr, 2) / uSolNorm2;
        pout() << "Level " << m_level << " u: inf, 1, and 2 norms / solNorms = "
               << unorm0 << ", " << unorm1 << ", " << unorm2 << endl;
        if (procID() == 0) {
            std::cout << "Level " << m_level << " u: inf, 1, and 2 norms / solNorms = "
                      << unorm0 << ", " << unorm1 << ", " << unorm2 << endl;
        }

        Real pSolNorm0 = computeMappedNorm(pSolRef, NULL, *m_levGeoPtr, 0);
        Real pSolNorm1 = computeMappedNorm(pSolRef, NULL, *m_levGeoPtr, 1);
        Real pSolNorm2 = computeMappedNorm(pSolRef, NULL, *m_levGeoPtr, 2);

        Real pnorm0 = computeMappedNorm(pError, NULL, *m_levGeoPtr, 0) / pSolNorm0;
        Real pnorm1 = computeMappedNorm(pError, NULL, *m_levGeoPtr, 1) / pSolNorm1;
        Real pnorm2 = computeMappedNorm(pError, NULL, *m_levGeoPtr, 2) / pSolNorm2;
        pout() << "Level " << m_level << " p: inf, 1, and 2 norms / solNorms = "
               << pnorm0 << ", " << pnorm1 << ", " << pnorm2 << endl;

        // Increment plot counter
        plot_data_counter += SpaceDim;  // velocity
        ++plot_data_counter;            // pressure
        ++plot_data_counter;            // FofT
    } // end if TaylorGreenBCUtil
}
