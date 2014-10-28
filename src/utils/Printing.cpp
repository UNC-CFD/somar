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
#include <fstream>
#include <sstream>
#include <iomanip>
#include "AMRIO.H"
#include "Printing.H"
#include "Constants.H"
#ifdef CH_USE_HDF5
#   include "CH_HDF5.H"
#endif


template <class T>
int writeLevel(HDF5Handle& a_handle,
               const int&  a_level,
               const T& a_data,
               const RealVect& a_dx, // dx for each direction
               const Real& a_dt,
               const Real& a_time,
               const Box&  a_domain,
               const IntVect&  a_refRatios, // ref ratio for each direction
               const IntVect& outputGhost,
               const Interval& comps)
{
    int error;
    char levelName[10];
    std::string currentGroup = a_handle.getGroup();
    sprintf(levelName, "/level_%i", a_level);
    error = a_handle.setGroup(currentGroup + levelName);
    if (error != 0) return 1;

    HDF5HeaderData meta;
    meta.m_realvect["vec_dx"] = a_dx;
    meta.m_real["dt"] = a_dt;
    meta.m_real["time"] = a_time;
    meta.m_box["prob_domain"] = a_domain;
    meta.m_intvect["vec_ref_ratio"] = a_refRatios;

    error = meta.writeToFile(a_handle);
    if (error != 0) return 2;

    error = write(a_handle, a_data.boxLayout());
    if (error != 0) return 3;

    error = write(a_handle, a_data, "data", outputGhost, comps);
    if (error != 0) return 4;

    a_handle.setGroup(currentGroup);

    return 0;
}

template <class T>
int readLevel(HDF5Handle&   a_handle,
              const int&    a_level,
              LevelData<T>& a_data,
              RealVect& a_dx,
              Real& a_dt,
              Real& a_time,
              Box&  a_domain,
              IntVect&  a_refRatio,
              const Interval&   a_comps,
              bool  setGhost)
{
    HDF5HeaderData header;
    header.readFromFile(a_handle);
    //unused
    // int nComp = header.m_int["num_components"];

    int error;
    char levelName[10];
    std::string currentGroup = a_handle.getGroup();
    sprintf(levelName, "/level_%i", a_level);
    error = a_handle.setGroup(currentGroup + levelName);
    if (error != 0) return 1;

    HDF5HeaderData meta;
    error = meta.readFromFile(a_handle);
    if (error != 0) return 2;
    a_dx       = meta.m_realvect["vec_dx"];
    a_dt       = meta.m_real["dt"];
    a_time     = meta.m_real["time"];
    a_domain   = meta.m_box["prob_domain"];
    a_refRatio = meta.m_intvect["vec_ref_ratio"];
    Vector<Box> boxes;
    error = read(a_handle, boxes);
    Vector<int> procIDs;
    LoadBalance(procIDs, boxes);

    DisjointBoxLayout layout(boxes, procIDs, a_domain);

    layout.close();
    if (error != 0) return 3;

    error = read(a_handle, a_data, "data", layout, a_comps, true);

    if (error != 0) return 4;

    a_handle.setGroup(currentGroup);

    return 0;
}


// -----------------------------------------------------------------------------
// Create a banner
// -----------------------------------------------------------------------------
std::string banner(const std::string& a_str)
{
    if(a_str.length() == 0) {
        return "\n" + std::string(80, '-') + "\n";
    }
    int leftdashlen  = (80 - a_str.length()) / 2 - 1;
    int rightdashlen = 78 - leftdashlen - a_str.length();
    return "\n" + std::string(leftdashlen, '-') + " " + a_str + " " + std::string(rightdashlen, '-') + "\n";
}


// -----------------------------------------------------------------------------
// Converts a time in seconds (double) into the more readable days, hours, mins,
// and seconds (std::string).
// -----------------------------------------------------------------------------
std::string formatTime (const double a_seconds)
{
    std::ostringstream retStr;
    retStr << setiosflags(ios::fixed);

    if (a_seconds < 60) {
        retStr << a_seconds << " secs";
    } else {
        const long longSecs = (long)a_seconds;
        const int seconds = longSecs % 60;
        const int minutes = (longSecs / 60) % 60;
        const int hours = (longSecs / (60 * 60)) % 24;
        const int days = longSecs / (24 * 60 * 60);

        if (days > 0)
            retStr << days << " days, ";
        if (days > 0 || hours > 0)
            retStr << hours << " hrs, ";
        if (days > 0 || hours > 0 || minutes > 0)
            retStr << minutes << " mins, ";
        retStr << seconds << " secs";
    }

    return retStr.str();
}


// -----------------------------------------------------------------------------
// Write a FArrayBox to a text file
// -----------------------------------------------------------------------------
void writeTextFile( const FArrayBox&    a_data,
                    const char*         a_filename,
                    const Real&         a_time)
{
    char fn[50];
    sprintf(fn, "%s.%d", a_filename, procID());
    fstream fstr;
    fstr.open(fn, fstream::out);

    if(a_time != -1) {
        fstr << "Time = " << a_time << "\n\n";
    }

    for (BoxIterator bit(a_data.box()); bit.ok(); ++bit) {
        fstr << bit();
        for (int comp = 0; comp < a_data.nComp(); ++comp) {
            fstr << "\t" << a_data(bit(), comp);
        }
        fstr << "\n";
    }

    fstr.close();
}


// -----------------------------------------------------------------------------
// Write a LevelData to a text file
// -----------------------------------------------------------------------------
void writeTextFile( const LevelData<FArrayBox>& a_data,
                    const char*                 a_filename,
                    const Real&                 a_time,
                    const bool                  a_writeGhosts)
{
    const DisjointBoxLayout& grids = a_data.getBoxes();

    char fn[50];
    sprintf(fn, "%s.%d", a_filename, procID());
    fstream fstr;
    fstr.open(fn, fstream::out);

    if(a_time != -1) {
        fstr << "Time = " << a_time << "\n\n"
             << "domain = " << grids.physDomain() << "\n\n"
             << "dbl = " << grids << std::endl;
    }

    for (DataIterator dit(grids); dit.ok(); ++dit) {
        Box region = a_data[dit].box();
        if (!a_writeGhosts) {
            region &= grids[dit];
        }
        for (BoxIterator bit(region); bit.ok(); ++bit) {
            fstr << bit();
            for (int comp = 0; comp < a_data.nComp(); ++comp) {
                fstr << "\t" << a_data[dit](bit(), comp);
            }
            fstr << "\n";
        }
        fstr << std::endl;
    }
    fstr.close();
}


// -----------------------------------------------------------------------------
void writeTextFile( const FluxBox&  a_data,
                    const char*     a_filename,
                    const Real&     a_time)
{
    char fn[50];
    sprintf(fn, "%s.%d", a_filename, procID());
    fstream fstr;
    fstr.open(fn, fstream::out);

    if(a_time != -1) {
        fstr << "Time = " << a_time << "\n\n";
    }

    for (int dir = 0; dir < SpaceDim; ++dir) {
        fstr << " ------------ dir = " << dir << " ------------------" << std::endl;
        for (BoxIterator bit(a_data[dir].box()); bit.ok(); ++bit) {
            fstr << bit();
            for (int comp = 0; comp < a_data.nComp(); ++comp) {
                fstr << "\t" << a_data[dir](bit(), comp);
            }
            fstr << "\n";
        }
    }

    fstr.close();
}


// -----------------------------------------------------------------------------
void writeTextFile( const LevelData<FluxBox>&   a_data,
                    const char*                 a_filename,
                    const Real&                 a_time)
{
    char fn[50];
    sprintf(fn, "%s.%d", a_filename, procID());
    fstream fstr;
    fstr.open(fn, fstream::out);

    if(a_time != -1) {
        fstr << "Time = " << a_time << "\n\n"
             << "domain = " << a_data.disjointBoxLayout().physDomain() << "\n\n"
             << "dbl = " << a_data.disjointBoxLayout() << std::endl;
    }

    for (int dir = 0; dir < SpaceDim; ++dir) {
        fstr << " ------------ dir = " << dir << " ------------------" << std::endl;
        for (DataIterator dit(a_data.disjointBoxLayout()); dit.ok(); ++dit) {
            for (BoxIterator bit(a_data[dit()].box()); bit.ok(); ++bit) {
                fstr << bit();
                for (int comp = 0; comp < a_data.nComp(); ++comp) {
                    fstr << "\t" << a_data[dit()][dir](bit(), comp);
                }
                fstr << "\n";
            }
            fstr << std::endl;
        }
    }
    fstr.close();
}


#include "AMRLESMeta.H"
// -----------------------------------------------------------------------------
// Write a LevelData to HDF5
// -----------------------------------------------------------------------------
void _writeLevelHDF5 (const LevelData<FArrayBox>& a_data,
                      const char*                 a_filename,
                      Real                        a_time,
                      bool                        a_oneGhost)
{
    {
        // Figure out the centering of the data
        IntVect dataType = IntVect::Zero;
        DataIterator dit = a_data.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            const Box thisBox = a_data[dit].box();
            if (thisBox.isEmpty()) continue;
            dataType = thisBox.type();
            break;
        }

#if CH_MPI
        // It's possible this rank had no grids. In that case, we must get the
        // dataType from the other ranks.
        for (int dir = 0; dir < SpaceDim; ++dir) {
            int localDataTypeDir = dataType[dir];
            int dataTypeDir;
            int ierr = MPI_Allreduce(&localDataTypeDir, &dataTypeDir, 1, MPI_INT, MPI_SUM, AMRLESMeta::amrComm);
            if (ierr != MPI_SUCCESS) {
                std::ostringstream errmsg;
                errmsg << "MPI_Allreduce failed. Error " << ierr << std::endl;
                MayDay::Error(errmsg.str().c_str());
            }
            dataType[dir] = ((dataTypeDir == 0)? 0: 1);
        }
#endif

        if (dataType == IntVect::Zero) {
            // Do nothing special
        } else if (dataType.sum() == 1) {
            // Call the FluxBox version
            LevelData<FluxBox> fluxData(a_data.getBoxes(), a_data.nComp(), a_data.ghostVect());
            D_TERM(int FCdir = 0;,
                   if (dataType[1] == 1) FCdir = 1;,
                   if (dataType[2] == 1) FCdir = 2;)

            for (dit.reset(); dit.ok(); ++dit) {
                fluxData[dit].setVal(0.0);
                fluxData[dit][FCdir].copy(a_data[dit]);
            }

            _writeLevelHDF5(fluxData, a_filename, a_time, a_oneGhost);
            return;

        } else {
            // Throw an error
            pout() << "dataType = " << dataType << endl;
            MayDay::Error("_writeLevelHDF5 only works with CC or FC data and needs comm code for dataType");
        }
    }

    Vector<LevelData<FArrayBox>*> vData(1);
    if (a_oneGhost) {
        const DisjointBoxLayout& grids = a_data.getBoxes();
        vData[0] = new LevelData<FArrayBox>(grids, a_data.nComp(), IntVect::Unit);

        DataIterator dit = grids.dataIterator();
        for (dit.reset(); dit.ok(); ++dit) {
            (*vData[0])[dit].copy(a_data[dit]);
        }
        // Copier oneGhostCopier(grids, grids, grids.physDomain(), IntVect::Unit, false);
        // a_data.copyTo(a_data.interval(), *vData[0], vData[0]->interval(), oneGhostCopier);
    } else {
        vData[0] = (LevelData<FArrayBox>*)(&a_data);
    }

    Vector<DisjointBoxLayout> vGrids(1);
    vGrids[0] = a_data.getBoxes();
    const Box& domainBox = vGrids[0].physDomain().domainBox();

    Vector<string> vNames(vData[0]->nComp());
    for (int comp = 0; comp < vData[0]->nComp(); ++comp) {
        ostringstream compName;
        compName << "comp " << comp;
        vNames[comp] = compName.str();
    }

    Vector<int> refRatios(1,1);
    int numLevels = 1;
    Real dx = 1.0;
    Real dt = 1.0;

    WriteAMRHierarchyHDF5(string(a_filename),
                          vGrids, vData, vNames,
                          domainBox, dx, dt, a_time,
                          refRatios, numLevels);

    if (a_oneGhost) {
        delete vData[0];
    }
}


// -----------------------------------------------------------------------------
void _writeLevelHDF5 (const LevelData<FluxBox>& a_data,
                      const char*               a_filename,
                      Real                      a_time,
                      bool                      a_oneGhost)
{
    const DisjointBoxLayout& grids = a_data.getBoxes();
    const Box& domainBox = grids.physDomain().domainBox();
    const int nOldComps = a_data.nComp();
    const int nNewComps = nOldComps * 2 * SpaceDim;
    DataIterator dit = a_data.dataIterator();
    const IntVect& ghostVect = a_oneGhost? IntVect::Unit: a_data.ghostVect();

    Vector<LevelData<FArrayBox>*> vData(1);
    vData[0] = new LevelData<FArrayBox>(grids, nNewComps, ghostVect);

    Vector<DisjointBoxLayout> vGrids(1);
    vGrids[0] = grids;

    Vector<string> vNames(nNewComps);
    for (int newComp = 0; newComp < nNewComps; newComp += 2*SpaceDim) {
        int oldComp = newComp / (2*SpaceDim);
        D_TERM(
        {
            ostringstream newCompName;
            newCompName << "left x comp " << oldComp;
            vNames[newComp] = newCompName.str();
        }
        {
            ostringstream newCompName;
            newCompName << "right x comp " << oldComp;
            vNames[newComp+1] = newCompName.str();
        },
        {
            ostringstream newCompName;
            newCompName << "left y comp " << oldComp;
            vNames[newComp+2] = newCompName.str();
        }
        {
            ostringstream newCompName;
            newCompName << "right y comp " << oldComp;
            vNames[newComp+3] = newCompName.str();
        },
        {
            ostringstream newCompName;
            newCompName << "left z comp " << oldComp;
            vNames[newComp+4] = newCompName.str();
        }
        {
            ostringstream newCompName;
            newCompName << "right z comp " << oldComp;
            vNames[newComp+5] = newCompName.str();
        })

        for (int FCdir = 0; FCdir < SpaceDim; ++FCdir) {
            for (dit.reset(); dit.ok(); ++dit) {
                FArrayBox& oldFAB = (FArrayBox&)(a_data[dit][FCdir]);
                FArrayBox& newFAB = (*vData[0])[dit];

                newFAB.setVal(quietNAN, newComp + 2*FCdir);
                newFAB.setVal(quietNAN, newComp + 2*FCdir + 1);

                oldFAB.shiftHalf(FCdir,1);
                newFAB.copy(oldFAB, oldComp, newComp + 2*FCdir, 1);

                oldFAB.shiftHalf(FCdir,-2);
                newFAB.copy(oldFAB, oldComp, newComp + 2*FCdir + 1, 1);

                oldFAB.shiftHalf(FCdir,1);
            }
        }
    }

    Vector<int> refRatios(1,1);
    int numLevels = 1;
    Real dx = 1.0;
    Real dt = 1.0;

    WriteAMRHierarchyHDF5(string(a_filename),
                          vGrids, vData, vNames,
                          domainBox, dx, dt, a_time,
                          refRatios, numLevels);

    delete vData[0];
}


// -----------------------------------------------------------------------------
void _writeHDF5 (const Vector<LevelData<FArrayBox>*>& a_data,
                 const LevelGeometry&                 a_levGeo,
                 int                                  a_lmin,
                 int                                  a_lmax,
                 const char*                          a_filename,
                 const Vector<std::string>&           a_names)
{
    const int finestLevel = a_data.size() - 1;
    if (a_lmax < 0) a_lmax = finestLevel;

    CH_assert(0 <= a_lmin);
    CH_assert(a_lmin <= a_lmax);
    CH_assert(a_lmax <= finestLevel);

    // Gather geometric data
    const Vector<const LevelGeometry*> vLevGeo = a_levGeo.getAMRLevGeos();
    const Vector<DisjointBoxLayout> vGrids = a_levGeo.getAMRGrids();
    const Vector<IntVect> vRefRatio = a_levGeo.getAMRRefRatios();
    const ProblemDomain& lev0Domain = vLevGeo[0]->getDomain();
    const RealVect lev0dx = vLevGeo[0]->getDx();

    // How many comps will we need to copy? How many ghosts are there?
    int numComps = -1;
    IntVect ghostVect = IntVect::Unit; // At least one for proper display of mapping.
    for (int lev = 0; lev <= finestLevel; ++lev) {
        if (a_data[lev] == NULL) continue;

        int thisNumComps = a_data[lev]->nComp();
        if (numComps == -1) numComps = thisNumComps;
        CH_assert(thisNumComps == numComps);

        IntVect thisGhostVect = a_data[lev]->ghostVect();
        D_TERM(ghostVect[0] = Max(ghostVect[0], thisGhostVect[0]);,
               ghostVect[1] = Max(ghostVect[1], thisGhostVect[1]);,
               ghostVect[2] = Max(ghostVect[2], thisGhostVect[2]);)
    }
    CH_assert(numComps > 0);

    // Package output into one CC holder
    Vector<LevelData<FArrayBox>*> vOutput(finestLevel+1, NULL);
    for (int lev = 0; lev <= finestLevel; ++lev) {
        const DisjointBoxLayout& grids = vGrids[lev];
        DataIterator dit = grids.dataIterator();

        // No matter what, initialize output and fill the displacement field.
        vOutput[lev] = new LevelData<FArrayBox>(grids, numComps+SpaceDim, ghostVect);
        Interval dispInt(numComps, numComps+SpaceDim-1);
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& outputFAB = (*vOutput[lev])[dit];
            outputFAB.setVal(quietNAN);

            FArrayBox dispFAB(dispInt, (*vOutput[lev])[dit]);
            vLevGeo[lev]->fill_displacement(dispFAB);
        }

        // Do we want this level's data?
        if (a_data[lev] == NULL) continue;
        if (lev < a_lmin) continue;
        if (a_lmax < lev) continue;

        // Copy the data
        for (dit.reset(); dit.ok(); ++dit) {
            FArrayBox& outputFAB = (*vOutput[lev])[dit];
            const FArrayBox& dataFAB = (*a_data[lev])[dit];

            outputFAB.copy(dataFAB, 0, 0, numComps);
        }
    }

    // Write to file
    {
        Vector<string> vNames(numComps + SpaceDim);
        int comp = 0;
        if (a_names.size() == 0) {
            if (numComps == SpaceDim) {
                // It's probably a vector.
                D_TERM(vNames[comp++] = "x_component";,
                       vNames[comp++] = "y_component";,
                       vNames[comp++] = "z_component";)
            } else {
                // It's probably not a vector.
                while (comp < numComps) {
                    char compName[40];
                    sprintf(compName, "comp_%d", comp);
                    vNames[comp++] = std::string(compName);
                }
            }
        } else {
            // The caller provided names.
            while (comp < numComps) {
                char compName[40];
                sprintf(compName, "comp_%d", comp);
                vNames[comp] = a_names[comp];
                ++comp;
            }
        }
        D_TERM(
        vNames[comp++] = "x_Displacement";,
        vNames[comp++] = "y_Displacement";,
        vNames[comp++] = "z_Displacement";)
        CH_assert(comp == numComps+SpaceDim);

        Real dt = 1.0;
        Real dummyTime = 0.0;
        WriteAnisotropicAMRHierarchyHDF5(string(a_filename),
                                         vGrids, vOutput, vNames,
                                         lev0Domain.domainBox(),
                                         lev0dx,
                                         dt, dummyTime,
                                         vRefRatio,
                                         vOutput.size());
    }

    // Free memory
    for (int lev = 0; lev <= finestLevel; ++lev) {
        delete vOutput[lev];
    }
}

void _writeHDF5 (const Vector<LevelData<FluxBox>*>& a_data,
                 const LevelGeometry&               a_levGeo,
                 int                                a_lmin,
                 int                                a_lmax,
                 const char*                        a_filename)
{
    const int finestLevel = a_data.size() - 1;
    if (a_lmax < 0) a_lmax = finestLevel;

    CH_assert(0 <= a_lmin);
    CH_assert(a_lmin <= a_lmax);
    CH_assert(a_lmax <= finestLevel);

    // How many comps will we need to copy? How many ghosts are there?
    int numComps = -1;
    IntVect ghostVect = IntVect::Unit; // At least one for proper display of mapping.
    for (int lev = 0; lev <= finestLevel; ++lev) {
        if (a_data[lev] == NULL) continue;

        int thisNumComps = a_data[lev]->nComp();
        if (numComps == -1) numComps = thisNumComps;
        CH_assert(thisNumComps == numComps);

        IntVect thisGhostVect = a_data[lev]->ghostVect();
        D_TERM(ghostVect[0] = Max(ghostVect[0], thisGhostVect[0]);,
               ghostVect[1] = Max(ghostVect[1], thisGhostVect[1]);,
               ghostVect[2] = Max(ghostVect[2], thisGhostVect[2]);)
    }
    CH_assert(numComps > 0);

    const int nNewComps = numComps * 2 * SpaceDim;
    Vector<string> vNames(nNewComps);
    for (int newComp = 0; newComp < nNewComps; newComp += 2*SpaceDim) {
        int oldComp = newComp / (2*SpaceDim);
        D_TERM(
        {
            ostringstream newCompName;
            newCompName << "left x comp " << oldComp;
            vNames[newComp] = newCompName.str();
        }
        {
            ostringstream newCompName;
            newCompName << "right x comp " << oldComp;
            vNames[newComp+1] = newCompName.str();
        },
        {
            ostringstream newCompName;
            newCompName << "left y comp " << oldComp;
            vNames[newComp+2] = newCompName.str();
        }
        {
            ostringstream newCompName;
            newCompName << "right y comp " << oldComp;
            vNames[newComp+3] = newCompName.str();
        },
        {
            ostringstream newCompName;
            newCompName << "left z comp " << oldComp;
            vNames[newComp+4] = newCompName.str();
        }
        {
            ostringstream newCompName;
            newCompName << "right z comp " << oldComp;
            vNames[newComp+5] = newCompName.str();
        })
    }

    Vector<LevelData<FArrayBox>*> ccData(finestLevel+1, NULL);
    for (int lev = 0; lev <= finestLevel; ++lev) {
        const DisjointBoxLayout& grids = a_data[lev]->getBoxes();
        DataIterator dit = grids.dataIterator();

        ccData[lev] = new LevelData<FArrayBox>(grids, nNewComps, ghostVect);
        for (dit.reset(); dit.ok(); ++dit) {
            (*ccData[lev])[dit].setVal(quietNAN);
        }

        if (a_data[lev] == NULL) continue;
        if (lev < a_lmin) continue;
        if (a_lmax < lev) continue;

        for (int newComp = 0; newComp < nNewComps; newComp += 2*SpaceDim) {
            int oldComp = newComp / (2*SpaceDim);

            for (int FCdir = 0; FCdir < SpaceDim; ++FCdir) {
                for (dit.reset(); dit.ok(); ++dit) {
                    FArrayBox& oldFAB = (FArrayBox&)((*a_data[lev])[dit][FCdir]);
                    FArrayBox& newFAB = (*ccData[lev])[dit];

                    newFAB.setVal(quietNAN, newComp + 2*FCdir);
                    newFAB.setVal(quietNAN, newComp + 2*FCdir + 1);

                    oldFAB.shiftHalf(FCdir,1);
                    newFAB.copy(oldFAB, oldComp, newComp + 2*FCdir, 1);

                    oldFAB.shiftHalf(FCdir,-2);
                    newFAB.copy(oldFAB, oldComp, newComp + 2*FCdir + 1, 1);

                    oldFAB.shiftHalf(FCdir,1);
                }
            }
        }
    }

    _writeHDF5(ccData, a_levGeo, a_lmin, a_lmax, a_filename, vNames);

    for (int lev = 0; lev <= finestLevel; ++lev) {
        delete ccData[lev];
    }
}


/*
\\ write out hierarchy of anisotropic amr data in HDF5 format
\\ filename,  == file to output to
\\ a_vectData == data at each level
\\ a_vectNames== names of variables
\\ a_domain == domain at coarsest level
\\ a_dx     == grid spacing in each direction at coarsest level
\\ a_dt     == time step at coarsest level
\\ a_time     == time
\\ a_vectRatio == refinement ratio in each direction at all levels
\\ (ith entry is refinement ratio between levels i and i + 1)
\\ a_numLevels == number of levels to output
*/
void
WriteAnisotropicAMRHierarchyHDF5(
    const string& filename,
    const Vector<DisjointBoxLayout>& a_vectGrids,
    const Vector<LevelData<FArrayBox>* >& a_vectData,
    const Vector<string>& a_vectNames,
    const Box& a_domain,
    const RealVect& a_dx,
    const Real& a_dt,
    const Real& a_time,
    const Vector<IntVect>& a_refRatios,
    const int& a_numLevels)
{
    HDF5Handle handle(filename.c_str(),  HDF5Handle::CREATE);

    WriteAnisotropicAMRHierarchyHDF5(
        handle, a_vectGrids, a_vectData, a_vectNames,
        a_domain, a_dx, a_dt, a_time, a_refRatios, a_numLevels);

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif

    handle.close();
}

void
WriteAnisotropicAMRHierarchyHDF5(
    HDF5Handle& handle,
    const Vector<DisjointBoxLayout>& a_vectGrids,
    const Vector<LevelData<FArrayBox>* >& a_vectData,
    const Vector<string>& a_vectNames,
    const Box& a_domain,
    const RealVect& a_dx, // Grid spacing in each direction
    const Real& a_dt,
    const Real& a_time,
    const Vector<IntVect>& a_refRatios, // for each level, in each direction
    const int& a_numLevels)
{
    CH_assert(a_numLevels > 0);
    CH_assert(a_vectData.size()  >= a_numLevels);
    CH_assert(a_refRatios.size() >= a_numLevels - 1);

    HDF5HeaderData header;
    int nComp = a_vectNames.size();

    string filedescriptor("VanillaAMRFileType");
    header.m_string ["filetype"]      = filedescriptor;
    header.m_int ["num_levels"]       = a_numLevels;
    header.m_int ["num_components"]    = nComp;

    for (int ivar = 0; ivar < nComp; ivar++) {
        char labelChSt[100];
        sprintf(labelChSt, "component_%d", ivar);
        string label(labelChSt);
        header.m_string[label] = a_vectNames[ivar];
    }
    header.writeToFile(handle);

    Box domainLevel = a_domain;
    Real dtLevel = a_dt;
    RealVect dxLevel = a_dx;
    for (int ilev = 0; ilev < a_numLevels; ilev++) {
        IntVect refLevel = IntVect::Unit;
        if (ilev != a_numLevels - 1) {
            refLevel = a_refRatios[ilev];
        }
        if (ilev != 0) {
            domainLevel.refine(a_refRatios[ilev - 1]);
            dtLevel /= a_refRatios[ilev - 1][0]; // HACK - just use 0 dir ref ratio
            dxLevel /= a_refRatios[ilev - 1];
        }
        CH_assert(a_vectData[ilev] != NULL);
        const LevelData<FArrayBox>& dataLevel = *a_vectData[ilev];
        CH_assert(dataLevel.nComp() == nComp);
        Interval comps(0, nComp - 1);
        IntVect ghostVect = a_vectData[0]->ghostVect();
        int eek = writeLevel(handle, ilev, dataLevel,
                             dxLevel, dtLevel, a_time,
                             domainLevel, refLevel, ghostVect, comps);
        if (eek != 0) {
            MayDay::Error("WriteAMRHierarchyHDF5: Error in writeLevel");
        }
    }
}

//
/*
\\ Read in hierarchy of amr data in ANISOTROPIC HDF5 format
\\ filename,  == file to output to
\\ a_vectData == data at each level
\\ a_vectNames== names of variables
\\ a_domain == domain at coarsest level
\\ a_dx     == grid spacing at coarsest level
\\ a_dt     == time step at coarsest level
\\ a_time     == time
\\ a_vectRatio == refinement ratio at all levels
\\ (ith entry is refinement ratio between levels i and i + 1)
\\ a_numLevels == number of levels to output

return values:
0: success
-1: bogus number of levels
-2: bogus number of components
-3: error in readlevel
-4: file open failed
*/
int
ReadAnisotropicAMRHierarchyHDF5(const string& filename,
                                Vector<DisjointBoxLayout>& a_vectGrids,
                                Vector<LevelData<FArrayBox>* >& a_vectData,
                                Vector<string>& a_vectNames,
                                Box& a_domain,
                                RealVect& a_dx,
                                Real& a_dt,
                                Real& a_time,
                                Vector<IntVect>& a_refRatio,
                                int& a_numLevels)
{
    HDF5Handle handle;
    int err = handle.open(filename.c_str(),  HDF5Handle::OPEN_RDONLY);
    if ( err < 0) {
        return -4;
    }
    int eekflag = ReadAnisotropicAMRHierarchyHDF5(handle, a_vectGrids, a_vectData,
                  a_vectNames, a_domain, a_dx, a_dt,
                  a_time, a_refRatio, a_numLevels);

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
    handle.close();

    return (eekflag);
}

int
ReadAnisotropicAMRHierarchyHDF5(HDF5Handle& handle,
                                Vector<DisjointBoxLayout>& a_vectGrids,
                                Vector<LevelData<FArrayBox>* >& a_vectData,
                                Vector<string>& a_vectNames,
                                Box& a_domain,
                                RealVect& a_dx,
                                Real& a_dt,
                                Real& a_time,
                                Vector<IntVect>& a_refRatio,
                                int& a_numLevels)
{

    HDF5HeaderData header;
    header.readFromFile(handle);

    a_numLevels = header.m_int["num_levels"];
    if (a_numLevels <= 0) {
        MayDay::Warning("ReadAnisotropicAMRHierarchyHDF5: Bogus number of levels");
        return (-1);
    }
    a_vectData.resize(a_numLevels);
    a_refRatio.resize(a_numLevels);
    a_vectGrids.resize(a_numLevels);

    int nComp = header.m_int["num_components"];
    if (nComp <= 0) {
        MayDay::Warning("ReadAnisotropicAMRHierarchyHDF5: Bogus number of Components");
        return (-2);
    }
    a_vectNames.resize(nComp);

    for (int ivar = 0; ivar < nComp; ivar++) {
        char labelChSt[100];
        sprintf(labelChSt, "component_%d", ivar);
        string label(labelChSt);
        a_vectNames[ivar] = header.m_string[label];
    }
    for (int ilev = 0; ilev < a_numLevels; ilev++) {
        IntVect refLevel = IntVect::Zero;
        Box domainLevel;
        Real dtLevel;
        RealVect dxLevel;
        a_vectData[ilev] = new LevelData<FArrayBox>();
        int eek = readLevel(handle, ilev, *(a_vectData[ilev]),
                            dxLevel, dtLevel,  a_time,
                            domainLevel, refLevel, Interval(), true);
        if (eek != 0) {
            MayDay::Warning("ReadAnisotropicAMRHierarchyHDF5: readLevel failed");
            return (-3);
        }

        const DisjointBoxLayout& dbl = a_vectData[ilev]->getBoxes();
        a_vectGrids[ilev] = dbl;

        if (ilev == 0) {
            a_domain = domainLevel;
            a_dt = dtLevel;
            a_dx = dxLevel;
        }
        a_refRatio[ilev] = refLevel;
    }
    return (0);
}

int
ReadAnisotropicAMRHierarchyHDF5(const string& filename,
                                Vector<DisjointBoxLayout>& a_vectGrids,
                                Vector<LevelData<FArrayBox>* >& a_vectData,
                                Box& a_domain,
                                Vector<IntVect>& a_refRatio,
                                int& a_numLevels)
{
    HDF5Handle handle;
    int err = handle.open(filename.c_str(),  HDF5Handle::OPEN_RDONLY);
    if ( err < 0) {
        return -4;
    }

    int eekflag = ReadAnisotropicAMRHierarchyHDF5(handle, a_vectGrids, a_vectData,
                  a_domain, a_refRatio, a_numLevels);

#ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
#endif
    handle.close();
    return (eekflag);
}

int
ReadAnisotropicAMRHierarchyHDF5(HDF5Handle& handle,
                                Vector<DisjointBoxLayout>& a_vectGrids,
                                Vector<LevelData<FArrayBox>* >& a_vectData,
                                Box& a_domain,
                                Vector<IntVect>& a_refRatio,
                                int& a_numLevels)
{
    HDF5HeaderData header;
    header.readFromFile(handle);

    a_numLevels = header.m_int["num_levels"];
    if (a_numLevels <= 0) {
        MayDay::Warning("ReadAnisotropicAMRHierarchyHDF5: Bogus number of levels");
        return (-1);
    }
    a_vectData.resize(a_numLevels);
    a_refRatio.resize(a_numLevels);
    a_vectGrids.resize(a_numLevels);

    //  int nComp = header.m_int["num_components"];
    for (int ilev = 0; ilev < a_numLevels; ilev++) {
        IntVect refLevel = IntVect::Zero;
        Box domainLevel;
        Real dtLevel;
        RealVect dxLevel;
        Real time;
        a_vectData[ilev] = new LevelData<FArrayBox>();
        int eek = readLevel(handle, ilev, *(a_vectData[ilev]),
                            dxLevel, dtLevel,  time,
                            domainLevel, refLevel, Interval(), true);
        if (eek != 0) {
            MayDay::Warning("ReadAnisotropicAMRHierarchyHDF5: readLevel failed");
            return (-3);
        }

        const DisjointBoxLayout& dbl = a_vectData[ilev]->getBoxes();
        a_vectGrids[ilev] = dbl;

        if (ilev == 0) {
            a_domain = domainLevel;
        }
        a_refRatio[ilev] = refLevel;
    }

    return (0);
}
