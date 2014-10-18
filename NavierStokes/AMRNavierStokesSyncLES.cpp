#include "AMRNavierStokes.H"
#include "Constants.H"
#include "AMRLESMeta.H"
#include "ExtrapolationUtils.H"
#include "MiscUtils.H"


// -----------------------------------------------------------------------------
// Provides syncing with LES code
// -----------------------------------------------------------------------------
void AMRNavierStokes::syncWithLES ()
{
    // This code is incomplete and not yet intended for public use.
    return;

#if 0
    CH_TIME("AMRNavierStokes::syncWithLES");
    TODO();

    // We only want to communicate with the LES when the entire hierarchy has
    // completed the timestep.
    if (m_level != 0) return;

    // Are we even using an LES?
    // if (AMRLESMeta::lesSize <= 0) return;

    // Compute phase, cycle number, eta of communication, etc.
    static const Real tidalPeriod = 2.0 * Pi / s_tidalOmega;
    Real timeIntoPhase = fmod(m_time, tidalPeriod);
    int periodNum = int(m_time / tidalPeriod);
    float phaseAngleOn2Pi = timeIntoPhase / tidalPeriod;

    Real eta = timeIntoPhase - (0.19101 * tidalPeriod); // Negative while approaching target
    // Real eta = timeIntoPhase - (0.25 * tidalPeriod); // Negative while approaching target

    // TEMPORARY!!! This should be in pout.
    // Write the phase and eta to stdout
    if (procID() == 0) {
        ostringstream msg;
        msg << "phase angle/(2*Pi) = " << phaseAngleOn2Pi << ", timeIntoPhase = " << timeIntoPhase << ", eta = " << eta << endl;
        cout << msg.str().c_str() << flush;
    }

    // Did we already sync this cycle?
    static bool calledForThisPhase = false;
    if (eta < -TIME_EPS) {
        calledForThisPhase = false;
        return;
    } else if (-TIME_EPS <= eta && !calledForThisPhase) {
        calledForThisPhase = true;
    } else {
        return;
    }

    // Is this a cycle that we want to sync?
    if (periodNum < 1) return;

    // TEMPORARY!!! This should be in pout.
    // All systems go. Let the user know.
    if (procID() == 0) {
        ostringstream msg;
        msg << "Syncing with LES. periodNum " << periodNum << endl;
        cout << msg.str().c_str() << flush;
    }


    // 1. Determine where the interesting section is.
    const Real X0 = -730.0; //-650.0;
    const Real height = 300.0;
    IntVect loIV, hiIV;
    Box testBox;
    {
        loIV[0] = -int(round(Real(m_problem_domain.domainBox().size(0)) * 0.4));
        loIV[1] = -int(round(Real(m_problem_domain.domainBox().size(1)) * 0.4));
        loIV[SpaceDim-1] = 0;

        hiIV[0] = 0;
        hiIV[1] = 0;
        hiIV[SpaceDim-1] = int(round(Real(m_problem_domain.domainBox().size(SpaceDim-1)) * 1.0));

        testBox.define(loIV, hiIV);
        pout() << "testBox = " << testBox << endl;
    }

    // // Let's see it!
    // MayDay::Warning("Destroying newScal");
    // LevelData<FArrayBox>& bpert = newScal(0);
    // for (DataIterator dit(grids); dit.ok(); ++dit) {
    //     FArrayBox& bFAB = bpert[dit];

    //     Box overlap = bFAB.box();
    //     overlap &= testBox;
    //     bFAB.setVal(1.0, overlap, 0);
    // }

    const int srcProc = uniqueProc(SerialTask::compute);
    const int thisProc = procID();
    Box lesBox;
    Real point0;
    Real zmin;
    {
        DisjointBoxLayout tgrids;
        defineOneProcGrids(tgrids, m_problem_domain, testBox, srcProc);
        DataIterator tdit;
        tdit = tgrids.dataIterator();
        tdit.reset();
        CH_assert(tdit.ok() || thisProc != srcProc);

        // Get the velocity over tgrids.
        LevelData<FArrayBox> tVel(tgrids, SpaceDim);
        newVel().copyTo(tVel);
        m_levGeoPtr->sendToCartesianBasis(tVel, true);

        LevelData<FArrayBox> testCartPos(tgrids, SpaceDim);
        if (thisProc == srcProc) {
            m_levGeoPtr->fill_physCoor(testCartPos[tdit]);

            BoxIterator bit(testBox);
            Real loZ = 1e10;
            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& cc = bit();
                if (testCartPos[tdit](cc,SpaceDim-1) < loZ) {
                    loZ = testCartPos[tdit](cc,SpaceDim-1);
                }
            }

            // This is the old version
            bit.reset();
            int loI = bit()[0];
            int hiK = bit()[SpaceDim-1];
            Real curX = testCartPos[tdit](bit(),0);
            Real hiZ = testCartPos[tdit](bit(),SpaceDim-1);
            for (; bit.ok(); ++bit) {
                const IntVect& cc = bit();
                if (abs(testCartPos[tdit](cc,0)-X0) < abs(curX-X0)) {
                    loI = cc[0];
                    curX = testCartPos[tdit](cc,0);
                }
                if (abs(testCartPos[tdit](cc,SpaceDim-1)-loZ-height) < abs(hiZ-loZ-height)) {
                    hiK = cc[SpaceDim-1];
                    hiZ = testCartPos[tdit](cc,SpaceDim-1);
                }
            }

            // // This version scans test box for the largest |u| profile.
            // bit.reset();
            // int loI = bit()[0];
            // int hiK = bit()[SpaceDim-1];
            // Real hiZ = testCartPos[tdit](bit(),SpaceDim-1);
            // Real maxU = 0.0;
            // for (bit.reset(); bit.ok(); ++bit) {
            //     const IntVect& cc = bit();

            //     if (Abs(tVel[tdit](cc,0)) > maxU) {
            //         maxU = Abs(tVel[tdit](cc,0));
            //         loI = cc[0];
            //     }
            // }
            // for (bit.reset(); bit.ok(); ++bit) {
            //     const IntVect& cc = bit();
            //     if (cc[0] != loI) continue;

            //     if (abs(testCartPos[tdit](cc,SpaceDim-1)-loZ-height) < abs(hiZ-loZ-height)) {
            //         hiK = cc[SpaceDim-1];
            //         hiZ = testCartPos[tdit](cc,SpaceDim-1);
            //     }
            // }
            // pout() << "maxU = " << maxU << endl;

            loIV[0] = loI;
            loIV[SpaceDim-1] = 0;
            hiIV[0] = loI;
            hiIV[SpaceDim-1] = hiK;
            point0 = (loZ + hiZ) * 0.5;
            lesBox.define(loIV, hiIV);

            Box bottomBox = bdryLo(lesBox, SpaceDim-1, 1);
            FArrayBox bottomCartPosFAB(bottomBox, SpaceDim);
            m_levGeoPtr->fill_physCoor(bottomCartPosFAB);
            zmin = bottomCartPosFAB(bottomBox.smallEnd(), SpaceDim-1);

            pout() << "point0 = " << point0 << endl;
            pout() << "lesBox = " << lesBox << endl;
            pout() << "cartesian box = (("
                   << testCartPos[tdit](loIV,0) << ", " << testCartPos[tdit](loIV,SpaceDim-1)
                   << "), ("
                   << testCartPos[tdit](hiIV,0) << ", " << testCartPos[tdit](hiIV,SpaceDim-1)
                   << "))" << endl;
            pout() << "zmin = " << zmin << endl;
        }
    }

    const DisjointBoxLayout& grids = m_levGeoPtr->getBoxes();

    // // Let's see it!
    // MayDay::Warning("Destroying newScal");
    // LevelData<FArrayBox>& bpert = newScal(0);
    // for (DataIterator dit(grids); dit.ok(); ++dit) {
    //     FArrayBox& bFAB = bpert[dit];

    //     Box overlap = bFAB.box();
    //     overlap &= lesBox;
    //     bFAB.setVal(1.0, overlap, 0);
    // }

    // 2.1 Create data holders on one proc covering lesBox.
    // const int srcProc = uniqueProc(SerialTask::compute);
    // const int thisProc = procID();
    LevelData<FArrayBox> lineCartPos;
    LevelData<FArrayBox> lineVel;
    LevelData<FArrayBox> lineBStar;
    LevelData<FArrayBox> lineBBar;
    DataIterator dit1;
    {
        broadcast(lesBox, srcProc);
        Vector<Box> boxArray(1, lesBox);
        Vector<int> procArray(1, srcProc);
        DisjointBoxLayout grids1;
        grids1.define(boxArray, procArray, m_problem_domain);

        dit1 = grids1.dataIterator();
        dit1.reset();
        CH_assert(dit1.ok() || thisProc != srcProc);

        lineCartPos.define(grids1, SpaceDim, BASISV(SpaceDim-1));
        lineVel.define(grids1, SpaceDim, BASISV(SpaceDim-1));
        lineBStar.define(grids1, 1, BASISV(SpaceDim-1));
        lineBBar.define(grids1, 1, BASISV(SpaceDim-1));
    }

    // 2.2 Fill everything.
    Copier cp(newVel().getBoxes(), lineVel.getBoxes(), m_problem_domain, IntVect::Zero, false);
    newVel().copyTo(lineVel, cp);
    newScal(0).copyTo(lineBStar, cp);

    if (thisProc == srcProc) {
        m_levGeoPtr->fill_physCoor(lineCartPos[dit1]);
        m_physBCPtr->setBackgroundScalar(lineBBar, 0, m_time, *m_levGeoPtr);

        ExtrapolateBCNoEV(lineVel[dit1], lesBox, lesBox, 2, BASISV(SpaceDim-1), BASISV(SpaceDim-1));
        ExtrapolateBCNoEV(lineBStar[dit1], lesBox, lesBox, 2, BASISV(SpaceDim-1), BASISV(SpaceDim-1));

        // Remember, we want the _Cartesian_ velocity components!
        m_levGeoPtr->sendToCartesianBasis(lineVel, true);
    }

    // 3. Transmit the data
    if (thisProc == srcProc) {
        FArrayBox& lineCartPosFAB = lineCartPos[dit1];
        FArrayBox& lineVelFAB = lineVel[dit1];
        FArrayBox& lineBStarFAB = lineBStar[dit1];
        FArrayBox& lineBBarFAB = lineBBar[dit1];

        lesBox.grow(SpaceDim-1, 1);

        int numFields = 5;
        int N = lesBox.size(SpaceDim-1);    // Num source points

        int sizeOfTransfer = 1  // time
                           + 1  // number of fields
                           + 1  // number of points
                           + 1  // Vertical Cartesian location at the hill
                           + numFields*N; // data
        double* linearData = new double[sizeOfTransfer];

        int idx = 0;
        linearData[idx++] = double(m_time);
        linearData[idx++] = double(numFields);
        linearData[idx++] = double(N);
        linearData[idx++] = double(zmin);

        IntVect cc = lesBox.smallEnd();
        for (int k = 0; k < N; ++k) {
            linearData[idx++] = lineCartPosFAB(cc,SpaceDim-1);
            ++cc[SpaceDim-1];
        }
        cc = lesBox.smallEnd();
        for (int k = 0; k < N; ++k) {
            linearData[idx++] = lineVelFAB(cc,0);
            ++cc[SpaceDim-1];
        }
        cc = lesBox.smallEnd();
        for (int k = 0; k < N; ++k) {
            linearData[idx++] = lineVelFAB(cc,SpaceDim-1);
            ++cc[SpaceDim-1];
        }
        cc = lesBox.smallEnd();
        for (int k = 0; k < N; ++k) {
            linearData[idx++] = lineBStarFAB(cc);
            ++cc[SpaceDim-1];
        }
        cc = lesBox.smallEnd();
        for (int k = 0; k < N; ++k) {
            linearData[idx++] = lineBBarFAB(cc);
            ++cc[SpaceDim-1];
        }

        { // TEMPORARY!!!
            pout() << "Broadcasted data:\n";
            pout() << "sizeOfTransfer (number of doubles) = " << sizeOfTransfer << endl;
            for (int idx = 0; idx < sizeOfTransfer; ++idx) {
                pout() << linearData[idx] << " ";
            }
            pout() << endl;
        }

        if (AMRLESMeta::lesSize > 0) {
            int ierr = MPI_Bcast(&sizeOfTransfer, 1, MPI_INT,
                                 AMRLESMeta::amr2lesLeader, AMRLESMeta::interComm);
            if (ierr != MPI_SUCCESS) {
                ostringstream msg;
                msg << "AMRNavierStokes::syncWithLES could not broadcast data from AMR to LES"
                    << " -- MPI_Bcast returned an error code = " << ierr;
                MayDay::Error(msg.str().c_str());
            }

            ierr = MPI_Bcast(&linearData[0], sizeOfTransfer, MPI_DOUBLE,
                             AMRLESMeta::amr2lesLeader, AMRLESMeta::interComm);
            if (ierr != MPI_SUCCESS) {
                ostringstream msg;
                msg << "AMRNavierStokes::syncWithLES could not broadcast data from AMR to LES"
                    << " -- MPI_Bcast returned an error code = " << ierr;
                MayDay::Error(msg.str().c_str());
            }
        }

        delete linearData;
        linearData = NULL;
    }
#endif
}
