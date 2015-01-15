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
#include "AMRLESTools.H"
#include "Constants.H"
#include "MayDay.H"
#include <sstream>
using std::ostringstream;


// -----------------------------------------------------------------------------
// Sends data from the AMR to the LES.
// -----------------------------------------------------------------------------
void AMRLESTools::sendDataToLES (const int* a_data, const int a_count)
{
    // Are we using an LES?
    if (AMRLESMeta::lesSize == 0) return;

    // Send the data.
    int ierr = MPI_Bcast((void*)a_data,             // data array pointer
                         a_count,                   // num array elements
                         MPI_INT,                   // MPI_Datatype
                         AMRLESMeta::amr2lesLeader, // root rank
                         AMRLESMeta::interComm);    // MPI_Comm

    // If there was an error, throw a message and halt.
    if (ierr != MPI_SUCCESS) {
        ostringstream msg;
        msg << "Could not broadcast " << a_count << " ints from AMR to LES"
            << " -- MPI_Bcast returned an error code = " << ierr;
        MayDay::Error(msg.str().c_str());
    }
}


// -----------------------------------------------------------------------------
// Sends data from the AMR to the LES.
// -----------------------------------------------------------------------------
void AMRLESTools::sendDataToLES (const double* a_data, const int a_count)
{
    // Are we using an LES?
    if (AMRLESMeta::lesSize == 0) return;

    // Send the data.
    int ierr = MPI_Bcast((void*)a_data,             // data array pointer
                         a_count,                   // num array elements
                         MPI_DOUBLE,                // MPI_Datatype
                         AMRLESMeta::amr2lesLeader, // root rank
                         AMRLESMeta::interComm);    // MPI_Comm

    // If there was an error, throw a message and halt.
    if (ierr != MPI_SUCCESS) {
        ostringstream msg;
        msg << "Could not broadcast " << a_count << " doubles from AMR to LES"
            << " -- MPI_Bcast returned an error code = " << ierr;
        MayDay::Error(msg.str().c_str());
    }
}


// -----------------------------------------------------------------------------
// Sends a FAB and metadata to the LES.
// -----------------------------------------------------------------------------
void AMRLESTools::sendDataToLES (const FArrayBox& a_dataFAB)
{
    // Collect the metadata.
    // We will send this first so that storage can be prepared.
    Vector<int> metadata(10);

    const IntVect& smallEnd = a_dataFAB.smallEnd();
    metadata[0] = smallEnd[0];
    metadata[1] = smallEnd[1];
    metadata[2] = ((SpaceDim > 2)? smallEnd[2]: 0);

    const IntVect& bigEnd = a_dataFAB.bigEnd();
    metadata[3] = bigEnd[0];
    metadata[4] = bigEnd[1];
    metadata[5] = ((SpaceDim > 2)? bigEnd[2]: 0);

    metadata[6] = a_dataFAB.nComp();

    const IntVect& boxType = a_dataFAB.box().type();
    metadata[7] = boxType[0];
    metadata[8] = boxType[1];
    metadata[9] = ((SpaceDim > 2)? boxType[2]: -1);

    // Now send everything.
    sendDataToLES(metadata);
    sendDataToLES(a_dataFAB.dataPtr(), a_dataFAB.box().numPts());
}
