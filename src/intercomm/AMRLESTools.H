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

#ifndef __AMRLESTools_H__INCLUDED__
#define __AMRLESTools_H__INCLUDED__

#include "AMRLESMeta.H"
#include "Vector.H"
#include "FArrayBox.H"


namespace AMRLESTools
{
    // Sends data from the AMR to the LES.
    void sendDataToLES (const int* a_data, const int a_count);
    void sendDataToLES (const double* a_data, const int a_count);

    // Allows vectors to be sent to the LES.
    template <typename T> inline void sendDataToLES (const std::vector<T>& a_data);
    template <typename T> inline void sendDataToLES (const Vector<T>& a_data);

    // Sends a FAB and metadata to the LES.
    void sendDataToLES (const FArrayBox& a_dataFAB);
};



// Inline function definitions...

// -----------------------------------------------------------------------------
// Allows std::vectors to be sent to the LES.
// -----------------------------------------------------------------------------
template <typename T>
void AMRLESTools::sendDataToLES (const std::vector<T>& a_data)
{
    sendDataToLES(
        // This pointer will be cast to void* anyway by the MPI_Bcast function,
        // so don't sweat the const_cast.
        const_cast <T*>  (&a_data[0]),
        // This cast is to avoid dealing function resolution issues. AMRLESTools
        // and MPI_Bcast want an int, while vector::size produces a size_t.
        static_cast<int> (a_data.size())
    );
}


// -----------------------------------------------------------------------------
// Allows Chombo Vectors to be sent to the LES.
// -----------------------------------------------------------------------------
template <typename T>
void AMRLESTools::sendDataToLES (const Vector<T>& a_data)
{
    sendDataToLES(
        // This pointer will be cast to void* anyway by the MPI_Bcast function,
        // so don't sweat the const_cast.
        const_cast <T*>  (&a_data[0]),
        // This cast is to avoid dealing function resolution issues. AMRLESTools
        // and MPI_Bcast want an int, while vector::size produces a size_t.
        static_cast<int> (a_data.size())
    );
}


#endif //!__AMRLESTools_H__INCLUDED__
