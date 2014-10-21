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
#ifndef __DataIndexCompare_H__INCLUDED__
#define __DataIndexCompare_H__INCLUDED__

#include "DataIndex.H"


// Used to lexicographically compare two DataIndices. This is needed
// when we use a DataIndex as the key in a std::map.
struct DataIndexCompare {
    bool operator() (const DataIndex& a_di1, const DataIndex& a_di2) const {
        return a_di1.datInd() < a_di2.datInd();
    }
};


#endif //!__DataIndexCompare_H__INCLUDED__
