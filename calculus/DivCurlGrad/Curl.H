/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Copyright (C) 2014 Edward Santilli & Alberto Scotti
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
#ifndef __CURL_H__INCLUDED__
#define __CURL_H__INCLUDED__

#include "LevelData.H"
class LevelGeometry;


#if   CH_SPACEDIM == 2
#   define CURL_NCOMP 1
#elif CH_SPACEDIM == 3
#   define CURL_NCOMP 3
#else
#   error Bad SPACEDIM
#endif


// -----------------------------------------------------------------------------
// Utility to compute the curl of a contravariant vector field.
// -----------------------------------------------------------------------------
class Curl
{
public:
    static void simpleCurlCC (FArrayBox&           a_curl,
                              const FArrayBox&     a_contraVel,
                              const DataIndex&     a_di,
                              const LevelGeometry& a_levGeo,
                              const bool           a_convertToVector = false);

    static void levelCurlCC (LevelData<FArrayBox>&       a_curl,
                             const LevelData<FArrayBox>& a_contraVel,
                             const LevelGeometry&        a_levGeo,
                             const bool                  a_convertToVector = false);
};



#endif //!__CURL_H__INCLUDED__
