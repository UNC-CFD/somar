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


// Utility to compute the curl of a contravariant vector field.
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
