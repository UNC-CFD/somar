/*******************************************************************************
 *    FILE: DataIndexCompare.H
 *    DESC: Allows a DataIndex to be used as the key in a std::map.
 *    DATE: Sun 10 Nov 2013 11:38:31 AM EST
 *    MAIL: santilli@physics.unc.edu
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
