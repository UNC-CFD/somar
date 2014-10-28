// This file was copied from Chombo/example/AMRINS/util.
// Chombo's copyright applies.
#ifndef _SETVALLEVEL_H_
#define _SETVALLEVEL_H_

#include "LevelData.H"
class FArrayBox;
class FluxBox;


// Set a_levels[a_min:a_max] to zero.
void setValLevels(Vector<LevelData<FArrayBox>*>& a_levels,
                  int                            a_min,
                  int                            a_max,
                  Real                           a_val);

// Set a_levels[a_min:a_max] to zero.
void setValLevels(Vector<BoxLayoutData<FArrayBox>*>& a_levels,
                  int                                a_min,
                  int                                a_max,
                  Real                               a_val);

// Set a_level to zero.
void setValLevel(BoxLayoutData<FArrayBox>& a_level,
                 Real                      a_val);


// Set a_level to zero.
void setValLevel(BoxLayoutData<FluxBox>& a_level,
                 Real                    a_val);


#endif
