// Source: Chombo/example/AMRINS/util

#include "SetValLevel.H"
#include "FluxBox.H"


/// -----------------------------------------------------------
void setValLevels(Vector<LevelData<FArrayBox>*>& a_levels,
                  int                            a_min,
                  int                            a_max,
                  Real                           a_val)
{
  for (int lev = a_min; lev <= a_max; lev++)
    {
      setValLevel(*a_levels[lev], a_val);
    }
}


/// -----------------------------------------------------------
void setValLevels(Vector<BoxLayoutData<FArrayBox>*>& a_levels,
                  int                                a_min,
                  int                                a_max,
                  Real                               a_val)
{
  for (int lev = a_min; lev <= a_max; lev++)
    {
      setValLevel(*a_levels[lev], a_val);
    }
}


/// -----------------------------------------------------------
void setValLevel(BoxLayoutData<FArrayBox>& a_level,
                 Real                      a_val)
{
  for (DataIterator dit = a_level.dataIterator(); dit.ok(); ++dit)
    {
      a_level[dit].setVal(a_val);
    }
}


/// -----------------------------------------------------------
void setValLevel(BoxLayoutData<FluxBox>& a_level,
                 Real                    a_val)
{
  for (DataIterator dit = a_level.dataIterator(); dit.ok(); ++dit)
    {
      a_level[dit].setVal(a_val);
    }
}
