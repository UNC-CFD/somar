/*******************************************************************************
 *    FILE: MappedAMRLevelFactory.H
 *    DESC: Almost identical to Chombo's AMRLevelFactory.
 *    DATE: Wed 14 Aug 2013 07:07:25 PM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __MappedAMRLevelFactory_H__INCLUDED__
#define __MappedAMRLevelFactory_H__INCLUDED__

#include "MappedAMRLevel.H"

/// Factory class to produce a MappedAMRLevel.
/**
   This is needed for construction of MappedAMRLevel-derived classes by
   the AMR class.
*/
class MappedAMRLevelFactory
{
public:
  ///
  /**
      Factory function.  Creates a new MappedAMRLevel.
      In your instantiation you should include defining
      all the stuff that is not part of the MappedAMRLevel API.
  */
  virtual MappedAMRLevel* new_amrlevel() const = 0;

  // to make gcc v4 happy
  virtual ~MappedAMRLevelFactory()
  {
  }
};


#endif //!__MappedAMRLevelFactory_H__INCLUDED__
