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

// BC interface for MappedAMRPoissonOp class. This is a rework of Chombo's
// BCFunc.H. Make sure these classes DO NOT own level-dependent data.

#ifndef __ELLIPTICBCINTERFACE_H__INCLUDED__
#define __ELLIPTICBCINTERFACE_H__INCLUDED__

#include "BCDescriptor.H"
#include "RefCountedPtr.H"
#include "Interval.H"
#include "CFRegion.H"
#include "DisjointBoxLayout.H"
class FArrayBox;
class FluxBox;
class DataIndex;
class RealVect;

extern const Real BOGUS_TIME;


// -----------------------------------------------------------------------------
// When creating a function that returns BC values, you should conform to this
// function typedef.
// a_pos array = position of BC (usually at center of cell face)
// a_dir int   = bdry direction, x being 0
// a_side int  = -1 for low, +1 = high
// a_value array = return values. (Usually, just a_value[0] is used.)
// -----------------------------------------------------------------------------
typedef void(*EllipticBCValueFunc) (Real*           a_pos,
                                    int*            a_dir,
                                    Side::LoHiSide* a_side,
                                    Real*           a_value,
                                    const RealVect& a_dx,
                                    Real            a_time);


// -----------------------------------------------------------------------------
// Similar to EllipticBCValueFunc, but this allows special data (such as the
// time, domain length, etc.) to be stored as members.
// -----------------------------------------------------------------------------
class EllipticBCValueClass
{
public:

    virtual ~EllipticBCValueClass ()
    {;}

    virtual void operator() (Real*           a_pos,
                             int*            a_dir,
                             Side::LoHiSide* a_side,
                             Real*           a_value,
                             const RealVect& a_dx,
                             Real            a_time) const = 0;
};


// There are 2 steps to create methods to set BCs. First, you must override
// the interfaces provided in either 1a or 1b. These methods will be called
// "boundary methods" from this point forward. (Not to be confused with "method"
// in the object-oriented sense!) They are what actually do the work when
// filling ghosts or setting fluxes. Once you have written your boundary
// methods, you can attach them to the utility class BCMethodHolder provided
// in step 2. This is the utility that MappedAMRPoissonOpFactory requests at its
// construction.


// -----------------------------------------------------------------------------
// Step #1a: Override these function prototypes to set BCs...
// -----------------------------------------------------------------------------

// Function interface to set BCs by filling ghost cells in MappedAMRPoissonOp.
// Typically used to set Diri BCs.
typedef void (*BCGhostFunc) (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,
                             const DataIndex&     a_index,
                             const FluxBox*       a_JgupPtr,
                             bool                 a_homogeneous,
                             Real                 a_time,
                             const Interval&      a_interval);

// Function interface to set BCs by filling bdry fluxes in MappedAMRPoissonOp.
// Typically used to set Neum BCs. a_state must be FC in a_dir. a_valid is CC.
typedef void (*BCFluxFunc) (FArrayBox&           a_state,
                            const FArrayBox*     a_extrapPtr,
                            const Box&           a_valid,
                            const ProblemDomain& a_domain,
                            const RealVect&      a_dx,
                            const DataIndex&     a_index,
                            const FluxBox*       a_JgupPtr,
                            int                  a_dir,
                            bool                 a_homogeneous,
                            Real                 a_time,
                            const Interval&      a_interval);


// -----------------------------------------------------------------------------
// Step #1b: Override these abstract base classes to set BCs if you need to
//           store additional data (time, domain length, etc.) as members.
// -----------------------------------------------------------------------------

// This serves the same purpose as the BCGhostFunc interface, but it
// allows you to store any additional data you may need (such as the time,
// domain length, etc.) as members.
class BCGhostClass
{
public:

    // Base class constructor. Called by all subclass constructors.
    // Make sure the constructor sets m_bcDescriptor!
    BCGhostClass ()
    {;}

    // Destructor.
    virtual ~BCGhostClass ()
    {;}

    // The pure virtual BCGhostFunc override.
    virtual void operator() (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,
                             const DataIndex&     a_index,
                             const FluxBox*       a_JgupPtr,
                             bool                 a_homogeneous,
                             Real                 a_time = BOGUS_TIME,
                             const Interval&      a_interval = Interval()) const = 0;

    virtual const BCDescriptor& getBCDescriptor() const
    {
        return m_bcDescriptor;
    }

protected:
    BCDescriptor m_bcDescriptor;

private:
    // Copy constructor and assignment operator are disallowed.
    BCGhostClass(const BCGhostClass&);
    BCGhostClass& operator=(const BCGhostClass&);
};


// This serves the same purpose as the BCFluxFunc interface, but it
// allows you to store any additional data you may need (such as the time,
// domain length, etc.) as members.
class BCFluxClass
{
public:

    // Base class constructor. Called by all subclass constructors.
    // Make sure this sets m_bcDescriptor!
    BCFluxClass ()
    {;}

    // Destructor.
    virtual ~BCFluxClass ()
    {;}

    // The pure virtual BCFluxFunc override.
    virtual void operator() (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,
                             const DataIndex&     a_index,
                             const FluxBox*       a_JgupPtr,
                             int                  a_dir,
                             bool                 a_homogeneous,
                             Real                 a_time = BOGUS_TIME,
                             const Interval&      a_interval = Interval()) const = 0;

    virtual const BCDescriptor& getBCDescriptor() const
    {
        return m_bcDescriptor;
    }

protected:
    BCDescriptor m_bcDescriptor;

private:
    // Copy constructor and assignment operator are disallowed.
    BCFluxClass(const BCFluxClass&);
    BCFluxClass& operator=(const BCFluxClass&);
};



// -----------------------------------------------------------------------------
// Step #2: Use this utility as a holder for your boundary methods...
// -----------------------------------------------------------------------------

// This is the main BC utility class that MappedAMRPoissonOpFactory accepts
// in its constructors. Both Diri and Neum BCs can be set using this interface.
// For pure Diri BCs, attach a function that overrides setBCGhostFunc and use
// a do-nothing setBCFluxFunc.
class BCMethodHolder
{
public:

    // Default constructor.
    // This object will not be useable until BCMethods are added.
    BCMethodHolder ();

    // Add a BCGhostFunc method.
    void addBCMethod(BCGhostFunc a_methodPtr);

    // Add a BCFluxFunc method.
    void addBCMethod(BCFluxFunc a_methodPtr);

    // Add a BCGhostClass method.
    void addBCMethod(RefCountedPtr<BCGhostClass>& a_methodPtr);

    // Add a BCFluxClass method.
    void addBCMethod(RefCountedPtr<BCFluxClass>& a_methodPtr);

    // Sets ghost cells
    const BCDescriptor& setGhosts (FArrayBox&           a_state,
                                   const FArrayBox*     a_extrapPtr,
                                   const Box&           a_valid,
                                   const ProblemDomain& a_domain,
                                   const RealVect&      a_dx,
                                   const DataIndex&     a_index,
                                   const FluxBox*       a_JgupPtr,
                                   bool                 a_homogeneous,
                                   Real                 a_time = BOGUS_TIME,
                                   const Interval&      a_interval = Interval()) const;

    // Sets fluxes
    const BCDirDescriptor& setFluxes (FArrayBox&           a_state,
                                      const FArrayBox*     a_extrapPtr,
                                      const Box&           a_valid,
                                      const ProblemDomain& a_domain,
                                      const RealVect&      a_dx,
                                      const DataIndex&     a_index,
                                      const FluxBox*       a_JgupPtr,
                                      int                  a_dir,
                                      bool                 a_homogeneous,
                                      Real                 a_time = BOGUS_TIME,
                                      const Interval&      a_interval = Interval()) const;

    // Sets fluxes in all dirs
    const BCDescriptor& setFluxes (FluxBox&             a_state,
                                   const FArrayBox*     a_extrapPtr,
                                   const Box&           a_valid,
                                   const ProblemDomain& a_domain,
                                   const RealVect&      a_dx,
                                   const DataIndex&     a_index,
                                   const FluxBox*       a_JgupPtr,
                                   bool                 a_homogeneous,
                                   Real                 a_time = BOGUS_TIME,
                                   const Interval&      a_interval = Interval()) const;

    // Returns true if all BCs on a level are periodic or Neumann.
    bool hasNullSpace (const DisjointBoxLayout& a_grids,
                       const CFRegion&          a_cfregion,
                       const IntVect&           a_activeDirs) const;

    // Used to identify what will be done when setGhosts is called.
    inline const BCDescriptor& getGhostDescriptor () const
    {
        return m_ghostDescriptor;
    }

    // Used to identify what will be done when setFluxes is called.
    inline const BCDescriptor& getFluxDescriptor () const
    {
        return m_fluxDescriptor;
    }

    // Checks if enough BCMethods have been added. This does NOT check the
    // validity of the method pointers themselves!
    inline bool isDefined() const
    {
        int numMethods = m_BCGhostFuncPtrs.size() + m_BCFluxFuncPtrs.size()
                       + m_BCGhostObjPtrs.size() + m_BCFluxObjPtrs.size();
        return (numMethods > 0);
    }

private:
    // Returns true if all BCs on a level are periodic or Neumann.
    // This does not perform communication.
    bool hasNullSpaceNoComm (const DisjointBoxLayout& a_grids,
                             const CFRegion&          a_cfregion,
                             const IntVect&           a_activeDirs) const;

    Vector<BCGhostFunc> m_BCGhostFuncPtrs;
    Vector<BCFluxFunc> m_BCFluxFuncPtrs;

    Vector<RefCountedPtr<BCGhostClass> > m_BCGhostObjPtrs;
    Vector<RefCountedPtr<BCFluxClass> > m_BCFluxObjPtrs;

    BCDescriptor m_ghostDescriptor;
    BCDescriptor m_fluxDescriptor;
};


#endif //!__ELLIPTICBCINTERFACE_H__INCLUDED__

