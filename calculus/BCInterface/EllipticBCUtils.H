/*******************************************************************************
 *	FILE: EllipticBCUtils.hpp
 *	DESC: Conveinience functions and classes that handle BCs.
 *	DATE: Thu 24 May 2012 01:01:52 PM EDT
 *	MAIL: santilli@physics.unc.edu
 *  NOTE: This is a rework of Chombo's BCFunc.H.
 *        Make sure these classes DO NOT own level-dependent data.
 *        These functions have the unfortunate "elliptic" prefix to avoid a
 *        clash with BCFunc.H.
 ******************************************************************************/
#ifndef __ELLIPTICBCUTILS_HPP__INCLUDED__
#define __ELLIPTICBCUTILS_HPP__INCLUDED__

#include "EllipticBCInterface.H"
#include "Tuple.H"
#include "RealVect.H"
#include "LevelData.H"
#include "BoundaryData.H"


// ***************************** Single-side BCs *******************************

// -----------------------------------------------------------------------------
// Sets Dirichlet BCs on a side
// For a CC state, this fills ghosts.
// For a FC state, this sets faces at the physical boundary.
// a_state and a_valid must have the same centering.
// -----------------------------------------------------------------------------
void setSideDiriBC (FArrayBox&           a_state,
                    const Box&           a_valid,
                    const ProblemDomain& a_domain,
                    const Real           a_value,
                    const int            a_dir,
                    const Side::LoHiSide a_side,
                    const bool           a_homogeneous,
                    const int            a_order);


// -----------------------------------------------------------------------------
// Sets the outward pointing normal covariant derivative of a_state on a side.
// For a CC state, this fills ghosts.
// For a FC state, this throws an error.
// a_state and a_valid must have the same centering.
// -----------------------------------------------------------------------------
void setSideNeumBC (FArrayBox&           a_state,
                    const Box&           a_valid,
                    const ProblemDomain& a_domain,
                    const Real           a_value,
                    const int            a_dir,
                    const Side::LoHiSide a_side,
                    const bool           a_homogeneous,
                    const FArrayBox&     a_Jgupi,
                    const FArrayBox*     a_extrapPtr,
                    const RealVect&      a_dx,
                    const DataIndex&     a_index);


// -----------------------------------------------------------------------------
// Sets extrapolation BCs on a side
// For a CC state, this fills ghosts.
// For a FC state, this sets faces at the physical boundary.
// a_state and a_valid must have the same centering.
// This function does NOT check for periodicity!
// -----------------------------------------------------------------------------
void setSideExtrapBC (FArrayBox&           a_state,
                      const Box&           a_valid,
                      const ProblemDomain& a_domain,
                      const int            a_dir,
                      const Side::LoHiSide a_side,
                      const int            a_order);


// -----------------------------------------------------------------------------
// Sets extrapolation BCs on a side.
// This version extrapolates CC data to faces and does not work in place.
// a_dest must be FC in a_dir.
// a_src and a_valid must be CC a_dir.
// a_dest, a_src, and a_valid must have the same centering in other directions.
// This function does NOT check for periodicity!
// -----------------------------------------------------------------------------
void setSideExtrapBC (FArrayBox&           a_dest,
                      const FArrayBox&     a_src,
                      const Box&           a_valid,
                      const ProblemDomain& a_domain,
                      const int            a_dir,
                      const Side::LoHiSide a_side,
                      const int            a_order);


// *************************** Complete BC methods *****************************

// -----------------------------------------------------------------------------
// A simple override of BCGhostClass to set constant Diri BCs.
// For CC dirs, fills ghosts. For FC dirs, sets bdry data directly.
// -----------------------------------------------------------------------------
class EllipticConstDiriBCGhostClass: public BCGhostClass
{
public:
    EllipticConstDiriBCGhostClass (const RealVect& a_loSideValue,
                                   const RealVect& a_hiSideValue,
                                   const IntVect   a_doLoSides = IntVect::Unit,
                                   const IntVect   a_doHiSides = IntVect::Unit)
    {
        m_bcDescriptor.set(BCType::Diri);

        m_BCValue[0] = a_loSideValue;
        m_BCValue[1] = a_hiSideValue;

        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    virtual ~EllipticConstDiriBCGhostClass ()
    {;}

    virtual void operator() (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,  // Just a dummy
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,         // Just a dummy
                             const DataIndex&     a_index,      // Just a dummy
                             const FluxBox*       a_JgupPtr,    // Just a dummy
                             bool                 a_homogeneous,
                             Real                 a_time = BOGUS_TIME,
                             const Interval&      a_interval = Interval()) const;

protected:
    Tuple<RealVect, 2> m_BCValue;
    Tuple<IntVect, 2>  m_doSide;

private:
    // Copy constructor and assignment operator are disallowed.
    EllipticConstDiriBCGhostClass(const EllipticConstDiriBCGhostClass&);
    EllipticConstDiriBCGhostClass& operator=(const EllipticConstDiriBCGhostClass&);
};


// -----------------------------------------------------------------------------
// A simple override of BCGhostClass that fills ghosts using constant Neum BCs.
// NOTE: It only makes sense to use this on CC data!
// -----------------------------------------------------------------------------
class EllipticConstNeumBCGhostClass: public BCGhostClass
{
public:
    EllipticConstNeumBCGhostClass (const RealVect& a_loSideValue,
                                   const RealVect& a_hiSideValue,
                                   const IntVect   a_doLoSides = IntVect::Unit,
                                   const IntVect   a_doHiSides = IntVect::Unit)
    {
        m_bcDescriptor.set(BCType::Neum);

        m_BCValue[0] = a_loSideValue;
        m_BCValue[1] = a_hiSideValue;

        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    virtual ~EllipticConstNeumBCGhostClass ()
    {;}

    virtual void operator() (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,
                             const DataIndex&     a_index,      // Just a dummy
                             const FluxBox*       a_JgupPtr,
                             bool                 a_homogeneous,
                             Real                 a_time = BOGUS_TIME,
                             const Interval&      a_interval = Interval()) const;

protected:
    Tuple<RealVect, 2> m_BCValue;
    Tuple<IntVect, 2>  m_doSide;

private:
    // Copy constructor and assignment operator are disallowed.
    EllipticConstNeumBCGhostClass(const EllipticConstNeumBCGhostClass&);
    EllipticConstNeumBCGhostClass& operator=(const EllipticConstNeumBCGhostClass&);
};


// -----------------------------------------------------------------------------
// A simple override of BCFluxClass that sets fluxes using constant Neum BCs.
// This is used to set BCs after taking the gradient of a field.
//
// NOTE: This function was not really made to apply zero-normal BCs to a
// velocity. For that, consider using a zero EllipticConstDiriBCGhostClass.
//
// NOTE: This can only set face values, not ghosts!
// -----------------------------------------------------------------------------
class EllipticConstNeumBCFluxClass: public BCFluxClass
{
public:
    EllipticConstNeumBCFluxClass (const RealVect& a_loSideValue,
                                  const RealVect& a_hiSideValue,
                                  const IntVect   a_doLoSides = IntVect::Unit,
                                  const IntVect   a_doHiSides = IntVect::Unit)
    {
        m_bcDescriptor.set(BCType::Neum);

        m_BCValue[0] = a_loSideValue;
        m_BCValue[1] = a_hiSideValue;

        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    virtual ~EllipticConstNeumBCFluxClass ()
    {;}

    virtual void operator() (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,  // Just a dummy
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,         // Just a dummy
                             const DataIndex&     a_index,      // Just a dummy
                             const FluxBox*       a_Jgup,       // Just a dummy
                             int                  a_dir,
                             bool                 a_homogeneous,
                             Real                 a_time = BOGUS_TIME,
                             const Interval&      a_interval = Interval()) const;

protected:
    Tuple<RealVect, 2> m_BCValue;
    Tuple<IntVect, 2>  m_doSide;

private:
    // Copy constructor and assignment operator are disallowed.
    EllipticConstNeumBCFluxClass(const EllipticConstNeumBCFluxClass&);
    EllipticConstNeumBCFluxClass& operator=(const EllipticConstNeumBCFluxClass&);
};


// -----------------------------------------------------------------------------
// A simple override of BCGhostClass that fills ghosts using general Diri BCs.
// -----------------------------------------------------------------------------
class EllipticDiriBCGhostClass: public BCGhostClass
{
public:
    EllipticDiriBCGhostClass (EllipticBCValueFunc a_funcPtr,
                              const IntVect       a_doLoSides = IntVect::Unit,
                              const IntVect       a_doHiSides = IntVect::Unit)
    : m_funcPtr(a_funcPtr),
      m_objPtr(RefCountedPtr<EllipticBCValueClass>(NULL))
    {
        m_bcDescriptor.set(BCType::Diri);

        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    EllipticDiriBCGhostClass (RefCountedPtr<EllipticBCValueClass> a_objPtr,
                              const IntVect                       a_doLoSides = IntVect::Unit,
                              const IntVect                       a_doHiSides = IntVect::Unit)
    : m_funcPtr(NULL),
      m_objPtr(a_objPtr)
    {
        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    virtual ~EllipticDiriBCGhostClass ()
    {;}

    virtual void operator() (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,      // Just a dummy
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,
                             const DataIndex&     a_index,          // Just a dummy
                             const FluxBox*       a_JgupPtr,        // Just a dummy
                             bool                 a_homogeneous,
                             Real                 a_time = BOGUS_TIME,
                             const Interval&      a_interval = Interval()) const;

protected:
    EllipticBCValueFunc                 m_funcPtr;
    RefCountedPtr<EllipticBCValueClass> m_objPtr;
    Tuple<IntVect, 2>                   m_doSide;

private:
    // Copy constructor and assignment operator are disallowed.
    EllipticDiriBCGhostClass(const EllipticDiriBCGhostClass&);
    EllipticDiriBCGhostClass& operator=(const EllipticDiriBCGhostClass&);
};


// -----------------------------------------------------------------------------
// A simple override of BCGhostClass that fills ghosts using general Neum BCs.
// -----------------------------------------------------------------------------
class EllipticNeumBCGhostClass: public BCGhostClass
{
public:
    EllipticNeumBCGhostClass (EllipticBCValueFunc a_funcPtr,
                              const IntVect       a_doLoSides = IntVect::Unit,
                              const IntVect       a_doHiSides = IntVect::Unit)
    : m_funcPtr(a_funcPtr),
      m_objPtr(RefCountedPtr<EllipticBCValueClass>(NULL)),
      m_bdryDataPtr(NULL)
    {
        m_bcDescriptor.set(BCType::Neum);

        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    EllipticNeumBCGhostClass (RefCountedPtr<EllipticBCValueClass> a_objPtr,
                              const IntVect                       a_doLoSides = IntVect::Unit,
                              const IntVect                       a_doHiSides = IntVect::Unit)
    : m_funcPtr(NULL),
      m_objPtr(a_objPtr),
      m_bdryDataPtr(NULL)
    {
        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    EllipticNeumBCGhostClass (RefCountedPtr<BoundaryData<Real> > a_bdryDataPtr,
                              const IntVect                      a_doLoSides = IntVect::Unit,
                              const IntVect                      a_doHiSides = IntVect::Unit)
    : m_funcPtr(NULL),
      m_objPtr(RefCountedPtr<EllipticBCValueClass>(NULL)),
      m_bdryDataPtr(a_bdryDataPtr)
    {
        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    virtual ~EllipticNeumBCGhostClass ()
    {;}

    virtual void operator() (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,
                             const DataIndex&     a_index,
                             const FluxBox*       a_JgupPtr,
                             bool                 a_homogeneous,
                             Real                 a_time = BOGUS_TIME,
                             const Interval&      a_interval = Interval()) const;

protected:
    EllipticBCValueFunc                 m_funcPtr;
    RefCountedPtr<EllipticBCValueClass> m_objPtr;
    RefCountedPtr<BoundaryData<Real> >  m_bdryDataPtr;
    Tuple<IntVect, 2>                   m_doSide;

private:
    // Copy constructor and assignment operator are disallowed.
    EllipticNeumBCGhostClass(const EllipticNeumBCGhostClass&);
    EllipticNeumBCGhostClass& operator=(const EllipticNeumBCGhostClass&);
};


// -----------------------------------------------------------------------------
// A simple override of BCFluxClass that sets fluxes using general Neum BCs.
// -----------------------------------------------------------------------------
class EllipticNeumBCFluxClass: public BCFluxClass
{
public:
    EllipticNeumBCFluxClass (EllipticBCValueFunc a_funcPtr,
                             const IntVect       a_doLoSides = IntVect::Unit,
                             const IntVect       a_doHiSides = IntVect::Unit)
    : m_funcPtr(a_funcPtr),
      m_objPtr(RefCountedPtr<EllipticBCValueClass>(NULL)),
      m_fluxPtr(NULL),
      m_bdryDataPtr(NULL)
    {
        m_bcDescriptor.set(BCType::Neum);

        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    EllipticNeumBCFluxClass (RefCountedPtr<EllipticBCValueClass> a_objPtr,
                             const IntVect                       a_doLoSides = IntVect::Unit,
                             const IntVect                       a_doHiSides = IntVect::Unit)
    : m_funcPtr(NULL),
      m_objPtr(a_objPtr),
      m_fluxPtr(NULL),
      m_bdryDataPtr(NULL)
    {
        // TODO: set the BC descriptor???
        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    EllipticNeumBCFluxClass (RefCountedPtr<LevelData<FluxBox> > a_fieldPtr,
                             const IntVect                      a_doLoSides = IntVect::Unit,
                             const IntVect                      a_doHiSides = IntVect::Unit)
    : m_funcPtr(NULL),
      m_objPtr(RefCountedPtr<EllipticBCValueClass>(NULL)),
      m_fluxPtr(a_fieldPtr),
      m_bdryDataPtr(NULL)
    {
        CH_assert(!m_fluxPtr.isNull());
        CH_assert(m_fluxPtr->nComp() == 1);

        m_bcDescriptor.set(BCType::Neum);

        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    EllipticNeumBCFluxClass (RefCountedPtr<BoundaryData<Real> > a_bdryDataPtr,
                             const IntVect                      a_doLoSides = IntVect::Unit,
                             const IntVect                      a_doHiSides = IntVect::Unit)
    : m_funcPtr(NULL),
      m_objPtr(RefCountedPtr<EllipticBCValueClass>(NULL)),
      m_fluxPtr(NULL),
      m_bdryDataPtr(a_bdryDataPtr)
    {
        m_bcDescriptor.set(BCType::Neum);

        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    virtual ~EllipticNeumBCFluxClass ()
    {;}

    virtual void operator() (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,      // Just a dummy
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,
                             const DataIndex&     a_index,          // Just a dummy
                             const FluxBox*       a_JgupPtr,        // Just a dummy
                             int                  a_dir,
                             bool                 a_homogeneous,
                             Real                 a_time = BOGUS_TIME,
                             const Interval&      a_interval = Interval()) const;

protected:
    EllipticBCValueFunc                 m_funcPtr;
    RefCountedPtr<EllipticBCValueClass> m_objPtr;
    RefCountedPtr<LevelData<FluxBox> >  m_fluxPtr;
    RefCountedPtr<BoundaryData<Real> >  m_bdryDataPtr;
    Tuple<IntVect, 2>                   m_doSide;

private:
    // Copy constructor and assignment operator are disallowed.
    EllipticNeumBCFluxClass(const EllipticNeumBCFluxClass&);
    EllipticNeumBCFluxClass& operator=(const EllipticNeumBCFluxClass&);
};


// -----------------------------------------------------------------------------
// A simple override of BCFluxClass that sets fluxes to zero if a_homogeneous.
// Does nothing if !a_homogeneous.
//
// NOTE: This function was designed to be used in an iterative solver to
// eliminate errors introduced by taking the gradient over extrapolated ghosts.
// -----------------------------------------------------------------------------
class EllipticDoNothingUnlessHomogNeumBCFluxClass: public BCFluxClass
{
public:
    EllipticDoNothingUnlessHomogNeumBCFluxClass (const IntVect a_doLoSides = IntVect::Unit,
                                                 const IntVect a_doHiSides = IntVect::Unit)
    {
        m_bcDescriptor.set(BCType::Neum);

        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    virtual ~EllipticDoNothingUnlessHomogNeumBCFluxClass ()
    {;}

    virtual void operator() (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,      // Just a dummy
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,             // Just a dummy
                             const DataIndex&     a_index,          // Just a dummy
                             const FluxBox*       a_JgupPtr,        // Just a dummy
                             int                  a_dir,
                             bool                 a_homogeneous,    // Just a dummy
                             Real                 a_time = BOGUS_TIME,
                             const Interval&      a_interval = Interval()) const;

protected:
    Tuple<IntVect, 2> m_doSide;

private:
    // Copy constructor and assignment operator are disallowed.
    EllipticDoNothingUnlessHomogNeumBCFluxClass(const EllipticDoNothingUnlessHomogNeumBCFluxClass&);
    EllipticDoNothingUnlessHomogNeumBCFluxClass& operator=(const EllipticDoNothingUnlessHomogNeumBCFluxClass&);
};


// -----------------------------------------------------------------------------
// A simple override of BCGhostClass that fills ghosts using extrapolation BCs.
// -----------------------------------------------------------------------------
class EllipticExtrapBCGhostClass: public BCGhostClass
{
public:
    EllipticExtrapBCGhostClass (int           a_order,
                                const IntVect a_doLoSides = IntVect::Unit,
                                const IntVect a_doHiSides = IntVect::Unit)
    : m_order(a_order)
    {
        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    virtual ~EllipticExtrapBCGhostClass ()
    {;}

    virtual void operator() (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,      // Just a dummy
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,             // Just a dummy
                             const DataIndex&     a_index,          // Just a dummy
                             const FluxBox*       a_JgupPtr,        // Just a dummy
                             bool                 a_homogeneous,    // Just a dummy
                             Real                 a_time = BOGUS_TIME,
                             const Interval&      a_interval = Interval()) const;
protected:
    int               m_order;
    Tuple<IntVect, 2> m_doSide;

private:
    // Copy constructor and assignment operator are disallowed.
    EllipticExtrapBCGhostClass(const EllipticExtrapBCGhostClass&);
    EllipticExtrapBCGhostClass& operator=(const EllipticExtrapBCGhostClass&);
};


// -----------------------------------------------------------------------------
// A simple override of BCFluxClass that sets fluxes using extrapolation BCs.
// -----------------------------------------------------------------------------
class EllipticExtrapBCFluxClass: public BCFluxClass
{
public:
    EllipticExtrapBCFluxClass (int           a_order,
                               const IntVect a_doLoSides = IntVect::Unit,
                               const IntVect a_doHiSides = IntVect::Unit,
                               const bool    a_useExtrapPtrAsSource = false)
    : m_order(a_order),
      m_useExtrapPtrAsSource(a_useExtrapPtrAsSource)
    {
        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    virtual ~EllipticExtrapBCFluxClass ()
    {;}

    virtual void operator() (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,  // Can be a CC state!
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,         // Just a dummy
                             const DataIndex&     a_index,      // Just a dummy
                             const FluxBox*       a_Jgup,       // Just a dummy
                             int                  a_dir,
                             bool                 a_homogeneous,
                             Real                 a_time = BOGUS_TIME,
                             const Interval&      a_interval = Interval()) const;

protected:
    int m_order;
    Tuple<IntVect, 2> m_doSide;
    bool m_useExtrapPtrAsSource;

private:
    // Copy constructor and assignment operator are disallowed.
    EllipticExtrapBCFluxClass(const EllipticExtrapBCFluxClass&);
    EllipticExtrapBCFluxClass& operator=(const EllipticExtrapBCFluxClass&);
};


// -----------------------------------------------------------------------------
// A simple override of BCGhostClass that fills ghosts using extrapolation
// in m_normDir and zero Diri BCs in all others.
// -----------------------------------------------------------------------------
class EllipticInviscidBCGhostClass: public BCGhostClass
{
public:
    EllipticInviscidBCGhostClass (int           a_normDir   = -1,
                                  const IntVect a_doLoSides = IntVect::Unit,
                                  const IntVect a_doHiSides = IntVect::Unit)
    : m_normDir(a_normDir)
    {
        CH_assert(0 <= m_normDir);
        CH_assert(m_normDir < SpaceDim);

        D_TERM(m_bcDescriptor.set(m_normDir, Side::Lo, BCType::Diri);
               m_bcDescriptor.set(m_normDir, Side::Hi, BCType::Diri);,
               m_bcDescriptor.set((m_normDir + 1) % SpaceDim, Side::Lo, BCType::Undefined);
               m_bcDescriptor.set((m_normDir + 1) % SpaceDim, Side::Hi, BCType::Undefined);,
               m_bcDescriptor.set((m_normDir + 2) % SpaceDim, Side::Lo, BCType::Undefined);
               m_bcDescriptor.set((m_normDir + 2) % SpaceDim, Side::Hi, BCType::Undefined);)

        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    virtual ~EllipticInviscidBCGhostClass ()
    {;}

    virtual void operator() (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,  // Just a dummy
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,         // Just a dummy
                             const DataIndex&     a_index,      // Just a dummy
                             const FluxBox*       a_JgupPtr,    // Just a dummy
                             bool                 a_homogeneous,
                             Real                 a_time = BOGUS_TIME,
                             const Interval&      a_interval = Interval()) const;

protected:
    int               m_normDir;
    Tuple<IntVect, 2> m_doSide;

private:
    // Copy constructor and assignment operator are disallowed.
    EllipticInviscidBCGhostClass(const EllipticInviscidBCGhostClass&);
    EllipticInviscidBCGhostClass& operator=(const EllipticInviscidBCGhostClass&);
};


// -----------------------------------------------------------------------------
// An override of BCGhostClass that sets inflow/outflow BCs on the velocity.
// -----------------------------------------------------------------------------
class BasicVelocityBCGhostClass: public BCGhostClass
{
public:
    typedef Real(*InflowVelFuncPtr)(Real);

    // Constructor -- constant inflow vel
    BasicVelocityBCGhostClass (const Real           a_inflowVel,
                               const int            a_inflowDir,
                               const Side::LoHiSide a_inflowSide,
                               const int            a_outflowDir,
                               const Side::LoHiSide a_outflowSide,
                               const int            a_velComp,
                               const bool           a_isViscous,
                               const IntVect        a_doLoSides = IntVect::Unit,
                               const IntVect        a_doHiSides = IntVect::Unit)
    : m_inflowVelFuncPtr(NULL)
    {
        this->define(a_inflowVel,
                     a_inflowDir,
                     a_inflowSide,
                     a_outflowDir,
                     a_outflowSide,
                     a_velComp,
                     a_isViscous,
                     a_doLoSides,
                     a_doHiSides);
    }

    // Constructor -- time varying inflow vel.
    BasicVelocityBCGhostClass (InflowVelFuncPtr     a_inflowVelFuncPtr,
                               const Real           a_inflowVelScale,
                               const int            a_inflowDir,
                               const Side::LoHiSide a_inflowSide,
                               const int            a_outflowDir,
                               const Side::LoHiSide a_outflowSide,
                               const int            a_velComp,
                               const bool           a_isViscous,
                               const IntVect        a_doLoSides = IntVect::Unit,
                               const IntVect        a_doHiSides = IntVect::Unit)
    : m_inflowVelFuncPtr(a_inflowVelFuncPtr),
      m_inflowVelScale(a_inflowVelScale)
    {
        CH_assert(m_inflowVelFuncPtr != NULL);

        this->define(1.0e300,
                     a_inflowDir,
                     a_inflowSide,
                     a_outflowDir,
                     a_outflowSide,
                     a_velComp,
                     a_isViscous,
                     a_doLoSides,
                     a_doHiSides);
    }

    virtual ~BasicVelocityBCGhostClass ()
    {;}

    // Define constructor
    virtual void define (const Real           a_inflowVel,
                         const int            a_inflowDir,
                         const Side::LoHiSide a_inflowSide,
                         const int            a_outflowDir,
                         const Side::LoHiSide a_outflowSide,
                         const int            a_velComp,
                         const bool           a_isViscous,
                         const IntVect        a_doLoSides = IntVect::Unit,
                         const IntVect        a_doHiSides = IntVect::Unit)
    {
        // Set the BC descriptors
        // First, set the descriptors ignoring inflow/outflow.
        if (a_isViscous) {
            m_bcDescriptor.set(BCType::Diri);
        } else {
            m_bcDescriptor.set(BCType::Undefined);
            m_bcDescriptor.set(a_velComp, Side::Lo, BCType::Diri);
            m_bcDescriptor.set(a_velComp, Side::Hi, BCType::Diri);
        }
        // Then, correct the inflow/outflow sides if they exist.
        if (0 <= a_outflowDir && a_outflowDir < SpaceDim) {
            m_bcDescriptor.set(a_outflowDir, a_outflowSide, BCType::Undefined); // Should be Neum
        }

        // Sanity checks
        CH_assert(0 <= a_velComp);
        CH_assert(a_velComp < SpaceDim);

        // Set the member variables
        m_inflowVel = a_inflowVel;
        m_inflowDir = a_inflowDir;
        m_inflowSide = a_inflowSide;
        m_outflowDir = a_outflowDir;
        m_outflowSide = a_outflowSide;

        m_velComp = a_velComp;
        m_isViscous = a_isViscous;

        m_doSide[0] = a_doLoSides;
        m_doSide[1] = a_doHiSides;
    }

    virtual void operator() (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,  // Just a dummy
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,         // Just a dummy
                             const DataIndex&     a_index,      // Just a dummy
                             const FluxBox*       a_JgupPtr,    // Just a dummy
                             bool                 a_homogeneous,
                             Real                 a_time = BOGUS_TIME,
                             const Interval&      a_interval = Interval()) const;

protected:
    InflowVelFuncPtr  m_inflowVelFuncPtr;
    Real              m_inflowVelScale;
    Real              m_inflowVel;
    int               m_inflowDir;
    Side::LoHiSide    m_inflowSide;
    int               m_outflowDir;
    Side::LoHiSide    m_outflowSide;
    int               m_velComp;
    bool              m_isViscous;
    Tuple<IntVect, 2> m_doSide;

private:
    // Copy constructor and assignment operator are disallowed.
    BasicVelocityBCGhostClass(const BasicVelocityBCGhostClass&);
    BasicVelocityBCGhostClass& operator=(const BasicVelocityBCGhostClass&);
};


#endif //!__ELLIPTICBCUTILS_HPP__INCLUDED__
