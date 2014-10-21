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
#ifndef __TAYLORGREENBCUTIL_H__INCLUDED__
#define __TAYLORGREENBCUTIL_H__INCLUDED__

#include "PhysBCUtil.H"


// -----------------------------------------------------------------------------
// IBCs for the Taylor-Green vortex test.
// -----------------------------------------------------------------------------
class TaylorGreenBCUtil: public PhysBCUtil
{
public:

    // Constructor
    TaylorGreenBCUtil ();

    // Destructor
    virtual ~TaylorGreenBCUtil ();

    // Define constructor
    virtual void define ();

    // Factory
    virtual PhysBCUtil* newPhysBCUtil () const
    {
        TaylorGreenBCUtil* newTGPtr = new TaylorGreenBCUtil;
        // newTGPtr->s_nu = this->s_nu;
        return static_cast<PhysBCUtil*>(newTGPtr);
    }

    // Calculates F(t)
    virtual Real FofT (const Real           a_time,
                       const LevelGeometry& a_levGeo) const;

    // Fills a FAB with the velocity solution
    // Results are in a Cartesian basis at mapped locations.
    virtual void fillVelSoln (FArrayBox&           a_velFAB,
                              const LevelGeometry& a_levGeo,
                              const DataIndex&     a_di,
                              const Real           a_time) const;

    // Fills a FluxBox with the velocity solution
    // Results are in a Cartesian basis at mapped locations.
    virtual void fillVelSoln (FluxBox&             a_velFB,
                              const LevelGeometry& a_levGeo,
                              const DataIndex&     a_di,
                              const Real           a_time) const;

    // Fills a FAB with the pressure solution
    virtual void fillPressureSoln (FArrayBox&           a_pressureFAB,
                                   const LevelGeometry& a_levGeo,
                                   const DataIndex&     a_di,
                                   const Real           a_time) const;

    // Computes the velocity error in an FArrayBox.
    virtual void computeVelError (FArrayBox&           a_errFAB,
                                  const FArrayBox&     a_velFAB,
                                  const LevelGeometry& a_levGeo,
                                  const DataIndex&     a_di,
                                  const Real           a_time,
                                  const bool           a_mapped = false,
                                  const bool           a_multByJ = false) const;

    // Fills a FAB with the initial velocity.
    // Locations are in mapped space, but components are Cartesian.
    virtual void setVelIC (FArrayBox&           a_velFAB,
                           const int            a_velComp,
                           const LevelGeometry& a_levGeo,
                           const DataIndex&     a_di) const;

    // Fills a FAB with the initial scalars
    virtual void setScalarIC (FArrayBox&           a_scalarFAB,
                              const int            a_scalarComp,
                              const LevelGeometry& a_levGeo,
                              const DataIndex&     a_di) const;

protected:
    // Sets physical BCs on velocities
    virtual BCMethodHolder basicVelFuncBC (int a_veldir, bool a_isViscous, Real a_scale) const;


    // Member variables and utility functions...

    // Have we read params from the input file yet?
    static bool s_staticParamsSet;

    static int      s_udir;
    static int      s_vdir;
    static Real     s_kx;
    static Real     s_ky;
    static RealVect s_L;
    static Real     s_nu;
};


// -----------------------------------------------------------------------------
// TaylorGreenBCGhostClass
// u = sin(k*x) * cos(l*y) * F(t)      <--- velDir = 0
// v = -cos(k*x) * sin(l*y) * F(t)     <--- velDir = 1
// w = 0                               <--- velDir = 2
// F(t) = exp(2*nu*t)
// -----------------------------------------------------------------------------
class TaylorGreenBCGhostClass: public BCGhostClass
{
public:
    //
    TaylorGreenBCGhostClass (int a_velDir, Real a_nu, Real a_scale)
    : m_velDir(a_velDir),
      s_nu(a_nu),
      m_scale(a_scale)
    {;}

    //
    virtual ~TaylorGreenBCGhostClass ()
    {;}

    //
    virtual void operator() (FArrayBox&           a_state,
                             const FArrayBox*     a_extrapPtr,  // Just a dummy
                             const Box&           a_valid,
                             const ProblemDomain& a_domain,
                             const RealVect&      a_dx,
                             const DataIndex&     a_index,      // Just a dummy
                             const FluxBox*       a_JgupPtr,    // Just a dummy
                             bool                 a_homogeneous,
                             Real                 a_time = BOGUS_TIME,
                             const Interval&      a_interval = Interval()) const;

protected:
    int  m_velDir;
    Real s_nu;
    Real m_scale;

private:
    // Copy constructor and assignment operator are disallowed.
    TaylorGreenBCGhostClass(const TaylorGreenBCGhostClass&);
    TaylorGreenBCGhostClass& operator=(const TaylorGreenBCGhostClass&);
};


#endif //!__TAYLORGREENBCUTIL_H__INCLUDED__

