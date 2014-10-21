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
#ifndef __TWISTEDMAP_H__INCLUDED__
#define __TWISTEDMAP_H__INCLUDED__

#include "GeoSourceInterface.H"


// -----------------------------------------------------------------------------
// The sinusoidally perturbed Cartesian map of JCOMP 2952-2976.
// -----------------------------------------------------------------------------
class TwistedMap: public GeoSourceInterface
{
public:
    // Constructor
    TwistedMap ();

    // Destructor
    virtual ~TwistedMap ();


    // Mandatory overrides...

    // 1. Must return the name of the coordinate mapping
    virtual const char* getCoorMapName () const;

    // 2. Must return whether or not this metric is diagonal
    virtual bool isDiagonal () const;

    // 3. Must return whether or not this metric is uniform
    virtual bool isUniform () const;

    // 4. Fills a mapped box with Cartesian locations.
    virtual void fill_physCoor (FArrayBox&      a_dest,
                                const int       a_destComp,
                                const int       a_mu,
                                const RealVect& a_dXi) const;


    // Optional overrides...

    // Fills a mapped box with Cartesian locations (a_dest must have SpaceDim comps)
    virtual void fill_physCoor (FArrayBox&      a_dest,
                                const RealVect& a_dXi,
                                const RealVect  a_scale = RealVect::Unit) const;

    // Fills an FArrayBox with the Jacobian matrix elements d[x^mu] / d[xi^nu].
    // Major speedup!!!
    virtual void fill_dxdXi (FArrayBox&      a_dest,
                             const int       a_destComp,
                             const int       a_mu,
                             const int       a_nu,
                             const RealVect& a_dXi,
                             const Real      a_scale = 1.0) const;

    // Fills an FArrayBox with J = det[Jacobian]
    // Major speedup!!!
    // Without the analytic version of this function, machine-epsilon sized
    // noise is introduced.
    virtual void fill_J (FArrayBox&      a_dest,
                         const int       a_destComp,
                         const RealVect& a_dXi,
                         const Real      a_scale = 1.0) const;

protected:
    RealVect m_L;
    RealVect m_pert;
    int      m_twistType;
};


#endif //!__TWISTEDMAP_H__INCLUDED__

