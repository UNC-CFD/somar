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
#ifndef __AlteredMetric_H__INCLUDED__
#define __AlteredMetric_H__INCLUDED__

#include "GeoSourceInterface.H"
class PhysBCUtil;


// -----------------------------------------------------------------------------
// This produces the metric used by the implicit BV projector.
// -----------------------------------------------------------------------------
class AlteredMetric: public FillJgupInterface
{
public:
    // Constructor
    AlteredMetric ();

    // Full constructor
    AlteredMetric (const GeoSourceInterface* a_geoSourcePtr,
                   const PhysBCUtil*         a_physBCPtr,
                   const Real                a_dtTheta,
                   const Real                a_coriolisF);

    // Destructor
    virtual ~AlteredMetric ();

    // Full constructor
    virtual void define (const GeoSourceInterface* a_geoSourcePtr,
                         const PhysBCUtil*         a_physBCPtr,
                         const Real                a_dtTheta,
                         const Real                a_coriolisF);

    // Set dest = J*(gup^{ij} - (omega^2 / (1 + omega^2)) * dXi^i/dz * dXi^j/dz)
    // where omega = (dt*N*theta)^2
    virtual void fill_Jgup (FArrayBox&       a_dest,
                            const int        a_destComp,
                            const int        a_mu,
                            const int        a_nu,
                            const RealVect&  a_dXi,
                            const Real       a_scale = 1.0) const;

protected:
    const GeoSourceInterface* m_geoSourcePtr;
    const PhysBCUtil*         m_physBCPtr;
    Real                      m_theta;
    Real                      m_coriolisF;
};



#endif //!__AlteredMetric_H__INCLUDED__
