/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Copyright (C) 2014 Edward Santilli & Alberto Scotti
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
#ifndef __StressMetric_H__INCLUDED__
#define __StressMetric_H__INCLUDED__


// To use a custom stress tensor, we override Jg^{i,j} in the viscous solver's
// elliptic operators. By defining USE_STRESSMETRIC, all of the needed code
// will be produced and you can rewrite the fill_Jgup function to include any
// subgrid model you like. If your subgrid model needs the state variables, the
// define function will need to be altered as well.
//
// If USE_STRESSMETRIC is undefined, then the StressMetric class will not exist
// anywhere in the (preprocessed) code. This will save some memory since the
// viscous solver's elliptic operators will use the Jg^{i,j} that has been
// cached by LevelGeometry.
//
// #define USE_STRESSMETRIC
#ifdef USE_STRESSMETRIC


#include "GeoSourceInterface.H"

// -----------------------------------------------------------------------------
// Creates J*T^{i j} for viscous solves.
// -----------------------------------------------------------------------------
class StressMetric: public FillJgupInterface
{
public:
    // Constructor
    StressMetric ();

    // Full constructor
    StressMetric (const GeoSourceInterface* a_geoSourcePtr);

    // Destructor
    virtual ~StressMetric ();

    // Full constructor
    virtual void define (const GeoSourceInterface* a_geoSourcePtr);

    // Is this defined?
    virtual inline bool isDefined () const;

    // Set dest = (1+turbVisc/dynVisc)*J*gup^{ij}
    // Remember: The solvers will multiply by dynVisc, that's why it's factored out.
    virtual void fill_Jgup (FArrayBox&       a_dest,
                            const int        a_destComp,
                            const int        a_mu,
                            const int        a_nu,
                            const RealVect&  a_dXi,
                            const Real       a_scale = 1.0) const;

protected:
    const GeoSourceInterface* m_geoSourcePtr;
    Real                      m_dynVisc;
};



// -----------------------------------------------------------------------------
// Is this defined?
// -----------------------------------------------------------------------------
bool StressMetric::isDefined () const
{
    return (m_geoSourcePtr != NULL);
}


#endif //USE_STRESSMETRIC
#endif //!__StressMetric_H__INCLUDED__
