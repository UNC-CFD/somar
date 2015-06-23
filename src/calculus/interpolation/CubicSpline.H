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
#ifndef CubicSpline_H__INCLUDED__
#define CubicSpline_H__INCLUDED__

#include "Vector.H"
#include "REAL.H"


// -----------------------------------------------------------------------------
// Utility to compute and store spline data.
// -----------------------------------------------------------------------------
class CubicSpline
{
public:
    // Default constructor
    CubicSpline ();

    // Destructor
    ~CubicSpline ();

    // Free memory and make object unusable.
    void clear ();

    // Computes the *natural* spline coefficients from the nodal data, (a_x, a_f).
    void solve (const Vector<Real>& a_f,
                const Vector<Real>& a_x);

    // Fills a_f with interpolated values at positions indicated by a_x.
    void interp (Vector<Real>&       a_f,
                 const Vector<Real>& a_x) const;

    // Returns the interpolated value at the position indicated by a_x.
    Real interp (const Real& a_x) const;
    // Fills a_df with the interpolated first derivatives at positions
    // indicated by a_x.
    // These inputs must have the same centering.
    void interpFirstDeriv (Vector<Real>&       a_df,
                           const Vector<Real>& a_x) const;

    // Fills a_d2f with the interpolated second derivatives at positions
    // indicated by a_x.
    void interpSecondDeriv (Vector<Real>&       a_d2f,
                            const Vector<Real>& a_x) const;

    // Uses a pre-computed set of data.
    void useSolution (const Vector<Real>& a_x,
                      const Vector<Real>& a_f,
                      const Vector<Real>& a_d2f);

protected:
    // The nodal coordinates
    Vector<Real> m_x;

    // The nodal function values
    Vector<Real> m_f;

    // The nodal second derivates needed to create the splines.
    Vector<Real> m_d2f;
};


#endif //!CubicSpline_H__INCLUDED__
