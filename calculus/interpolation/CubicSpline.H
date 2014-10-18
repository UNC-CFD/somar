/*******************************************************************************
 *    FILE: CubicSpline.H
 *    DESC: Utility to compute and store spline data.
 *    DATE: Sat 24 May 2014 09:59:42 AM EDT  (UTC-0400)
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef CubicSpline_H__INCLUDED__
#define CubicSpline_H__INCLUDED__

#include "Vector.H"
#include "REAL.H"


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
