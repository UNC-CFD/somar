#include "CubicSpline.H"
#include "CubicSplineF_F.H"


// -----------------------------------------------------------------------------
// Default constructor
// -----------------------------------------------------------------------------
CubicSpline::CubicSpline ()
{;}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
CubicSpline::~CubicSpline ()
{
    this->clear();
}


// -----------------------------------------------------------------------------
// Free memory and make object unusable.
// -----------------------------------------------------------------------------
void CubicSpline::clear ()
{
    m_x.clear();
    m_f.clear();
    m_d2f.clear();
}


// -----------------------------------------------------------------------------
// Computes the *natural* spline coefficients from the nodal data, (a_x, a_f).
// -----------------------------------------------------------------------------
void CubicSpline::solve (const Vector<Real>& a_f,
                         const Vector<Real>& a_x)
{
    // Start clean.
    this->clear();

    // Make sure input vectors are the same size
    CH_assert(a_x.size() == a_f.size());

    // Create our own copies of the input data.
    m_x = a_x;
    m_f = a_f;

    // Set Neumann BCs. The fortran function will create natural splines if
    // the BCs are chosen to be larger than 1e30.
    Real lofbc = 1.0e50;
    Real hifbc = 1.0e50;

    // Create space for the second derivatives
    m_d2f.resize(a_f.size());

    // Create workspace for the solver
    Vector<Real> workspace(a_f.size());

    // Solve for the second derivatives
    FORT_CUBICSPLINE_SOLVESECONDDERIV (
        CHF_VR(m_d2f),
        CHF_CONST_VR(m_x),
        CHF_CONST_VR(m_f),
        CHF_CONST_REAL(lofbc),
        CHF_CONST_REAL(hifbc),
        CHF_VR(workspace));
}


// -----------------------------------------------------------------------------
// Fills a_f with interpolated values at positions indicated by a_x.
// -----------------------------------------------------------------------------
void CubicSpline::interp (Vector<Real>&       a_f,
                          const Vector<Real>& a_x) const
{
    // Make sure the solve has been called.
    CH_assert(m_d2f.size() > 0);

    // Make sure the output vector has room for the results.
    const int Nax = a_x.size();
    a_f.resize(Nax);

    // Loop over a_x and compute f(x).
    Real x, xlo, xhi, dx;
    Real A, B, C, D;
    int klo, khi, k;
    const int Nmx = m_x.size();

    for (int i = 0; i < Nax; ++i) {
        // Search for x(i) in m_x.
        x = a_x[i];
        klo = 0;
        khi = Nmx-1;

        while (khi-klo > 1) {
            k = (khi + klo) >> 1;
            if (m_x[k] > x) {
                khi = k;
            } else {
                klo = k;
            }
        }

        // Compute cell extents and width.
        xlo = m_x[klo];
        xhi = m_x[khi];
        dx = xhi - xlo;
        CH_assert(dx != 0.0);

        // Predict value using linear interpolant.
        A = (xhi - x) / dx;
        B = (x - xlo) / dx;
        a_f[i] = A*m_f[klo] + B*m_f[khi];

        // Apply cubic correction.
        C = A*(A*A-1.0);
        D = B*(B*B-1.0);
        a_f[i] += (C*m_d2f[klo] + D*m_d2f[khi]) * dx*dx/6.0;
    }
}


// -----------------------------------------------------------------------------
// Fills a_df with the interpolated first derivatives at positions
// indicated by a_x.
// -----------------------------------------------------------------------------
void CubicSpline::interpFirstDeriv (Vector<Real>&       a_df,
                                    const Vector<Real>& a_x) const
{
    // Make sure the solve has been called.
    CH_assert(m_d2f.size() > 0);

    // Make sure the output vector has room for the results.
    const int Nax = a_x.size();
    a_df.resize(Nax);

    // Loop over a_x and compute dx/dx.
    Real x, xlo, xhi, dx, df;
    Real A, B;
    int klo, khi, k;
    const int Nmx = m_x.size();

    for (int i = 0; i < Nax; ++i) {
        // Search for x(i) in m_x.
        x = a_x[i];
        klo = 0;
        khi = Nmx-1;

        while (khi-klo > 1) {
            k = (khi + klo) >> 1;
            if (m_x[k] > x) {
                khi = k;
            } else {
                klo = k;
            }
        }

        // Compute cell extents and width.
        xlo = m_x[klo];
        xhi = m_x[khi];
        dx = xhi - xlo;
        CH_assert(dx != 0.0);

        // Compute finite difference approximation to the derivative.
        df = m_f[khi] - m_f[klo];
        a_df[i] = df/dx;

        // Apply higher-order corrections
        A = (xhi - x) / dx;
        B = (x - xlo) / dx;
        a_df[i] += (- (3.0*A*A-1.0) * m_d2f[klo]
                    + (3.0*B*B-1.0) * m_d2f[khi]) * (dx / 6.0);
    }
}


// -----------------------------------------------------------------------------
// Fills a_d2f with the interpolated second derivatives at positions
// indicated by a_x.
// -----------------------------------------------------------------------------
void CubicSpline::interpSecondDeriv (Vector<Real>&       a_d2f,
                                     const Vector<Real>& a_x) const
{
    // Make sure the solve has been called.
    CH_assert(m_d2f.size() > 0);

    // Make sure the output vector has room for the results.
    const int Nax = a_x.size();
    a_d2f.resize(Nax);

    // Loop over a_x and compute second derivative.
    Real x, xlo, xhi, dx;
    Real A, B;
    int klo, khi, k;
    const int Nmx = m_x.size();

    for (int i = 0; i < Nax; ++i) {
        // Search for x(i) in m_x.
        x = a_x[i];
        klo = 0;
        khi = Nmx-1;

        while (khi-klo > 1) {
            k = (khi + klo) >> 1;
            if (m_x[k] > x) {
                khi = k;
            } else {
                klo = k;
            }
        }

        // Compute cell extents and width.
        xlo = m_x[klo];
        xhi = m_x[khi];
        dx = xhi - xlo;
        CH_assert(dx != 0.0);

        // Second derivative is just a linear combination of values at nodes.
        A = (xhi - x) / dx;
        B = (x - xlo) / dx;
        a_d2f[i] = A * m_d2f[klo] + B * m_d2f[khi];
    }
}


// -----------------------------------------------------------------------------
// Uses a pre-computed set of data.
// -----------------------------------------------------------------------------
void CubicSpline::useSolution (const Vector<Real>& a_x,
                               const Vector<Real>& a_f,
                               const Vector<Real>& a_d2f)
{
    const int numNodes = a_x.size();
    CH_assert(a_f.size() == numNodes);
    CH_assert(a_d2f.size() == numNodes);

    m_x   = a_x;
    m_f   = a_f;
    m_d2f = a_d2f;
}
