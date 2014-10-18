/*******************************************************************************
 *    FILE: AMREllipticSolver.H
 *    DESC: A unified interface for all AMR elliptic solvers.
 *    DATE: Mon 19 May 2014 03:32:09 PM EDT  (UTC-0400)
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef AMREllipticSolver_H__INCLUDED__
#define AMREllipticSolver_H__INCLUDED__

#include "Vector.H"


template<class T>
class AMREllipticSolver
{
public:
    virtual ~AMREllipticSolver ()
    {;}

    virtual void solve (Vector<T*>&       a_phi,
                        const Vector<T*>& a_rhs,
                        const int         l_max,
                        const int         l_base,
                        const bool        a_zeroPhi,
                        const bool        a_forceHomogeneous) = 0;
};


#endif //!AMREllipticSolver_H__INCLUDED__
