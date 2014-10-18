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
#ifndef AMREllipticSolver_H__INCLUDED__
#define AMREllipticSolver_H__INCLUDED__

#include "Vector.H"


// -----------------------------------------------------------------------------
// A unified interface for all AMR elliptic solvers.
// -----------------------------------------------------------------------------
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
