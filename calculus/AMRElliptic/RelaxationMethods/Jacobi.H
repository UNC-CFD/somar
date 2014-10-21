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
#ifndef __Jacobi_H__INCLUDED__
#define __Jacobi_H__INCLUDED__

#include "RelaxationMethod.H"
#include "LinearSolver.H"


// -----------------------------------------------------------------------------
// The Jacobi relaxation method.
// -----------------------------------------------------------------------------
class Jacobi: public RelaxationMethod
{
    typedef LinearOp<LevelData<FArrayBox> > OperatorType;
public:
    virtual ~Jacobi ();

    Jacobi (OperatorType* a_opPtr,
            const Real    a_alpha,
            const Real    a_beta);

    virtual void relax (LevelData<FArrayBox>&       a_phi,
                        const LevelData<FArrayBox>& a_rhs);

protected:
    OperatorType* m_opPtr;
    Real          m_alpha;
    Real          m_beta;
};


#endif //!__Jacobi_H__INCLUDED__
