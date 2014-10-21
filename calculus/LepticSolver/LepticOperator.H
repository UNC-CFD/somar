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
#ifndef __LepticOperator_H__INCLUDED__
#define __LepticOperator_H__INCLUDED__

#include "RealVect.H"
#include "LevelData.H"
#include "EllipticBCInterface.H"


// To use an operator in a LepticSolver, it must inherit from LinearOp
// and this pure abstract base class.
class LepticOperator {
public:
    virtual ~LepticOperator () {;}

    virtual const RealVect& getDx () const = 0;
    virtual const RealVect& getDxCrse () const = 0; // AMR, not MG
    virtual void setDxCrse (const RealVect& a_dxCrse) = 0; // AMR, not MG
    virtual const RefCountedPtr<LevelData<FluxBox> > getFCJgup () const = 0;
    virtual const RefCountedPtr<LevelData<FArrayBox> > getCCJinv () const = 0;
    virtual BCMethodHolder getBCs () const = 0;
    virtual int smallestGridSize () const = 0;
};


#endif //!__LepticOperator_H__INCLUDED__
