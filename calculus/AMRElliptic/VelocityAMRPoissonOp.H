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

// This is a rework of code found in the Chombo library. While these functions
// have been reproduced for ease of modification, some of these functions may
// not have been altered at all. To that end, you can find Chombo's copyright
// file at somar/Chombo/Copyright.txt.

#ifndef __VELOCITYAMRPOISSONOP_H__INCLUDED__
#define __VELOCITYAMRPOISSONOP_H__INCLUDED__

#include "VelBCHolder.H"
#include "MappedAMRPoissonOp.H"


// The reason for this class -- want to be able to apply operator
// with different BCs for each component in a_phi.
class VelocityAMRPoissonOp: public MappedAMRPoissonOp
{
public:
    // Default constructor
    VelocityAMRPoissonOp()
    {;}

    // Destructor
    virtual ~VelocityAMRPoissonOp()
    {;}

    // Defines this op to simply calculate Laplacians and residuals.
    void define(const LevelGeometry&                        a_levGeo,
                const Real                                  a_alpha,
                const Real                                  a_beta,
                const RefCountedPtr<LevelData<FArrayBox> >& a_lapDiagsPtr = RefCountedPtr<LevelData<FArrayBox> >(NULL));

    //
    virtual void applyOpI(LevelData<FArrayBox>&        a_lhs,
                          const LevelData<FArrayBox>&  a_phi,
                          bool                         a_homogeneous,
                          VelBCHolder*                 a_velBCPtr);

    ///
    virtual void AMROperatorNF(LevelData<FArrayBox>&        a_LofPhi,
                               const LevelData<FArrayBox>&  a_phi,
                               const LevelData<FArrayBox>&  a_phiCoarse,
                               bool                         a_homogeneousPhysBC,
                               VelBCHolder*                 a_velBCPtr);

    ///
    virtual void applyOp(LevelData<FArrayBox>&        a_lhs,
                         const LevelData<FArrayBox>&  a_phi,
                         const LevelData<FArrayBox>*  a_phiCoarsePtr,
                         bool                         a_homogeneous,
                         VelBCHolder*                 a_velBCPtr);
};

#endif //!__VELOCITYAMRPOISSONOP_H__INCLUDED__

