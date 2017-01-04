/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2016
 *      Philadelphia University &
 *      University of North Carolina at Chapel Hill
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

#ifndef ___BGScalarProfiles_H__INCLUDED___
#define ___BGScalarProfiles_H__INCLUDED___

#include "EllipticBCInterface.H"


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class LinearBGScalarProfile: public EllipticBCValueClass
{
public:
    LinearBGScalarProfile (Real a_Nsq)
    : m_Nsq(a_Nsq)
    {;}

    virtual ~LinearBGScalarProfile ()
    {;}

    virtual void operator() (Real*           a_pos,
                             int*            a_dir,
                             Side::LoHiSide* a_side,
                             Real*           a_value,
                             const RealVect& a_dx,
                             Real            a_time) const
    {
        a_value[0] = -m_Nsq * a_pos[CH_SPACEDIM-1];
    }

protected:
    Real m_Nsq;
};


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class QuadraticBGScalarProfile: public EllipticBCValueClass
{
public:
    QuadraticBGScalarProfile (Real a_NsqBottom,
                              Real a_NsqTop)
    : m_NsqBottom(a_NsqBottom),
      m_NsqTop(a_NsqTop)
    {
        LevelGeometry::staticDefine();
        m_H = LevelGeometry::getDomainLength(CH_SPACEDIM-1);
        m_zCoeff = (m_NsqBottom - m_NsqTop) / m_H;
    }

    virtual ~QuadraticBGScalarProfile ()
    {;}

    virtual void operator() (Real*           a_pos,
                             int*            a_dir,
                             Side::LoHiSide* a_side,
                             Real*           a_value,
                             const RealVect& a_dx,
                             Real            a_time) const
    {
        Real z = a_pos[CH_SPACEDIM-1];
        Real Nsq = m_NsqBottom - m_zCoeff * z;
        a_value[0] = -Nsq * z;
    }

protected:
    Real m_NsqBottom;
    Real m_NsqTop;
    Real m_H;
    Real m_zCoeff;
};


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class TanhBGScalarProfile: public EllipticBCValueClass
{
public:
    TanhBGScalarProfile (Real a_z0,
                         Real a_delta,
                         Real a_rho0,
                         Real a_drho)
    : m_z0(a_z0),
      m_delta(a_delta),
      m_rho0(a_rho0),
      m_drho(a_drho)
    {;}

    virtual ~TanhBGScalarProfile ()
    {;}

    virtual void operator() (Real*           a_pos,
                             int*            a_dir,
                             Side::LoHiSide* a_side,
                             Real*           a_value,
                             const RealVect& a_dx,
                             Real            a_time) const
    {
        Real arg = (a_pos[CH_SPACEDIM-1] - m_z0) / m_delta;
        a_value[0] = m_rho0 - m_drho*tanh(arg);
    }

protected:
    Real m_z0;      // Center location
    Real m_delta;   // Thickness
    Real m_rho0;    // Median buoyancy value
    Real m_drho;    // Max buoyancy difference
};


#endif //!___BGScalarProfiles_H__INCLUDED__

