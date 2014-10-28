#ifdef CH_LANG_CC
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

// This file contains a set of classes that help the elliptic operators identify
// what BCs are being applyed in an attempt to minimize numerical noise.
#endif // CH_LANG_CC

#ifndef __BCDescriptor_H__INCLUDED__
#define __BCDescriptor_H__INCLUDED__

#ifdef CH_LANG_FORT
#   define BCType_Undefined -2
#   define BCType_None      -1
#   define BCType_Neum       0
#   define BCType_Diri       1
#   define BCType_Periodic   2
#   define BCType_CF         3
#endif //CH_LANG_FORT


#ifdef CH_LANG_CC
// The rest of this file contains C++ code

#include "ProblemDomain.H"
#include "CFRegion.H"


// The standard BC types
class BCType
{
public:
    enum {
        Undefined = -2,
        None      = -1,
        Neum      =  0,
        Diri      =  1,
        Periodic  =  2,
        CF        =  3
    };
};

// Used as a return value by BC-setting functions.
class BCDirDescriptor
{
public:
    BCDirDescriptor ()
    {
        m_bctype[0] = BCType::None;
        m_bctype[1] = BCType::None;
    }

    BCDirDescriptor (const BCDirDescriptor& a_src)
    {
        if (&a_src != this) {
            m_bctype[0] = a_src.m_bctype[0];
            m_bctype[1] = a_src.m_bctype[1];
        }
    }

    ~BCDirDescriptor ()
    {;}

    int& operator[] (Side::LoHiSide a_side)
    {
        return m_bctype[a_side];
    }

    const int& operator[] (Side::LoHiSide a_side) const
    {
        return m_bctype[a_side];
    }

    BCDirDescriptor& operator= (const BCDirDescriptor& a_src)
    {
        if (&a_src != this) {
            m_bctype[0] = a_src.m_bctype[0];
            m_bctype[1] = a_src.m_bctype[1];
        }
        return *this;
    }

private:
    int m_bctype[2];
};

class BCDescriptor
{
public:
    BCDescriptor ()
    {;}

    BCDescriptor (const BCDescriptor& a_src)
    {
        if (&a_src != this) {
            for (int dir = 0; dir < SpaceDim; ++dir) {
                m_bcdir[dir] = a_src.m_bcdir[dir];
            }
        }
    }

    BCDescriptor (int a_bcType)
    {
        this->set(a_bcType);
    }

    BCDescriptor (int a_dir, int a_bcType)
    {
        this->set(a_dir, Side::Lo, a_bcType);
        this->set(a_dir, Side::Hi, a_bcType);
    }

    BCDescriptor (int a_dir, Side::LoHiSide a_side, int a_bcType)
    {
        this->set(a_dir, a_side, a_bcType);
    }

    ~BCDescriptor ()
    {;}

    void set (int a_dir, Side::LoHiSide a_side, int a_bcType)
    {
        CH_assert(0 <= a_dir && a_dir < SpaceDim);
        m_bcdir[a_dir][a_side] = a_bcType;
    }

    void set (int a_bcType)
    {
        for (int dir = 0; dir < SpaceDim; ++dir) {
            m_bcdir[dir][Side::Lo] = a_bcType;
            m_bcdir[dir][Side::Hi] = a_bcType;
        }
    }

    BCDirDescriptor& operator[] (int a_dir)
    {
        CH_assert(0 <= a_dir && a_dir < SpaceDim);
        return m_bcdir[a_dir];
    }

    const BCDirDescriptor& operator[] (int a_dir) const
    {
        CH_assert(0 <= a_dir && a_dir < SpaceDim);
        return m_bcdir[a_dir];
    }

    BCDescriptor& operator= (const BCDescriptor& a_src)
    {
        if (&a_src != this) {
            for (int dir = 0; dir < SpaceDim; ++dir) {
                m_bcdir[dir] = a_src.m_bcdir[dir];
            }
        }
        return *this;
    }

    BCDescriptor& operator|= (const BCDescriptor& a_src)
    {
        for (int idir = 0; idir < SpaceDim; ++idir) {
            SideIterator sit;
            for (sit.reset(); sit.ok(); ++sit) {
                Side::LoHiSide iside = sit();
                int& dest = m_bcdir[idir][iside];
                const int& src = a_src[idir][iside];

                if ((dest == BCType::Undefined || src == BCType::Undefined) || (dest != BCType::None && src != BCType::None))
                    dest = BCType::Undefined;
                else if (dest == BCType::None)
                    dest = src;
            }
        }
        return *this;
    }

    BCDescriptor operator| (const BCDescriptor& a_src)
    {
        BCDescriptor ret(*this);
        ret |= a_src;
        return ret;
    }


    // This version needs to be fast. It will be used in GSRB::boundaryGSRB()
    // and really only needs to detect Neumann BCs, which can only occur at the
    // domain boundary.
    inline int stencil (const Box&           a_region,
                        const ProblemDomain& a_domain,
                        int                  a_dir,
                        Side::LoHiSide       a_side) const
    {
        if ((!a_domain.isPeriodic(a_dir)) &&
            (a_region.sideEnd(a_side)[a_dir] == a_domain.domainBox().sideEnd(a_side)[a_dir]))
            return m_bcdir[a_dir][a_side];
        return BCType::None;
    }

    // This version also tests for Diri BCs due to a CF interface.
    // This is intended for use by the Leptic solver with a_dir set to the vertical.
    // It will fail on general grids.
    int stencil (const Box&           a_region,
                 const ProblemDomain& a_domain,
                 int                  a_dir,
                 Side::LoHiSide       a_side,
                 CFRegion&            a_cfregion,
                 const DataIndex&     a_index) const;

private:
    BCDirDescriptor m_bcdir[CH_SPACEDIM];
};


#endif //CH_LANG_CC
#endif //!__BCDescriptor_H__INCLUDED__
