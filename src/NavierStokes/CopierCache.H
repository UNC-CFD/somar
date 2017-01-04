/*******************************************************************************
 *  SOMAR - Stratified Ocean Model with Adaptive Refinement
 *  Developed by Ed Santilli & Alberto Scotti
 *  Copyright (C) 2015
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
#ifndef ___CopierCache_H__INCLUDED___
#define ___CopierCache_H__INCLUDED___

#include "DisjointBoxLayout.H"
#include "IntVect.H"
#include "Copier.H"
#include "CornerCopier.H"


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
class CopierCache
{
private:
    struct CopierEntry
    {
        Copier            cp;
        DisjointBoxLayout grids;
        IntVect           ghostVect;
        bool              isExchange;

        void define (const DisjointBoxLayout& a_grids,
                     const IntVect&           a_ghostVect,
                     const bool               a_isExchange);
    };

    struct CornerCopierEntry
    {
        CornerCopier      cp;
        DisjointBoxLayout grids;
        IntVect           ghostVect;
        bool              isExchange;

        void define (const DisjointBoxLayout& a_grids,
                     const IntVect&           a_ghostVect,
                     const bool               a_isExchange);
    };

public:
    // Constructor.
    CopierCache (const IntVect& a_advectGhostVect);

    // Destructor.
    ~CopierCache ();

    // Accessors. These check the grids and redefine if needed.
    const Copier& getTracingExCopier (const DisjointBoxLayout& a_grids);
    const CornerCopier& getTracingExCornerCopier (const DisjointBoxLayout& a_grids);
    const Copier& getOneGhostExCopier (const DisjointBoxLayout& a_grids);

    inline const IntVect& getTracingGhosts () const {return m_tracingGhosts;}

    // Defines the copiers with new grids.
    // This function re-defines the copiers even if they are already defined
    // over the same set of grids.
    void regrid (const DisjointBoxLayout& a_grids);

private:
    IntVect           m_tracingGhosts;
    CopierEntry       m_tracingExCopierEntry;
    CornerCopierEntry m_tracingExCornerCopierEntry;
    CopierEntry       m_oneGhostExCopierEntry;
};


#endif //!___CopierCache_H__INCLUDED___
