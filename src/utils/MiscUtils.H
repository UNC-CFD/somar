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

// This file contains random tools that have no home.
#ifndef __MiscUtils_H__INCLUDED__
#define __MiscUtils_H__INCLUDED__

#include "Copier.H"
#include "BoxLayout.H"
#include "BaseFab.H"
#include "LevelData.H"
#include "LevelGeometry.H"


// Takes a DBL and returns a new one whose boxes are merged in a_dir.
void mergeLayout (DisjointBoxLayout&       a_newLayout,
                  const DisjointBoxLayout& a_origLayout,
                  const int                a_dir);

// Splits a box into a load balanced set of boxes that are not split in a_dir.
void lineLayout (DisjointBoxLayout&   a_newLayout,
                 const ProblemDomain& a_domain,
                 const int            a_dir);

// Define a dbl with just one box.
// This function must be called on ALL procs, but a_box only needs to be
// defined on a_srcProc.
void defineOneProcGrids (DisjointBoxLayout&   a_grids,
                         const ProblemDomain& a_domain,
                         Box                  a_box,
                         const int            a_srcProc = uniqueProc(SerialTask::compute));


// This will define a copier that does not care about valid vs invalid data -
// it will just copy everything.
void defineImpartialCopier (Copier&                  a_copier,
                            const DisjointBoxLayout& a_srcGrids,
                            const BoxLayout&         a_destGrids,
                            const IntVect&           a_ghostVect,
                            const IntVect&           a_shift = IntVect::Zero);


// Shifts all boxes in a BoxLayout while preserving proc assignments, etc.
class ShiftTransform: public BaseTransform
{
public:
    ShiftTransform (const IntVect& a_shift)
    : m_shift(a_shift)
    {;}

    ShiftTransform (const int a_dir, const int a_shift)
    : m_shift(a_shift * BASISV(a_dir))
    {;}

    virtual ~ShiftTransform ()
    {;}

    virtual Box operator() (const Box& a_inputBox)
    {
        return Box(a_inputBox).shift(m_shift);
    }

protected:
    IntVect m_shift;
};


// Slices all boxes in a BoxLayout while preserving proc assignments, etc.
class SliceTransform: public BaseTransform
{
public:
    SliceTransform (const int a_dir, const int a_pos)
    : m_dir(a_dir), m_pos(a_pos)
    {;}

    virtual ~SliceTransform ()
    {;}

    virtual Box operator() (const Box& a_inputBox)
    {
        Box retBox = a_inputBox;
        retBox.setBig(m_dir, a_inputBox.smallEnd(m_dir));
        retBox.shift(m_dir, m_pos - a_inputBox.smallEnd(m_dir));
        CH_assert(!retBox.isEmpty());
        return retBox;
    }

protected:
    int m_dir, m_pos;
};


// This shifts shiftFAB so that its lower vertical end coincides with the lower
// vertical end of target*. The centering will also be changed. You MUST call
// restore() before this object goes out of scope or an error will be thrown.
template<typename T>
class VertShifter
{
    typedef BaseFab<T> DataType;

public:
    VertShifter (DataType&       a_shiftFAB,
                 const DataType& a_targetFAB)
    : m_dataPtr(NULL)
    {
        this->shift(a_shiftFAB, a_targetFAB.box());
    }

    VertShifter (DataType&  a_shiftFAB,
                 const Box& a_targetBox)
    : m_dataPtr(NULL)
    {
        this->shift(a_shiftFAB, a_targetBox);
    }

    ~VertShifter ()
    {
#ifndef NDEBUG
        if (m_dataPtr != NULL) {
            MayDay::Error("You forgot to call VertShifter<T>::restore");
        }
#endif
    }

    void restore ()
    {
        CH_assert(m_dataPtr != NULL);
        m_dataPtr->shift(SpaceDim-1, -m_totalShift);
        if (m_typeIsChanged) {
            m_dataPtr->shiftHalf(SpaceDim-1, -1);
        }
        m_dataPtr = NULL;
    }

private:
    // This is private so that we can only use this object once.
    void shift (DataType&  a_shiftFAB,
                const Box& a_targetBox)
    {
        CH_assert(m_dataPtr == NULL);
        m_dataPtr = &a_shiftFAB;

        // Change type
        if (a_targetBox.type(SpaceDim-1) != a_shiftFAB.box().type(SpaceDim-1)) {
            a_shiftFAB.shiftHalf(SpaceDim-1, 1);
            m_typeIsChanged = true;
        } else {
            m_typeIsChanged = false;
        }

        // Shift
        m_totalShift = a_targetBox.smallEnd(SpaceDim-1)
                     - a_shiftFAB.box().smallEnd(SpaceDim-1);
        a_shiftFAB.shift(SpaceDim-1, m_totalShift);
    }

    int  m_totalShift;
    int  m_typeIsChanged;
    DataType* m_dataPtr;
};


#endif //!__MiscUtils_H__INCLUDED__
