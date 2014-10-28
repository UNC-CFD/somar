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
#ifndef __BoundaryData_H__INCLUDED__
#define __BoundaryData_H__INCLUDED__

#include "FArrayBox.H"
#include "LevelData.H"
#include "RealVect.H"
#include "DisjointBoxLayout.H"
#include "DataIterator.H"
#include "RefCountedPtr.H"
#include "DataIndexCompare.H"
#include "Subspace.H"
#include "EllipticBCInterface.H"


// -----------------------------------------------------------------------------
// The BoundaryData template.
// Caches data at the physical boundary faces.
// -----------------------------------------------------------------------------
template <class T>
class BoundaryData
{
public:
    // Data types --------------------------------------------------------------

    // Pointer to a data holder
    typedef RefCountedPtr<BaseFab<T> > DataPtrType;

    // A DataIndex -> BaseFab<T> map.
    typedef std::map<DataIndex, DataPtrType, DataIndexCompare> BoundaryDataMapType;

    // The BoundaryDataMapType's iterator
    typedef typename BoundaryDataMapType::iterator       iterator;
    typedef typename BoundaryDataMapType::const_iterator const_iterator;


    // Construction / destruction ----------------------------------------------

    // Default constructor (leaves unusable)
    BoundaryData ();

    // Full constructor (calls define)
    BoundaryData (const DisjointBoxLayout& a_grids,
                  const ProblemDomain&     a_domain,
                  const int                a_ncomp = 1);

    // Copy constructor
    BoundaryData (const BoundaryData<T>& a_src);

    // Destructor
    virtual ~BoundaryData ();

    // Full define
    virtual void define (const DisjointBoxLayout& a_grids,
                         const ProblemDomain&     a_domain,
                         const int                a_ncomp = 1);

    // Copy operator
    virtual BoundaryData<T>& operator= (const BoundaryData<T>& a_src);

    // Copies a_src entirely without sharing pointers.
    // NOTE: m_grids will still be a shallow copy.
    virtual void deepCopy (const BoundaryData<T>& a_src);

    // Clears the BoundaryRegion and leaves this object unusable.
    virtual void clear ();


    // Inquiries ---------------------------------------------------------------

    // Checks if a DataIndex is compatible with the grids that defined this object.
    virtual bool check (const DataIndex& a_di) const;

    // Checks if the grids are compatible and data lies on the exact same boxes.
    virtual bool sameLayout (const BoundaryData<T>& a_src) const;

    // Checks if a_src is a flattened version of this object.
    virtual bool flattenedLayout (const BoundaryData<T>& a_flatSrc) const;

    // Checks if this object can be flattened
    virtual bool flattenable () const;

    // Checks if this object is vertically flat.
    virtual bool isFlat () const;


    // Iterators ---------------------------------------------------------------

    // Returns an iterator at the start of the appropriate map. (const version)
    virtual const_iterator begin (const int             a_dir,
                                  const Side::LoHiSide& a_side) const;

    // Returns an iterator at the start of the appropriate map.
    virtual iterator begin (const int             a_dir,
                            const Side::LoHiSide& a_side);

    // Returns an iterator at the end of the appropriate map. (const version)
    virtual const_iterator end (const int             a_dir,
                                const Side::LoHiSide& a_side) const;

    // Returns an iterator at the end of the appropriate map.
    virtual iterator end (const int             a_dir,
                          const Side::LoHiSide& a_side);


    // Data accessors ----------------------------------------------------------

    // Returns the appropriate FC boundary data holder. (const version)
    virtual const BaseFab<T>& getData (const DataIndex&      a_di,
                                       const int             a_dir,
                                       const Side::LoHiSide& a_side) const;

    // Returns the appropriate FC boundary data holder.
    virtual BaseFab<T>& getData (const DataIndex&      a_di,
                                 const int             a_dir,
                                 const Side::LoHiSide& a_side);


    // Data modification -------------------------------------------------------

    // Sets all data to a uniform value.
    virtual void setVal (T a_val);

protected:
    // Member variables --------------------------------------------------------

    // This is used to call check() on the input DataIndex.
    DisjointBoxLayout m_grids;

    // The map of boundary faces.
    BoundaryDataMapType m_bdryFaceMap[CH_SPACEDIM][2];

    // Has this object been vertically averaged or set up to hold flat data?
    bool m_isFlat;

    // An empty data holder
    static BaseFab<T> s_emptyFab;
};



// -----------------------------------------------------------------------------
// The BoundaryData template specialization: T = FArrayBox.
// -----------------------------------------------------------------------------
template <>
class BoundaryData<Real>
{
    // Creates a coarsened BoundaryData object.
    friend void coarsen (BoundaryData<Real>&       a_crse,
                         const BoundaryData<Real>& a_fine,
                         const IntVect&            a_ref);

public:
    // Data types --------------------------------------------------------------

    // Pointer to a data holder
    typedef RefCountedPtr<FArrayBox> DataPtrType;

    // A DataIndex -> BaseFab<T> map.
    typedef std::map<DataIndex, DataPtrType, DataIndexCompare> BoundaryDataMapType;

    // The BoundaryDataMapType's iterator
    typedef typename BoundaryDataMapType::iterator iterator;
    typedef typename BoundaryDataMapType::const_iterator const_iterator;


    // Construction / destruction ----------------------------------------------

    // Default constructor (leaves unusable)
    BoundaryData ();

    // Full constructor (calls define)
    BoundaryData (const DisjointBoxLayout& a_grids,
                  const ProblemDomain&     a_domain,
                  const int                a_ncomp = 1);

    // Full constructor (calls define)
    BoundaryData (const LevelData<FluxBox>& a_src);

    // Copy constructor
    BoundaryData (const BoundaryData<Real>& a_src);

    // Destructor
    virtual ~BoundaryData ();

    // Full define
    virtual void define (const DisjointBoxLayout& a_grids,
                         const ProblemDomain&     a_domain,
                         const int                a_ncomp = 1);

    // Full define and copy
    virtual void define (const LevelData<FluxBox>& a_src);

    // Copy operator
    virtual BoundaryData<Real>& operator= (const BoundaryData<Real>& a_src);

    // Copies a_src entirely without sharing pointers.
    // NOTE: m_grids will still be a shallow copy.
    virtual void deepCopy (const BoundaryData<Real>& a_src);

    // Clears the BoundaryRegion and leaves this object unusable.
    virtual void clear ();


    // Inquiries ---------------------------------------------------------------

    // Checks if a DataIndex is compatible with the grids that defined this object.
    virtual bool check (const DataIndex& a_di) const;

    // Checks if the grids are compatible and data lies on the exact same boxes.
    virtual bool sameLayout (const BoundaryData<Real>& a_src) const;

    // Checks if a_src is a flattened version of this object.
    virtual bool flattenedLayout (const BoundaryData<Real>& a_flatSrc) const;

    // Checks if this object can be flattened
    virtual bool flattenable () const;

    // Checks if this object is vertically flat.
    virtual bool isFlat () const;

    // Performs a simple integral of the boundary fluxes, assuming the data actually
    // represents fluxes (scaled by J). The bottom boundaries will be scaled by
    // a_loMult. Set this to -1.0 (default) to compute the net flux into the volume.
    virtual RealVect integrate (const RealVect& a_dx,
                                const Real      a_loMult = -1.0) const;

    // Computes Integral[J.rhs.dV] - Integral[bdry fluxes.dA]
    virtual Real consistencyCheck (const LevelData<FArrayBox>& a_Jrhs,
                                   const RealVect&             a_dx,
                                   const Real                  a_loMult = -1.0) const;


    // Iterators ---------------------------------------------------------------

    // Returns an iterator at the start of the appropriate map. (const version)
    virtual const_iterator begin (const int             a_dir,
                                  const Side::LoHiSide& a_side) const;

    // Returns an iterator at the start of the appropriate map.
    virtual iterator begin (const int             a_dir,
                            const Side::LoHiSide& a_side);

    // Returns an iterator at the end of the appropriate map. (const version)
    virtual const_iterator end (const int             a_dir,
                                const Side::LoHiSide& a_side) const;

    // Returns an iterator at the end of the appropriate map.
    virtual iterator end (const int             a_dir,
                          const Side::LoHiSide& a_side);


    // Data accessors ----------------------------------------------------------

    // Returns the appropriate FC boundary data holder. (const version)
    virtual const FArrayBox& getData (const DataIndex&      a_di,
                                      const int             a_dir,
                                      const Side::LoHiSide& a_side) const;

    // Returns the appropriate FC boundary data holder.
    virtual FArrayBox& getData (const DataIndex&      a_di,
                                const int             a_dir,
                                const Side::LoHiSide& a_side);

    // Writes boundary data to a holder that can be sent to HDF5.
    virtual void copyTo (LevelData<FluxBox>& a_dest) const;


    // Data modification -------------------------------------------------------

    // Copies boundary values. Grids must be compatible.
    virtual void copy (const LevelData<FluxBox>& a_src);

    // Sets all data to a uniform value.
    virtual void setVal (Real a_val);

    // Gathers boundary fluxes from a BCFluxClass.
    virtual void collectFluxes (const BCFluxClass& a_bc,
                                const Real         a_time);

    // Accumulates values from another BoundaryData object.
    virtual const BoundaryData<Real>& plus (const BoundaryData<Real>& a_src,
                                            const Real                a_scale);

    // Accumulates values from a flat BoundaryData object.
    // This only operates on the horizontal BCs.
    virtual const BoundaryData<Real>& plusFlat (const BoundaryData<Real>& a_flatSrc,
                                                const Real                a_scale);

    // Accumulates values from a flat LevelData.
    // This only operates on the vertical BCs.
    virtual const BoundaryData<Real>& vertPlus (const BoxLayoutData<FArrayBox>& a_flatSrc,
                                                const Real                      a_scale,
                                                const Side::LoHiSide            a_side);

    // Converts the horizontal boundary data to a flattened vertical average.
    // The vertical boundary data is deallocated.
    virtual void vertAvg ();


protected:
    // Member variables --------------------------------------------------------

    // This is used to call check() on the input DataIndex.
    DisjointBoxLayout m_grids;

    // The map of boundary faces.
    BoundaryDataMapType m_bdryFaceMap[CH_SPACEDIM][2];

    // Has this object been vertically averaged or set up to hold flat data?
    bool m_isFlat;

    // An empty data holder.
    static FArrayBox s_emptyFab;
};



// Templated definitions
#include "BoundaryDataI.H"

#endif //!__BoundaryData_H__INCLUDED__

