#include "BoundaryData.H"
#include "BoxIterator.H"
#include "computeMappedSum.H"
#include "SetValLevel.H"
#include "AnisotropicRefinementTools.H"
#include "MiscUtils.H"
#include "Debug.H"


// An empty data holder
FArrayBox BoundaryData<Real>::s_emptyFab;


// -----------------------------------------------------------------------------
// Default constructor (leaves unusable)
// -----------------------------------------------------------------------------
BoundaryData<Real>::BoundaryData ()
{
    this->clear();
}


// -----------------------------------------------------------------------------
// Full constructor (calls define)
// -----------------------------------------------------------------------------
BoundaryData<Real>::BoundaryData (const DisjointBoxLayout& a_grids,
                                  const ProblemDomain&     a_domain,
                                  const int                a_ncomp)
{
    this->define(a_grids, a_domain, a_ncomp);
}


// -----------------------------------------------------------------------------
// Full constructor (calls define)
// -----------------------------------------------------------------------------
BoundaryData<Real>::BoundaryData (const LevelData<FluxBox>& a_src)
{
    this->define(a_src);
}


// -----------------------------------------------------------------------------
// Copy constructor
// -----------------------------------------------------------------------------
BoundaryData<Real>::BoundaryData (const BoundaryData<Real>& a_src)
{
    this->deepCopy(a_src);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
BoundaryData<Real>::~BoundaryData ()
{
    this->clear();
}


// -----------------------------------------------------------------------------
// Full define
// -----------------------------------------------------------------------------
void BoundaryData<Real>::define (const DisjointBoxLayout& a_grids,
                                 const ProblemDomain&     a_domain,
                                 const int                a_ncomp)
{
    CH_assert(a_ncomp >= 1);

    // Remember the source layout so we can call check() later.
    m_grids = a_grids;

    // This object will not be flat until we call a function that flattens it.
    m_isFlat = false;

    const Box& domBox = a_domain.domainBox();
    DataIterator dit = a_grids.dataIterator();

    for (int dir = 0; dir < SpaceDim; ++dir) {
        // Only process non-periodic boundaries
        if (a_domain.isPeriodic()) continue;

        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide& iside = sit();
            const int s = int(iside);
            const Box domBdryFaces = bdryBox(domBox, dir, iside, 1);
            BoundaryDataMapType& thisMap = m_bdryFaceMap[dir][s];

            for (dit.reset(); dit.ok(); ++dit) {
                // Compute valid faces at the physical boundary.
                const Box& valid = a_grids[dit];
                Box validBdryFaces = bdryBox(valid, dir, iside, 1);
                validBdryFaces &= domBdryFaces;

                // If there are no cells at the boundary, just move on.
                if (validBdryFaces.isEmpty()) continue;

                // If we get here, then we have faces to store in the map.
                FArrayBox* newDataPtr = new FArrayBox(validBdryFaces, a_ncomp);
                thisMap[dit()] = DataPtrType(newDataPtr);

            } // end loop over grids (dit)
        } // end loop over sides (sit)
    } // end loop over directions (dir)
}


// -----------------------------------------------------------------------------
// Full define and copy
// -----------------------------------------------------------------------------
void BoundaryData<Real>::define (const LevelData<FluxBox>& a_src)
{
    const DisjointBoxLayout& grids = a_src.getBoxes();
    this->define(grids, grids.physDomain(), a_src.nComp());
    this->copy(a_src);
}


// -----------------------------------------------------------------------------
// Copy operator
// -----------------------------------------------------------------------------
BoundaryData<Real>& BoundaryData<Real>::operator= (const BoundaryData<Real>& a_src)
{
    this->deepCopy(a_src);
    return *this;
}


// -----------------------------------------------------------------------------
// Copies a_src entirely without sharing pointers.
// NOTE: m_grids will still be a shallow copy.
// -----------------------------------------------------------------------------
void BoundaryData<Real>::deepCopy (const BoundaryData<Real>& a_src)
{
    this->clear();

    m_grids = a_src.m_grids;
    m_isFlat = a_src.m_isFlat;

    const ProblemDomain& domain = m_grids.physDomain();
    const Box& domBox = domain.domainBox();
    DataIterator dit = m_grids.dataIterator();

    for (int dir = 0; dir < SpaceDim; ++dir) {
        // Only process non-periodic boundaries
        if (domain.isPeriodic()) continue;

        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide& iside = sit();
            const int s = int(iside);

            const_iterator cit = a_src.begin(dir, iside);
            const_iterator citEnd = a_src.end(dir, iside);
            for (; cit != citEnd; ++cit) {
                const DataIndex srcDi = cit->first;
                const RefCountedPtr<FArrayBox>& srcFABPtr = cit->second;
                const Box& srcBox = srcFABPtr->box();
                const int srcComps = srcFABPtr->nComp();

                RefCountedPtr<FArrayBox> newFABPtr(new FArrayBox(srcBox, srcComps));
                newFABPtr->copy(*srcFABPtr);
                m_bdryFaceMap[dir][s][srcDi] = newFABPtr;
            }
        } // end loop over sides (sit)
    } // end loop over directions (dir)
}


// -----------------------------------------------------------------------------
// Clears the std::maps and leaves this object unusable.
// -----------------------------------------------------------------------------
void BoundaryData<Real>::clear ()
{
    // We no longer need the grids for index checking.
    m_grids = DisjointBoxLayout();

    m_isFlat = false;

    // Loop over directions and sides, freeing memory.
    for (int dir = 0; dir < SpaceDim; ++dir) {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const int s = int(sit());
            m_bdryFaceMap[dir][s].clear();
        }
    }
}


// -----------------------------------------------------------------------------
// Checks if a DataIndex is compatible with the grids that defined this object.
// -----------------------------------------------------------------------------
bool BoundaryData<Real>::check (const DataIndex& a_di) const
{
    if (!m_grids.isClosed()) return false;
    return (m_grids.check(a_di));
}


// -----------------------------------------------------------------------------
// Checks if the grids are compatible and data lies on the exact same boxes.
// -----------------------------------------------------------------------------
bool BoundaryData<Real>::sameLayout (const BoundaryData<Real>& a_src) const
{
    // Both grids need to be well-defined.
    if (!m_grids.isClosed()) return false;
    if (!a_src.m_grids.isClosed()) return false;

    // Are the grids compatible?
    if (!m_grids.compatible(a_src.m_grids)) return false;

    // Are they both flat / not flat?
    if (isFlat() != a_src.isFlat()) return false;

    // Loop over directions and side. Check if maps are identical.
    for (int dir = 0; dir < SpaceDim; ++dir) {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const int s = int(sit());
            const BoundaryDataMapType& thisMap = m_bdryFaceMap[dir][s];
            const BoundaryDataMapType& srcMap = a_src.m_bdryFaceMap[dir][s];

            if (thisMap.size() != srcMap.size()) return false;

            const_iterator thisIt = thisMap.begin();
            const_iterator srcIt = srcMap.begin();
            const const_iterator endIt = thisMap.end();

            while (thisIt != endIt) {
                const DataIndex& thisDit = thisIt->first;
                const DataPtrType& thisFABPtr = thisIt->second;

                const DataIndex& srcDit = srcIt->first;
                const DataPtrType& srcFABPtr = srcIt->second;

                if (thisDit != srcDit) return false;
                if (thisFABPtr->box() != srcFABPtr->box()) return false;

                ++thisIt;
                ++srcIt;
            }
        } // end loop over sides (sit)
    } // end loop over directions (dir)

    // If we get here, the layouts are the same.
    return true;
}


// -----------------------------------------------------------------------------
// Checks if a_src is a flattened version of this object.
// -----------------------------------------------------------------------------
bool BoundaryData<Real>::flattenedLayout (const BoundaryData<Real>& a_flatSrc) const
{
    // If this object is flat or the src is not flat, why would we be here?!
    CH_assert(!isFlat());
    CH_assert(a_flatSrc.isFlat());

    // Both grids need to be well-defined.
    if (!m_grids.isClosed()) return false;
    if (!a_flatSrc.m_grids.isClosed()) return false;

    // Are the grids compatible?
    if (!m_grids.compatible(a_flatSrc.m_grids)) return false;

    // Loop over directions and side. Check if maps are identical.
    // NOTE: This does not check the vertical boundaries.
    for (int dir = 0; dir < SpaceDim-1; ++dir) {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const int s = int(sit());
            const BoundaryDataMapType& thisMap = m_bdryFaceMap[dir][s];
            const BoundaryDataMapType& srcMap = a_flatSrc.m_bdryFaceMap[dir][s];

            if (thisMap.size() != srcMap.size()) return false;

            const_iterator thisIt = thisMap.begin();
            const_iterator srcIt = srcMap.begin();
            const const_iterator endIt = thisMap.end();

            while (thisIt != endIt) {
                const DataIndex& thisDit = thisIt->first;
                const DataIndex& srcDit = srcIt->first;
                if (thisDit != srcDit) return false;

                const Box& thisFlatBox = flattenBox(thisIt->second->box(), SpaceDim-1);
                const Box& srcBox = srcIt->second->box();
                if (thisFlatBox != srcBox) return false;

                ++thisIt;
                ++srcIt;
            }
        } // end loop over sides (sit)
    } // end loop over directions (dir)

    // If we get here, the layouts are the same.
    return true;
}


// -----------------------------------------------------------------------------
// Checks if this object can be flattened
// -----------------------------------------------------------------------------
bool BoundaryData<Real>::flattenable () const
{
    CH_assert(m_grids.isClosed());
    CH_assert(!this->isFlat());

    const Box& domBox = m_grids.physDomain().domainBox();
    const int vertSize = domBox.size(SpaceDim-1);
    DataIterator dit = m_grids.dataIterator();

    for (dit.reset(); dit.ok(); ++dit) {
        if (m_grids[dit].size(SpaceDim-1) != vertSize) return false;
    }
    return true;
}


// -----------------------------------------------------------------------------
// Checks if this object is vertically flat.
// -----------------------------------------------------------------------------
bool BoundaryData<Real>::isFlat () const
{
    return m_isFlat;
}


// -----------------------------------------------------------------------------
// Performs a simple integral of the boundary fluxes, assuming the data actually
// represents fluxes (scaled by J). The bottom boundaries will be scaled by
// a_loMult. Set this to -1.0 (default) to compute the net flux into the volume.
// -----------------------------------------------------------------------------
RealVect BoundaryData<Real>::integrate (const RealVect& a_dx,
                                        const Real      a_loMult) const
{
    const bool horizOnly = this->isFlat();

    // Compute the scaling of each integral.
    RealVect localDx = a_dx;
    if (horizOnly) {
        localDx[SpaceDim-1] = 1.0;
    }
    const RealVect dA = localDx.product() / localDx;

    // Loop over boundary faces and integrate.
    RealVect localSum = RealVect::Zero;

    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (dir == SpaceDim-1 && horizOnly) continue;

        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide iside = sit();
            Real scale = dA[dir] * ((iside == Side::Lo)? a_loMult: 1.0);

            const_iterator it;
            for (it = begin(dir,iside); it != end(dir,iside); ++it) {
                localSum[dir] += scale * it->second->sum(0);
            }
        } // end of loop over sides (sit)
    } // end of loop over directions (dir)

    // Compute global sum (this is where the MPI communication happens)
#ifdef CH_MPI
    RealVect globalSum = RealVect::Zero;
    int result = MPI_Allreduce(localSum.dataPtr(), globalSum.dataPtr(), SpaceDim,
                               MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);

    if (result != MPI_SUCCESS) {
        MayDay::Error("Sorry, but I had a communication error in computeMappedSum");
    }

#else
    Real globalSum = localSum;
#endif

    // If this is a flat object, fix the z-comp.
    if (horizOnly) {
        globalSum[SpaceDim-1] = 0.0;
    }

    return globalSum;
}


// -----------------------------------------------------------------------------
// Computes Integral[J*rhs*dV] - Integral[bdry fluxes.dA]
// -----------------------------------------------------------------------------
Real BoundaryData<Real>::consistencyCheck (const LevelData<FArrayBox>& a_Jrhs,
                                           const RealVect&             a_dx,
                                           const Real                  a_loMult) const
{
    CH_assert(a_Jrhs.nComp() == 1);

    const bool horizOnly = this->isFlat();
    RealVect localDx = a_dx;
    if (horizOnly) {
        localDx[SpaceDim-1] = 1.0;
    }

    LevelData<FArrayBox> Jinv(a_Jrhs.getBoxes(), 1);
    setValLevel(Jinv, 1.0);

    Real vol = 0.0;
    Real intRhs = computeUnmappedSum(vol, a_Jrhs, NULL, IntVect::Unit, localDx, Jinv, 0);
    Real intFlux = this->integrate(a_dx).sum();

    return (intRhs - intFlux);
}


// -----------------------------------------------------------------------------
// Returns an iterator at the start of the appropriate map. (const version)
// -----------------------------------------------------------------------------
BoundaryData<Real>::const_iterator
BoundaryData<Real>::begin (const int             a_dir,
                           const Side::LoHiSide& a_side) const
{
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    CH_assert(!isFlat() || a_dir != SpaceDim-1);

    const int s = int(a_side);
    const_iterator it = m_bdryFaceMap[a_dir][s].begin();
    return it;
}


// -----------------------------------------------------------------------------
// Returns an iterator at the start of the appropriate map.
// -----------------------------------------------------------------------------
BoundaryData<Real>::iterator
BoundaryData<Real>::begin (const int             a_dir,
                           const Side::LoHiSide& a_side)
{
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    CH_assert(!isFlat() || a_dir != SpaceDim-1);

    const int s = int(a_side);
    iterator it = m_bdryFaceMap[a_dir][s].begin();
    return it;
}


// -----------------------------------------------------------------------------
// Returns an iterator at the end of the appropriate map. (const version)
// -----------------------------------------------------------------------------
BoundaryData<Real>::const_iterator
BoundaryData<Real>::end (const int             a_dir,
                         const Side::LoHiSide& a_side) const
{
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    CH_assert(!isFlat() || a_dir != SpaceDim-1);

    const int s = int(a_side);
    const_iterator it = m_bdryFaceMap[a_dir][s].end();
    return it;
}


// -----------------------------------------------------------------------------
// Returns an iterator at the end of the appropriate map.
// -----------------------------------------------------------------------------
BoundaryData<Real>::iterator
BoundaryData<Real>::end (const int             a_dir,
                         const Side::LoHiSide& a_side)
{
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    CH_assert(!isFlat() || a_dir != SpaceDim-1);

    const int s = int(a_side);
    iterator it = m_bdryFaceMap[a_dir][s].end();
    return it;
}


// -----------------------------------------------------------------------------
// Returns the appropriate FC boundary data holder. (const version)
// -----------------------------------------------------------------------------
const FArrayBox& BoundaryData<Real>::getData (const DataIndex&      a_di,
                                              const int             a_dir,
                                              const Side::LoHiSide& a_side) const
{
    // Sanity checks
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    CH_assert(!isFlat() || a_dir != SpaceDim-1);
    CH_assert(this->check(a_di));

    // Grab a reference to the appropriate map.
    const int s = int(a_side);
    const BoundaryDataMapType& thisMap = m_bdryFaceMap[a_dir][s];

    // Find the DataIndex in the map.
    const_iterator it;
    it = thisMap.find(a_di);

    // If we didn't find anything, return the default empty holder.
    if (it == thisMap.end()) return s_emptyFab;

    // We found a holder at the boundary. Return it.
    const DataPtrType thisDataPtr = it->second;
    CH_assert(!thisDataPtr.isNull());
    return *thisDataPtr;
}


// -----------------------------------------------------------------------------
// Returns the appropriate FC boundary data holder.
// -----------------------------------------------------------------------------
FArrayBox& BoundaryData<Real>::getData (const DataIndex&      a_di,
                                        const int             a_dir,
                                        const Side::LoHiSide& a_side)
{
    // Sanity checks
    CH_assert(0 <= a_dir && a_dir < SpaceDim);
    CH_assert(!isFlat() || a_dir != SpaceDim-1);
    CH_assert(this->check(a_di));

    // Grab a reference to the appropriate map.
    const int s = int(a_side);
    BoundaryDataMapType& thisMap = m_bdryFaceMap[a_dir][s];

    // Find the DataIndex in the map.
    iterator it;
    it = thisMap.find(a_di);

    // If we didn't find anything, return the default empty holder.
    if (it == thisMap.end()) return s_emptyFab;

    // We found a holder at the boundary. Return it.
    DataPtrType thisDataPtr = it->second;
    CH_assert(!thisDataPtr.isNull());
    return *thisDataPtr;
}


// -----------------------------------------------------------------------------
// Writes boundary data to a holder that can be sent to HDF5.
// -----------------------------------------------------------------------------
void BoundaryData<Real>::copyTo (LevelData<FluxBox>& a_dest) const
{
    // Sanity checks
    CH_assert(a_dest.isDefined());
    CH_assert(a_dest.getBoxes().compatible(m_grids));
    CH_assert(a_dest.getBoxes().physDomain() == m_grids.physDomain());

    // Loop over directions, sides, and grids.
    for (int dir = 0; dir < SpaceDim; ++dir) {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const_iterator it;
            for (it = begin(dir,sit()); it != end(dir,sit()); ++it) {
                const DataIndex di = it->first;
                const FArrayBox& srcFAB = *(it->second);
                FArrayBox& destFAB = a_dest[di][dir];

                CH_assert(destFAB.box().type() == srcFAB.box().type());
                const Box overlap = destFAB.box() & srcFAB.box();

                destFAB.copy(srcFAB, overlap);

            } // end loop over grids (it)
        } // end loop over sides (sit)
    } // end loop over directions (dir)
}


// -----------------------------------------------------------------------------
// Copies boundary values. Grids must be compatible.
// -----------------------------------------------------------------------------
void BoundaryData<Real>::copy (const LevelData<FluxBox>& a_src)
{
    // Sanity checks
    CH_assert(a_src.isDefined());
    CH_assert(a_src.getBoxes().compatible(m_grids));
    CH_assert(a_src.getBoxes().physDomain() == m_grids.physDomain());

    // Loop over grids, directions, and sides.
    DataIterator dit = a_src.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        for (int dir = 0; dir < SpaceDim; ++dir) {
            SideIterator sit;
            for (sit.reset(); sit.ok(); ++sit) {
                const int s = int(sit());
                BoundaryDataMapType& thisMap = m_bdryFaceMap[dir][s];

                // Find the DataIndex in the map.
                iterator it;
                it = thisMap.find(dit());

                // If we didn't find anything, move on.
                if (it == thisMap.end()) continue;

                // Looks like we are at a physical boundary. Copy what we can.
                DataPtrType& destFABPtr = it->second;
                const FArrayBox& srcFAB = a_src[dit][dir];

                const Box& region = destFABPtr->box();
                CH_assert(srcFAB.box().contains(region));

                destFABPtr->copy(srcFAB, region);

            } // end loop over sides (sit)
        } // end loop over directions (dir)
    } // end loop over grids (dit)
}


// -----------------------------------------------------------------------------
// Sets all data to a uniform value.
// -----------------------------------------------------------------------------
void BoundaryData<Real>::setVal (Real a_val)
{
    int numDirs = SpaceDim;
    if (isFlat()) --numDirs;

    for (int dir = 0; dir < numDirs; ++dir) {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const int s = int(sit());
            BoundaryDataMapType& thisMap = m_bdryFaceMap[dir][s];

            iterator it = thisMap.begin();
            const iterator endIt = thisMap.end();

            while (it != endIt) {
                it->second->setVal(a_val);
                ++it;
            }
        } // end loop over sides (sit)
    } // end loop over directions (dir)
}


// -----------------------------------------------------------------------------
// Gathers boundary fluxes from a BCFluxClass.
// -----------------------------------------------------------------------------
void BoundaryData<Real>::collectFluxes (const BCFluxClass& a_bc,
                                        const Real         a_time)
{
    const ProblemDomain& domain = m_grids.physDomain();
    const RealVect dx = RealVect::Zero; // Just a dummy.

    for (int dir = 0; dir < SpaceDim; ++dir) {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const int s = int(sit());
            BoundaryDataMapType& thisMap = m_bdryFaceMap[dir][s];

            iterator it = thisMap.begin();
            for (; it != thisMap.end(); ++it) {
                FArrayBox& thisFAB = *(it->second);
                const DataIndex& di = it->first;
                const Box& valid = m_grids[di];

                a_bc(thisFAB, NULL, valid, domain, dx,
                     di, NULL, dir, false, a_time);

            } // end loop over map elements (it)
        } // end loop over sides (sit)
    } // end loop over directions (dir)
}


// -----------------------------------------------------------------------------
// Accumulates values from another BoundaryData object.
// -----------------------------------------------------------------------------
const BoundaryData<Real>& BoundaryData<Real>::plus (const BoundaryData<Real>& a_src,
                                                    const Real                a_scale)
{
    // Sanity checks
    CH_assert(this->sameLayout(a_src));

    for (int dir = 0; dir < SpaceDim; ++dir) {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide iside = sit();

            iterator thisIt = this->begin(dir, iside);
            const_iterator srcIt = a_src.begin(dir, iside);
            const iterator endIt = this->end(dir, iside);

            while (thisIt != endIt) {
                DataPtrType& thisFABPtr = thisIt->second;
                const DataPtrType& srcFABPtr = srcIt->second;

                CH_assert(thisFABPtr->box() == srcFABPtr->box());
                thisFABPtr->plus(*srcFABPtr, a_scale);

                ++thisIt;
                ++srcIt;
            }
        } // end loop over sides (sit)
    } // end loop over directions (dir)

    return *this;
}


// -----------------------------------------------------------------------------
// Accumulates values from a flat BoundaryData object.
// This only operates on the horizontal BCs.
// -----------------------------------------------------------------------------
const BoundaryData<Real>& BoundaryData<Real>::plusFlat (const BoundaryData<Real>& a_flatSrc,
                                                        const Real                a_scale)
{
    // Sanity checks
    CH_assert(!this->isFlat());
    CH_assert(a_flatSrc.isFlat());
    CH_assert(this->flattenedLayout(a_flatSrc));

    for (int dir = 0; dir < SpaceDim-1; ++dir) {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const Side::LoHiSide iside = sit();

            iterator thisIt = this->begin(dir, iside);
            const_iterator srcIt = a_flatSrc.begin(dir, iside);
            const iterator endIt = this->end(dir, iside);

            while (thisIt != endIt) {
                // Get this object's info (this will be modified)
                DataPtrType& thisFABPtr = thisIt->second;
                const Box& thisBox = thisFABPtr->box();
                const int loZ = thisBox.smallEnd(SpaceDim-1);
                const int hiZ = thisBox.bigEnd(SpaceDim-1);

                // Get the src info (Not modified. This is the increment.)
                const DataPtrType& srcFABPtr = srcIt->second;
                const Box& srcBox = srcFABPtr->box();
                CH_assert(srcBox.smallEnd(SpaceDim-1) == 0);
                CH_assert(srcBox.bigEnd(SpaceDim-1) == 0);

                // Loop over each element of the flattened box, get the increment,
                // and apply it to entire vertical lines.
                BoxIterator bit(srcBox);
                for (bit.reset(); bit.ok(); ++bit) {
                    const IntVect& flatIV = bit();
                    Real srcVal = (*srcFABPtr)(flatIV);

                    Box destBox(flatIV, flatIV, thisBox.type());
                    destBox.shift(SpaceDim-1, loZ);
                    destBox.setBig(SpaceDim-1, hiZ);
                    CH_assert(thisBox.sameType(destBox));
                    CH_assert(thisBox.contains(destBox));

                    thisFABPtr->plus(srcVal*a_scale, destBox);
                }

                ++thisIt;
                ++srcIt;
            }
        } // end loop over sides (sit)
    } // end loop over directions (dir)

    return *this;
}


// -----------------------------------------------------------------------------
// Accumulates values from a flat BoxLayoutData.
// This only operates on the vertical BCs.
// -----------------------------------------------------------------------------
const BoundaryData<Real>& BoundaryData<Real>::vertPlus (const BoxLayoutData<FArrayBox>& a_flatSrc,
                                                        const Real                      a_scale,
                                                        const Side::LoHiSide            a_side)
{
    CH_assert(!this->isFlat());
    CH_assert(m_grids.compatible(a_flatSrc.boxLayout()));

    DataIterator dit = m_grids.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& bdryFAB = this->getData(dit(), SpaceDim-1, a_side);
        const FArrayBox& srcFAB = a_flatSrc[dit];

        CH_assert(srcFAB.box().size(SpaceDim-1) == 1);

        VertShifter<Real> shifty(bdryFAB, srcFAB);
        bdryFAB.plus(srcFAB, a_scale);
        shifty.restore();
    }

    return *this;
}


// -----------------------------------------------------------------------------
// Converts the horizontal boundary data to a flattened vertical average.
// The vertical boundary data is deallocated.
// -----------------------------------------------------------------------------
void BoundaryData<Real>::vertAvg ()
{
    // Sanity check
    CH_assert(!this->isFlat());
    CH_assert(this->flattenable());

    // Deallocate the vertical boundary data.
    m_bdryFaceMap[SpaceDim-1][0].clear();
    m_bdryFaceMap[SpaceDim-1][1].clear();

    // Create a new map.
    BoundaryDataMapType newBdryFaceMap[CH_SPACEDIM][2];

    for (int dir = 0; dir < SpaceDim-1; ++dir) {
        SideIterator sit;
        for (sit.reset(); sit.ok(); ++sit) {
            const int s = int(sit());
            BoundaryDataMapType& newMap = newBdryFaceMap[dir][s];
            BoundaryDataMapType& oldMap = m_bdryFaceMap[dir][s];

            const_iterator oldIt = oldMap.begin();
            for (; oldIt != oldMap.end(); ++oldIt) {
                const DataIndex& oldDi = oldIt->first;
                const FArrayBox& oldFAB = *(oldIt->second);
                const Box& oldBox = oldFAB.box();
                const int loZ = oldBox.smallEnd(SpaceDim-1);
                const int hiZ = oldBox.bigEnd(SpaceDim-1);
                const int vertSize = oldBox.size(SpaceDim-1);

                Box newBox = flattenBox(oldBox, SpaceDim-1);
                FArrayBox* newFABPtr = new FArrayBox(newBox, 1);
                newMap[oldDi] = RefCountedPtr<FArrayBox>(newFABPtr);

                BoxIterator bit(newBox);
                for (bit.reset(); bit.ok(); ++bit) {
                    const IntVect& flatIV = bit();

                    Box lineBox(flatIV, flatIV, oldBox.type());
                    lineBox.shift(SpaceDim-1, loZ);
                    lineBox.setBig(SpaceDim-1, hiZ);
                    CH_assert(oldBox.contains(lineBox));

                    (*newFABPtr)(flatIV) = oldFAB.sum(lineBox, 0, 1) / Real(vertSize);
                }
            } // end loop over old FABs (oldIt)

            // Deallocate old map and install new one.
            oldMap = newMap;

        } // end loop over sides (sit)
    } // end loop over directions (dir)

    m_isFlat = true;
}
