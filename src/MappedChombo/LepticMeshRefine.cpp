#include "LepticMeshRefine.H"
#include "AnisotropicRefinementTools.H"
#include "RealVect.H"
#include "MergeBoxesOnLines.H"


// -----------------------------------------------------------------------------
// Default constructor -- leaves object in an unusable state
// -----------------------------------------------------------------------------
LepticMeshRefine::LepticMeshRefine ()
{
    // Do nothing special.
    // The base class default constructors will be called automatically.
}


// -----------------------------------------------------------------------------
// Full constructor -- leaves object in usable state
// -----------------------------------------------------------------------------
LepticMeshRefine::LepticMeshRefine (const Box&             a_baseDomain,    // Level 0 domain
                                    const Vector<IntVect>& a_refRatios,     // Refinement ratios -- refRatio[0] is btwn levels 0 and 1
                                    const Real             a_fillRatio,     // Measure of how efficiently tagged cells will be covered
                                    const int              a_blockFactor,   // Amount by which grids are guaranteed to be coarsenable
                                    const int              a_bufferSize,    // Proper nesting buffer amount
                                    const IntVect&         a_maxSize,       // Maximum grid length in any direction -- 0 means no limit.
                                    const IntVect&         a_spanDirs)      // Set to 1 for new boxes to span the dim.
{
    const ProblemDomain crseDom(a_baseDomain);

    this->define(crseDom,
                 a_refRatios,
                 a_fillRatio,
                 a_blockFactor,
                 a_bufferSize,
                 a_maxSize,
                 a_spanDirs);
}


// -----------------------------------------------------------------------------
// Full constructor -- leaves object in usable state
// -----------------------------------------------------------------------------
LepticMeshRefine::LepticMeshRefine (const ProblemDomain&   a_baseDomain,    // Level 0 domain
                                    const Vector<IntVect>& a_refRatios,     // Refinement ratios -- refRatio[0] is btwn levels 0 and 1
                                    const Real             a_fillRatio,     // Measure of how efficiently tagged cells will be covered
                                    const int              a_blockFactor,   // Amount by which grids are guaranteed to be coarsenable
                                    const int              a_bufferSize,    // Proper nesting buffer amount
                                    const IntVect&         a_maxSize,       // Maximum grid length in any direction -- 0 means no limit.
                                    const IntVect&         a_spanDirs)      // Set to 1 for new boxes to span the dim.
{
    this->define(a_baseDomain,
                 a_refRatios,
                 a_fillRatio,
                 a_blockFactor,
                 a_bufferSize,
                 a_maxSize,
                 a_spanDirs);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
LepticMeshRefine::~LepticMeshRefine ()
{
    // Do nothing special.
    // The base class destructors will be called automatically.
}


// -----------------------------------------------------------------------------
// Define function -- size of RefRatios will define maximum number of levels
// -----------------------------------------------------------------------------
void LepticMeshRefine::define (const Box&             a_baseDomain,     // Level 0 domain
                               const Vector<IntVect>& a_refRatios,      // Refinement ratios -- refRatio[0] is btwn levels 0 and 1
                               const Real             a_fillRatio,      // Measure of how efficiently tagged cells will be covered
                               const int              a_blockFactor,    // Amount by which grids are guaranteed to be coarsenable
                               const int              a_bufferSize,     // Proper nesting buffer amount
                               const IntVect&         a_maxSize,        // Maximum grid length in any direction -- 0 means no limit.
                               const IntVect&         a_spanDirs)       // Set to 1 for new boxes to span the dim.
{
    const ProblemDomain crseDom(a_baseDomain);

    this->define(crseDom,
                 a_refRatios,
                 a_fillRatio,
                 a_blockFactor,
                 a_bufferSize,
                 a_maxSize,
                 a_spanDirs);
}


// -----------------------------------------------------------------------------
// Define function -- size of RefRatios will define maximum number of levels
// This is the function that all other full constructors call.
// -----------------------------------------------------------------------------
void LepticMeshRefine::define (const ProblemDomain&   a_baseDomain,     // Level 0 domain
                               const Vector<IntVect>& a_refRatios,      // Refinement ratios -- refRatio[0] is btwn levels 0 and 1
                               const Real             a_fillRatio,      // Measure of how efficiently tagged cells will be covered
                               const int              a_blockFactor,    // Amount by which grids are guaranteed to be coarsenable
                               const int              a_bufferSize,     // Proper nesting buffer amount
                               const IntVect&         a_maxSize,        // Maximum grid length in any direction -- 0 means no limit.
                               const IntVect&         a_spanDirs)       // Set to 1 for new boxes to span the dim.
{
    // Clobber MeshRefine variables that we are overriding.
    MeshRefine::m_maxSize = -1;
    MeshRefine::m_nRefVect.resize(0);
    MeshRefine::m_level_blockfactors.resize(0);

    // Change default PND mode from 0 to 1 (ndk 8.4.2008)
    //m_PNDMode = 0; // old, tested version of using nestingRegion
    m_PNDMode = 1; // New version which does not use expensive IntVectSet nestingRegion
    int maxLevel = a_refRatios.size();
    m_nRefVect.resize(maxLevel);
    m_vectDomains.resize(maxLevel+1);
    m_pnds.resize(maxLevel);
    //m_lastBase = maxLevel+1;
    //m_lastTop  = 0;
    //m_lastBuffer = 0;
    m_level_blockfactors.resize(maxLevel);

    m_vectDomains[0] = a_baseDomain;
    for (int lev = 0; lev < maxLevel; ++lev) {
        m_nRefVect[lev] = a_refRatios[lev];
        m_vectDomains[lev+1] = refine(m_vectDomains[lev], a_refRatios[lev]);
    }

    // do some quick sanity checks
    CH_assert( a_blockFactor >= 1 );
    CH_assert( a_bufferSize >= 0 );
    D_TERM(CH_assert( (a_maxSize[0] >= 2*a_bufferSize) || (a_maxSize[0] == 0) );,
           CH_assert( (a_maxSize[1] >= 2*a_bufferSize) || (a_maxSize[1] == 0) );,
           CH_assert( (a_maxSize[2] >= 2*a_bufferSize) || (a_maxSize[2] == 0) );)
    D_TERM(CH_assert( (a_blockFactor <= a_maxSize[0]) || (a_maxSize[0] == 0) );,
           CH_assert( (a_blockFactor <= a_maxSize[1]) || (a_maxSize[1] == 0) );,
           CH_assert( (a_blockFactor <= a_maxSize[2]) || (a_maxSize[2] == 0) );)
    CH_assert( a_fillRatio > 0.0 || a_fillRatio <= 1.0 );

    D_TERM(CH_assert(a_spanDirs[0] == 0 || a_spanDirs[0] == 1);,
           CH_assert(a_spanDirs[1] == 0 || a_spanDirs[1] == 1);,
           CH_assert(a_spanDirs[2] == 0 || a_spanDirs[2] == 1);)

    m_fillRatio = a_fillRatio;
    m_blockFactor = a_blockFactor;
    m_bufferSize = a_bufferSize;
    m_maxSize = a_maxSize * (IntVect::Unit - a_spanDirs); // If spanDirs is set, maxSize should be zero.
    m_spanDirs = a_spanDirs;

    computeLocalBlockFactors();

    m_isDefined = true;
}


// -----------------------------------------------------------------------------
// Create hierarchy of grids from a single level of tags.
//  This function creates a hierarchy of grids from a single level of
//  tags on BaseLevel. If tags exist, then all levels will have grids.
//  Returns the new finest level of grids.
// -----------------------------------------------------------------------------
int LepticMeshRefine::regrid (Vector<Vector<Box> >&       a_newmeshes,  // new set of grids at every level
                              const IntVectSet&           a_tags,       // tagged cells on baseLevel
                              const int                   a_baseLevel,  // index of base mesh level (finest unchanged level)
                              const int                   a_topLevel,   // top level to refine (one less than finest possible level)
                              const Vector<Vector<Box> >& a_oldMeshes)  // existing grids (if no previous grids, set to domains)
{
    MayDay::Error("LepticMeshRefine::regrid 1 not yet written");
    return -1; // Just to shut the compiler up.
}


// -----------------------------------------------------------------------------
// Create hierarchy of grids from tags at all levels.
//  This function creates a hierarchy of grids from tags at all
//  refinement levels.  It is possible that not all levels will
//  return with grids, since there may not be tags at all levels.
//  Returns the new finest level of grids.
// -----------------------------------------------------------------------------
int LepticMeshRefine::regrid (Vector<Vector<Box> >&       a_newmeshes,  // new set of grids at every level
                              Vector<IntVectSet>&         a_tags,       // tagged cells on each existing level
                              const int                   a_baseLevel,  // index of base mesh level (finest unchanged level)
                              const int                   a_topLevel,   // top level to refine (one less than finest possible level)
                              const Vector<Vector<Box> >& a_OldMeshes)  // existing grids (if no previous grids, set to domains)
{
    CH_TIME("LepticMeshRefine::regrid");
    //
    // Validate arguments and handle special cases
    //
    CH_assert( isDefined());
    CH_assert( a_topLevel >= 0 );
    CH_assert( a_baseLevel < (a_topLevel+1) && a_baseLevel >= 0 );
    CH_assert( a_OldMeshes.size() >= a_topLevel + 1 );
    CH_assert( m_vectDomains.size() >= a_topLevel + 1 );
    CH_assert( m_nRefVect.size() >= a_topLevel + 1 );       // TODO: Relax this restriction.
    CH_assert( a_tags.size() >= a_topLevel+1 );

    // Span whatever directions we need to.
    const Vector<Vector<Box> >* oldMeshesPtr = &a_OldMeshes;
    Vector<Vector<Box> > stretchedOldMeshes(0);
    if (m_spanDirs.sum() > 0) {
        stretchedOldMeshes = a_OldMeshes;
        for (int lev = a_baseLevel; lev <= a_topLevel; ++lev) {
            // pout() << "Old level " << lev << " flags:\n";
            bool recheck = LepticMeshRefine::spanBoxes(stretchedOldMeshes[lev],
                                                       m_vectDomains[lev].domainBox(),
                                                       m_spanDirs);
            if (recheck) --lev;
        }

        oldMeshesPtr = &stretchedOldMeshes;
    }

    const Vector<Vector<Box> >& oldMeshesRef = *oldMeshesPtr;
    Vector<IntVectSet>& modifiedTags = a_tags;

    // set the top level to be the finest level which actually has tags
    int TopLevel = a_topLevel;
    int new_finest_level;
    int isize = a_tags.size();
    int level;
    for (level = Min(a_topLevel, isize-1); level >= a_baseLevel; --level) {
        if (!a_tags[level].isEmpty()) break;
    }

#ifdef CH_MPI
    int mlevel;
    MPI_Allreduce(&level, &mlevel, 1, MPI_INT, MPI_MAX, Chombo_MPI::comm);
    level=mlevel;
#endif

    if (level >= a_baseLevel) {
        TopLevel = level;
        Vector<IntVect> totalBufferSize(a_tags.size(), IntVect::Zero);

        // reinitialize the output array .... (noel: why +2 ?)
        a_newmeshes.resize( TopLevel+2 );

        //
        // Generate new meshes if requested.
        if ( TopLevel+1 > a_baseLevel ) {
            IntVect vectBF;
            D_TERM(vectBF[0] = (m_spanDirs[0] == 1)? 1: m_blockFactor;,
                   vectBF[1] = (m_spanDirs[1] == 1)? 1: m_blockFactor;,
                   vectBF[2] = (m_spanDirs[2] == 1)? 1: m_blockFactor;)

            Box domaint = m_vectDomains[a_baseLevel].domainBox();
            Box testdom = coarsen(domaint, vectBF);
            testdom.refine(vectBF);
            if (domaint != testdom)
                MayDay::Error("LepticMeshRefine:domain and Blocking Factor incompatible");
            // [NOTE: the following validations are expensive so only do them when
            //        debugging.]

#if ! defined( NDEBUG )
            // all existing meshes should be within the boundaries
            for ( int i = a_baseLevel ; i <= TopLevel ; i++ ) {
                if (oldMeshesRef[i].size() > 0) {
                    Box minbox = oldMeshesRef[i][0];
                    for (int ibox = 0; ibox < oldMeshesRef[i].size(); ++ibox) {
                        minbox.minBox(oldMeshesRef[i][ibox]);
                    }
                    CH_assert(m_vectDomains[i].contains(minbox));
                }
            }
#endif
            // all tagged cells must be within the existing meshes
            // make sure that each level of \var{Domains} is consistent with the
            // previous level and the refinement ratio.
            for ( int i = a_baseLevel ; i < TopLevel ; i++ ) {
                if (m_vectDomains[i+1].domainBox() != refine(m_vectDomains[i].domainBox(), m_nRefVect[i]))
                    MayDay::Error("LepticMeshRefine:domains and refratios incompatible");
            }

            //
            // coarsen old meshes, tags, problem domains, and buffers by the
            // appropriate blocking factor
            Vector<Vector<Box> > OldMeshes = oldMeshesRef;
            Vector<ProblemDomain> Domains = m_vectDomains;
            Vector<IntVect> blocked_BufferSize(TopLevel+1);

            for (int level = a_baseLevel; level <= TopLevel; ++level) {
                // first intersect tags with domains (note that in periodic case, this can
                // result in wrapping of tags that are outside the base domain box)
                // int nBefore = modifiedTags[level].numPts();
                // IntVectSet tmp(modifiedTags[level]);
                //tmp &= m_vectDomains[level].domainBox();
                //int nTmp = tmp.numPts();
                modifiedTags[level] &= m_vectDomains[level];
                //int nAfter = modifiedTags[level].numPts();
                //pout()<<"nBefore:"<<nBefore<<" nAfter:"<<nAfter<<" nTmp:"<<nTmp
                //            <<" diff:"<<nBefore-nAfter<<" newDiff:"<<nTmp-nAfter <<"\n";
                // coarsen by blocking factor, rounding upwards the goal of this is to
                // coarsen everything down to a level which is m_blockFactor coarser than
                // the new fine level we're going to generate.  By generating grids at
                // this fake level, we can then refine up to the new level, and the blocking
                // factor will be automatically enforced.

                // TODO: This is a bottleneck!
                // CH_assert(   m_level_blockfactors[level][0] == m_level_blockfactors[level][1]
                //           && m_level_blockfactors[level][0] == m_level_blockfactors[level][2]   );
{
CH_TIME("Coarsen_modifiedTags");
coarsen(modifiedTags[level], m_level_blockfactors[level]);
// modifiedTags[level].coarsen(m_level_blockfactors[level][0]);

// int lvl = level;
// pout() << "m_vectDomains[" << lvl << "] = " << flush;
// pout() << m_vectDomains[lvl] << endl;

// pout() << "m_level_blockfactors[" << lvl << "] = " << flush;
// pout() << m_level_blockfactors[lvl] << endl;

// pout() << "modifiedTags[" << lvl << "].numPts() = " << flush;
// pout() << modifiedTags[lvl].numPts() << endl;
}

                // Same block-factor coarsening for the domains
                Box newDomBox(Domains[level].domainBox());
                newDomBox.coarsen(m_level_blockfactors[level]);
                bool isPeriodic[CH_SPACEDIM];
                D_TERM(isPeriodic[0] = Domains[level].isPeriodic(0);,
                       isPeriodic[1] = Domains[level].isPeriodic(1);,
                       isPeriodic[2] = Domains[level].isPeriodic(2);)
                Domains[level].define(newDomBox, isPeriodic);

                // Same block-factor coarsening for the buffers
                blocked_BufferSize[level] =
                    (  m_level_blockfactors[level]
                     + IntVect(D_DECL(m_bufferSize-1,m_bufferSize-1,m_bufferSize-1))  )
                    / m_level_blockfactors[level];
            }
            // We only need the base mesh coarsened by the blocking factor
            Vector<Box> OldBaseMesh = oldMeshesRef[a_baseLevel];
            {
                const IntVect crFactor = m_level_blockfactors[a_baseLevel];
                for (int i = 0; i != OldBaseMesh.size(); ++i) {
                    OldBaseMesh[i].coarsen(crFactor);
                }
            }

            for (int i=a_baseLevel; i<=TopLevel; i++) {
                m_pnds[i].makeEmpty();
            }
            this->makePNDs(m_pnds, totalBufferSize, a_baseLevel, TopLevel,
                           Domains, OldBaseMesh, blocked_BufferSize);

            // Clip the tagged cells to the proper nesting domains by
            // intersecting the two sets.
            // Note: if m_PNDMode = 1, this only ensures the tags are in the base
            // PND!  Remaining adherence to each level's PND is performed in the
            // box generation algorithm
            for ( int lvl = a_baseLevel ; lvl <= TopLevel ; lvl++ ) {
                modifiedTags[lvl] &= m_pnds[lvl] ;
            }

            //
            // Generate new meshes.
            //

            // At each level, starting at the top, generate boxes that cover the tag
            // cells using makeBoxes(), refine these boxes and save them in the
            // output \var{a_newmeshes}.  Take the unrefined boxes, add a buffer zone
            // around each box, coarsen, and union with the tag cells on the next
            // coarser level.  This modifies the tags variable.  To handle
            // \var{BlockFactor}, coarsen everything before making the new
            // meshes and then refine the resulting mesh boxes.
            Vector<Box> lvlboxes ;  // new boxes on this level
            for ( int lvl = TopLevel ; lvl >= a_baseLevel ; lvl-- ) {
                // make a new mesh at the same level as the tags

                const int dest_proc = uniqueProc(SerialTask::compute);

                Vector<IntVectSet> all_tags;
                gather(all_tags, modifiedTags[lvl], dest_proc);

                if (procID() == dest_proc) {
                    for (int i = 0; i < all_tags.size(); ++i) {
                        //                     modifiedTags[lvl] |= all_tags[i];
                        //**FIXME -- revert to above line when IVS is fixed.
                        //**The following works around a bug in IVS that appears if
                        //**the above line is used.  This bug is observed when there
                        //**is a coarsening of an IVS containing only IntVect::Zero
                        //**followed by an IVS |= IVS.
                        for (IVSIterator ivsit(all_tags[i]); ivsit.ok(); ++ivsit) {
                            modifiedTags[lvl] |= ivsit();
                        }
                        //**FIXME -- end
                        // Regain memory used (BVS,NDK 6/30/2008)
                        all_tags[i].makeEmpty();
                    }
                }

                broadcast( modifiedTags[lvl] , dest_proc);

                // Move this union _after_ the above gather/broadcast to
                // reduce memory -- shouldn't have other effects. (BVS,NDK 6/30/2008)
                // Union the meshes from the previous level with the tags on this
                // level to guarantee that proper nesting is satisfied.  On the
                // first iteration this is a no-op because \var{lvlboxes} is empty.
                // [NOTE: for every iteration after the first, \var{lvlboxes} will
                //        already be coarsened by \var{BlockFactor} so it will be
                //        at the same refinement level as \var{tags[lvl]}, which
                //        has also been coarsened]
                // this is simple in the non-periodic case, more complicated
                // in the periodic case
                ProblemDomain lvldomain = Domains[lvl]; // domain of this level
                if (lvldomain.isPeriodic() ) {
                    const Box domainBox = lvldomain.domainBox();
                    ShiftIterator shiftIt = lvldomain.shiftIterator();
                    IntVect shiftMult(domainBox.size());
                    for (int i=0; i<lvlboxes.size(); i++) {
                        Box localBox(lvlboxes[i]);
                        // will handle periodic wraparound through shifting and
                        // adding shifted image to tags, which will all remain
                        // within the domainBox
                        localBox &= domainBox;
                        modifiedTags[lvl] |= localBox;
                        // now do shifts to capture periodic images necessary to
                        // enforce proper nesting only do this if original box was
                        // not contained by domainBox
                        if (localBox != lvlboxes[i]) {
                            for (shiftIt.begin(); shiftIt.ok(); ++shiftIt) {
                                IntVect shiftVect(shiftIt()*shiftMult);
                                Box localShiftedBox(lvlboxes[i]);
                                localShiftedBox.shift(shiftVect);
                                localShiftedBox &= domainBox;
                                if (!localShiftedBox.isEmpty()) {
                                    modifiedTags[lvl] |= localShiftedBox;
                                }
                            } // end loop over shift directions
                        } // end whether periodic checking was needed
                    } // end loop over finer-level boxes to enforce nesting
                } else {
                    // non periodic case is simple
                    for ( int i = 0 ; i < lvlboxes.size() ; i++ ) {
                        modifiedTags[lvl] |= lvlboxes[i] ;
                    }
                }

                // this is the maximum allowable box size at this resolution
                // which will result in satisfying the maxSize restriction when
                // everything is refined up to the new level
                const IntVect maxBoxSizeLevel = m_maxSize/(m_level_blockfactors[lvl]*m_nRefVect[lvl]);
                this->makeBoxes(lvlboxes, modifiedTags[lvl], m_pnds[lvl],
                                lvldomain, maxBoxSizeLevel, totalBufferSize[lvl]);

                // This ensures the m_spanDirs requirements.
                for (int d = 0; d < SpaceDim; ++d) {
                    if (m_spanDirs[d] == 0) continue;

                    for (int i = 0; i < lvlboxes.size(); ++i) {
                        Box& b = lvlboxes[i];
                        b.shift(d, lvldomain.domainBox().smallEnd(d) - b.smallEnd(d));
                        b.setBig(d, lvldomain.domainBox().bigEnd(d));
                    }

                    MergeBoxesOnLines().mergeBoxes(lvlboxes, d);
                }

                // After change to reduce memory, this may now be needed.
                // Previously, there were a_tags.makeEmpty() calls in BRMesh.cpp, and now,
                // if there are a few tags leftover here, they will get added onto the mix -- which is not
                // the behavior of the code before these changes (ndk) (leaving commented for now)
                modifiedTags[lvl].makeEmpty(); // <-- why is this necessary now?

                // refine the new mesh and save it
                //[NOTE: we have to undo the coarsening by \var{BlockFactor} as well
                //       as refine to the next level.]
                //[NOTE: refine() operates in-place so copy first then refine()
                //       because the unrefined mesh will be needed later.]
                a_newmeshes[lvl+1] = lvlboxes ;
                for (int ibox = 0; ibox < a_newmeshes[lvl+1].size(); ++ibox) {
                    a_newmeshes[lvl+1][ibox].refine(m_level_blockfactors[lvl]*m_nRefVect[lvl]);
                }
                // Make the boxes ready for the next iteration.
                // Don't have to do this for the last iteration.
                if ( lvl > a_baseLevel ) {
                    // coarsen the unrefined new mesh so it matches the tags
                    // at that level, then add the buffer cells, clip at the domain
                    // boundaries so we can union them in the next iteration
                    // allInOne_nRef does:
                    //   a) refine  by BF_ref[lvl];
                    //   b) coarsen by  n_ref[lvl-1];
                    //   c) coarsen by BF_ref[lvl-1];
                    IntVect allInOne_nRef = m_nRefVect[lvl-1];
                    allInOne_nRef *= m_level_blockfactors[lvl-1];
                    allInOne_nRef /= m_level_blockfactors[lvl];

                    for (int ibox = 0 ; ibox < lvlboxes.size() ; ++ibox) {
                        lvlboxes[ibox].grow(blocked_BufferSize[lvl]); // Growing in vertical should have no consequence...
                        lvlboxes[ibox] &= Domains[lvl];               // ...because of this line.
                        lvlboxes[ibox].coarsen(allInOne_nRef);
                    }
                }
            } // end loop over levels
        } // end if TopLevel+1 > baseLevel

        //
        // Finally, copy the old mesh levels that didn't change.
        //
        for ( int i = 0 ; i <= a_baseLevel ; i++ ) {
            a_newmeshes[i] = oldMeshesRef[i] ;
        }

        //
        // Done generating grids
        //

        // set new finest level
        new_finest_level = TopLevel+1;
        // this is designed to catch the pathological but possible case
        // where there were tags on the TopLevel, but no grids were generated
        // (possibly if all tags were outside pnds)
        while (a_newmeshes[new_finest_level].size() == 0
               && new_finest_level > a_baseLevel) {
            new_finest_level -= 1;
        }

    } else {
        // if no tags on any level, just return
        new_finest_level = a_baseLevel;
        a_newmeshes.resize(oldMeshesRef.size());
        for ( int i = 0 ; i <= a_baseLevel ; i++ ) {
            a_newmeshes[i] = oldMeshesRef[i] ;
        }
    }

    // Span whatever directions we need to.
    if (m_spanDirs.sum() > 0) {
        for (int dir = 0; dir < SpaceDim; ++dir) {
            if (m_spanDirs[dir] != 0) {
                for (int lev = a_baseLevel + 1; lev <= a_topLevel; ++lev) {
                    MergeBoxesOnLines().mergeBoxes(a_newmeshes[lev], dir);
                }
            }
        }
        // for (int lev = a_baseLevel; lev <= a_topLevel; ++lev) {
        //     // pout() << "Level " << lev << " flags:\n";
        //     bool recheck = LepticMeshRefine::spanBoxes(a_newmeshes[lev],
        //                                                m_vectDomains[lev].domainBox(),
        //                                                m_spanDirs);
        //     if (recheck) --lev;
        // }
    }

    return new_finest_level;
}


// -----------------------------------------------------------------------------
// This streches boxes to conform to m_spanDir and fixes any overlaps
// that arise. This function returns true if boxes were removed or fixed to
// avoid an overlap.
// -----------------------------------------------------------------------------
bool LepticMeshRefine::spanBoxes (Vector<Box>&   a_boxes,
                                  const Box&     a_domBox,
                                  const IntVect& a_spanDirs)
{
    CH_TIME("LepticMeshRefine::spanBoxes");

    // Sanity checks
    CH_assert(!a_domBox.isEmpty());
    D_TERM(CH_assert(a_spanDirs[0] == 0 || a_spanDirs[0] == 1);,
           CH_assert(a_spanDirs[1] == 0 || a_spanDirs[1] == 1);,
           CH_assert(a_spanDirs[2] == 0 || a_spanDirs[2] == 1);)

    // Are there any directions that need stretching?
    if (a_spanDirs.sum() == 0) return false;

    // Is this level empty?
    int numBoxes = a_boxes.size();
    if (numBoxes == 0) return false;


    // // Stretch all boxes
    // Vector<Box> oldBoxes = a_boxes;
    // for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
    //     // If we aren't stretching boxes in this dir, move on.
    //     if (a_spanDirs[dir] == 0) continue;

    //     const int loEnd = a_domBox.smallEnd(dir);
    //     const int hiEnd = a_domBox.bigEnd(dir);

    //     // Strectch every box at this level
    //     for (int idx = 0; idx < numBoxes; ++idx) {
    //         Box& b = a_boxes[idx];
    //         b.shift(dir, loEnd-b.smallEnd(dir));
    //         b.setBig  (dir, hiEnd);
    //     } // end loop over stretched boxes (idx)

    //     MergeBoxesOnLines().mergeBoxes(a_boxes, dir);

    // } // end loop over stretching directions (dir)

    // // Has anything changed?
    // if (a_boxes.size() != oldBoxes.size()) {
    //     return true;
    // }
    // for (int idx = 0; idx < a_boxes.size(); ++idx) {
    //     if (a_boxes[idx] != oldBoxes[idx]) {
    //         return true;
    //     }
    // }
    // return false;


    // Stretch all boxes
    for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
        // If we aren't stretching boxes in this dir, move on.
        if (a_spanDirs[dir] == 0) continue;

        const int loEnd = a_domBox.smallEnd(dir);
        const int hiEnd = a_domBox.bigEnd(dir);

        // Strectch every box at this level
        for (int idx = 0; idx < numBoxes; ++idx) {
            Box& b = a_boxes[idx];
            b.shift(dir, loEnd-b.smallEnd(dir));
            b.setBig  (dir, hiEnd);
        } // end loop over stretched boxes (idx)
    } // end loop over stretching directions (dir)

    // Search for boxes that can be removed.
    enum {KEEP = 0, REMOVE, FIX};
    Vector<int> flags(numBoxes, KEEP);
    int needsRemoval = 0;

    for (int i1 = 0; i1 < numBoxes; ++i1) {
        const Box& b1 = a_boxes[i1];
        for (int i2 = i1 + 1; i2 < numBoxes; ++i2) {
            const Box& b2 = a_boxes[i2];

            if (!b1.intersects(b2)) continue;

            if (b1.contains(b2)) {
                if (flags[i2] != REMOVE) {
                    flags[i2] = REMOVE;
                    ++needsRemoval;
                }
                continue;
            }

            if (b2.contains(b1)) {
                if (flags[i1] != REMOVE) {
                    flags[i1] = REMOVE;
                    ++needsRemoval;
                }
                continue;
            }
        } // end loop over i2
    } // end loop over i1

    // for (int i1 = 0; i1 < numBoxes; ++i1) {
    //     pout() << "\tBox[" << i1 << "]: ";
    //     if (flags[i1] == KEEP) {
    //         pout() << "KEEP\n";
    //     } else if (flags[i1] == REMOVE) {
    //         pout() << "REMOVE\n";
    //     } else if (flags[i1] == FIX) {
    //         pout() << "FIX\n";
    //     } else {
    //         MayDay::Error("A flag was improperly set");
    //     }
    // } // end loop over i1
    // pout() << std::flush;

    // Perform the removal.
    if (needsRemoval > 0) {
        CH_assert(numBoxes - needsRemoval > 0);
        Vector<Box> newBoxes(numBoxes - needsRemoval);
        int newIdx = 0;
        for (int oldIdx = 0; oldIdx < numBoxes; ++oldIdx) {
            if (flags[oldIdx] == REMOVE) continue;
            newBoxes[newIdx] = a_boxes[oldIdx];
            ++newIdx;
        }
        a_boxes = newBoxes;
        numBoxes = a_boxes.size();
        flags.resize(numBoxes, KEEP);
    }

    // The remaining boxes may still overlap. Let's check.
    int needsFixing = 0;

    for (int i1 = 0; i1 < numBoxes; ++i1) {
        const Box& b1 = a_boxes[i1];
        for (int i2 = i1 + 1; i2 < numBoxes; ++i2) {
            const Box& b2 = a_boxes[i2];

            if (!b1.intersects(b2)) continue;

            const IntVect b1SmallEnd = b1.smallEnd() * a_spanDirs;
            const IntVect b2SmallEnd = b2.smallEnd() * a_spanDirs;
            const IntVect domSmallEnd = a_domBox.smallEnd() * a_spanDirs;

            const IntVect b1BigEnd = b1.bigEnd() * a_spanDirs;
            const IntVect b2BigEnd = b2.bigEnd() * a_spanDirs;
            const IntVect domBigEnd = a_domBox.bigEnd() * a_spanDirs;

            if (b1SmallEnd != domSmallEnd || b1BigEnd != domBigEnd) {
                if (flags[i1] != FIX) {
                    flags[i1] = FIX;
                    ++needsFixing;
                }
            }

            if (b2SmallEnd != domSmallEnd || b2BigEnd != domBigEnd) {
                if (flags[i2] != FIX) {
                    flags[i2] = FIX;
                    ++needsFixing;
                }
            }

            if (flags[i1] == KEEP && flags[i2] == KEEP) {
                if (flags[i2] != FIX) {
                    flags[i2] = FIX;
                    ++needsFixing;
                }
            }
        } // end loop over i2
    } // end loop over i1

    // Perform the fixing.
    if (needsFixing > 0) {
        // Do we ever actually need to deal with this?
        MayDay::Error("Looks like you will need to fix some boxes");
    }

    // If boxes were removed or fixed, let the caller know.
    if (needsRemoval || needsFixing) {
        return true;
    }

    // No boxes were removed or fixed.
    return false;
}


// -----------------------------------------------------------------------------
// Static utility
// Splits domain into vector of disjoint boxes with max size maxsize.
// This version does not split the domain in planes perpendicular to the
// vertical. This means the resulting grids will be suitable for leptic solves.
// If a_maxBoxSize[i] == 0, the domain will not be split in the ith direction.
// -----------------------------------------------------------------------------
void LepticMeshRefine::domainSplit (const ProblemDomain& a_domain,
                                    Vector<Box>&         a_vbox,
                                    const IntVect&       a_maxBoxSize,
                                    int                  a_blockFactor)
{
    const Box& domBox = a_domain.domainBox();
    LepticMeshRefine::domainSplit(domBox, a_vbox, a_maxBoxSize, a_blockFactor);
}


// -----------------------------------------------------------------------------
// Static utility
// Splits domain into vector of disjoint boxes with max size maxsize.
// This version does not split the domain in planes perpendicular to the
// vertical. This means the resulting grids will be suitable for leptic solves.
// If a_maxBoxSize[i] == 0, the domain will not be split in the ith direction.
// -----------------------------------------------------------------------------
void LepticMeshRefine::domainSplit (const Box&     a_domain,
                                    Vector<Box>&   a_vbox,
                                    const IntVect& a_maxBoxSize,
                                    int            a_blockFactor)
{
    // Sanity checks
    D_TERM(CH_assert(a_maxBoxSize[0] >= 0);,
           CH_assert(a_maxBoxSize[1] >= 0);,
           CH_assert(a_maxBoxSize[2] >= 0);)

    // Convert the block factor into a vector.
    IntVect maskedBF(D_DECL(a_blockFactor, a_blockFactor, a_blockFactor));
    D_TERM(if (a_maxBoxSize[0] == 0) maskedBF[0] = 1;,
           if (a_maxBoxSize[1] == 0) maskedBF[1] = 1;,
           if (a_maxBoxSize[2] == 0) maskedBF[2] = 1;)

    // Initialize workspace. Note that we work with a coarsened box to ensure
    // our results satisfy the blockFactor requirement.
    a_vbox.resize(0);
    Box d(a_domain);
    d.coarsen(maskedBF);

    // Test if the original box could not be coarsened by the block factor.
    // If so, something went wrong.
    if (refine(d, maskedBF) != a_domain) {
        MayDay::Error("LepticMeshRefine::domainSplit: a_domain not coarsenable by blockingFactor");
    }

    // Save this box and check if it needs to be broken up along any directions.
    a_vbox.push_back(d);
    for (int i = 0; i < CH_SPACEDIM; ++i) {
        if (a_maxBoxSize[i] == 0) continue;

        const int coarsenedMaxBoxSize = a_maxBoxSize[i] / a_blockFactor;
        LepticMeshRefine::breakBoxes(a_vbox, coarsenedMaxBoxSize, i);
    }

    // Finally, refine our results back to the original domain size.
    for (int i = 0; i < a_vbox.size(); ++i) {
        a_vbox[i].refine(maskedBF);
    }
}


// -----------------------------------------------------------------------------
//  Given a set of tagged cells defined on a single level of an AMR grid,
//  construct a BoxArray that covers all these cells that minimizes the
//  number of boxes and maximizes the ratio of tagged cells in each box.
//  This is part of the process of refining a mesh hierarchy based on error
//  estimates.
//
//  This gets tricky when the minbox around a set of tags has a non-zero offset.
//  Have to be careful to remember this when looking at box sizes.
// -----------------------------------------------------------------------------
void LepticMeshRefine::makeBoxes (Vector<Box>&         a_mesh,
                                  const IntVectSet&    a_tags,
                                  const IntVectSet&    a_pnd,
                                  const ProblemDomain& a_domain,
                                  const IntVect&       a_maxBoxSize,
                                  const IntVect&       a_totalBufferSize) const
{
    CH_TIME("LepticMeshRefine::makeBoxes");
    std::list<Box> boxes;

#ifdef CH_MPI
    int size;
    MPI_Comm_size (Chombo_MPI::comm, &size );
    Interval interval(0, size - 1);
    this->makeBoxesParallel(boxes, (IntVectSet&)a_tags, a_pnd, a_domain,
                            a_maxBoxSize, 0, a_totalBufferSize,
                            100, interval);
#else
    this->makeBoxes(boxes, (IntVectSet&)a_tags, a_pnd, a_domain,
                    a_maxBoxSize, 0, a_totalBufferSize);
#endif

    //boxes.sort();
    a_mesh.resize(boxes.size());
    std::list<Box>::iterator it = boxes.begin();
    for (int i = 0; i < a_mesh.size(); ++i, ++it) a_mesh[i] = *it;
}


// -----------------------------------------------------------------------------
// Constructs a set of boxes which covers a set of tagged cells by using the
// Berger-Rigoutsos algorithm.  Everything should be on the same level, and
// blocking factor is not applied. Boxes will be on the same refinement level as
// the tags. This would normally be a protected function, but it can be useful
// to call it on it's own, so it has been left public.
// -----------------------------------------------------------------------------
void LepticMeshRefine::makeBoxes (std::list<Box>&      a_mesh,
                                  IntVectSet&          a_tags,
                                  const IntVectSet&    a_pnd,
                                  const ProblemDomain& a_domain,
                                  const IntVect&       a_maxBoxSize,
                                  const int            a_depth,
                                  const IntVect&       a_totalBufferSize) const
{
    long long int Ntags  = a_tags.numPts();
    // Box minbx ;                                  //min box around Tags
    std::list<Box> mesh_hi ;                     //boxes from recursion
    IntVectSet& tags_lo = a_tags;
    IntVectSet  tags_hi ;                //tags for recursion
    a_mesh.clear();

    //
    // Handle special cases of no tags
    //
    if ( a_tags.isEmpty() ) {
        // return null box
        return;
    }

    //
    // Validate inputs
    //

    // The number of tagged cells in a box cannot exceed the number
    // of cells in the box so enforce an upper bound on \var{FillRatio}.
    //[NOTE: 0 or negative values are allowed -- they mean any box
    //       is acceptable.  This will probably be a mistake by the
    //       caller, but there is no obvious lower valid value for
    //       this variable (1e-6 is just as likely a mistake as 0).]
    CH_assert ( m_fillRatio <= 1.0 );

    //
    // Handle the case of all tags fitting in one box.
    // This always happens at the bottom of the recursion.
    //
    Box minbox = a_tags.minBox() ;
    bool nested = properlyNested(minbox, a_domain, a_pnd, a_totalBufferSize);

    if (!nested && Ntags == 1) {
        CH_assert(minbox.numPts() == 1);
        return; // IntVect wasn't in PND after all (bvs)
    }

    // If minbox has enough tagged cells and is properly nested, want
    // to add this minbox to the mesh.
    // If not, continue with splitting it and recursing on the pieces.
    if ( (Ntags >= (minbox.numPts() * m_fillRatio)) && nested ) {
        // no IntVects in the box are outside the PND so this box can
        // be accepted.  If it is larger than the maximum box size, split
        // it into two boxes
        //[NOTE: one of the boxes may violate the FillRatio requirement, but
        //       we will ignore this.]
        //  pout()<< "split Box "<<minbox<<"  maxsize "<<a_maxBoxSize<<std::endl;
        a_mesh.push_front(minbox) ;
        if (a_maxBoxSize.sum() > 0) {
            for (std::list<Box>::iterator it = a_mesh.begin(); it != a_mesh.end(); ++it) {
                this->splitBox( a_mesh, it, a_maxBoxSize ) ;
            }
        } // end if we are enforcing a max box size
    } else {
        // if efficiency criterion not met or box not properly nested...

        // Note tags_lo contains a_tags going in
        this->splitTagsInBestDimension(tags_lo, tags_hi, a_maxBoxSize);

        if (!(tags_lo.isEmpty())) {
            // low interval
            this->makeBoxes( a_mesh, tags_lo, a_pnd, a_domain, a_maxBoxSize, a_depth + 1, a_totalBufferSize);
        }
        if (!(tags_hi.isEmpty())) {
            // high interval
            this->makeBoxes( mesh_hi, tags_hi, a_pnd, a_domain, a_maxBoxSize, a_depth + 1, a_totalBufferSize);
        }

        // combine the results into a single mesh
        a_mesh.splice(a_mesh.begin(), mesh_hi);
    } // done if we need to split the box

    // Done
} //end of makeBoxes


// -----------------------------------------------------------------------------
// Does the same thing as makeBoxes, but across multiple processors.
// -----------------------------------------------------------------------------
void LepticMeshRefine::makeBoxesParallel (std::list<Box>&      a_mesh,
                                          IntVectSet&          a_tags,
                                          const IntVectSet&    a_pnd,
                                          const ProblemDomain& a_domain,
                                          const IntVect&       a_maxBoxSize,
                                          const int            a_depth,
                                          const IntVect&       a_totalBufferSize,
                                          const int            a_minSize,
                                          const Interval&      a_procInterval) const
{
    if (a_procInterval.size() == 1) {
        this->makeBoxes(a_mesh, a_tags, a_pnd, a_domain,
                        a_maxBoxSize, a_depth, a_totalBufferSize);
        return;
    }

    long long int Ntags  = a_tags.numPts();
    std::list<Box> mesh_hi ;             //boxes from recursion
    IntVectSet& tags_lo = a_tags;
    IntVectSet  tags_hi ;                //tags for recursion
    a_mesh.clear() ;

    //
    // Handle special cases of no tags
    //
    if ( a_tags.isEmpty() ) {
        //return null box
        return;
    }

    //
    // Validate inputs
    //

    // The number of tagged cells in a box cannot exceed the number
    // of cells in the box so enforce an upper bound on \var{FillRatio}.
    //[NOTE: 0 or negative values are allowed -- they mean any box
    //       is acceptable.  This will probably be a mistake by the
    //       caller, but there is no obvious lower valid value for
    //       this variable (1e-6 is just as likely a mistake as 0).]
    CH_assert ( m_fillRatio <= 1.0 );

    //
    // Handle the case of all tags fitting in one box.
    // This always happens at the bottom of the recursion.
    //
    Box minbox = a_tags.minBox() ;
    bool nested = properlyNested(minbox, a_domain, a_pnd, a_totalBufferSize);

    if (!nested && Ntags == 1) {
        CH_assert(minbox.numPts() == 1);
        return; // IntVect was not in PND after all
    }

    // If minbox has enough tagged cells and is properly nested, want
    // to add this minbox to the mesh.
    // If not, continue with splitting it and recursing on the pieces.
    if ( (Ntags >= minbox.numPts() * m_fillRatio) && nested) {
        // no IntVects in the box are outside the PND so this box can
        // be accepted.  If it is larger than the maximum box size, split
        // it into two boxes
        //[NOTE: one of the boxes may violate the FillRatio requirement, but
        //       we will ignore this.]
        //  pout()<< "split Box "<<minbox<<"  maxsize "<<a_maxBoxSize<<std::endl;
        a_mesh.push_front(minbox) ;
        if (a_maxBoxSize.sum() > 0) {
            for (std::list<Box>::iterator it = a_mesh.begin(); it != a_mesh.end(); ++it) {
                this->splitBox(a_mesh, it, a_maxBoxSize);
            }
        } // end if we are enforcing a max box size
    } else
    {
        // if efficiency criterion not met or box not properly nested...
        //

        // Note tags_lo = a_tags going in
        this->splitTagsInBestDimension(tags_lo, tags_hi, a_maxBoxSize);
        //a_tags.makeEmpty();

        if (a_procInterval.size() == 1 || Ntags <= a_minSize  ) { // do the regular algorithm for BRMeshRefine
            // Recurse on the two halves of the Tags
            if ( !tags_lo.isEmpty() ) {
                this->makeBoxes( a_mesh, tags_lo, a_pnd, a_domain, a_maxBoxSize, a_depth + 1, a_totalBufferSize);
            }
            if ( !tags_hi.isEmpty() ) {
                this->makeBoxes( mesh_hi, tags_hi, a_pnd, a_domain, a_maxBoxSize, a_depth + 1, a_totalBufferSize);
            }
        } else // do makeBoxes in Parallel
        {
            //pout()<<"depth, interval "<<a_depth<<" "<<a_procInterval.begin()
            //      <<a_procInterval.end()<<std::endl;

            // first, split interval in two
            Interval lo_interval(a_procInterval.begin(),
            (a_procInterval.end() + a_procInterval.begin() - 1) / 2);
            Interval hi_interval(lo_interval.end() + 1, a_procInterval.end());
            // pout()<<"lo "<<lo_interval.begin()<<lo_interval.end()
            //       <<"\nhi "<<hi_interval.begin()<<hi_interval.end()<<std::endl;
            if (lo_interval.contains(procID()) && !tags_lo.isEmpty()) {

                this->makeBoxesParallel(a_mesh, tags_lo, a_pnd, a_domain, a_maxBoxSize,
                                        a_depth + 1, a_totalBufferSize, a_minSize, lo_interval);
                sendBoxesParallel(a_mesh, a_depth);
            }
            if (hi_interval.contains(procID()) &&  !tags_hi.isEmpty()    ) {

                this->makeBoxesParallel(mesh_hi, tags_hi, a_pnd, a_domain, a_maxBoxSize,
                                        a_depth + 1, a_totalBufferSize, a_minSize, hi_interval);
                sendBoxesParallel(mesh_hi, a_depth);
            }
            if (hi_interval.contains(procID()) &&  !tags_lo.isEmpty()    ) {
                receiveBoxesParallel(hi_interval, lo_interval, a_mesh, a_depth);
            }
            if (lo_interval.contains(procID()) &&  !tags_hi.isEmpty()    ) {
                receiveBoxesParallel(lo_interval, hi_interval, mesh_hi, a_depth);
            }
        }
        // combine the results into a single mesh
        a_mesh.splice(a_mesh.begin(), mesh_hi);
    }
} //end of makeBoxesParallel


// -----------------------------------------------------------------------------
// Simply checks if a_pnd contains points that are properly nested in a_box.
// -----------------------------------------------------------------------------
bool LepticMeshRefine::properlyNested (const Box&           a_box,
                                       const ProblemDomain& a_domain,
                                       const IntVectSet&    a_pnd,
                                       const IntVect&       a_totalBufferSize) const
{
    if (m_PNDMode == 0) {
        MayDay::Error("LepticMeshRefine::properlyNested: m_PNDMode should not be zero");
        return a_pnd.contains(a_box);
    } else {
        Box growBox(a_box);
        growBox.grow(a_totalBufferSize);
        growBox &= a_domain.domainBox();
        if (!a_pnd.contains(growBox)) return false; //typical case
        if (a_domain.isPeriodic()) {
            Box growPeriodic(a_box);
            growPeriodic.grow(a_totalBufferSize);
            growPeriodic &= a_domain; // now intersect with (possibly) periodic domain
            if (growPeriodic != growBox) {
                //has periodic images
                //  for (BoxIterator bit(growPeriodic); bit.ok(); ++bit)
                //        {
                //          if (!growBox.contains(bit())){
                //            IntVect image=bit();
                //            bool isImage = a_domain.image(image);
                //            CH_assert(isImage);
                //            if (!a_pnd.contains(image)) return false;
                //          }
                //        }
                static ImageIterator images(a_domain);
                images.checkDefine(a_domain);
                for (images.begin(growPeriodic); images.ok(); ++images) {
                    if (!images.box().isEmpty()) {
                        if (!a_pnd.contains(images.box())) {
                            return false;
                        }
                    }
                }
            }
        }
    }
    return true;
}


// -----------------------------------------------------------------------------
// Compute the traces (signatures) of the minbox in each direction, and
// find a hole in the trace (zero value) and an inflection point (zero
// Laplacian) in each direction; keep the best of each.
//[NOTE: both \func{find*} functions return -1 if nothing was found.]
//[NOTE: \var{infl_val} is modified by \func{findMaxInflectionPoint}.]
//[NOTE: \var{trace} indexes from 0, so indices into \var{trace} must
//       be shifted to be used as indices into the \var{a_tags}.]
// a_tags_inout_lo holds the input tags coming in, and the lo tags going out.
// maxBoxSize is the maxboxsize and is used in determining where
// to begin searching for a split index.
// -----------------------------------------------------------------------------
void LepticMeshRefine::splitTagsInBestDimension(IntVectSet&    a_tags_inout_lo,
                                                IntVectSet&    a_tags_hi,
                                                const IntVect& a_maxBoxSize) const
{
    int hole_indx[SpaceDim], best_hole_dim;     //holes in traces
    int infl_indx[SpaceDim], best_infl_dim;     //inflection points in traces
    int infl_val [SpaceDim] ;                   //magnitudes of infl.points
    Vector<int> traces[SpaceDim];

    //  makeTraces( a_tags_inout_lo, traces) ;
    Box minbox = a_tags_inout_lo.minBox() ;
    IntVect offset = minbox.smallEnd() ;
    const IntVect& size = minbox.size();
    D_TERM6(traces[0].resize(size[0], 0); ,
    traces[1].resize(size[1], 0); ,
    traces[2].resize(size[2], 0); ,
    traces[3].resize(size[3], 0); ,
    traces[4].resize(size[4], 0); ,
    traces[5].resize(size[5], 0););
    IntVect iv;
    IVSIterator i(a_tags_inout_lo);
    for (i.begin() ; i.ok() ; ++i ) {
        iv = i() - offset;
        D_TERM6(traces[0][iv[0]]++,; traces[1][iv[1]]++,; traces[2][iv[2]]++,;
        traces[3][iv[3]]++,; traces[4][iv[4]]++,; traces[5][iv[5]]++);
    }

    for ( int idim = 0 ; idim < SpaceDim ; idim++ ) {
        if (a_maxBoxSize[idim] == 0) {
            hole_indx[idim] = -1;
            infl_indx[idim] = -1;
            continue;
        }

        //hole_indx[idim] = findSplit( traces[idim] ) ;
        //infl_indx[idim] = findMaxInflectionPoint(traces[idim], infl_val[idim] ) ;
        // The following two functions, with the a_maxBoxSize argument,
        //  help balance the tag splitting by changing where to begin
        //  searching for a split or inflection index.
        hole_indx[idim] = findSplit( traces[idim], a_maxBoxSize[idim]);
        infl_indx[idim] = findMaxInflectionPoint(traces[idim], infl_val[idim], a_maxBoxSize[idim]) ;
    }
    // Take the highest index as the best one because we want to take as large
    // a box as possible  (fewer large boxes are better than many small ones)
    best_hole_dim = maxloc( hole_indx, SpaceDim ) ;
    best_infl_dim = maxloc( infl_indx, SpaceDim ) ;

    //
    // Split the Tag set at a hole in one of the traces, if there is one, or an
    // inflection point in the Laplacian of the traces.  Failing that, split
    // at the middle of the longest dimension of the enclosing box.
    int split_dim, split_index;
    if ( hole_indx[best_hole_dim] >= 0 ) {
        // split at a hole in the trace, adjusting the trace index for the
        // offset into \var{a_tags_inout_lo} indices
        split_dim = best_hole_dim;
        split_index = hole_indx[best_hole_dim] + minbox.smallEnd(best_hole_dim);
    } else if ( infl_indx[best_infl_dim] >= 0 ) {
        // split at an inflection point in the trace, adjusting the trace
        // index for the offset into \var{a_tags_inout_lo} indices
        split_dim = best_infl_dim;
        split_index = infl_indx[best_infl_dim] + minbox.smallEnd(best_infl_dim);
    } else
    {
        // split on the midpoint of the longest side of \var{minbox}, rounding up,
        // allowing for \var{minbox} to have a non-zero offset
        minbox.longside(split_dim); //[NOTE: split_dim is set by \func(longside)]
        split_index = (minbox.smallEnd(split_dim) + minbox.bigEnd(split_dim) + 1) / 2;
    }

    splitTagsInPlace( split_dim, split_index, a_tags_inout_lo, a_tags_hi );
}


// -----------------------------------------------------------------------------
// Takes an element of the a_boxes list and splits it as needed.
// -----------------------------------------------------------------------------
void LepticMeshRefine::splitBox (std::list<Box>&                 a_boxes,
                                 const std::list<Box>::iterator& a_box,
                                 const IntVect&                  a_maxBoxSize) const
{
    // Find the direction that needs to be split before others. We can do this
    // by looking for the largest boxSize / maxBoxSize that is > 1.0.
    Box& b = *a_box;
    const IntVect boxSize = b.size();
    const RealVect ratios(D_DECL(
        (a_maxBoxSize[0] > 0)? Real(boxSize[0]) / Real(a_maxBoxSize[0]): 0.0,
        (a_maxBoxSize[1] > 0)? Real(boxSize[1]) / Real(a_maxBoxSize[1]): 0.0,
        (a_maxBoxSize[2] > 0)? Real(boxSize[2]) / Real(a_maxBoxSize[2]): 0.0));
    const int dir = ratios.maxDir(false); // Do not take abs vals.

    // If we don't have anything to do, just leave, ending the recursion.
    if (ratios[dir] <= 1.0) return;

    // Break up the box.
    int midpt = (b.smallEnd(dir) + b.bigEnd(dir) + 1) / 2;
    std::list<Box>::iterator chppt = a_boxes.insert(a_box, b.chop(dir, midpt));

    // Recurse on the two boxes
    this->splitBox(a_boxes, a_box, a_maxBoxSize);
    this->splitBox(a_boxes, chppt, a_maxBoxSize);
}


// -----------------------------------------------------------------------------
// Computes amount that tags on a level l must be coarsened so that grids on
// next finer level l+1 will be guaranteed to satisfy the blocking factor.
//
// The easiest way to do this is to coarsen everything down to a level which is
// m_blockFactor coarser than the new fine level which will be generated. So,
// m_level_blockfactor[lvl] is the amount by which the tags, etc at level lvl
// must be coarsened so that they are m_blockFactor coarser than the (lvl+1)
// grids which will be generated.  If m_BlockFactor is less than the refinement
// ratio between levels (lvl) and (lvl+1), then no coarsening needs to be done
// so we default to one, in that case.
// -----------------------------------------------------------------------------
void LepticMeshRefine::computeLocalBlockFactors ()
{
    for (int lev = 0; lev < (m_level_blockfactors.size()); ++lev) {
        // This is simply ceil(m_blockFactor / m_nRefVect[lev]).
        m_level_blockfactors[lev] =
            (  m_nRefVect[lev]
             + IntVect(D_DECL(m_blockFactor-1,m_blockFactor-1,m_blockFactor-1))  )
            / m_nRefVect[lev];

        // for (int dir = 0; dir < SpaceDim; ++dir) {
        //     if (m_spanDirs[dir] != 0) {
        //         m_level_blockfactors[lev][dir] = m_vectDomains[lev].size(dir);
        //     }
        // }

        // The coarsening ratio due to blocking needs to be a power of 2 otherwise
        // IntVectSet will complain.  Catch here to better describe the error.
        for (int dir = 0; dir < CH_SPACEDIM; ++dir) {
            if (!LepticMeshRefine::isPower2(m_level_blockfactors[lev][dir])) {
                pout() << "Unable to implement blocking for level " << lev + 1 << ".  "
                       "Blocking requires ceil(blockFactor/nRef) to be a power of 2 but "
                       "for nRef[" << lev << "] = " << m_nRefVect[lev]
                       << " and blockFactor = " << m_blockFactor << ", this is "
                       << m_level_blockfactors[lev] << '.' << endl;
                MayDay::Error("aborting LepticMeshRefine::regrid");
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
//
// Function:
// ---------
//  int makePNDs()
//
// Short Description:
// ------------------
//  Compute the proper nesting domains for a subset of the levels in a
//  multilevel adaptive mesh.
//
// Usage:
// ------
//  This is called by \func{regrid}.  It wasn't designed to be called
//  directly so if you do so you're on your own.  The output is a
//  \var{Vector<IntVectSet>}.  This function is called to compute the proper
//  nesting domains needed to do mesh refinement.  The function that calls this
//  one must have the mesh and problem domain data structures already defined
//  for all the levels and have allocated a vector<IntVectSet> for the proper
//  nesting domains.  This function computes the proper nesting domains of the
//  levels to be refined (i.e., the levels \var{a_baseLevel} and above) and
//  modifies the \var{pnds} array.  The output vector is undefined for
//  indices corresponding to mesh levels below \var{a_baseLevel}
//  (i.e. [0:\var{a_baseLevel}-1]).
//
//  NOTE: in MeshRefine::regrid, this function is called _after_ everthing
//        has been coarsened to enforce the blocking factor, so this behaves
//        differently than one might expect -- call outside of the regrid fn
//        at your own risk! (DFM 8/24/01)
//
// Arguments:
// ----------
//  \var{pnds}                     type: Vector<IntVectSet>
//      (output, index[\var{a_baseLevel}:\var{TopLevel}]) proper nesting domains
//
//  \var{a_baseLevel}                type: int
//      (input, range[0:\var{TopLevel}]) index of the mesh level to use as the
//      source for the proper nesting domains
//
//  \var{TopLevel}                 type: int
//      (input, range[\var{a_baseLevel}:\var{OldMeshes.size()-1]) index of the
//      finest mesh level to compute a PND for
//
//  \var{BaseMesh}                 type: Vector<Box>
//      (input) the boxes at level \var{a_baseLevel} in the mesh hierarchy.
//
//  \var{Domains}                  type: Vector<Box>
//      [same as for \func{MeshRefine}.]
//
//  \var{BufferSize}              type: int
//      [same as for \func{MeshRefine}.]
//
// Returns:
// --------
//  An int giving an exception number is returned.  The argument
//  \var{pnds} is modified.
//
// References:
// -----------
//  See ??? for a description of Proper Nesting Domains.
//  meshRefine -- calls this function
//
// Numerical Algorithm:
// --------------------
//  The proper nesting domain (PND) of level L+1 is defined as the interior of
//  the PND for level L (the next coarser) except at the boundaries of the
//  problem domain (ie, the level 0 box), in which case the PNDs of levels L
//  and L+1 both contain the boundary.  The PND of the base level is defined as
//  the interior of the union of all the grid boxes in the base level.  Given
//  the PND of level L (denoted by PND[L]) this function computes PND[L+1] by
//  iterating over the coordinate directions and computing the cells to remove
//  from the proper nesting domain due to the boundaries in that direction.
//  For each direction the domain is shifted (first up, then down) and only the
//  part that is shifted _out_ of the PND but remains _inside_ the problem
//  boundary is kept and then shifted back into the PND and subtracted
//  (removed) from it.  So if the shift moves part of the PND outside the
//  problem boundary then no cells will be removed.  If the cells that get
//  shifted out of the PND remain inside the boundary then they will be removed
//  from the PND when they get shifted back.  When all the directions have been
//  handled, the PND is refined to make it correct for level L+1.  The PND of
//  the base level is the union of the mesh boxes in the base level and serves
//  as the starting point.  In set notation the operation is:
//
//            D - ( ( ( (D << d) - D ) * B ) >> d )
//
//  where:
//    d is the current direction (+i,-i,+j,etc)
//    B is the boundary of the problem domain at the current level
//    - is the subtraction (removal) operator
//    * is the intersection operator
//   << is the down-shift operator (shift in negative direction)
//   >> is the up-shift operator (shift in positive direction)
//
//
// Implementation Method:
//  The PNDs are implemented as a vector of \type{IntVectSet}s.  A scratch
//  IntVectSet is used to store the PND of the current level while it is
//  being computed.  When the PND for a level is complete, it is refined
//  and used for the next level.  The base PND is computed as the union of
//  the mesh boxes in the base level mesh.  The PNDs for each level are
//  stored in the appropriate elements of the pnds vector as they are
//  computed.  Levels below a_baseLevel are not modified.  The code loops
//  through each level starting from Baselevel.  Loop through each coordinate
//  direction applying the operation defined in the "Numerical Algorithm"
//  section above.  Copy the final result into the PNDs vector.  If an
//  exception occurs, do not change PNDs on the level that took the
//   exception.
//
// Implementation Notes:
//  This assumes cell-centers are being manipulated, not vertices.  The whole
//  vector of pnds is passed even though only some are accessed
//  because it is more convenient since these variables are likely to exist
//  already.
//
///////////////////////////////////////////////////////////////////////////////
void LepticMeshRefine::makePNDs (Vector<IntVectSet>&          a_pnds,
                                 Vector<IntVect>&             a_totalBufferSize,
                                 const int                    a_baseLevel,
                                 const int                    a_topLevel,
                                 const Vector<ProblemDomain>& a_domains,
                                 const IntVectSet&            a_baseMesh,
                                 const Vector<IntVect>&       a_bufferSize) const
{
    // Validate inputs
    CH_assert( a_baseLevel <= a_topLevel && a_baseLevel >= 0 );
    CH_assert( a_domains.size() >= a_topLevel + 1 );
    CH_assert( m_nRefVect.size() >= a_topLevel ); // ES: Changed from toplevel+1
    // all existing boxes on the base level must be aligned to BuffSize

    a_pnds[a_baseLevel] = a_baseMesh;
    a_totalBufferSize[a_baseLevel] = IntVect::Zero;

    for ( int lvl = a_baseLevel ; lvl <= a_topLevel ; lvl++) {
        IntVectSet& pnd = a_pnds[lvl];

        if (m_PNDMode == 0) {
            MayDay::Error("LepticMeshRefine::m_PNDMode should not be zero");
            // pnd.nestingRegion(a_bufferSize[lvl], a_domains[lvl], m_granularity);
        }
        a_totalBufferSize[lvl] += a_bufferSize[lvl];

        if ( (a_topLevel - lvl) > 0) {
            a_pnds[lvl + 1] = pnd;

            // This does:
            //   a) refine  by BF_ref[lvl];
            //   b) refine  by  n_ref[lvl];
            //   c) coarsen by BF_ref[lvl+1];
            IntVect allInOne_nRef = m_nRefVect[lvl];
            allInOne_nRef *= m_level_blockfactors[lvl];
            allInOne_nRef /= m_level_blockfactors[lvl + 1];

            a_totalBufferSize[lvl + 1] = a_totalBufferSize[lvl] * allInOne_nRef;

            refine(a_pnds[lvl + 1], allInOne_nRef);
        }
    }
}


void LepticMeshRefine::makePNDs (Vector<IntVectSet>&          a_pnds,
                                 Vector<IntVect>&             a_totalBufferSize,
                                 const int                    a_baseLevel,
                                 const int                    a_topLevel,
                                 const Vector<ProblemDomain>& a_domains,
                                 const Vector<Box>&           a_oldMeshes,
                                 const Vector<IntVect>&       a_bufferSize) const
{
    IntVectSet mesh;
    for (int box = 0; box < a_oldMeshes.size(); ++box) {
        const Box& b = a_oldMeshes[box];
        mesh |= b;
    }
    this->makePNDs(a_pnds, a_totalBufferSize, a_baseLevel, a_topLevel,
                   a_domains, mesh, a_bufferSize);
}


// -----------------------------------------------------------------------------
// Static utility
// Recursive function to enforce max size of boxes in a given direction.
// -----------------------------------------------------------------------------
void LepticMeshRefine::breakBoxes (Vector<Box>& a_vboxin,
                                   const int&   a_maxBoxSize,
                                   const int&   a_idir)
{
    if(a_maxBoxSize == 0) return;

    int nboxes = a_vboxin.size();
    //need to use STL vector for bools.
    using std::vector;
    vector<bool> splitDec(nboxes);
    bool anyBoxesToSplit = false;
    //find out which boxes need to be chopped and in what direction
    for (int ibox = 0; ibox < nboxes; ibox++) {
        if ( a_vboxin[ibox].size(a_idir ) > a_maxBoxSize ) {
            splitDec[ibox] = true;
            anyBoxesToSplit = true;
        } else {
            splitDec[ibox] = false;
        }
    }
    //if there are no boxes to split, just return
    //otherwise, split all the boxes that need to be
    //split ONCE and then call function recursively
    // and set the return vector to the temporary after
    // the recursion
    if (anyBoxesToSplit) {
        Vector<Box> vboxtemp;
        for (int ibox = 0; ibox < nboxes; ibox++) {
            Box boxone = a_vboxin[ibox];
            if (splitDec[ibox]) {
                // Chombo does this in different ways. This seems to be the
                // way that works.
                int mid = boxone.smallEnd(a_idir) + a_maxBoxSize;

                Box boxtwo = boxone.chop(a_idir, mid);
                vboxtemp.push_back(boxone);
                vboxtemp.push_back(boxtwo);
            } else {
                vboxtemp.push_back(boxone);
            }
        } // end loop over boxes
        breakBoxes(vboxtemp, a_maxBoxSize, a_idir);
        a_vboxin = vboxtemp;
    }
    return;
}


// -----------------------------------------------------------------------------
// All of the code below this line was added to avoid the BRMeshRefine memory
// leak. It adds no new functionality.
// -----------------------------------------------------------------------------
#define MAXBOXES     130000
#define P_BUFFERSIZE (130000 * 8 * CH_SPACEDIM)
extern std::list<int*> sendBuffers;
extern std::list<int>  ch_count;

int* LepticMeshRefine::s_recBuffer = NULL;

void LepticMeshRefine::deleteBuffer ()
{
#ifdef CH_MPI
    if (s_recBuffer != NULL) {
        freeMT(s_recBuffer);
    }
    s_recBuffer = NULL;
#endif
}


void LepticMeshRefine::receiveBoxesParallel (const Interval& a_from,
                                             const Interval& a_to,
                                             std::list<Box>& a_mesh,
                                             int             tag) const
{
#ifdef CH_MPI

    // pout()<<"from "<<a_from.begin()<<a_from.end()<<"\n"
    //      <<"to   "<<a_to.begin()<<a_to.end()<<"\n";
    // NOTE: This static malloc is currently not being freed anywhere (ndk)
    // static int* recBuffer = (int*)mallocMT(P_BUFFERSIZE);
    if (s_recBuffer == NULL) {
        s_recBuffer = (int*)mallocMT(P_BUFFERSIZE);
    }

    int* next = s_recBuffer;
    MPI_Status status;
    const int boxSize = 2 * CH_SPACEDIM;

    int source, dest;

    source = procID() - a_from.begin() + a_to.begin();
    dest = source;

    bool hang = false;
    if (a_from.size() > a_to.size() && procID() == a_from.end()) hang = true;
    if (a_to.size() > a_from.size() && procID() == a_to.end()) hang = true;

    if (a_from.size() == a_to.size() || !hang) {
        //pout()<<"expecting boxes from "<<source<<"\n";
        //pout()<<"sending "<<ch_count.back()/boxSize<<" boxes to "<<dest<<std::endl;
        MPI_Sendrecv( sendBuffers.back(), ch_count.back(), MPI_INT,
                      dest, tag,
                      s_recBuffer, MAXBOXES * boxSize, MPI_INT,
                      source, tag, Chombo_MPI::comm, &status );
    }
    // OK, need to pick up the oddball stragler here.
    if (a_from.size() < a_to.size() && procID() == a_from.end()) {
        dest = a_to.end();
        //pout()<<"SEnding "<<ch_count.back()/boxSize<<" boxes to "<<dest<<std::endl;
        MPI_Send(sendBuffers.back(), ch_count.back(), MPI_INT,
                 a_to.end(), tag, Chombo_MPI::comm);
    }
    if (a_from.size() > a_to.size() && procID() == a_from.end()) {
        source = a_to.end();
        //pout()<<"EXpecting boxes from "<<source<<"\n";
        MPI_Recv(s_recBuffer, MAXBOXES * boxSize, MPI_INT,
                 source, tag, Chombo_MPI::comm, &status );
    }

    delete[] sendBuffers.back();
    sendBuffers.pop_back();
    ch_count.pop_back();

    int numInt;
    MPI_Get_count(&status, MPI_INT, &numInt);

    int* end = next + numInt;
    Box b;
    IntVect& lo = (IntVect&)(b.smallEnd());
    IntVect& hi = (IntVect&)(b.bigEnd());
    for (; next < end; next += (2 * CH_SPACEDIM)) {
        D_TERM6(lo[0] = next[0]; hi[0] = next[1]; ,
                lo[1] = next[2]; hi[1] = next[3]; ,
                lo[2] = next[4]; hi[2] = next[5]; ,
                lo[3] = next[6]; hi[3] = next[7]; ,
                lo[4] = next[8]; hi[4] = next[9]; ,
                lo[5] = next[10]; hi[5] = next[11];);

        b.computeBoxLenNotEmpty();
        a_mesh.push_back(b);
    }
    //pout()<<"received "<<a_mesh.size()<<" boxes from "<<source<<std::endl;
#endif
}

