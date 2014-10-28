#include "HeaderOverrides/DisjointBoxLayout.H"
#include "HeaderOverrides/Copier.H"
#include "HeaderOverrides/CFRegion.H"
#include "HeaderOverrides/CFIVS.H"
#include "HeaderOverrides/TreeIntVectSet.H"
#include "HeaderOverrides/IntVectSet.H"
#include "AnisotropicRefinementTools.H"


// -----------------------------------------------------------------------------
// Returns true if a_box can be coarsened by a_refRatio and
// return back to the original Box when refined by a_refRatio.
// -----------------------------------------------------------------------------
bool coarsenable (const Box&     a_box,
                  const IntVect& a_refRatio)
{
    Box b = a_box;
    b.coarsen(a_refRatio);
    b.refine(a_refRatio);
    return (b == a_box);
}


// -----------------------------------------------------------------------------
// A very simple utility to calculate how much refinement is needed to send
// data from a_src to a_dest. Returns a negative value if coarsening is
// needed and zero if the two are incompatible.
// -----------------------------------------------------------------------------
IntVect calculateRefinementRatio (const Box& a_src,
                                  const Box& a_dest)
{
    Box crseBox = a_src;
    Box fineBox = a_dest;

    CH_assert(!crseBox.isEmpty());
    CH_assert(!fineBox.isEmpty());

    // Check if refinement or coarsening is needed.
    IntVect refRatio = IntVect::Zero;
    for (int dir = 0; dir < SpaceDim; ++dir) {
        if (fineBox.size(dir) == crseBox.size(dir)) {
            refRatio[dir] = 1;
        } else if (fineBox.size(dir) > crseBox.size(dir)) {
            if (fineBox.size(dir) % crseBox.size(dir) == 0) {
                refRatio[dir] = fineBox.size(dir) / crseBox.size(dir);
            }
        } else {
            if (crseBox.size(dir) % fineBox.size(dir) == 0) {
                refRatio[dir] = -crseBox.size(dir) / fineBox.size(dir);
            }
        }
    }

    return refRatio;
}


// -----------------------------------------------------------------------------
// ProblemDomain coarsen.
// Use this instead of Chombo's coarsen.
// -----------------------------------------------------------------------------
void coarsen (ProblemDomain&       a_crseDomain,
              const ProblemDomain& a_fineDomain,
              const IntVect&       a_refRatio)
{
    Box newDomBox(a_fineDomain.domainBox());
    newDomBox.coarsen(a_refRatio);
    bool isPeriodic[CH_SPACEDIM];
    D_TERM(isPeriodic[0] = a_fineDomain.isPeriodic(0);,
           isPeriodic[1] = a_fineDomain.isPeriodic(1);,
           isPeriodic[2] = a_fineDomain.isPeriodic(2);)
    a_crseDomain.define(newDomBox, isPeriodic);
}


// -----------------------------------------------------------------------------
// ProblemDomain refine.
// Use this instead of Chombo's refine.
// -----------------------------------------------------------------------------
void refine (ProblemDomain&        a_fineDomain,
             const ProblemDomain& a_crseDomain,
             const IntVect&       a_refRatio)
{
    Box newDomBox(a_crseDomain.domainBox());
    newDomBox.refine(a_refRatio);
    bool isPeriodic[CH_SPACEDIM];
    D_TERM(isPeriodic[0] = a_crseDomain.isPeriodic(0);,
           isPeriodic[1] = a_crseDomain.isPeriodic(1);,
           isPeriodic[2] = a_crseDomain.isPeriodic(2);)
    a_fineDomain.define(newDomBox, isPeriodic);
}


// -----------------------------------------------------------------------------
// Returns true if every Box in the BoxLayout can be coarsened by a_refRatio and
// return back to the original Box when refined by a_refRatio.
// -----------------------------------------------------------------------------
bool coarsenable (const BoxLayout& a_layout,
                  const IntVect&   a_refRatio)
{
    Vector<Box> boxes = a_layout.boxArray();

    for (int idx = 0; idx < boxes.size(); ++idx) {
        Box b = boxes[idx];
        b.coarsen(a_refRatio);
        b.refine(a_refRatio);
        if (b != boxes[idx]) return false;
    }
    return true;
}


// -----------------------------------------------------------------------------
// BoxLayout coarsen
// -----------------------------------------------------------------------------
void coarsen (BoxLayout&       a_output,
              const BoxLayout& a_input,
              const IntVect&   a_refinement)
{
    if (!a_input.isClosed()) {
        MayDay::Error("input to coarsen must be called with closed BoxLayout");
    }
    if (a_output.isClosed()) {
        MayDay::Error("output of coarsen must be called on open BoxLayout");
    }
    //a_output.deepCopy(a_input);
    a_output.m_boxes      = RefCountedPtr<Vector<Entry> >(new Vector<Entry>(*(a_input.m_boxes)));
    a_output.m_layout     = a_input.m_layout;
#ifdef CH_MPI
    a_output.m_dataIndex  = a_input.m_dataIndex;
#endif

    for (int ivec = 0; ivec < a_output.m_boxes->size(); ivec++) {
        (*a_output.m_boxes)[ivec].box.coarsen(a_refinement);
    }
    a_output.close();
}


// -----------------------------------------------------------------------------
// BoxLayout refine
// -----------------------------------------------------------------------------
void refine (BoxLayout&       a_output,
             const BoxLayout& a_input,
             const IntVect&   a_refinement)
{
    if (!a_input.isClosed()) {
        MayDay::Error("input to refine must be called with closed BoxLayout");
    }
    if (a_output.isClosed()) {
        MayDay::Error("output of refine must be called on open BoxLayout");
    }
    a_output.deepCopy(a_input);

    for (int ivec = 0; ivec < a_output.m_boxes->size(); ivec++) {
        (*a_output.m_boxes)[ivec].box.refine(a_refinement);
    }
    a_output.close();
}


// -----------------------------------------------------------------------------
// DisjointBoxLayout coarsen
// -----------------------------------------------------------------------------
void coarsen (DisjointBoxLayout&       a_output,
              const DisjointBoxLayout& a_input,
              const IntVect&           a_refinement)
{
    if (a_input.size() == 0) {
        a_output = a_input;
        return;
    }
    CH_assert(coarsenable(a_input, a_refinement));
    if (!a_input.isClosed()) {
        MayDay::Error("input to coarsen must be called with closed BoxLayout");
    }
    if (a_output.isClosed()) {
        MayDay::Error("output of coarsen must be called on open BoxLayout");
    }

    // copy first, then coarsen everything
    // a_output.deepCopy(a_input);
    a_output.m_boxes      = RefCountedPtr<Vector<Entry> >(new Vector<Entry>(*(a_input.m_boxes)));
    a_output.m_layout     = a_input.m_layout;
#ifdef CH_MPI
    a_output.m_dataIndex  = a_input.m_dataIndex;
#endif
    // now coarsen the physDomain
    a_output.m_physDomain = coarsen(a_input.m_physDomain, a_refinement);

    Vector<Entry>& boxes = *(a_output.m_boxes);
    int j = 0;
    for (int i = 0 ; i <= (int)boxes.size() - 4; i += 4) {
        boxes[i].box.coarsen(a_refinement);
        boxes[i + 1].box.coarsen(a_refinement);
        boxes[i + 2].box.coarsen(a_refinement);
        boxes[i + 3].box.coarsen(a_refinement);
        j += 4;
    }
    for (; j < boxes.size(); j++) boxes[j].box.coarsen(a_refinement);

    //    for (LayoutIterator it(a_input.layoutIterator()); it.ok(); ++it)
    //      {
    //        a_output.ref(it()).coarsen(a_refinement);
    //     }
    a_output.closeN(a_input.m_neighbors);
}


// -----------------------------------------------------------------------------
// DisjointBoxLayout refine
// -----------------------------------------------------------------------------
void refine (DisjointBoxLayout&       a_output,
             const DisjointBoxLayout& a_input,
             const IntVect&           a_refinement)
{
    if (!a_input.isClosed()) {
        MayDay::Error("input to refine must be called with closed BoxLayout");
    }
    if (a_output.isClosed()) {
        MayDay::Error("output of refine must be called on open BoxLayout");
    }

    // first copy, then refine everything
    a_output.deepCopy(a_input);

    // start by refining the physDomain
    a_output.m_physDomain = refine(a_input.m_physDomain, a_refinement);

    for (int ivec = 0; ivec < a_output.m_boxes->size(); ivec++) {
        (*a_output.m_boxes)[ivec].box.refine(a_refinement);
    }
    a_output.closeN(a_input.m_neighbors);
}


// -----------------------------------------------------------------------------
// Copier coarsen
// -----------------------------------------------------------------------------
void coarsen (Copier&        a_copier,
              const IntVect& a_refRatio)
{
    for (int i = 0; i < a_copier.m_localMotionPlan.size(); ++i) {
        a_copier.m_localMotionPlan[i]->fromRegion.coarsen(a_refRatio);
        a_copier.m_localMotionPlan[i]->toRegion.coarsen(a_refRatio);
    }
    for (int i = 0; i < a_copier.m_fromMotionPlan.size(); ++i) {
        a_copier.m_fromMotionPlan[i]->fromRegion.coarsen(a_refRatio);
        a_copier.m_fromMotionPlan[i]->toRegion.coarsen(a_refRatio);
    }
    for (int i = 0; i < a_copier.m_toMotionPlan.size(); ++i) {
        a_copier.m_toMotionPlan[i]->fromRegion.coarsen(a_refRatio);
        a_copier.m_toMotionPlan[i]->toRegion.coarsen(a_refRatio);
    }
}






#include "BoxIterator.H"
#include "BRMeshRefine.H"
// -----------------------------------------------------------------------------
// IntVectSet coarsen
// -----------------------------------------------------------------------------
void coarsen (IntVectSet&    a_ivs,
              const IntVect& a_ref)
{
    // If a_ref is isotropic, just run Chombo's fast coarsen function.
    if (D_TERM(true, && a_ref[0] == a_ref[1], && a_ref[0] == a_ref[2])) {
        a_ivs.coarsen(a_ref[0]);
        return;
    }

    if (a_ivs.isDense()) {
        CH_TIME("coarsenIVS_dense_branch");

        a_ivs.compact(); // TODO: Do we want this?
        const Vector<Box> boxes = a_ivs.boxes();
        const int numBoxes = boxes.size();

        DenseIntVectSet denseIVS(coarsen(a_ivs.minBox(), a_ref),
                                 false); // Start with no tags set.
        for (int idx = 0; idx < numBoxes; ++idx) {
            denseIVS |= coarsen(boxes[idx], a_ref);
        }
        a_ivs.define(denseIVS);

    } else {
        CH_TIME("coarsenIVS_tree_branch");

// #define DO_PARALLEL_VERSION                  // Parallel version is extremely slow!
#ifdef DO_PARALLEL_VERSION
        // Get the coarsened min box.
        a_ivs.compact();
        Box crseMinBox = a_ivs.minBox();
        crseMinBox.coarsen(a_ref);

        // Distribute this box over the procs
        Vector<Box> vBox;
        domainSplit(crseMinBox, vBox, crseMinBox.size(0) / numProc());
        DisjointBoxLayout dbl;
        dbl.defineAndLoadBalance(vBox, NULL);
        DataIterator dit = dbl.dataIterator();

        // Loop over the IVs in the box and find cells to remove.
        IntVectSet crseIVS;
        for (dit.reset(); dit.ok(); ++dit) {
            // Begin with all cells in the coarsened minBox tagged.
            crseIVS |= dbl[dit];

            // Iterator through the coarsened minBox, searching for
            // IntVects to remove.
            BoxIterator bit(dbl[dit]);
            for (bit.reset(); bit.ok(); ++bit) {
                const IntVect& crseIV = bit();

                // This coarse IV corresponds to a set of points, fineTestBox,
                // that need to be checked. If any of those points are in the
                // IVS, then we cannot remove this coarse IV.
                Box fineTestBox(crseIV, crseIV);
                fineTestBox.refine(a_ref);

                bool doRemove = true;
                BoxIterator testBit(fineTestBox);
                for (testBit.reset(); testBit.ok(); ++testBit) {
                    const IntVect& testIV = testBit();
                    if (a_ivs.contains(testIV)) {
                        doRemove = false;
                        break;
                    }
                }

                // Remove the coarse IV in nothing was found in the test box.
                if (doRemove) crseIVS -= crseIV;
            }
        }

        // Gather data into vector on src proc
        const int srcProc = uniqueProc(SerialTask::compute);
        Vector<IntVectSet> vdata;
        gather(vdata, crseIVS, srcProc);

        // Put gathered data into dest data holder on src proc
        crseIVS.makeEmpty();
        if (procID() == srcProc) {
            for (int idx = 0; idx < vdata.size(); ++idx) {
                crseIVS |= vdata[idx];
            }
        }

        // Distribute data among all procs
        broadcast(crseIVS, srcProc);

        // Copy the result to the output holder.
        crseIVS.compact();
        a_ivs.define(crseIVS);
#else
        // Begin with all cells in the coarsened minBox tagged.
        a_ivs.compact();
        Box crseMinBox = a_ivs.minBox();
        crseMinBox.coarsen(a_ref);
        DenseIntVectSet crseIVS(crseMinBox, true); // Should this be dense?
        // pout() << "crseMinBox = " << crseMinBox << endl;

        // Iterator through the coarsened minBox, searching for
        // IntVects to remove.
        BoxIterator bit(crseMinBox);
        for (bit.reset(); bit.ok(); ++bit) {
            const IntVect& crseIV = bit();

            // This coarse IV corresponds to a set of points, fineTestBox,
            // that need to be checked. If any of those points are in the
            // IVS, then we cannot remove this coarse IV.
            Box fineTestBox(crseIV, crseIV);
            fineTestBox.refine(a_ref);

            bool doRemove = true;
            BoxIterator testBit(fineTestBox);
            for (testBit.reset(); testBit.ok(); ++testBit) {
                const IntVect& testIV = testBit();
                if (a_ivs.contains(testIV)) {
                    doRemove = false;
                    break;
                }
            }

            // Remove the coarse IV in nothing was found in the test box.
            if (doRemove) crseIVS -= crseIV;
        }

        // Copy the result to the output holder.
        crseIVS.compact();
        a_ivs.define(crseIVS);
#endif

// #ifndef CH_MPI
//         a_ivs.makeEmpty();
//         for (int idx = 0; idx < numBoxes; ++idx) {
//             a_ivs |= coarsen(boxes[idx], a_ref);
//         }
// #else
//         const int thisProc = procID();
//         const int numProcs = numProc();
//         const int serialProc = uniqueProc(SerialTask::compute);

//         // Each processor takes a box, coarsens it and throws it on a local IVS.
//         // This shares the burden of the two most expensive operations, calling
//         // the IntVectSet::operator|= and Box::coarsen functions.
//         IntVectSet localCrseIVS;
//         {
//             // TODO: This is a bottleneck. See if you can speed this up.
//             CH_TIME("create_localCrseIVS");
//             for (int sidx = thisProc; sidx < numBoxes; sidx += numProcs) {
//                 localCrseIVS |= coarsen(boxes[sidx], a_ref);
//             }

//             // TODO: Does this slow things down?
//             localCrseIVS.compact();
//         }

//         // Gather the coarsened IVSs onto the serial task proc.
//         Vector<IntVectSet> allIVS(numProcs);
//         {
//             CH_TIME("gather_localCrseIVS");
//             gather(allIVS, localCrseIVS, serialProc);
//         }

//         // Free up resources.
//         localCrseIVS.makeEmpty();
//         a_ivs.makeEmpty();

//         // Join the IVSs on the serial task proc.
//         {
//             CH_TIME("serial_work");
//             if (thisProc == serialProc) {
//                 for (int i = 0; i < allIVS.size(); ++i) {
//                     a_ivs |= allIVS[i];
//                 }
//             }
//         }

//         // Free up resources.
//         allIVS.resize(0);

//         // Broadcast the results.
//         {
//             CH_TIME("broadcast_ivs");
//             broadcast (a_ivs, serialProc);
//         }
// #endif
    } // end if a_ivs is not dense
}


// -----------------------------------------------------------------------------
// IntVectSet refine
// -----------------------------------------------------------------------------
void refine (IntVectSet&    a_ivs,
             const IntVect& a_ref)
{
    Vector<Box> boxes = a_ivs.boxes();

    if (a_ivs.isDense()) {
        CH_TIME("refineIVS_dense_branch");

        DenseIntVectSet denseIVS(refine(a_ivs.minBox(), a_ref),
                                 false); // Start with no tags set.
        for (int idx = 0; idx < boxes.size(); ++idx) {
            denseIVS |= refine(boxes[idx], a_ref);
        }
        a_ivs.define(denseIVS);

    } else {
        CH_TIME("refineIVS_tree_branch");

        IntVectSet treeIVS;
        for (int idx = 0; idx < boxes.size(); ++idx) {
            treeIVS |= refine(boxes[idx], a_ref);
        }
        a_ivs.define(treeIVS);
    }

    // // TODO: Check if this is a bottleneck.
    // a_ivs.compact();
}


// -----------------------------------------------------------------------------
// CFRegion coarsen
// -----------------------------------------------------------------------------
void coarsen (CFRegion&      a_cfRegion,
              const IntVect& a_refRatio)
{
    for (DataIterator dit = a_cfRegion.m_loCFIVS[0].dataIterator(); dit.ok(); ++dit) {
        for (int i = 0; i < CH_SPACEDIM; ++i) {
            coarsen((a_cfRegion.m_loCFIVS[i])[dit], a_refRatio);
            coarsen((a_cfRegion.m_hiCFIVS[i])[dit], a_refRatio);
        }
    }
}


// -----------------------------------------------------------------------------
// CFIVS coarsen
// -----------------------------------------------------------------------------
void coarsen (CFIVS&         a_cfivs,
              const IntVect& a_refRatio)
{
    CH_assert(a_cfivs.m_defined);
    if (!a_cfivs.m_empty) {
        coarsen(a_cfivs.m_IVS, a_refRatio);
        a_cfivs.m_packedBox.coarsen(a_refRatio);
    }
}
