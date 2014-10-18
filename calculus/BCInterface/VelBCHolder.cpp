#include "VelBCHolder.H"
#include "FluxBox.H"


// -------------------------------------------------------------
VelBCHolder::VelBCHolder()
{;}


// --------------------------------------------------------------
VelBCHolder::~VelBCHolder()
{;}


// --------------------------------------------------------------
VelBCHolder::VelBCHolder(const Tuple<BCMethodHolder, SpaceDim>& a_componentBC)
{
    for (int idir = 0; idir < SpaceDim; ++idir) {
        m_componentBC[idir] = a_componentBC[idir];
    }
}


// --------------------------------------------------------------
void VelBCHolder::setGhosts(FArrayBox&           a_state,
                            const FArrayBox*     a_extrapPtr,
                            const Box&           a_valid,
                            const ProblemDomain& a_domain,
                            const RealVect&      a_dx,
                            const DataIndex&     a_index,
                            const FluxBox*       a_JgupPtr,
                            bool                 a_homogeneous,
                            Real                 a_time)
{
    CH_assert(a_state.nComp() == SpaceDim);

    for (int idir = 0; idir < SpaceDim; ++idir) {
        m_componentBC[idir].setGhosts(a_state,
                                      a_extrapPtr,
                                      a_valid,
                                      a_domain,
                                      a_dx,
                                      a_index,
                                      a_JgupPtr,
                                      a_homogeneous,
                                      a_time,
                                      Interval(idir,idir));
    }
}

// --------------------------------------------------------------
void VelBCHolder::setGhosts(LevelData<FArrayBox>&       a_state,
                            const LevelData<FArrayBox>* a_extrapPtr,
                            const RealVect&             a_dx,
                            const LevelData<FluxBox>*   a_JgupPtr,
                            bool                        a_homogeneous,
                            Real                        a_time)
{
    CH_assert(a_state.nComp() == SpaceDim);

    const DisjointBoxLayout& grids = a_state.getBoxes();
    const ProblemDomain& domain = grids.physDomain();
    DataIterator dit = a_state.dataIterator();

#ifndef NDEBUG
    if (a_JgupPtr != NULL) {
        if (!a_JgupPtr->getBoxes().compatible(grids)) {
            pout() << "ERROR: These DBLs are incompatible.\n"
                   << "-- a_JgupPtr:\n"
                   << "domain = " << a_JgupPtr->getBoxes().physDomain() << "\n"
                   << "dbl = " << a_JgupPtr->getBoxes()
                   << "-- a_state:\n"
                   << "domain = " << a_state.getBoxes().physDomain() << "\n"
                   << "dbl = " << a_state.getBoxes() << endl;
            CH_assert(a_JgupPtr->getBoxes().compatible(grids));
        }
    }
#endif

    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& stateFAB = a_state[dit];
        const FArrayBox* extrapFABPtr = (a_extrapPtr == NULL)? NULL: &((*a_extrapPtr)[dit]);
        const FluxBox* JgupFBPtr = (a_JgupPtr == NULL)? NULL: &((*a_JgupPtr)[dit]);

        // We need to fill boundary ghosts adjacent to exchange ghosts too!
        const Box valid = a_state[dit].box() & domain;

        this->setGhosts(stateFAB,
                        extrapFABPtr,
                        valid,
                        domain,
                        a_dx,
                        dit(),
                        JgupFBPtr,
                        a_homogeneous,
                        a_time);
    }
}


// --------------------------------------------------------------
void VelBCHolder::setFluxes(FluxBox&             a_state,
                            const FluxBox*       a_extrapPtr,
                            const Box&           a_validCC,
                            const ProblemDomain& a_domain,
                            const RealVect&      a_dx,
                            const DataIndex&     a_index,
                            const FluxBox*       a_JgupPtr,
                            bool                 a_homogeneous,
                            Real                 a_time)
{

    for (int edgedir = 0; edgedir < SpaceDim; ++edgedir) {
        const FArrayBox* extrapFABPtr = (a_extrapPtr == NULL)? NULL: &((*a_extrapPtr)[edgedir]);
        Interval interv(edgedir, edgedir);  // Old way
        // const Interval& interv = a_state[edgedir].interval();  // Needed for ZeroFluxBC to work

        for (int velcomp = 0; velcomp < SpaceDim; ++velcomp) {
            m_componentBC[edgedir].setFluxes(a_state[edgedir],
                                             extrapFABPtr,
                                             a_validCC,
                                             a_domain,
                                             a_dx,
                                             a_index,
                                             a_JgupPtr,
                                             velcomp,
                                             a_homogeneous,
                                             a_time,
                                             interv);   // TODO: I think this should be Interval(velcomp,velcomp)
        }
    }
}


// -------------------------------------------------------------
EdgeVelBCHolder::EdgeVelBCHolder()
{
}

// --------------------------------------------------------------
EdgeVelBCHolder::~EdgeVelBCHolder()
{
}

// --------------------------------------------------------------
EdgeVelBCHolder::EdgeVelBCHolder(const Tuple<BCMethodHolder, SpaceDim>& a_componentBC)
{
    for (int idir = 0; idir < SpaceDim; ++idir) {
        m_componentBC[idir] = a_componentBC[idir];
    }
}

// --------------------------------------------------------------
void EdgeVelBCHolder::setGhosts(LevelData<FluxBox>&       a_state,
                                const LevelData<FluxBox>* a_extrapPtr,
                                const RealVect&           a_dx,
                                const LevelData<FluxBox>* a_JgupPtr,
                                bool                      a_homogeneous,
                                Real                      a_time)
{

    const DisjointBoxLayout& grids = a_state.getBoxes();
    const ProblemDomain& domain = grids.physDomain();
    DataIterator dit = a_state.dataIterator();

    if (a_state.nComp() == 1) {

        for (dit.reset(); dit.ok(); ++dit) {
            FluxBox& stateFB = a_state[dit];
            const FluxBox* JgupFBPtr = (a_JgupPtr == NULL)?
                                       NULL:
                                       &((*a_JgupPtr)[dit]);

            for (int FCdir = 0; FCdir < SpaceDim; ++FCdir) {
                const FArrayBox* extrapFABPtr = (a_extrapPtr == NULL)?
                                                NULL:
                                                &((*a_extrapPtr)[dit][FCdir]);

                Box bx = a_state[dit][FCdir].box();
                bx &= surroundingNodes(domain.domainBox(), FCdir);

                Interval interv(0,0);
                m_componentBC[FCdir].setGhosts(stateFB[FCdir],
                                               extrapFABPtr,
                                               bx,
                                               domain,
                                               a_dx,
                                               dit(),
                                               JgupFBPtr,
                                               a_homogeneous,
                                               a_time,
                                               interv);
            } // end loop over FCdir
        } // end loop over grids (dit)

    } else if (a_state.nComp() == SpaceDim) {

        for (dit.reset(); dit.ok(); ++dit) {
            FluxBox& stateFB = a_state[dit];
            const FluxBox* JgupFBPtr = (a_JgupPtr == NULL)?
                                       NULL:
                                       &((*a_JgupPtr)[dit]);

            for (int FCdir = 0; FCdir < SpaceDim; ++FCdir) {
                const FArrayBox* extrapFABPtr = (a_extrapPtr == NULL)?
                                                NULL:
                                                &((*a_extrapPtr)[dit][FCdir]);

                Box bx = a_state[dit][FCdir].box();
                bx &= surroundingNodes(domain.domainBox(), FCdir);

                for (int comp = 0; comp < stateFB[FCdir].nComp(); ++comp) {
                    Interval interv(comp,comp);
                    m_componentBC[comp].setGhosts(stateFB[FCdir],
                                                  extrapFABPtr,
                                                  bx,
                                                  domain,
                                                  a_dx,
                                                  dit(),
                                                  JgupFBPtr,
                                                  a_homogeneous,
                                                  a_time,
                                                  interv);
                } // end loop over state comps
            } // end loop over FCdir
        } // end loop over grids (dit)
    } else {
        MayDay::Error("a_state must have 1 or SpaceDim comps");
    }
}
