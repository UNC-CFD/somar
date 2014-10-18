#include "MappedQuadCFInterp.H"
#include "MappedQuadCFInterpF_F.H"
#include "MappedCFStencil.H"
#include "AnisotropicRefinementTools.H"
#include "CFIVS.H"
#include "TensorCFInterp.H" // Only needed for gradIndex static utility.

bool MappedQuadCFInterp::newCFInterMode = true;


// -----------------------------------------------------------------------------
// Default constructor leaves QCFI undefined.
// -----------------------------------------------------------------------------
MappedQuadCFInterp::MappedQuadCFInterp ()
{
    this->clear();
}


// -----------------------------------------------------------------------------
// Full Constructor. Makes all coarse-fine information and sets internal
// variables calls full define.
// -----------------------------------------------------------------------------
MappedQuadCFInterp::MappedQuadCFInterp (const DisjointBoxLayout& a_fineBoxes,
                                        const DisjointBoxLayout* a_coarBoxes,
                                        const RealVect&          a_dxFine,
                                        const IntVect&           a_refRatio,
                                        const int                a_nComp,
                                        const Box&               a_domf)
{
    ProblemDomain fineProbDomain(a_domf);
    this->define(a_fineBoxes, a_coarBoxes,
                 a_dxFine, a_refRatio,
                 a_nComp, fineProbDomain);
}


// -----------------------------------------------------------------------------
// Full Constructor. Makes all coarse-fine information and sets internal
// variables calls full define.
// -----------------------------------------------------------------------------
MappedQuadCFInterp::MappedQuadCFInterp (const DisjointBoxLayout& a_fineBoxes,
                                        const DisjointBoxLayout* a_coarBoxes,
                                        const RealVect&          a_dxFine,
                                        const IntVect&           a_refRatio,
                                        const int                a_nComp,
                                        const ProblemDomain&     a_domf)
{
    this->define(a_fineBoxes, a_coarBoxes,
                 a_dxFine, a_refRatio,
                 a_nComp, a_domf);
}


// -----------------------------------------------------------------------------
// Destructor
// -----------------------------------------------------------------------------
MappedQuadCFInterp::~MappedQuadCFInterp ()
{
    this->clear();
}


// -----------------------------------------------------------------------------
// Full define function. Makes all coarse-fine information and sets
// internal variables.
// -----------------------------------------------------------------------------
void MappedQuadCFInterp::define (const DisjointBoxLayout& a_fineBoxes,
                                 const DisjointBoxLayout* a_coarBoxesPtr,
                                 const RealVect&          a_dxLevel,
                                 const IntVect&           a_refRatio,
                                 const int                a_nComp,
                                 const ProblemDomain&     a_domf)
{
    CH_TIME("MappedQuadCFInterp::define");

    this->clear();
    m_isDefined = true;

    CH_assert(a_nComp > 0);
    CH_assert(!a_domf.isEmpty());
    CH_assert(a_fineBoxes.checkPeriodic(a_domf));

    m_domainFine = a_domf;
    m_dxFine = a_dxLevel;
    m_refRatio = a_refRatio;
    m_nComp = a_nComp;
    m_inputFineLayout = a_fineBoxes;

    int dims = SpaceDim;
    if (a_refRatio[SpaceDim-1] == 1 && a_domf.size(SpaceDim-1) == 1) {
        m_isFlat = true;
        --dims;
    } else {
        m_isFlat = false;
    }

    if (SpaceDim < 3 && m_isFlat) {
        MayDay::Error("MappedQuadCFInterp and the CFStencils have not been tested when SpaceDim < 3 and m_isFlat! Just wanted to let you know");
    }

    bool fineCoversCoarse = false;
    if (a_coarBoxesPtr != NULL) {
        int factor = a_refRatio.product();
        long long numPts = a_fineBoxes.numCells() / factor;
        numPts -= a_coarBoxesPtr->numCells();
        if (numPts == 0) fineCoversCoarse = true;
    }
    m_fineCoversCoarse = fineCoversCoarse;

    // This is whey m_level is "fake." It is always 0 or 1.
    m_level = (a_coarBoxesPtr == NULL || fineCoversCoarse)? 0: 1;

    if (m_level > 0) {
        // (DFM) only check for valid refRatio if a coarser level exists
        D_TERM(CH_assert(a_refRatio[0] >= 1);,
               CH_assert(a_refRatio[1] >= 1);,
               CH_assert(a_refRatio[2] >= 1);)

        const DisjointBoxLayout& coarBoxes = *a_coarBoxesPtr;
        m_inputCoarLayout = coarBoxes;

#ifndef NDEBUG
        {
            ProblemDomain crseDomF;
            coarsen(crseDomF, a_domf, a_refRatio);
            CH_assert(coarBoxes.checkPeriodic(crseDomF));
        }
#endif

        for (int idir = 0; idir < SpaceDim; ++idir) {
            m_loQCFS[idir].define(a_fineBoxes);
            m_hiQCFS[idir].define(a_fineBoxes);
        }

        //locoarboxes and hicoarboxes are now open
        //and have same processor mapping as a_fineboxes
        //make boxes for coarse buffers
        m_coarBoxes = DisjointBoxLayout();
        coarsen(m_coarBoxes, a_fineBoxes, m_refRatio);
        if (!m_isFlat) {
            m_coarBoxes.grow(2);
        } else {
            m_coarBoxes.grow(2*(IntVect::Unit-BASISV(CH_SPACEDIM-1)));
        }

        m_coarBoxes.close();
        m_coarBuffer.define(m_coarBoxes, m_nComp);
        m_copier.define(coarBoxes, m_coarBoxes);

        if (!newCFInterMode) { //old n^2 algorithm (bvs)
            //make cfstencils and boxes for coarse buffers
            DataIterator dit = a_fineBoxes.dataIterator();
            for (dit.begin(); dit.ok(); ++dit) {
                const Box& fineBox = a_fineBoxes[dit()];
                for (int idir = 0; idir < dims; ++idir) {
                    //low side cfstencil
                    m_loQCFS[idir][dit()].define(a_domf,
                                                 fineBox,
                                                 a_fineBoxes,
                                                 coarBoxes,
                                                 a_refRatio,
                                                 idir,
                                                 Side::Lo,
                                                 m_isFlat);
                    //high side cfstencil
                    m_hiQCFS[idir][dit()].define(a_domf,
                                                 fineBox,
                                                 a_fineBoxes,
                                                 coarBoxes,
                                                 a_refRatio,
                                                 idir,
                                                 Side::Hi,
                                                 m_isFlat);

                }
            }
        } else {
            //new "moving window" version of CF stencil building
            Vector<Box> periodicFine;
            CFStencil::buildPeriodicVector(periodicFine, a_domf, a_fineBoxes);
            Vector<Box> coarsenedFine(periodicFine);
            for (int i = 0; i < coarsenedFine.size(); ++i) {
                coarsenedFine[i].coarsen(a_refRatio);
            }
            DataIterator dit = a_fineBoxes.dataIterator();
            for (dit.begin(); dit.ok(); ++dit) {
                const Box& fineBox = a_fineBoxes[dit()];
                for (int idir = 0; idir < dims; ++idir) {
                    //low side cfstencil
                    m_loQCFS[idir][dit()].define(a_domf,
                                                 fineBox,
                                                 periodicFine,
                                                 coarsenedFine,
                                                 coarBoxes,
                                                 a_refRatio,
                                                 idir,
                                                 Side::Lo,
                                                 m_isFlat);
                    //high side cfstencil
                    m_hiQCFS[idir][dit()].define(a_domf,
                                                 fineBox,
                                                 periodicFine,
                                                 coarsenedFine,
                                                 coarBoxes,
                                                 a_refRatio,
                                                 idir,
                                                 Side::Hi,
                                                 m_isFlat);
                }
            }
        }
    }
}


// -----------------------------------------------------------------------------
// Has full define function been called?
// -----------------------------------------------------------------------------
bool MappedQuadCFInterp::isDefined () const
{
    return m_isDefined;
}


// -----------------------------------------------------------------------------
// Return QCFI to undefined state
// -----------------------------------------------------------------------------
void MappedQuadCFInterp::clear ()
{
    m_isDefined = false;
    m_level = -1;
    m_dxFine = RealVect(D_DECL(-1.0,-1.0,-1.0));
}


// -----------------------------------------------------------------------------
// Apply coarse-fine boundary conditions -- assume that phi grids
// are grown by one.
// -----------------------------------------------------------------------------
void MappedQuadCFInterp::coarseFineInterp (BaseFab<Real>&             a_phif,
                                           const BaseFab<Real>&       a_phic,
                                           const MappedQuadCFStencil& a_qcfs,
                                           const Side::LoHiSide       a_hiorlo,
                                           const int                  a_idir,
                                           const Interval&            a_variables) const
{
    CH_TIME("MappedQuadCFInterp::coarseFineInterp(BaseFab<Real> & a_phif,...)");
    CH_assert(isDefined());
    CH_assert(!m_isFlat || a_idir < CH_SPACEDIM-1);

    //nothing happens if m_level == 0
    if (m_level > 0) {
        if (!a_qcfs.isEmpty()) {
            BaseFab<Real>  phistar;

            //first find extended value phistar
            //includes finding slopes of coarse solution bar
            getPhiStar(phistar, a_phic, a_qcfs, a_hiorlo, a_idir, a_variables);

            //given phistar, interpolate on fine ivs
            interpOnIVS(a_phif, phistar, a_qcfs, a_hiorlo, a_idir, a_variables);
        }
    }
}


// -----------------------------------------------------------------------------
// Get extended phi (lives next to interpivs)
// -----------------------------------------------------------------------------
void MappedQuadCFInterp::getPhiStar (BaseFab<Real>&             a_phistar,
                                     const BaseFab<Real>&       a_phic,
                                     const MappedQuadCFStencil& a_qcfs,
                                     const Side::LoHiSide       a_hiorlo,
                                     const int                  a_idir,
                                     const Interval&            a_variables) const
{
    CH_TIMERS("MappedQuadCFInterp::getPhiStar");
    //CH_TIMER("MappedQuadCFInterp::computeFirstDerivative",  t1st);
    //CH_TIMER("MappedQuadCFInterp::computesecondDerivative", t2nd);
    //CH_TIMER("MappedQuadCFInterp::computemixedDerivative",  tmixed);
    CH_TIMER("MappedQuadCFInterp::slopes", tslopes);
    CH_TIMER("MappedQuadCFInterp::notPacked", tnp);
    CH_TIMER("MappedQuadCFInterp::preamble", tpreamble);
    CH_assert(isDefined());
    CH_assert(a_qcfs.isDefined());
#if (CH_SPACEDIM > 1)
    const RealVect dxf = m_dxFine;
    const RealVect dxc = m_dxFine * RealVect(m_refRatio);
#endif

    // For flat 3D grids.
    int dims = SpaceDim;
    if (m_isFlat) --dims;

    // if we think of a_idir as the "me" direction, then
    // the other directions can be "you1" and "you2"
#if (CH_SPACEDIM == 3)
    int you1, you2;
    if (!m_isFlat) {
        // Full 3D...
        if (a_idir == 0) {
            you1 = 1;
            you2 = 2;
        } else if (a_idir == 1) {
            you1 = 0;
            you2 = 2;
        } else
        {
            you1 = 0;
            you2 = 1;
        }
    } else {
        // Flat 3D...
        you1 = (a_idir + 1) % dims;
        you2 = -1;
    }
#else // (CH_SPACEDIM == 2)
    int you1;
    if (a_idir == 0) {
        you1 = 1;
    } else
    {
        you1 = 0;
    }
#endif

    //if cfsten is empty, nothing to interpolate.
    if (!a_qcfs.isEmpty()) {
        CH_START(tpreamble);
        CH_assert(m_level > 0);
        const IntVectSet& interp_ivs = a_qcfs.getFineIVS();
        const IntVectSet& coarsl_ivs = a_qcfs.getCoarIVS();
        if (!coarsl_ivs.isDense()) {
            MayDay::Error("What the hell?? TreeIntVectSet ???");
        }
        if (!interp_ivs.isDense()) {
            MayDay::Error("What the hell?? TreeIntVectSet ???");
        }

        Box coarinterpbox = coarsl_ivs.minBox();
        int ncomp = a_phic.nComp();
        CH_assert(ncomp == m_nComp);
        CH_assert(a_phic.box().contains((coarinterpbox)));

        // allocate phistar here
        int ihilo = sign(a_hiorlo);
        Box phistarbox = interp_ivs.minBox();
        phistarbox.shift(a_idir, ihilo);
        a_phistar.resize(phistarbox, ncomp);
        CH_STOP(tpreamble);
        for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {

            CH_START(tslopes);
            //phi = phino + slope*x + half*x*x*curvature
            BaseFab<Real>  coarslope(coarinterpbox, dims);
            BaseFab<Real>  coarcurva(coarinterpbox, dims);
#if (CH_SPACEDIM == 3)
            BaseFab<Real>  coarmixed;
            if (!m_isFlat) {
                coarmixed.define(coarinterpbox, 1);
            }
#endif

            // coarslope.setVal(0.);

            //first find extended value phistar. get slopes of coarse solution
            IVSIterator coar_ivsit(coarsl_ivs);

            if (SpaceDim == 3 && m_isFlat) {
                // Flat 3D...
#if (CH_SPACEDIM > 1)
                for (coar_ivsit.begin(); coar_ivsit.ok(); ++coar_ivsit) {
                    // this isn't relevant for 1D
                    const IntVect& coariv = coar_ivsit();

                    coarslope(coariv, you1) =
                        a_qcfs.computeFirstDerivative (a_phic, you1, ivar, coariv, dxc[you1]);
                    coarcurva(coariv, you1) =
                        a_qcfs.computeSecondDerivative(a_phic, you1, ivar, coariv, dxc[you1]);
                } //end loop over coarse intvects
#endif
            } else {
                // Full ND...
                for (coar_ivsit.begin(); coar_ivsit.ok(); ++coar_ivsit) {
                    // this isn't relevant for 1D
#if (CH_SPACEDIM > 1)
                    const IntVect& coariv = coar_ivsit();

                    // coarslope(coariv, a_idir) = 0.0;
                    // coarcurva(coariv, a_idir) = 0.0;

                    coarslope(coariv, you1) =
                        a_qcfs.computeFirstDerivative (a_phic, you1, ivar, coariv, dxc[you1]);
                    coarcurva(coariv, you1) =
                        a_qcfs.computeSecondDerivative(a_phic, you1, ivar, coariv, dxc[you1]);

#endif

#if (CH_SPACEDIM == 3)
                    coarslope(coariv, you2) =
                        a_qcfs.computeFirstDerivative (a_phic, you2, ivar, coariv, dxc[you2]);
                    coarcurva(coariv, you2) =
                        a_qcfs.computeSecondDerivative(a_phic, you2, ivar, coariv, dxc[you2]);
                    coarmixed(coariv) =
                        a_qcfs.computeMixedDerivative(a_phic, ivar, coariv, dxc);
#endif
                } //end loop over coarse intvects
            }
            CH_STOP(tslopes);

            if (a_qcfs.finePacked() && CH_SPACEDIM == 3 && !m_isFlat) {
                const IntVect& iv = phistarbox.smallEnd();
                IntVect civ(iv);
                civ.coarsen(m_refRatio);
                Box region = a_qcfs.packedBox();
#if (CH_SPACEDIM == 3)
                FORT_MAPPEDPHISTAR(CHF_FRA_SHIFT(a_phistar, iv),
                                   CHF_BOX_SHIFT(region, iv),
                                   CHF_CONST_FRA_SHIFT(a_phic, civ),
                                   CHF_FRA_SHIFT(coarslope, civ),
                                   CHF_FRA_SHIFT(coarcurva, civ),
                                   CHF_FRA_SHIFT(coarmixed, civ),
                                   CHF_CONST_REALVECT(dxf),
                                   CHF_CONST_INT(ivar),
                                   CHF_CONST_INT(a_idir),
                                   CHF_CONST_INT(ihilo),
                                   CHF_CONST_INTVECT(m_refRatio));
#endif
            } else {
                CH_START(tnp);
                IntVect ivf, ivc, ivstar;
                // ifdef is here to prevent unused variable warnings in 1D
#if (CH_SPACEDIM > 1)
                int jf, jc;
                Real xf, xc, x1;
#endif
                Real pc, update1 = 0, update2 = 0, update3 = 0;
                IVSIterator fine_ivsit(interp_ivs);

                if (SpaceDim == 3 && m_isFlat) {
                    // Flat 3D...
                    for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit) {
                        ivf = fine_ivsit();
                        ivc = coarsen(ivf, m_refRatio);
                        ivstar = ivf;
                        ivstar.shift(a_idir, ihilo);
                        pc = a_phic(ivc, ivar);

                        // for 1D, none of this is necessary -- just copy
                        // coarse value into phiStar
#if (CH_SPACEDIM > 1)
                        jf = ivf[you1];
                        jc = ivc[you1];
                        xf = (jf + 0.5) * dxf[you1];
                        xc = (jc + 0.5) * dxc[you1];
                        x1 = xf - xc;
                        update1 = x1 * coarslope(ivc, you1) + 0.5 * x1 * x1 * coarcurva(ivc, you1);
#endif
                        a_phistar(ivstar, ivar) = pc + update1;
                    } //end loop over fine intvects
                } else {
                    // Full ND...
                    for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit) {
                        ivf = fine_ivsit();
                        ivc = coarsen(ivf, m_refRatio);
                        ivstar = ivf;
                        ivstar.shift(a_idir, ihilo);
                        pc = a_phic(ivc, ivar);

                        // for 1D, none of this is necessary -- just copy
                        // coarse value into phiStar
#if (CH_SPACEDIM > 1)
                        if ( !(SpaceDim == 2 && m_isFlat) ) {
                            jf = ivf[you1];
                            jc = ivc[you1];
                            xf = (jf + 0.5) * dxf[you1];
                            xc = (jc + 0.5) * dxc[you1];
                            x1 = xf - xc;
                            update1 = x1 * coarslope(ivc, you1) + 0.5 * x1 * x1 * coarcurva(ivc, you1);
                        }
#endif

#if (CH_SPACEDIM==3)
                        Real x2;
                        jf = ivf[you2];
                        jc = ivc[you2];
                        xf = (jf + 0.5) * dxf[you2];
                        xc = (jc + 0.5) * dxc[you2];
                        x2 = xf - xc;
                        update2 =  x2 * coarslope(ivc, you2) + 0.5 * x2 * x2 * coarcurva(ivc, you2);

                        //add in mixed derivative component
                        update3 =  x1 * x2 * coarmixed(ivc);
#endif
                        a_phistar(ivstar, ivar) = pc + update1 + update2 + update3;
                    } //end loop over fine intvects
                }

                CH_STOP(tnp);
            } // end if for not packed optimization
        }//end loop over variables
    } //end if (level>0 && !hocfs.isempty())
} //end function getphistar


// -----------------------------------------------------------------------------
// Interpolate over correct intvectset
// -----------------------------------------------------------------------------
void MappedQuadCFInterp::interpOnIVS (BaseFab<Real>&             a_phif,
                                      const BaseFab<Real>&       a_phistar,
                                      const MappedQuadCFStencil& a_qcfs,
                                      const Side::LoHiSide       a_hiorlo,
                                      const int                  a_idir,
                                      const Interval&            a_variables) const
{
    CH_TIME("MappedQuadCFInterp::interpOnIVS");
    CH_assert(isDefined());
    CH_assert(a_qcfs.isDefined());
    CH_assert(!m_isFlat || a_idir < CH_SPACEDIM-1);

    //if cfsten is empty, nothing to interpolate.
    if (!a_qcfs.isEmpty()) {
        //if there IS something to interpolate, the level ident
        //had better be greater than zero.  Otherwise a null
        //was sent in as coarse grids on construction
        CH_assert(m_level > 0);
        const IntVectSet& interp_ivs = a_qcfs.getFineIVS();

        int ihilo = sign(a_hiorlo);
        int nref = m_refRatio[a_idir];  // TODO: Is this right?
        Real dxf = m_dxFine[a_idir];    // TODO: Is this right?

        if (!a_qcfs.finePacked()) {
            IVSIterator fine_ivsit(interp_ivs);
            CH_assert(a_phistar.nComp() == a_phif.nComp());
            CH_assert(a_phistar.nComp() == m_nComp);
            const IntVect e = ihilo * BASISV(a_idir);

            for (fine_ivsit.begin(); fine_ivsit.ok(); ++fine_ivsit) {
                IntVect ivf = fine_ivsit();
                // quadratic interpolation
                for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {
                    Real pa =    a_phif(ivf - 2*e, ivar);
                    Real pb =    a_phif(ivf -   e, ivar);
                    Real ps = a_phistar(ivf +   e, ivar);
                    //phi = ax**2 + bx + c, x = 0 at pa
                    Real h = dxf;
                    Real a = (2. / h / h) * (2.*ps + pa * (nref + 1.0) - pb * (nref + 3.0)) /
                             (nref * nref + 4 * nref + 3.0);
                    Real b = (pb - pa) / h - a * h;
                    Real c = pa;
                    Real x = 2.*h;
                    a_phif (ivf, ivar) = a * x * x + b * x + c;
                } //end loop over components
            } //end loop over fine intvects
        } else {
            // data is packed, just call Fortran
            int b = a_variables.begin();
            int e = a_variables.end();
            FORT_MAPPEDQUADINTERP(CHF_FRA(a_phif),
                                  CHF_CONST_FRA(a_phistar),
                                  CHF_BOX(a_qcfs.packedBox()),
                                  CHF_CONST_INT(ihilo),
                                  CHF_CONST_REAL(dxf),
                                  CHF_CONST_INT(a_idir),
                                  CHF_CONST_INT(b),
                                  CHF_CONST_INT(e),
                                  CHF_CONST_INT(nref));
        }
    } //end if (level>0 && !oscfs.isempty())
} //end function interponivs


// -----------------------------------------------------------------------------
// Apply coarse-fine boundary conditions -- assume that phi grids
// are grown by one.
// -----------------------------------------------------------------------------
void MappedQuadCFInterp::coarseFineInterp (LevelData<FArrayBox>&       a_phif,
                                           const LevelData<FArrayBox>& a_phic) const
{
    CH_TIME("MappedQuadCFInterp::coarseFineInterp");
    CH_assert(isDefined());

    Interval variables = a_phic.interval();

    if (m_level > 0) {
        CH_assert(a_phif.nComp() == m_nComp);
        CH_assert(a_phic.nComp() == m_nComp);
        CH_assert(a_phic.boxLayout() == m_inputCoarLayout);
        CH_assert(a_phif.boxLayout() == m_inputFineLayout);

        a_phic.copyTo(a_phic.interval(),
                      m_coarBuffer,
                      m_coarBuffer.interval(),
                      m_copier);

        for (int idir = 0; idir < CH_SPACEDIM; ++idir) {
            if (a_phif.ghostVect()[idir] == 0) continue;

            DataIterator ditFine = a_phif.dataIterator();
            for (ditFine.begin(); ditFine.ok(); ++ditFine) {
                DataIndex datIndGlo = ditFine();
                BaseFab<Real>& phif = a_phif[datIndGlo];
                const BaseFab<Real>& phiC = m_coarBuffer[datIndGlo];
                //lo side cfinterp
                //recall that buffers have fine processor mapping
                {
                    const MappedQuadCFStencil& loQCFS = m_loQCFS[idir][datIndGlo];
                    coarseFineInterp(phif, phiC, loQCFS, Side::Lo, idir, variables);

                }
                //hi side cfinterp
                {
                    const MappedQuadCFStencil& hiQCFS = m_hiQCFS[idir][datIndGlo];
                    coarseFineInterp(phif, phiC, hiQCFS, Side::Hi, idir, variables);
                }
            }//end iteration over boxes in fine grid
        } //end iteration over directions
    }
}

