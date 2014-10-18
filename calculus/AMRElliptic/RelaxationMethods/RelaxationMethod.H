/*******************************************************************************
 *    FILE: RelaxationMethod.H
 *    DESC: Interface for all relaxation methods. This does not conform to the
 *          Chombo LinearOp standard. Children of this class are meant to be
 *          used in MappedAMRPoissonOp.
 *    DATE: Fri 22 Nov 2013 03:56:44 PM EST
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __RelaxationMethod_H__INCLUDED__
#define __RelaxationMethod_H__INCLUDED__

#include "Copier.H"
#include "CFRegion.H"
#include "LevelData.H"
#include "FArrayBox.H"
#include "FluxBox.H"
#include "IntVect.H"
#include "RealVect.H"
#include "RefCountedPtr.H"
#include "EllipticBCInterface.H"


class RelaxationMethod
{
public:
    // The relaxation methods. If you add a method,
    // you should add it to this list.
    enum {
        NORELAX     = -1,
        JACOBI      =  0,
        LEVEL_GSRB  =  1,
        LOOSE_GSRB  =  2,
        LINE_GSRB   =  3,
        NUM_RELAXATION_METHODS
    };

    // Constructor -- leaves object unusable.
    RelaxationMethod ();

    // Destructor
    virtual ~RelaxationMethod ();

    // Define constructor.
    virtual void define (const Real                                  a_alpha,
                         const Real                                  a_beta,
                         const RealVect&                             a_dx,
                         const RealVect&                             a_crseDx,
                         const RefCountedPtr<LevelData<FluxBox> >&   a_FCJgup,
                         const RefCountedPtr<LevelData<FArrayBox> >& a_CCJinv,
                         const RefCountedPtr<LevelData<FArrayBox> >& a_lapDiag,
                         BCMethodHolder&                             a_bc,
                         const Copier&                               a_exchangeCopier,
                         const CFRegion&                             a_cfRegion,
                         const bool                                  a_isDiagonal,
                         const IntVect&                              a_activeDirs = IntVect::Unit);

    // Is this object ready for use?
    virtual inline bool isDefined () const;

    // The relaxation method.
    virtual void relax (LevelData<FArrayBox>&       a_phi,
                        const LevelData<FArrayBox>& a_rhs) = 0;

    // This allows the TGA solver to change alpha and beta.
    virtual inline void setAlphaAndBeta (const Real a_alpha,
                                         const Real a_beta);

protected:
    // This collects the boundary box data needed by (simple)boundaryGSRB.
    virtual void collectBoundaryData ();

    // As advertised, this will fill the ghosts of a_phi at the physical
    // boundary using the BCs and at the CF interface using a homogeneous
    // coarser level. This function will also extrapolate a_phi's valid data to
    // help compute mixed second derivatives when the metric is not diagonal.
    // This function is not blocking and the exchange is left to the caller.
    virtual void fillGhostsAndExtrapolate (LevelData<FArrayBox>& a_phi,
                                           const bool            a_doCFInterp = true,
                                           const bool            a_doExtrapAndCopy = true,
                                           const bool            a_doBCs = true);

protected:
    bool            m_isDefined;

    BCMethodHolder  m_bc;
    Copier          m_exchangeCopier;
    CFRegion        m_cfRegion;
    IntVect         m_activeDirs;
    bool            m_isDiagonal;

    Real            m_alpha;
    Real            m_beta;

    RealVect        m_dx;
    RealVect        m_crseDx;

    RefCountedPtr<LevelData<FluxBox> >   m_FCJgup;
    RefCountedPtr<LevelData<FArrayBox> > m_CCJinv;
    RefCountedPtr<LevelData<FArrayBox> > m_lapDiag;

    LevelData<FArrayBox> m_extrap;

    typedef struct {
        DataIndex index;
        Box       validBdry;
        Box       valid;
        int       stencil[CH_SPACEDIM][2];
    } BoundaryBoxData;

    Vector<BoundaryBoxData> m_simpleBoundaryBoxData;
    Vector<BoundaryBoxData> m_boundaryBoxData;
};


// -----------------------------------------------------------------------------
// Is this object ready for use?
// -----------------------------------------------------------------------------
bool RelaxationMethod::isDefined () const
{
    return m_isDefined;
}


// -----------------------------------------------------------------------------
// This allows the TGA solver to change alpha and beta.
// -----------------------------------------------------------------------------
void RelaxationMethod::setAlphaAndBeta (const Real a_alpha,
                                        const Real a_beta)
{
    m_alpha = a_alpha;
    m_beta = a_beta;
}


#endif //!__RelaxationMethod_H__INCLUDED__
