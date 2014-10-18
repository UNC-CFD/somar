/*******************************************************************************
 *    FILE: LevelMACProjector.H
 *    DESC: Projects a single-level, FC velocity field.
 *    DATE: Mon 14 Apr 2014 09:26:59 AM EDT  (UTC-0400)
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef LevelMACProjector_H__INCLUDED__
#define LevelMACProjector_H__INCLUDED__

#include "BaseProjector.H"
#include "FluxBox.H"
#include "MappedQuadCFInterp.H"
#include "EllipticBCInterface.H"
#include "PhysBCUtil.H"
#include "Divergence.H"
#include "Gradient.H"


class LevelMACProjector: public BaseProjector<FluxBox>
{
public:
    // Default constructor
    // This sets the solver parameters, but leaves object unusable.
    LevelMACProjector ();

    // Destructor
    virtual ~LevelMACProjector ();

    // Allocates memory and leaves object useable.
    virtual void define (LevelData<FArrayBox>*       a_phiPtr,
                         const LevelData<FArrayBox>* a_crsePhiPtr,
                         const PhysBCUtil&           a_physBCUtil,
                         const LevelGeometry&        a_levGeo,
                         const FillJgupInterface*    a_customFillJgupPtr = NULL);

    // Frees memory and leaves object unuseable.
    virtual void undefine ();

    // Is this object in a useable state?
    virtual inline bool isDefined() const;

protected:
    // From BaseProjector:
    //  Computes J^{-1}\partial_i(J u^i) over an AMR hierarchy.
    //  This must be overriden or an error will be thrown.
    virtual void computeDiv (Vector<LevelData<FArrayBox>*>&     a_div,
                             const Vector<LevelData<FluxBox>*>& a_flux,
                             const int                          a_lmin,
                             const int                          a_lmax) const;

    // From BaseProjector:
    //  Computes Jg^{i,j}\partial_j(phi) over an AMR hierarchy.
    //  This must be overriden or an error will be thrown.
    virtual void computeGrad (Vector<LevelData<FluxBox>*>&         a_flux,
                              const Vector<LevelData<FArrayBox>*>& a_phi,
                              const int                            a_lmin,
                              const int                            a_lmax) const;

    // From BaseProjector:
    //  Applies vel = vel - dt*corr.
    virtual void applyCorrection (Vector<LevelData<FluxBox>*>&       a_amrVel,
                                  const Vector<LevelData<FluxBox>*>& a_amrCorr,
                                  const Real                         a_dt,
                                  const int                          a_lmin,
                                  const int                          a_lmax) const;

    // Is this object in a useable state?
    bool m_isDefined;

    const LevelGeometry* m_levGeoPtr;
    MappedQuadCFInterp   m_cfInterp;

    BCMethodHolder       m_solverBC;
    Divergence::BC_type  m_divBC;
    Gradient::BC_type    m_gradBC;

private:
    // Copy and assignment not allowed
    LevelMACProjector (const LevelMACProjector&);
    void operator= (const LevelMACProjector&);
};


// -----------------------------------------------------------------------------
// Is this object in a useable state?
// -----------------------------------------------------------------------------
inline bool LevelMACProjector::isDefined () const
{
    return m_isDefined;
}


#endif //!LevelMACProjector_H__INCLUDED__
