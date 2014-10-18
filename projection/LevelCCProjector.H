/*******************************************************************************
 *    FILE: LevelCCProjector.H
 *    DESC: Approximately projects a single-level, CC velocity field.
 *    DATE: Fri 02 May 2014 02:49:13 PM EDT  (UTC-0400)
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef LevelCCProjector_H__INCLUDED__
#define LevelCCProjector_H__INCLUDED__

#include "BaseProjector.H"
#include "MappedQuadCFInterp.H"
#include "EllipticBCInterface.H"
#include "PhysBCUtil.H"
#include "Divergence.H"
#include "Gradient.H"


class LevelCCProjector: public BaseProjector<FArrayBox>
{
public:
    // Default constructor
    // This sets the solver parameters, but leaves object unusable.
    LevelCCProjector ();

    // Destructor
    virtual ~LevelCCProjector ();

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
    virtual void computeDiv (Vector<LevelData<FArrayBox>*>&       a_div,
                             const Vector<LevelData<FArrayBox>*>& a_flux,
                             const int                            a_lmin,
                             const int                            a_lmax) const;

    // From BaseProjector:
    //  Computes Jg^{i,j}\partial_j(phi) over an AMR hierarchy.
    virtual void computeGrad (Vector<LevelData<FArrayBox>*>&       a_flux,
                              const Vector<LevelData<FArrayBox>*>& a_phi,
                              const int                            a_lmin,
                              const int                            a_lmax) const;

    // From BaseProjector:
    //  Applies vel = vel - dt*corr.
    virtual void applyCorrection (Vector<LevelData<FArrayBox>*>&       a_amrVel,
                                  const Vector<LevelData<FArrayBox>*>& a_amrCorr,
                                  const Real                           a_dt,
                                  const int                            a_lmin,
                                  const int                            a_lmax) const;

    // Is this object in a useable state?
    bool m_isDefined;

    const LevelGeometry* m_levGeoPtr;
    MappedQuadCFInterp   m_pressureCFInterp;
    MappedQuadCFInterp   m_velCFInterp;

    BCMethodHolder       m_solverBC;
    Divergence::BC_type  m_divBC;
    Gradient::BC_type    m_gradBC;

private:
    // Copy and assignment not allowed
    LevelCCProjector (const LevelCCProjector&);
    void operator= (const LevelCCProjector&);
};


// -----------------------------------------------------------------------------
// Is this object in a useable state?
// -----------------------------------------------------------------------------
inline bool LevelCCProjector::isDefined () const
{
    return m_isDefined;
}


#endif //!LevelCCProjector_H__INCLUDED__
