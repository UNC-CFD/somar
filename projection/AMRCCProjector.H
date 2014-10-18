/*******************************************************************************
 *    FILE: AMRCCProjector.H
 *    DESC: Performs the composite projections during initialization and
 *          synchronization.
 *    DATE: Sun 11 May 2014 11:36:32 PM EDT  (UTC-0400)
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef AMRCCProjector_H__INCLUDED__
#define AMRCCProjector_H__INCLUDED__

#include "BaseProjector.H"
#include "MappedQuadCFInterp.H"
#include "MappedCoarseAverage.H"
#include "PhysBCUtil.H"
#include "Divergence.H"
#include "Gradient.H"


class AMRCCProjector: public BaseProjector<FArrayBox>
{
public:
    // Default constructor
    // This sets the solver parameters, but leaves object unusable.
    AMRCCProjector ();

    // Destructor
    virtual ~AMRCCProjector ();

    // Allocates memory and leaves object useable.
    virtual void define (Vector<LevelData<FArrayBox>*>& a_amrPressure,
                         const PhysBCUtil&              a_physBCUtil,
                         const LevelGeometry&           a_levGeo,
                         const FillJgupInterface*       a_customFillJgupPtr = NULL);

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

    // Applies vel = vel - dt*corr.
    virtual void applyCorrection (Vector<LevelData<FArrayBox>*>&       a_amrVel,
                                  const Vector<LevelData<FArrayBox>*>& a_amrCorr,
                                  const Real                           a_dt,
                                  const int                            a_lmin,
                                  const int                            a_lmax) const;

    // Is this object in a useable state?
    bool m_isDefined;

    Vector<const LevelGeometry*>               m_amrLevGeo;
    Vector<RefCountedPtr<MappedQuadCFInterp> > m_pressureCFInterp;
    Vector<RefCountedPtr<MappedQuadCFInterp> > m_velCFInterp;

    BCMethodHolder       m_solverBC;
    Divergence::BC_type  m_divBC;
    Gradient::BC_type    m_gradBC;

    bool m_applySyncProjection;

private:
    // Copy and assignment not allowed
    AMRCCProjector (const AMRCCProjector&);
    void operator= (const AMRCCProjector&);
};



// -----------------------------------------------------------------------------
// Is this object in a useable state?
// -----------------------------------------------------------------------------
inline bool AMRCCProjector::isDefined () const
{
    return m_isDefined;
}


#endif //!AMRCCProjector_H__INCLUDED__
