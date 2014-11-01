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
#ifndef AMRCCProjector_H__INCLUDED__
#define AMRCCProjector_H__INCLUDED__

#include "BaseProjector.H"
#include "MappedQuadCFInterp.H"
#include "MappedCoarseAverage.H"
#include "PhysBCUtil.H"
#include "Divergence.H"
#include "Gradient.H"


// -----------------------------------------------------------------------------
// Performs the composite projections during initialization and synchronization.
// -----------------------------------------------------------------------------
class AMRCCProjector: public BaseProjector<FArrayBox>
{
public:
    // Default constructor
    // This sets the solver parameters, but leaves object unusable.
    AMRCCProjector ();

    // Destructor
    virtual ~AMRCCProjector ();

    // Allocates memory and leaves object usable.
    virtual void define (Vector<LevelData<FArrayBox>*>& a_amrPressure,
                         const PhysBCUtil&              a_physBCUtil,
                         const LevelGeometry&           a_levGeo,
                         const FillJgupInterface*       a_customFillJgupPtr = NULL);

    // Frees memory and leaves object unusable.
    virtual void undefine ();

    // Is this object in a usable state?
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

    // Is this object in a usable state?
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
// Is this object in a usable state?
// -----------------------------------------------------------------------------
inline bool AMRCCProjector::isDefined () const
{
    return m_isDefined;
}


#endif //!AMRCCProjector_H__INCLUDED__
