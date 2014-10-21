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
#ifndef LevelMACProjector_H__INCLUDED__
#define LevelMACProjector_H__INCLUDED__

#include "BaseProjector.H"
#include "FluxBox.H"
#include "MappedQuadCFInterp.H"
#include "EllipticBCInterface.H"
#include "PhysBCUtil.H"
#include "Divergence.H"
#include "Gradient.H"


// -----------------------------------------------------------------------------
// Projects a single-level, FC velocity field.
// -----------------------------------------------------------------------------
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
