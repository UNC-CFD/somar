/*******************************************************************************
 *    FILE: AlteredMetric.H
 *    DESC: This produces the metric used by the implicit BV projector.
 *    DATE: Fri 11 Oct 2013 08:37:34 AM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __AlteredMetric_H__INCLUDED__
#define __AlteredMetric_H__INCLUDED__

#include "GeoSourceInterface.H"
class PhysBCUtil;


class AlteredMetric: public FillJgupInterface
{
public:
    // Constructor
    AlteredMetric ();

    // Full constructor
    AlteredMetric (const GeoSourceInterface* a_geoSourcePtr,
                   const PhysBCUtil*         a_physBCPtr,
                   const Real                a_dtTheta);

    // Destructor
    virtual ~AlteredMetric ();

    // Full constructor
    virtual void define (const GeoSourceInterface* a_geoSourcePtr,
                         const PhysBCUtil*         a_physBCPtr,
                         const Real                a_dtTheta);

    // Set dest = J*(gup^{ij} - (omega^2 / (1 + omega^2)) * dXi^i/dz * dXi^j/dz)
    // where omega = (dt*N*theta)^2
    virtual void fill_Jgup (FArrayBox&       a_dest,
                            const int        a_destComp,
                            const int        a_mu,
                            const int        a_nu,
                            const RealVect&  a_dXi,
                            const Real       a_scale = 1.0) const;

protected:
    const GeoSourceInterface* m_geoSourcePtr;
    const PhysBCUtil*         m_physBCPtr;
    Real                      m_theta;
};



#endif //!__AlteredMetric_H__INCLUDED__
