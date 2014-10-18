/*******************************************************************************
 *    FILE: StressMetric.H
 *    DESC: Creates J*T^{i j} for viscous solves.
 *    DATE: Wed 6 Aug 2014 10:58 AM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __StressMetric_H__INCLUDED__
#define __StressMetric_H__INCLUDED__

// To use a custom stress tensor, we override Jg^{i,j} in the viscous solver's
// elliptic operators. By defining USE_STRESSMETRIC, all of the needed code
// will be produced and you can rewrite the fill_Jgup function to include any
// subgrid model you like. If your subgrid model needs the state variables, the
// define function will need to be altered as well.
//
// If USE_STRESSMETRIC is undefined, then the StressMetric class will not exist
// anywhere in the (preprocessed) code. This will save some memory since the
// viscous solver's elliptic operators will use the Jg^{i,j} that has been
// cached by LevelGeometry.
//
// #define USE_STRESSMETRIC
#ifdef USE_STRESSMETRIC

#include "GeoSourceInterface.H"


class StressMetric: public FillJgupInterface
{
public:
    // Constructor
    StressMetric ();

    // Full constructor
    StressMetric (const GeoSourceInterface* a_geoSourcePtr);

    // Destructor
    virtual ~StressMetric ();

    // Full constructor
    virtual void define (const GeoSourceInterface* a_geoSourcePtr);

    // Is this defined?
    virtual inline bool isDefined () const;

    // Set dest = (1+turbVisc/dynVisc)*J*gup^{ij}
    // Remember: The solvers will multiply by dynVisc, that's why it's factored out.
    virtual void fill_Jgup (FArrayBox&       a_dest,
                            const int        a_destComp,
                            const int        a_mu,
                            const int        a_nu,
                            const RealVect&  a_dXi,
                            const Real       a_scale = 1.0) const;

protected:
    const GeoSourceInterface* m_geoSourcePtr;
    Real                      m_dynVisc;
};



// -----------------------------------------------------------------------------
// Is this defined?
// -----------------------------------------------------------------------------
bool StressMetric::isDefined () const
{
    return (m_geoSourcePtr != NULL);
}


#endif //USE_STRESSMETRIC
#endif //!__StressMetric_H__INCLUDED__
