/*******************************************************************************
 *    FILE: LepticOperator.H
 *    DESC: To use an operator in a LepticSolver, it must inherit from LinearOp
 *          and this pure abstract base class.
 *    DATE: Thu 11 Jul 2013 04:02:29 PM EDT
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __LepticOperator_H__INCLUDED__
#define __LepticOperator_H__INCLUDED__

#include "RealVect.H"
#include "LevelData.H"
#include "EllipticBCInterface.H"


// To use an operator in a LepticSolver, it must inherit from LinearOp
// and this pure abstract base class.
class LepticOperator {
public:
    virtual ~LepticOperator () {;}

    virtual const RealVect& getDx () const = 0;
    virtual const RealVect& getDxCrse () const = 0; // AMR, not MG
    virtual void setDxCrse (const RealVect& a_dxCrse) = 0; // AMR, not MG
    virtual const RefCountedPtr<LevelData<FluxBox> > getFCJgup () const = 0;
    virtual const RefCountedPtr<LevelData<FArrayBox> > getCCJinv () const = 0;
    virtual BCMethodHolder getBCs () const = 0;
    virtual int smallestGridSize () const = 0;
};


#endif //!__LepticOperator_H__INCLUDED__
