/*******************************************************************************
 *    FILE: Jacobi.H
 *    DESC: Jacobi relaxation method.
 *    DATE: Tue 26 Nov 2013 02:48:44 PM EST
 *    MAIL: santilli@physics.unc.edu
 ******************************************************************************/

#ifndef __Jacobi_H__INCLUDED__
#define __Jacobi_H__INCLUDED__

#include "RelaxationMethod.H"
#include "LinearSolver.H"


// The Jacobi relaxation method.
class Jacobi: public RelaxationMethod
{
    typedef LinearOp<LevelData<FArrayBox> > OperatorType;
public:
    virtual ~Jacobi ();

    Jacobi (OperatorType* a_opPtr,
            const Real    a_alpha,
            const Real    a_beta);

    virtual void relax (LevelData<FArrayBox>&       a_phi,
                        const LevelData<FArrayBox>& a_rhs);

protected:
    OperatorType* m_opPtr;
    Real          m_alpha;
    Real          m_beta;
};


#endif //!__Jacobi_H__INCLUDED__
