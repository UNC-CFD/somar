#include "Jacobi.H"
#include "JacobiF_F.H"


// -----------------------------------------------------------------------------
// The Jacobi relaxation method.
// -----------------------------------------------------------------------------
Jacobi::Jacobi (OperatorType* a_opPtr,
                const Real    a_alpha,
                const Real    a_beta)
: m_opPtr(a_opPtr),
  m_alpha(a_alpha),
  m_beta(a_beta)
{
    CH_assert(m_opPtr != NULL);
}

Jacobi::~Jacobi ()
{
    m_opPtr = NULL;
}


void Jacobi::relax (LevelData<FArrayBox>&       a_phi,
                    const LevelData<FArrayBox>& a_rhs)
{
    CH_TIME("Jacobi::relax");

    // Sanity checks
    CH_assert(isDefined());
    CH_assert(a_phi.isDefined());
    CH_assert(a_rhs.isDefined());
    CH_assert(a_phi.ghostVect() >= m_activeDirs);
    CH_assert(a_phi.nComp() == a_rhs.nComp());
    CH_assert(a_phi.getBoxes().compatible(a_rhs.getBoxes()));
    CH_assert(a_phi.getBoxes().compatible(m_FCJgup->getBoxes()));

    LevelData<FArrayBox> resid;
    m_opPtr->create(resid, a_rhs);

    // Get the residual
    m_opPtr->residual(resid, a_phi, a_rhs, true);

    // Do the Jacobi relaxation
    DisjointBoxLayout grids = a_phi.getBoxes();
    DataIterator dit = a_phi.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& phiFAB = a_phi[dit];
        const FArrayBox& resFAB = resid[dit];
        const FArrayBox& lapDiagFAB = (*m_lapDiag)[dit];
        const Box& region = grids[dit];

        FORT_JACOBIITER(CHF_FRA(phiFAB),
                        CHF_CONST_FRA(resFAB),
                        CHF_CONST_FRA1(lapDiagFAB,0),
                        CHF_BOX(region),
                        CHF_CONST_REAL(m_alpha),
                        CHF_CONST_REAL(m_beta));
    }
}
