#include "Curl.H"
#include "DivCurlGradF_F.H"
#include "LevelGeometry.H"


// -----------------------------------------------------------------------------
// Compute the curl of a contravariant vector field.
// -----------------------------------------------------------------------------
void Curl::simpleCurlCC (FArrayBox&           a_curl,
                         const FArrayBox&     a_contraVector,
                         const DataIndex&     a_di,
                         const LevelGeometry& a_levGeo,
                         const bool           a_convertToVector)
{
    CH_TIME("Curl::simpleCurlCC");

    // Sanity checks
    const Box& valid = a_levGeo.getBoxes()[a_di];
    CH_assert(a_contraVector.box().contains(grow(valid,1)));
    CH_assert(a_contraVector.nComp() == SpaceDim);
    CH_assert(a_curl.nComp() == CURL_NCOMP);

    // Convert a_contraVector to covariant.
    FArrayBox coVel(a_contraVector.box(), a_contraVector.nComp());
    a_levGeo.makeCovariant(coVel, a_contraVector, a_contraVector.box(), a_di);

    // Compute d[coVel]
    const RealVect& dx = a_levGeo.getDx();
    FORT_DONEFORMCC(CHF_FRA(a_curl),
                    CHF_CONST_FRA(coVel),
                    CHF_BOX(valid),
                    CHF_CONST_REALVECT(dx));

    // a_curl is a two-form. Should we turn this into a vector?
    if (a_convertToVector) {
        FArrayBox twoForm(a_curl.box(), a_curl.nComp());
        twoForm.copy(a_curl);

        a_levGeo.divByJ(twoForm, a_di);
        a_levGeo.makeContravariant(a_curl, twoForm, valid, a_di);
    }
}


// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void Curl::levelCurlCC (LevelData<FArrayBox>&       a_curl,
                        const LevelData<FArrayBox>& a_contraVel,
                        const LevelGeometry&        a_levGeo,
                        const bool                  a_convertToVector)
{
    CH_TIME("Curl::levelCurlCC");

    // Sanity checks
    CH_assert(a_curl.nComp() == CURL_NCOMP);
    CH_assert(a_contraVel.nComp() == SpaceDim);

    // Calculate curl on each grid
    DataIterator dit = a_contraVel.dataIterator();
    for (dit.reset(); dit.ok(); ++dit) {
        FArrayBox& curlFAB = a_curl[dit];
        const FArrayBox& contraVelFAB = a_contraVel[dit];

        Curl::simpleCurlCC(curlFAB,
                           contraVelFAB,
                           dit(),
                           a_levGeo);
    }
}
