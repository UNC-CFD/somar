#include "HermiteInterp.H"
#include "HermiteInterpF_F.H"


// -----------------------------------------------------------------------------
// The nodes will be labeled as A, B, C, and D in this order:   C---D
//                                                              |   |
//                                                              A---B
// f_{ABCD}(u,v) = f_{AB}(u)*h1(v) + f_{CD}(u)*h2(v)
//               + f_{AC}(v)*h1(u) + f_{BD}(v)*h2(u)
//               - f_{A}*h1(u)*h1(v) - f_{B}*h2(u)*h1(v)
//               - f_{C}*h1(u)*h2(v) - f_{D}*h2(u)*h2(v)
// where f_{AB}(u) = f_{A}*h1(u) + f_{B}*h2(u) + f'_{A}*h3(u) + f'_{B}*h4(u),
// and so on...
//
// This is the same as
// f_{ABCD}(u,v) = f_{A}*h1(u)*h1(v)
//               + f_{B}*h2(u)*h1(v)
//               + f_{C}*h1(u)*h2(v)
//               + f_{D}*h2(u)*h2(v)
//               + fu'_{A}*h3(u)*h1(v)
//               + fu'_{B}*h4(u)*h1(v)
//               + fu'_{C}*h3(u)*h2(v)
//               + fu'_{D}*h4(u)*h2(v)
//               + fv'_{A}*h1(u)*h3(v)
//               + fv'_{B}*h2(u)*h3(v)
//               + fv'_{C}*h1(u)*h4(v)
//               + fv'_{D}*h2(u)*h4(v)
// -----------------------------------------------------------------------------
void HermiteInterp2D (FArrayBox&          a_fInterp,
                      const FArrayBox&    a_xInterp,
                      const FArrayBox&    a_yInterp,
                      const Box&          a_interpBox,
                      const int           a_xdir,
                      const int           a_ydir,
                      const Vector<Real>& a_x,
                      const Vector<Real>& a_y,
                      const FArrayBox&    a_f,
                      const FArrayBox&    a_dfdx,
                      const FArrayBox&    a_dfdy)
{
#ifndef NDEBUG
    {
        // Check centerings.
        CH_assert(a_fInterp.box().type() == a_interpBox.type());
        CH_assert(a_xInterp.box().type() == a_interpBox.type());
        CH_assert(a_yInterp.box().type() == a_interpBox.type());
        CH_assert(a_dfdx   .box().type() == a_f.box().type());
        CH_assert(a_dfdy   .box().type() == a_f.box().type());

        // Check FAB regions.
        CH_assert(a_fInterp.box().contains(a_interpBox));
        CH_assert(a_xInterp.box().contains(a_interpBox));
        CH_assert(a_yInterp.box().contains(a_interpBox));
        CH_assert(a_dfdx.box() == a_f.box());
        CH_assert(a_dfdy.box() == a_f.box());
        CH_assert(a_f.box().size(a_xdir) == a_x.size());
        CH_assert(a_f.box().size(a_ydir) == a_y.size());

        // Check number of comps.
        CH_assert(a_xInterp.nComp() == 1);
        CH_assert(a_yInterp.nComp() == 1);
        CH_assert(a_f      .nComp() == a_fInterp.nComp());
        CH_assert(a_dfdx   .nComp() == a_fInterp.nComp());
        CH_assert(a_dfdy   .nComp() == a_fInterp.nComp());

        // Check dirs
        CH_assert(0 <= a_xdir);
        CH_assert(a_xdir < SpaceDim);
        CH_assert(0 <= a_ydir);
        CH_assert(a_ydir < SpaceDim);
        CH_assert(a_xdir != a_ydir);
    }
#endif

    // The fortran function can only handle right-handed
    // permutations of the x, y, and z directions.
    if ((a_xdir == 0 && a_ydir == 1) ||
        (a_xdir == 1 && a_ydir == 2) ||
        (a_xdir == 2 && a_ydir == 0)) {

        // The fortran function requires the vector indices to
        // coincide with the constraint indices.
        IntVect shift = IntVect::Zero;
        shift[a_xdir] = a_f.box().smallEnd(a_xdir);
        shift[a_ydir] = a_f.box().smallEnd(a_ydir);

        // Interpolate
        FORT_HERMITEINTERP2DF (
            CHF_FRA(a_fInterp),
            CHF_CONST_FRA1(a_xInterp,0),
            CHF_CONST_FRA1(a_yInterp,0),
            CHF_BOX(a_interpBox),
            CHF_CONST_INT(a_xdir),
            CHF_CONST_INT(a_ydir),
            CHF_CONST_VR(a_x),
            CHF_CONST_VR(a_y),
            CHF_CONST_FRA_SHIFT(a_f,shift),
            CHF_CONST_FRA_SHIFT(a_dfdx,shift),
            CHF_CONST_FRA_SHIFT(a_dfdy,shift));
    } else {

        // The fortran function requires the vector indices to
        // coincide with the constraint indices.
        IntVect shift = IntVect::Zero;
        shift[a_ydir] = a_f.box().smallEnd(a_xdir);
        shift[a_xdir] = a_f.box().smallEnd(a_ydir);

        // Interpolate
        FORT_HERMITEINTERP2DF (
            CHF_FRA(a_fInterp),
            CHF_CONST_FRA1(a_yInterp,0),
            CHF_CONST_FRA1(a_xInterp,0),
            CHF_BOX(a_interpBox),
            CHF_CONST_INT(a_ydir),
            CHF_CONST_INT(a_xdir),
            CHF_CONST_VR(a_y),
            CHF_CONST_VR(a_x),
            CHF_CONST_FRA_SHIFT(a_f,shift),
            CHF_CONST_FRA_SHIFT(a_dfdy,shift),
            CHF_CONST_FRA_SHIFT(a_dfdx,shift));
    }
}
