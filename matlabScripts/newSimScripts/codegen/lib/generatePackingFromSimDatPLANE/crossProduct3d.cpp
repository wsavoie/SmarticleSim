//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: crossProduct3d.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "crossProduct3d.h"

// Function Definitions

//
// CROSSPRODUCT3D Vector cross product faster than inbuilt MATLAB cross.
//
//    C = crossProduct3d(A, B)
//    returns the cross product of the two 3D vectors A and B, that is:
//        C = A x B
//    A and B must be N-by-3 element vectors. If either A or B is a 1-by-3
//    row vector, the result C will have the size of the other input and will
//    be the  concatenation of each row's cross product.
//
//    Example
//      v1 = [2 0 0];
//      v2 = [0 3 0];
//      crossProduct3d(v1, v2)
//      ans =
//          0   0   6
//
//
//    Class support for inputs A,B:
//       float: double, single
//
//    See also DOT.
// Arguments    : const double a[3]
//                const double b[3]
//                double c[3]
// Return Type  : void
//
void crossProduct3d(const double a[3], const double b[3], double c[3])
{
  int k;
  static const signed char iv1[3] = { 1, 2, 0 };

  static const signed char iv2[3] = { 2, 0, 1 };

  //    Sven Holcombe
  //  HISTORY
  //  2017-11-24 rename from vectorCross3d to crossProduct3d
  //  size of inputs
  //  Initialise c to the size of a or b, whichever has more dimensions. If
  //  they have the same dimensions, initialise to the larger of the two
  for (k = 0; k < 3; k++) {
    c[k] = a[iv1[k]] * b[iv2[k]] - b[iv1[k]] * a[iv2[k]];
  }
}

//
// File trailer for crossProduct3d.cpp
//
// [EOF]
//
