//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: normalizeVector3d.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include <math.h>
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "normalizeVector3d.h"

// Function Definitions

//
// NORMALIZEVECTOR3D Normalize a 3D vector to have norm equal to 1
//
//    V2 = normalizeVector3d(V);
//    Returns the normalization of vector V, such that ||V|| = 1. Vector V is
//    given as a row vector.
//
//    If V is a N-by-3 array, normalization is performed for each row of the
//    input array.
//
//    See also:
//    vectors3d, vectorNorm3d
// Arguments    : const double v[3]
//                double vn[3]
// Return Type  : void
//
void normalizeVector3d(const double v[3], double vn[3])
{
  int k;
  double y;
  double z1[3];

  //    ---------
  //    author : David Legland
  //    INRA - TPV URPOI - BIA IMASTE
  //    created the 29/11/2004.
  //
  //  HISTORY
  //  2005-11-30 correct a bug
  //  2009-06-19 rename as normalizeVector3d
  //  2010-11-16 use bsxfun (Thanks to Sven Holcombe)
  for (k = 0; k < 3; k++) {
    z1[k] = v[k] * v[k];
  }

  y = z1[0];
  for (k = 0; k < 2; k++) {
    y += z1[k + 1];
  }

  y = sqrt(y);
  for (k = 0; k < 3; k++) {
    vn[k] = v[k] / y;
  }
}

//
// File trailer for normalizeVector3d.cpp
//
// [EOF]
//
