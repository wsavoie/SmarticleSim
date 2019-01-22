//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: vectorNorm3d.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include <math.h>
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "vectorNorm3d.h"

// Function Definitions

//
// VECTORNORM3D Norm of a 3D vector or of set of 3D vectors
//
//    N = vectorNorm3d(V);
//    Returns the norm of vector V.
//
//    When V is a N-by-3 array, compute norm for each vector of the array.
//    Vector are given as rows. Result is then a N-by-1 array.
//
//    NOTE: compute only euclidean norm.
//
//    See Also
//    vectors3d, normalizeVector3d, vectorAngle3d, hypot3
// Arguments    : const double v[3]
// Return Type  : double
//
double vectorNorm3d(const double v[3])
{
  int k;
  double y;
  double x[3];

  //    ---------
  //    author : David Legland
  //    INRA - TPV URPOI - BIA IMASTE
  //    created the 21/02/2005.
  //    HISTORY
  //    19/06/2009 rename as vectorNorm3d
  for (k = 0; k < 3; k++) {
    x[k] = v[k] * v[k];
  }

  y = x[0];
  for (k = 0; k < 2; k++) {
    y += x[k + 1];
  }

  return sqrt(y);
}

//
// File trailer for vectorNorm3d.cpp
//
// [EOF]
//
