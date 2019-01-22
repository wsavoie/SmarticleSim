//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: dot.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "dot.h"

// Function Definitions

//
// Arguments    : const double a[9]
//                const double b[9]
//                double c[3]
// Return Type  : void
//
void b_dot(const double a[9], const double b[9], double c[3])
{
  int ic;
  int i1;
  int j;
  double b_c;
  int ix;
  int iy;
  int k;
  ic = -1;
  i1 = -1;
  for (j = 0; j < 3; j++) {
    ic++;
    i1++;
    b_c = 0.0;
    ix = i1;
    iy = i1;
    for (k = 0; k < 3; k++) {
      b_c += a[ix] * b[iy];
      ix += 3;
      iy += 3;
    }

    c[ic] = b_c;
  }
}

//
// Arguments    : const double a[12]
//                const double b[12]
//                double c[4]
// Return Type  : void
//
void dot(const double a[12], const double b[12], double c[4])
{
  int ic;
  int i1;
  int j;
  double b_c;
  int ix;
  int iy;
  int k;
  ic = -1;
  i1 = -1;
  for (j = 0; j < 4; j++) {
    ic++;
    i1++;
    b_c = 0.0;
    ix = i1;
    iy = i1;
    for (k = 0; k < 3; k++) {
      b_c += a[ix] * b[iy];
      ix += 4;
      iy += 4;
    }

    c[ic] = b_c;
  }
}

//
// File trailer for dot.cpp
//
// [EOF]
//
