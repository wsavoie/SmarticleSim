//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: bsxfun.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "bsxfun.h"

// Function Definitions

//
// Arguments    : const double a[3]
//                const double b[9]
//                double c[9]
// Return Type  : void
//
void b_bsxfun(const double a[3], const double b[9], double c[9])
{
  int k;
  int b_k;
  for (k = 0; k < 3; k++) {
    for (b_k = 0; b_k < 3; b_k++) {
      c[b_k + 3 * k] = a[k] * b[b_k + 3 * k];
    }
  }
}

//
// Arguments    : const double a[3]
//                const double b[9]
//                double c[9]
// Return Type  : void
//
void bsxfun(const double a[3], const double b[9], double c[9])
{
  int k;
  int b_k;
  for (k = 0; k < 3; k++) {
    for (b_k = 0; b_k < 3; b_k++) {
      c[b_k + 3 * k] = a[k] - b[b_k + 3 * k];
    }
  }
}

//
// Arguments    : const double a[9]
//                const double b[9]
//                double c[9]
// Return Type  : void
//
void c_bsxfun(const double a[9], const double b[9], double c[9])
{
  int k;
  int b_k;
  for (k = 0; k < 3; k++) {
    for (b_k = 0; b_k < 3; b_k++) {
      c[b_k + 3 * k] = a[b_k + 3 * k] * b[b_k + 3 * k];
    }
  }
}

//
// Arguments    : const double a[9]
//                const double b[9]
//                double c[9]
// Return Type  : void
//
void d_bsxfun(const double a[9], const double b[9], double c[9])
{
  int k;
  int b_k;
  for (k = 0; k < 3; k++) {
    for (b_k = 0; b_k < 3; b_k++) {
      c[b_k + 3 * k] = a[b_k + 3 * k] + b[b_k + 3 * k];
    }
  }
}

//
// Arguments    : const double a[12]
//                const double b[3]
//                double c[12]
// Return Type  : void
//
void e_bsxfun(const double a[12], const double b[3], double c[12])
{
  int k;
  int b_k;
  for (k = 0; k < 3; k++) {
    for (b_k = 0; b_k < 4; b_k++) {
      c[b_k + (k << 2)] = a[b_k + (k << 2)] - b[k];
    }
  }
}

//
// Arguments    : const double a[9]
//                const double b[3]
//                double c[9]
// Return Type  : void
//
void f_bsxfun(const double a[9], const double b[3], double c[9])
{
  int k;
  int b_k;
  for (k = 0; k < 3; k++) {
    for (b_k = 0; b_k < 3; b_k++) {
      c[b_k + 3 * k] = a[b_k + 3 * k] - b[k];
    }
  }
}

//
// Arguments    : const double a_data[]
//                const int a_size[2]
//                const double b_data[]
//                const int b_size[2]
//                double c_data[]
//                int c_size[2]
// Return Type  : void
//
void g_bsxfun(const double a_data[], const int a_size[2], const double b_data[],
              const int b_size[2], double c_data[], int c_size[2])
{
  int acoef;
  int bcoef;
  int sck;
  int b_b_size;
  int k;
  int b_k;
  acoef = b_size[0];
  bcoef = a_size[0];
  if (acoef < bcoef) {
    bcoef = acoef;
  }

  if (b_size[0] == 1) {
    sck = a_size[0];
  } else if (a_size[0] == 1) {
    sck = b_size[0];
  } else if (a_size[0] == b_size[0]) {
    sck = a_size[0];
  } else {
    sck = bcoef;
  }

  acoef = b_size[0];
  bcoef = a_size[0];
  if (acoef < bcoef) {
    bcoef = acoef;
  }

  if (b_size[0] == 1) {
    c_size[0] = (signed char)a_size[0];
  } else if (a_size[0] == 1) {
    c_size[0] = (signed char)b_size[0];
  } else if (a_size[0] == b_size[0]) {
    c_size[0] = (signed char)a_size[0];
  } else {
    c_size[0] = (signed char)bcoef;
  }

  c_size[1] = 3;
  acoef = b_size[0];
  bcoef = a_size[0];
  if (acoef < bcoef) {
    bcoef = acoef;
  }

  if (b_size[0] == 1) {
    b_b_size = a_size[0];
  } else if (a_size[0] == 1) {
    b_b_size = b_size[0];
  } else if (a_size[0] == b_size[0]) {
    b_b_size = a_size[0];
  } else {
    b_b_size = bcoef;
  }

  if ((signed char)b_b_size != 0) {
    acoef = (a_size[0] != 1);
    bcoef = (b_size[0] != 1);
    for (k = 0; k < 3; k++) {
      for (b_k = 0; b_k < (signed char)sck; b_k++) {
        c_data[b_k + (signed char)sck * k] = a_data[acoef * b_k + a_size[0] * k]
          - b_data[bcoef * b_k + b_size[0] * k];
      }
    }
  }
}

//
// File trailer for bsxfun.cpp
//
// [EOF]
//
