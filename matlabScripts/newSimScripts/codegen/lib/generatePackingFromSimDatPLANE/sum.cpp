//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: sum.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "sum.h"

// Function Definitions

//
// Arguments    : const double x_data[]
//                const int x_size[2]
//                double y_data[]
//                int y_size[1]
// Return Type  : void
//
void b_sum(const double x_data[], const int x_size[2], double y_data[], int
           y_size[1])
{
  int vstride;
  int j;
  int k;
  int xoffset;
  if (x_size[0] == 0) {
    y_size[0] = 0;
  } else {
    vstride = x_size[0];
    y_size[0] = (signed char)x_size[0];
    for (j = 0; j < vstride; j++) {
      y_data[j] = x_data[j];
    }

    for (k = 0; k < 2; k++) {
      xoffset = (k + 1) * vstride;
      for (j = 0; j < vstride; j++) {
        y_data[j] += x_data[xoffset + j];
      }
    }
  }
}

//
// Arguments    : const double x[9]
//                double y[3]
// Return Type  : void
//
void sum(const double x[9], double y[3])
{
  int j;
  int k;
  int xoffset;
  for (j = 0; j < 3; j++) {
    y[j] = x[j];
  }

  for (k = 0; k < 2; k++) {
    xoffset = (k + 1) * 3;
    for (j = 0; j < 3; j++) {
      y[j] += x[xoffset + j];
    }
  }
}

//
// File trailer for sum.cpp
//
// [EOF]
//
