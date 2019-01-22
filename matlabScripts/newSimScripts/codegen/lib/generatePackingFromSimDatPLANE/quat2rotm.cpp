//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: quat2rotm.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include <math.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "quat2rotm.h"

// Function Definitions

//
// Arguments    : const double q[4]
//                double R[9]
// Return Type  : void
//
void quat2rotm(const double q[4], double R[9])
{
  int k;
  double y;
  double normRowMatrix[4];
  double tempR[9];
  int b_k;
  for (k = 0; k < 4; k++) {
    normRowMatrix[k] = q[k] * q[k];
  }

  y = normRowMatrix[0];
  for (k = 0; k < 3; k++) {
    y += normRowMatrix[k + 1];
  }

  y = 1.0 / sqrt(y);
  for (k = 0; k < 4; k++) {
    normRowMatrix[k] = q[k] * y;
  }

  tempR[0] = 1.0 - 2.0 * (normRowMatrix[2] * normRowMatrix[2] + normRowMatrix[3]
    * normRowMatrix[3]);
  tempR[1] = 2.0 * (normRowMatrix[1] * normRowMatrix[2] - normRowMatrix[0] *
                    normRowMatrix[3]);
  tempR[2] = 2.0 * (normRowMatrix[1] * normRowMatrix[3] + normRowMatrix[0] *
                    normRowMatrix[2]);
  tempR[3] = 2.0 * (normRowMatrix[1] * normRowMatrix[2] + normRowMatrix[0] *
                    normRowMatrix[3]);
  tempR[4] = 1.0 - 2.0 * (normRowMatrix[1] * normRowMatrix[1] + normRowMatrix[3]
    * normRowMatrix[3]);
  tempR[5] = 2.0 * (normRowMatrix[2] * normRowMatrix[3] - normRowMatrix[0] *
                    normRowMatrix[1]);
  tempR[6] = 2.0 * (normRowMatrix[1] * normRowMatrix[3] - normRowMatrix[0] *
                    normRowMatrix[2]);
  tempR[7] = 2.0 * (normRowMatrix[2] * normRowMatrix[3] + normRowMatrix[0] *
                    normRowMatrix[1]);
  tempR[8] = 1.0 - 2.0 * (normRowMatrix[1] * normRowMatrix[1] + normRowMatrix[2]
    * normRowMatrix[2]);
  memcpy(&R[0], &tempR[0], 9U * sizeof(double));
  for (k = 0; k < 3; k++) {
    for (b_k = 0; b_k < 3; b_k++) {
      R[k + 3 * b_k] = tempR[b_k + 3 * k];
    }
  }
}

//
// File trailer for quat2rotm.cpp
//
// [EOF]
//
