//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: pdist.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include <math.h>
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "pdist.h"

// Function Definitions

//
// Arguments    : const double Xin[6]
// Return Type  : double
//
double pdist(const double Xin[6])
{
  double Y;
  int kk;
  int ii;
  int i2;
  boolean_T logIndX[2];
  double X[6];
  boolean_T nanflag;
  int jj;
  double tempSum;
  boolean_T exitg1;
  for (kk = 0; kk < 2; kk++) {
    for (i2 = 0; i2 < 3; i2++) {
      X[i2 + 3 * kk] = Xin[kk + (i2 << 1)];
    }

    logIndX[kk] = true;
  }

#pragma omp parallel for \
 num_threads(omp_get_max_threads()) \
 private(nanflag,jj,exitg1)

  for (ii = 1; ii < 3; ii++) {
    nanflag = false;
    jj = 1;
    exitg1 = false;
    while ((!exitg1) && (jj < 4)) {
      if (rtIsNaN(X[(jj + 3 * (ii - 1)) - 1])) {
        nanflag = true;
        exitg1 = true;
      } else {
        jj++;
      }
    }

    if (nanflag) {
      logIndX[ii - 1] = false;
    }
  }

#pragma omp parallel for \
 num_threads(omp_get_max_threads()) \
 private(tempSum,jj)

  for (kk = 0; kk < 1; kk++) {
    tempSum = 0.0;
    if (logIndX[0] && logIndX[1]) {
      for (jj = 0; jj < 3; jj++) {
        tempSum += (X[3 + jj] - X[jj]) * (X[3 + jj] - X[jj]);
      }

      Y = sqrt(tempSum);
    }
  }

  return Y;
}

//
// File trailer for pdist.cpp
//
// [EOF]
//
