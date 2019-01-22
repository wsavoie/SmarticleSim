//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: sum.h
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//
#ifndef SUM_H
#define SUM_H

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "omp.h"
#include "generatePackingFromSimDatPLANE_types.h"

// Function Declarations
extern void b_sum(const double x_data[], const int x_size[2], double y_data[],
                  int y_size[1]);
extern void sum(const double x[9], double y[3]);

#endif

//
// File trailer for sum.h
//
// [EOF]
//
