//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: inpolygon.h
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//
#ifndef INPOLYGON_H
#define INPOLYGON_H

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "omp.h"
#include "generatePackingFromSimDatPLANE_types.h"

// Function Declarations
extern void b_inpolygon(const double x[3], const double y[3], const double xv[4],
  const double yv[4], boolean_T in[3]);
extern void inpolygon(const double x[3], const double y[3], const double
                      xv_data[], const int xv_size[1], const double yv_data[],
                      boolean_T in[3]);

#endif

//
// File trailer for inpolygon.h
//
// [EOF]
//
