//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: bsxfun.h
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//
#ifndef BSXFUN_H
#define BSXFUN_H

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "omp.h"
#include "generatePackingFromSimDatPLANE_types.h"

// Function Declarations
extern void b_bsxfun(const double a[3], const double b[9], double c[9]);
extern void bsxfun(const double a[3], const double b[9], double c[9]);
extern void c_bsxfun(const double a[9], const double b[9], double c[9]);
extern void d_bsxfun(const double a[9], const double b[9], double c[9]);
extern void e_bsxfun(const double a[12], const double b[3], double c[12]);
extern void f_bsxfun(const double a[9], const double b[3], double c[9]);
extern void g_bsxfun(const double a_data[], const int a_size[2], const double
                     b_data[], const int b_size[2], double c_data[], int c_size
                     [2]);

#endif

//
// File trailer for bsxfun.h
//
// [EOF]
//
