//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: projPointOnPlane.h
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//
#ifndef PROJPOINTONPLANE_H
#define PROJPOINTONPLANE_H

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "omp.h"
#include "generatePackingFromSimDatPLANE_types.h"

// Function Declarations
extern void b_projPointOnPlane(double point[9], const double plane[9]);
extern void projPointOnPlane(double point[12], const double plane[9]);

#endif

//
// File trailer for projPointOnPlane.h
//
// [EOF]
//
