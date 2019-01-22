//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: splitPolygons.h
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//
#ifndef SPLITPOLYGONS_H
#define SPLITPOLYGONS_H

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "omp.h"
#include "generatePackingFromSimDatPLANE_types.h"

// Function Declarations
extern void splitPolygons(const double polygon[8], cell_wrap_3 polygons_data[],
  int polygons_size[1]);

#endif

//
// File trailer for splitPolygons.h
//
// [EOF]
//
