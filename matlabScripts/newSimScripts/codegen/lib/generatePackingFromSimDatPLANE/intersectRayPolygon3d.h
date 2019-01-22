//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: intersectRayPolygon3d.h
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//
#ifndef INTERSECTRAYPOLYGON3D_H
#define INTERSECTRAYPOLYGON3D_H

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "omp.h"
#include "generatePackingFromSimDatPLANE_types.h"

// Function Declarations
extern void intersectRayPolygon3d(const double ray[18], const double poly[12],
  double inter[9], boolean_T inside[3]);

#endif

//
// File trailer for intersectRayPolygon3d.h
//
// [EOF]
//
