//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: splitPolygons.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include <string.h>
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "splitPolygons.h"

// Function Definitions

//
// SPLITPOLYGONS Convert a NaN separated polygon list to a cell array of polygons
//
//    POLYGONS = splitPolygons(POLYGON);
//    POLYGON is a N-by-2 array of points, possibly with pairs of NaN values.
//    The functions separates each component separated by NaN values, and
//    returns a cell array of polygons.
//
//    See also:
//    polygons2d
// Arguments    : const double polygon[8]
//                cell_wrap_3 polygons_data[]
//                int polygons_size[1]
// Return Type  : void
//
void splitPolygons(const double polygon[8], cell_wrap_3 polygons_data[], int
                   polygons_size[1])
{
  int k;
  double y;
  boolean_T b[8];
  cell_wrap_3 r0;
  int idx;
  signed char b_y[4];
  boolean_T exitg1;
  int loop_ub;
  int ii_data[4];
  int i_data[4];
  int inds_data[6];
  int i3;
  int b_loop_ub;
  int i4;

  //  ------
  //  Author: David Legland
  //  e-mail: david.legland@inra.fr
  //  Created: 2007-10-12,    using Matlab 7.4.0.287 (R2007a)
  //  Copyright 2007 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
  for (k = 0; k < 8; k++) {
    b[k] = rtIsNaN(polygon[k]);
  }

  y = b[0];
  for (k = 0; k < 7; k++) {
    y += (double)b[k + 1];
  }

  if (y == 0.0) {
    //  single polygon -> no break
    r0.f1.size[0] = 4;
    r0.f1.size[1] = 2;
    memcpy(&r0.f1.data[0], &polygon[0], sizeof(double) << 3);
    polygons_size[0] = 1;
    polygons_data[0] = r0;
  } else {
    //  find indices of NaN couples
    for (k = 0; k < 8; k++) {
      b[k] = rtIsNaN(polygon[k]);
    }

    for (k = 0; k < 4; k++) {
      b_y[k] = (signed char)(b[k] + b[k + 4]);
    }

    idx = 0;
    k = 1;
    exitg1 = false;
    while ((!exitg1) && (k < 5)) {
      if (b_y[k - 1] > 0) {
        idx++;
        ii_data[idx - 1] = k;
        if (idx >= 4) {
          exitg1 = true;
        } else {
          k++;
        }
      } else {
        k++;
      }
    }

    if (1 > idx) {
      loop_ub = 0;
      k = 0;
    } else {
      loop_ub = idx;
      k = idx;
    }

    if (0 <= loop_ub - 1) {
      memcpy(&i_data[0], &ii_data[0], (unsigned int)(loop_ub * (int)sizeof(int)));
    }

    //  number of polygons
    //  iterate over NaN-separated regions to create new polygon
    inds_data[0] = 0;
    if (0 <= k - 1) {
      memcpy(&inds_data[1], &i_data[0], (unsigned int)(k * (int)sizeof(int)));
    }

    inds_data[1 + k] = 5;
    polygons_size[0] = (signed char)(loop_ub + 1);
    for (idx = 0; idx <= loop_ub; idx++) {
      if ((double)inds_data[idx] + 1.0 > (double)inds_data[idx + 1] - 1.0) {
        k = 1;
        i3 = 0;
      } else {
        k = inds_data[idx] + 1;
        i3 = (int)((double)inds_data[idx + 1] - 1.0);
      }

      polygons_data[idx].f1.size[0] = (i3 - k) + 1;
      polygons_data[idx].f1.size[1] = 2;
      b_loop_ub = (i3 - k) + 1;
      for (i3 = 0; i3 < 2; i3++) {
        for (i4 = 0; i4 < b_loop_ub; i4++) {
          polygons_data[idx].f1.data[i4 + polygons_data[idx].f1.size[0] * i3] =
            polygon[((k + i4) + (i3 << 2)) - 1];
        }
      }
    }
  }
}

//
// File trailer for splitPolygons.cpp
//
// [EOF]
//
