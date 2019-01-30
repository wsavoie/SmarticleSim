//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: polygonArea.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include <string.h>
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "polygonArea.h"

// Function Definitions

//
// POLYGONAREA Compute the signed area of a polygon
//
//    A = polygonArea(POINTS);
//    Compute area of a polygon defined by POINTS. POINTS is a N-by-2 array
//    of double containing coordinates of vertices.
//
//    Vertices of the polygon are supposed to be oriented Counter-Clockwise
//    (CCW). In this case, the signed area is positive.
//    If vertices are oriented Clockwise (CW), the signed area is negative.
//
//    If polygon is self-crossing, the result is undefined.
//
//    Examples
//      % compute area of a simple shape
//      poly = [10 10;30 10;30 20;10 20];
//      area = polygonArea(poly)
//      area =
//          200
//
//      % compute area of CW polygon
//      area2 = polygonArea(poly(end:-1:1, :))
//      area2 =
//          -200
//
//      % Computes area of a paper hen
//      x = [0 10 20  0 -10 -20 -10 -10  0];
//      y = [0  0 10 10  20  10  10  0 -10];
//      poly = [x' y'];
//      area = polygonArea(poly)
//      area =
//         400
//
//      % Area of unit square with 25% hole
//      pccw = [0 0; 1 0; 1 1; 0 1];
//      pcw = pccw([1 4 3 2], :) * .5 + .25;
//      polygonArea ([pccw; nan(1,2); pcw])
//      ans =
//         0.75
//
//    References
//    algo adapted from P. Bourke web page
//    http://paulbourke.net/geometry/polygonmesh/
//
//    See also:
//    polygons2d, polygonCentroid, polygonSecondAreaMoments, triangleArea
// Arguments    : const cell_wrap_3 poly_data[]
//                const int poly_size[1]
// Return Type  : double
//
double b_polygonArea(const cell_wrap_3 poly_data[], const int poly_size[1])
{
  double area;
  int i;

  //    ---------
  //    author : David Legland
  //    INRA - TPV URPOI - BIA IMASTE
  //    created the 05/05/2004.
  //
  //    HISTORY
  //    25/04/2005: add support for multiple polygons
  //    12/10/2007: update doc
  //  Process special cases
  //  in case of polygon sets, computes the sum of polygon areas
  area = 0.0;
  for (i = 0; i < poly_size[0]; i++) {
    area += polygonArea(poly_data[i].f1.data, poly_data[i].f1.size);
  }

  return area;
}

//
// POLYGONAREA Compute the signed area of a polygon
//
//    A = polygonArea(POINTS);
//    Compute area of a polygon defined by POINTS. POINTS is a N-by-2 array
//    of double containing coordinates of vertices.
//
//    Vertices of the polygon are supposed to be oriented Counter-Clockwise
//    (CCW). In this case, the signed area is positive.
//    If vertices are oriented Clockwise (CW), the signed area is negative.
//
//    If polygon is self-crossing, the result is undefined.
//
//    Examples
//      % compute area of a simple shape
//      poly = [10 10;30 10;30 20;10 20];
//      area = polygonArea(poly)
//      area =
//          200
//
//      % compute area of CW polygon
//      area2 = polygonArea(poly(end:-1:1, :))
//      area2 =
//          -200
//
//      % Computes area of a paper hen
//      x = [0 10 20  0 -10 -20 -10 -10  0];
//      y = [0  0 10 10  20  10  10  0 -10];
//      poly = [x' y'];
//      area = polygonArea(poly)
//      area =
//         400
//
//      % Area of unit square with 25% hole
//      pccw = [0 0; 1 0; 1 1; 0 1];
//      pcw = pccw([1 4 3 2], :) * .5 + .25;
//      polygonArea ([pccw; nan(1,2); pcw])
//      ans =
//         0.75
//
//    References
//    algo adapted from P. Bourke web page
//    http://paulbourke.net/geometry/polygonmesh/
//
//    See also:
//    polygons2d, polygonCentroid, polygonSecondAreaMoments, triangleArea
// Arguments    : const double poly_data[]
//                const int poly_size[2]
// Return Type  : double
//
double polygonArea(const double poly_data[], const int poly_size[2])
{
  double area;
  int loop_ub;
  int i5;
  boolean_T b_data[8];
  int i2;
  boolean_T y[2];
  int iy;
  int i;
  boolean_T b_y;
  int i1;
  boolean_T exitg1;
  double px_data[4];
  double py_data[4];
  signed char y_data[3];
  signed char iNext_data[4];
  double c_y;
  boolean_T b_b_data[8];
  double x_data[4];
  cell_wrap_3 r1;
  emxArray_cell_wrap_3_5 polygons;
  boolean_T b_x_data[4];
  unsigned char ii_data[4];
  int inds_data[6];

  //    ---------
  //    author : David Legland
  //    INRA - TPV URPOI - BIA IMASTE
  //    created the 05/05/2004.
  //
  //    HISTORY
  //    25/04/2005: add support for multiple polygons
  //    12/10/2007: update doc
  //  Process special cases
  //  in case of polygon sets, computes the sum of polygon areas
  //  check there are enough points
  if (poly_size[0] < 2) {
    area = 0.0;
  } else {
    //  case of polygons with holes -> computes the sum of areas
    loop_ub = poly_size[0] * poly_size[1];
    for (i5 = 0; i5 < loop_ub; i5++) {
      b_data[i5] = rtIsNaN(poly_data[i5]);
    }

    for (i5 = 0; i5 < 2; i5++) {
      y[i5] = false;
    }

    i2 = 1;
    iy = -1;
    for (i = 0; i < 2; i++) {
      i1 = i2;
      i2 += poly_size[0];
      iy++;
      exitg1 = false;
      while ((!exitg1) && (i1 <= i2 - 1)) {
        if (b_data[i1 - 1]) {
          y[iy] = true;
          exitg1 = true;
        } else {
          i1++;
        }
      }
    }

    b_y = true;
    i2 = 1;
    exitg1 = false;
    while ((!exitg1) && (i2 < 3)) {
      if (!y[i2 - 1]) {
        b_y = false;
        exitg1 = true;
      } else {
        i2++;
      }
    }

    if (b_y) {
      // SPLITPOLYGONS Convert a NaN separated polygon list to a cell array of polygons 
      //
      //    POLYGONS = splitPolygons(POLYGON);
      //    POLYGON is a N-by-2 array of points, possibly with pairs of NaN values. 
      //    The functions separates each component separated by NaN values, and
      //    returns a cell array of polygons.
      //
      //    See also:
      //    polygons2d
      //
      //  ------
      //  Author: David Legland
      //  e-mail: david.legland@inra.fr
      //  Created: 2007-10-12,    using Matlab 7.4.0.287 (R2007a)
      //  Copyright 2007 INRA - BIA PV Nantes - MIAJ Jouy-en-Josas.
      iy = poly_size[0] << 1;
      loop_ub = poly_size[0] << 1;
      for (i5 = 0; i5 < loop_ub; i5++) {
        b_b_data[i5] = rtIsNaN(poly_data[i5]);
      }

      c_y = b_b_data[0];
      for (i2 = 2; i2 <= iy; i2++) {
        c_y += (double)b_b_data[i2 - 1];
      }

      if (c_y == 0.0) {
        //  single polygon -> no break
        r1.f1.size[0] = poly_size[0];
        r1.f1.size[1] = 2;
        loop_ub = poly_size[0] * poly_size[1];
        if (0 <= loop_ub - 1) {
          memcpy(&r1.f1.data[0], &poly_data[0], (unsigned int)(loop_ub * (int)
                  sizeof(double)));
        }

        polygons.size[0] = 1;
        polygons.data[0] = r1;
      } else {
        //  find indices of NaN couples
        loop_ub = poly_size[0] * poly_size[1];
        for (i5 = 0; i5 < loop_ub; i5++) {
          b_data[i5] = rtIsNaN(poly_data[i5]);
        }

        i2 = poly_size[0];
        for (iy = 0; iy < i2; iy++) {
          px_data[iy] = b_data[iy];
        }

        for (iy = 0; iy < i2; iy++) {
          px_data[iy] += (double)b_data[i2 + iy];
        }

        loop_ub = (signed char)poly_size[0];
        for (i5 = 0; i5 < loop_ub; i5++) {
          b_x_data[i5] = (px_data[i5] > 0.0);
        }

        i2 = (signed char)poly_size[0];
        i1 = 0;
        iy = 1;
        exitg1 = false;
        while ((!exitg1) && (iy <= i2)) {
          if (b_x_data[iy - 1]) {
            i1++;
            ii_data[i1 - 1] = (unsigned char)iy;
            if (i1 >= i2) {
              exitg1 = true;
            } else {
              iy++;
            }
          } else {
            iy++;
          }
        }

        if (1 > i1) {
          loop_ub = 0;
          i2 = 0;
        } else {
          loop_ub = i1;
          i2 = i1;
        }

        for (i5 = 0; i5 < loop_ub; i5++) {
          px_data[i5] = ii_data[i5];
        }

        //  number of polygons
        //  iterate over NaN-separated regions to create new polygon
        inds_data[0] = 0;
        for (i5 = 0; i5 < i2; i5++) {
          inds_data[i5 + 1] = (int)px_data[i5];
        }

        inds_data[1 + i2] = poly_size[0] + 1;
        polygons.size[0] = (signed char)(loop_ub + 1);
        for (i = 0; i <= loop_ub; i++) {
          if ((double)inds_data[i] + 1.0 > (double)inds_data[i + 1] - 1.0) {
            i5 = 1;
            i2 = 0;
          } else {
            i5 = inds_data[i] + 1;
            i2 = (int)((double)inds_data[i + 1] - 1.0);
          }

          polygons.data[i].f1.size[0] = (i2 - i5) + 1;
          polygons.data[i].f1.size[1] = 2;
          iy = (i2 - i5) + 1;
          for (i2 = 0; i2 < 2; i2++) {
            for (i1 = 0; i1 < iy; i1++) {
              polygons.data[i].f1.data[i1 + polygons.data[i].f1.size[0] * i2] =
                poly_data[((i5 + i1) + poly_size[0] * i2) - 1];
            }
          }
        }
      }

      area = b_polygonArea(polygons.data, polygons.size);
    } else {
      //  Process single polygons or single rings
      //  extract coordinates
      //  polygon given as N-by-2 array
      loop_ub = poly_size[0];
      memcpy(&px_data[0], &poly_data[0], (unsigned int)(loop_ub * (int)sizeof
              (double)));
      loop_ub = poly_size[0];
      for (i5 = 0; i5 < loop_ub; i5++) {
        py_data[i5] = poly_data[i5 + poly_size[0]];
      }

      //  indices of next vertices
      loop_ub = (signed char)poly_size[0] - 2;
      for (i5 = 0; i5 <= loop_ub; i5++) {
        y_data[i5] = (signed char)(2 + (signed char)i5);
      }

      loop_ub = (signed char)poly_size[0] - 1;
      if (0 <= loop_ub - 1) {
        memcpy(&iNext_data[0], &y_data[0], (unsigned int)(loop_ub * (int)sizeof
                (signed char)));
      }

      iNext_data[(signed char)poly_size[0] - 1] = 1;

      //  compute area (vectorized version)
      loop_ub = poly_size[0];
      for (i5 = 0; i5 < loop_ub; i5++) {
        x_data[i5] = px_data[i5] * py_data[iNext_data[i5] - 1] -
          px_data[iNext_data[i5] - 1] * py_data[i5];
      }

      c_y = x_data[0];
      for (i2 = 2; i2 <= poly_size[0]; i2++) {
        c_y += x_data[i2 - 1];
      }

      area = c_y / 2.0;
    }
  }

  return area;
}

//
// File trailer for polygonArea.cpp
//
// [EOF]
//
