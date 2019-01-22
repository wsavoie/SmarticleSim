//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: projPointOnPlane.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include <string.h>
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "projPointOnPlane.h"
#include "bsxfun.h"
#include "sum.h"
#include "crossProduct3d.h"

// Function Definitions

//
// PROJPOINTONPLANE Return the orthogonal projection of a point on a plane
//
//    PT2 = projPointOnPlane(PT1, PLANE);
//    Compute the (orthogonal) projection of point PT1 onto the plane PLANE,
//    given as [X0 Y0 Z0  VX1 VY1 VZ1  VX2 VY2 VZ2] (origin point, first
//    direction vector, second directionvector).
//
//    The function is fully vectorized, in that multiple points may be
//    projected onto multiple planes in a single call, returning multiple
//    points. With the exception of the second dimension (where
//    SIZE(PT1,2)==3, and SIZE(PLANE,2)==9), each dimension of PT1 and PLANE
//    must either be equal or one, similar to the requirements of BSXFUN. In
//    basic usage, point PT1 is a [N*3] array, and PLANE is a [N*9] array
//    (see createPlane for details). Result PT2 is a [N*3] array, containing
//    coordinates of orthogonal projections of PT1 onto planes PLANE. In
//    vectorised usage, PT1 is an [N*3*M*P...] matrix, and PLANE is an
//    [X*9*Y...] matrix, where (N,X), (M,Y), etc, are either equal pairs, or
//    one of the two is one.
//
//    See also:
//    planes3d, points3d, planePosition, intersectLinePlane
// Arguments    : double point[9]
//                const double plane[9]
// Return Type  : void
//
void b_projPointOnPlane(double point[9], const double plane[9])
{
  double normals[3];
  int k;
  double y;
  double z1[3];
  double c[9];
  double dv6[9];
  double t[3];
  double b_point[9];
  int b_k;

  //    ---------
  //    author : David Legland
  //    INRA - TPV URPOI - BIA IMASTE
  //    created the 18/02/2005.
  //
  //    HISTORY
  //    21/08/2006: debug support for multiple points or planes
  //    22/04/2013: uses bsxfun for mult. pts/planes in all dimensions (Sven H)
  //  Unpack the planes into origins and normals, keeping original shape
  crossProduct3d(*(double (*)[3])&plane[3], *(double (*)[3])&plane[6], normals);

  //  difference between origins of plane and point
  //  relative position of point on normal's line
  for (k = 0; k < 3; k++) {
    z1[k] = normals[k] * normals[k];
  }

  y = z1[0];
  for (k = 0; k < 2; k++) {
    y += z1[k + 1];
  }

  bsxfun(*(double (*)[3])&plane[0], point, c);
  b_bsxfun(normals, c, dv6);
  sum(dv6, z1);
  for (k = 0; k < 3; k++) {
    t[k] = z1[k] / y;
  }

  //  add relative difference to project point back to plane
  for (k = 0; k < 3; k++) {
    for (b_k = 0; b_k < 3; b_k++) {
      c[b_k + 3 * k] = t[b_k] * normals[k];
    }
  }

  memcpy(&b_point[0], &point[0], 9U * sizeof(double));
  d_bsxfun(b_point, c, point);
}

//
// PROJPOINTONPLANE Return the orthogonal projection of a point on a plane
//
//    PT2 = projPointOnPlane(PT1, PLANE);
//    Compute the (orthogonal) projection of point PT1 onto the plane PLANE,
//    given as [X0 Y0 Z0  VX1 VY1 VZ1  VX2 VY2 VZ2] (origin point, first
//    direction vector, second directionvector).
//
//    The function is fully vectorized, in that multiple points may be
//    projected onto multiple planes in a single call, returning multiple
//    points. With the exception of the second dimension (where
//    SIZE(PT1,2)==3, and SIZE(PLANE,2)==9), each dimension of PT1 and PLANE
//    must either be equal or one, similar to the requirements of BSXFUN. In
//    basic usage, point PT1 is a [N*3] array, and PLANE is a [N*9] array
//    (see createPlane for details). Result PT2 is a [N*3] array, containing
//    coordinates of orthogonal projections of PT1 onto planes PLANE. In
//    vectorised usage, PT1 is an [N*3*M*P...] matrix, and PLANE is an
//    [X*9*Y...] matrix, where (N,X), (M,Y), etc, are either equal pairs, or
//    one of the two is one.
//
//    See also:
//    planes3d, points3d, planePosition, intersectLinePlane
// Arguments    : double point[12]
//                const double plane[9]
// Return Type  : void
//
void projPointOnPlane(double point[12], const double plane[9])
{
  double normals[3];
  int k;
  int b_k;
  double y[4];
  double c[12];
  double dp[12];
  int xoffset;
  double b_y;
  double z1[3];
  double t[4];

  //    ---------
  //    author : David Legland
  //    INRA - TPV URPOI - BIA IMASTE
  //    created the 18/02/2005.
  //
  //    HISTORY
  //    21/08/2006: debug support for multiple points or planes
  //    22/04/2013: uses bsxfun for mult. pts/planes in all dimensions (Sven H)
  //  Unpack the planes into origins and normals, keeping original shape
  crossProduct3d(*(double (*)[3])&plane[3], *(double (*)[3])&plane[6], normals);

  //  difference between origins of plane and point
  //  relative position of point on normal's line
  for (k = 0; k < 3; k++) {
    for (b_k = 0; b_k < 4; b_k++) {
      dp[b_k + (k << 2)] = plane[k] - point[b_k + (k << 2)];
      c[b_k + (k << 2)] = normals[k] * dp[b_k + (k << 2)];
    }
  }

  for (b_k = 0; b_k < 4; b_k++) {
    y[b_k] = c[b_k];
  }

  for (k = 0; k < 2; k++) {
    xoffset = (k + 1) << 2;
    for (b_k = 0; b_k < 4; b_k++) {
      y[b_k] += c[xoffset + b_k];
    }
  }

  for (k = 0; k < 3; k++) {
    z1[k] = normals[k] * normals[k];
  }

  b_y = z1[0];
  for (k = 0; k < 2; k++) {
    b_y += z1[k + 1];
  }

  for (k = 0; k < 4; k++) {
    t[k] = y[k] / b_y;
  }

  //  add relative difference to project point back to plane
  for (k = 0; k < 3; k++) {
    for (b_k = 0; b_k < 4; b_k++) {
      c[b_k + (k << 2)] = t[b_k] * normals[k];
    }
  }

  memcpy(&dp[0], &point[0], 12U * sizeof(double));
  for (k = 0; k < 3; k++) {
    for (b_k = 0; b_k < 4; b_k++) {
      point[b_k + (k << 2)] = dp[b_k + (k << 2)] + c[b_k + (k << 2)];
    }
  }
}

//
// File trailer for projPointOnPlane.cpp
//
// [EOF]
//
