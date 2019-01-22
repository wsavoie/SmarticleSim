//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: distancePointEdge3d.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include <math.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "distancePointEdge3d.h"
#include "sum.h"
#include "bsxfun.h"

// Function Definitions

//
// DISTANCEPOINTEDGE3D Minimum distance between a 3D point and a 3D edge
//
//    DIST = distancePointEdge3d(POINT, EDGE);
//    Return the euclidean distance between edge EDGE and point POINT.
//    EDGE has the form: [x1 y1 z1 x2 y2 z2], and POINT is [x y z].
//
//    If EDGE is N-by-6 array, result is N-by-1 array computed for each edge.
//    If POINT is a N-by-3 array, the result is computed for each point.
//    If both POINT and EDGE are array, they must have the same number of
//    rows, and the result is computed for each couple point(i,:);edge(i,:).
//
//    [DIST POS] = distancePointEdge3d(POINT, EDGE);
//    Also returns the position of closest point on the edge. POS is
//    comprised between 0 (first point) and 1 (last point).
//
//    See also:
//    edges3d, points3d, distancePoints3d, distancePointLine3d
// Arguments    : const double point_data[]
//                const int point_size[2]
//                const double edge_data[]
//                const int edge_size[2]
//                double dist_data[]
//                int dist_size[1]
// Return Type  : void
//
void distancePointEdge3d(const double point_data[], const int point_size[2],
  const double edge_data[], const int edge_size[2], double dist_data[], int
  dist_size[1])
{
  int loop_ub;
  int k;
  int nx;
  double vl_data[9];
  int result;
  boolean_T empty_non_axis_sizes;
  int acoef;
  int bcoef;
  int b_loop_ub;
  double b_edge_data[9];
  double result_data[18];
  int result_size[2];
  double dp_data[9];
  int dp_size[2];
  int y_size[2];
  double b_vl_data[9];
  double denom_data[3];
  int denom_size[1];
  boolean_T invalidLine_data[3];
  double b_denom_data;
  int sck;
  int vl_size[2];
  double a_data[3];
  int a_size[1];
  signed char csz_idx_0;
  double t_data[3];
  int b_vl_size[2];
  int c_vl_size[2];

  //    ---------
  //    author : David Legland
  //    INRA - CEPIA URPOI - MIA MathCell
  //    created the 07/04/2004.
  //
  //    HISTORY
  //    2005-06-24 rename, and change arguments sequence
  //    2009-04-30 add possibility to return position of closest point
  //    2011-04-14 add checkup for degenerate edges, improve speed, update doc
  //  direction vector of each edge
  loop_ub = edge_size[0];
  for (k = 0; k < 3; k++) {
    for (nx = 0; nx < loop_ub; nx++) {
      vl_data[nx + loop_ub * k] = edge_data[nx + edge_size[0] * (3 + k)] -
        edge_data[nx + edge_size[0] * k];
    }
  }

  //  compute position of points projected on the supporting line
  //  (Size of t is the max number of edges or points)
  if (!(edge_size[0] == 0)) {
    result = edge_size[0];
  } else {
    result = 0;
  }

  empty_non_axis_sizes = (result == 0);
  if (empty_non_axis_sizes || (!(edge_size[0] == 0))) {
    acoef = 3;
  } else {
    acoef = 0;
  }

  if (empty_non_axis_sizes || (!(edge_size[0] == 0))) {
    bcoef = 3;
  } else {
    bcoef = 0;
  }

  b_loop_ub = edge_size[0];
  for (k = 0; k < 3; k++) {
    for (nx = 0; nx < b_loop_ub; nx++) {
      b_edge_data[nx + b_loop_ub * k] = edge_data[nx + edge_size[0] * k];
    }
  }

  for (k = 0; k < acoef; k++) {
    for (nx = 0; nx < result; nx++) {
      result_data[nx + result * k] = b_edge_data[nx + result * k];
    }
  }

  for (k = 0; k < bcoef; k++) {
    for (nx = 0; nx < result; nx++) {
      result_data[nx + result * (k + acoef)] = vl_data[nx + result * k];
    }
  }

  // LINEPOSITION3D Return the position of a 3D point projected on a 3D line
  //
  //    T = linePosition3d(POINT, LINE)
  //    Computes position of point POINT on the line LINE, relative to origin
  //    point and direction vector of the line.
  //    LINE has the form [x0 y0 z0 dx dy dy],
  //    POINT has the form [x y z], and is assumed to belong to line.
  //    The result T is the value such that POINT = LINE(1:3) + T * LINE(4:6).
  //    If POINT does not belong to LINE, the position of its orthogonal
  //    projection is computed instead.
  //
  //    T = linePosition3d(POINT, LINES)
  //    If LINES is an array of NL lines, return NL positions, corresponding to
  //    each line.
  //
  //    T = linePosition3d(POINTS, LINE)
  //    If POINTS is an array of NP points, return NP positions, corresponding
  //    to each point.
  //
  //    See also:
  //    lines3d, createLine3d, distancePointLine3d, projPointOnLine3d
  //
  //    ---------
  //    author : David Legland
  //    INRA - TPV URPOI - BIA IMASTE
  //    created the 17/02/2005.
  //
  //    HISTORY
  //    05/01/2007 update doc
  //    28/10/2010 change to bsxfun calculation for arbitrary input sizes
  //        (Thanks to Sven Holcombe)
  //  vector from line origin to point
  result_size[0] = result;
  result_size[1] = 3;
  for (k = 0; k < 3; k++) {
    for (nx = 0; nx < result; nx++) {
      b_edge_data[nx + result * k] = result_data[nx + result * k];
    }
  }

  g_bsxfun(point_data, point_size, b_edge_data, result_size, dp_data, dp_size);

  //  direction vector of the line
  for (k = 0; k < 3; k++) {
    for (nx = 0; nx < result; nx++) {
      b_vl_data[nx + result * k] = result_data[nx + result * (3 + k)];
    }
  }

  //  precompute and check validity of denominator
  y_size[0] = (signed char)result;
  y_size[1] = 3;
  nx = (signed char)result * 3;
  for (b_loop_ub = 0; b_loop_ub < nx; b_loop_ub++) {
    b_edge_data[b_loop_ub] = b_vl_data[b_loop_ub] * b_vl_data[b_loop_ub];
  }

  b_sum(b_edge_data, y_size, denom_data, denom_size);
  b_loop_ub = denom_size[0];
  for (k = 0; k < b_loop_ub; k++) {
    invalidLine_data[k] = (denom_data[k] < 2.2204460492503131E-16);
  }

  nx = denom_size[0];
  for (acoef = 0; acoef < nx; acoef++) {
    b_denom_data = denom_data[acoef];
    if (denom_data[acoef] < 2.2204460492503131E-16) {
      b_denom_data = 1.0;
    }

    denom_data[acoef] = b_denom_data;
  }

  //  compute position using dot product normalized with norm of line vector.
  if (result == 1) {
    sck = dp_size[0];
  } else if (dp_size[0] == 1) {
    sck = result;
  } else if (dp_size[0] == result) {
    sck = dp_size[0];
  } else if (result < dp_size[0]) {
    sck = result;
  } else {
    sck = dp_size[0];
  }

  vl_size[0] = (signed char)sck;
  vl_size[1] = 3;
  if ((signed char)sck != 0) {
    acoef = (dp_size[0] != 1);
    bcoef = (result != 1);
    for (b_loop_ub = 0; b_loop_ub < 3; b_loop_ub++) {
      for (k = 0; k < (signed char)sck; k++) {
        b_vl_data[k + (signed char)sck * b_loop_ub] = dp_data[acoef * k +
          dp_size[0] * b_loop_ub] * result_data[bcoef * k + result * (3 +
          b_loop_ub)];
      }
    }
  }

  b_sum(b_vl_data, vl_size, a_data, a_size);
  nx = denom_size[0];
  acoef = a_size[0];
  if (nx < acoef) {
    acoef = nx;
  }

  if (denom_size[0] == 1) {
    sck = a_size[0];
  } else if (a_size[0] == 1) {
    sck = denom_size[0];
  } else if (a_size[0] == denom_size[0]) {
    sck = a_size[0];
  } else {
    sck = acoef;
  }

  nx = denom_size[0];
  acoef = a_size[0];
  if (nx < acoef) {
    acoef = nx;
  }

  if (denom_size[0] == 1) {
    csz_idx_0 = (signed char)a_size[0];
  } else if (a_size[0] == 1) {
    csz_idx_0 = (signed char)denom_size[0];
  } else if (a_size[0] == denom_size[0]) {
    csz_idx_0 = (signed char)a_size[0];
  } else {
    csz_idx_0 = (signed char)acoef;
  }

  if (csz_idx_0 != 0) {
    acoef = (a_size[0] != 1);
    bcoef = (denom_size[0] != 1);
    for (b_loop_ub = 0; b_loop_ub < (signed char)sck; b_loop_ub++) {
      t_data[b_loop_ub] = a_data[acoef * b_loop_ub] / denom_data[bcoef *
        b_loop_ub];
    }
  }

  //  position on a degenerated line is set to 0
  nx = denom_size[0];
  for (acoef = 0; acoef < nx; acoef++) {
    if (invalidLine_data[acoef]) {
      t_data[acoef] = 0.0;
    }
  }

  //  change position to ensure projected point is located on the edge
  nx = csz_idx_0;
  for (acoef = 0; acoef < nx; acoef++) {
    b_denom_data = t_data[acoef];
    if (t_data[acoef] < 0.0) {
      b_denom_data = 0.0;
    }

    t_data[acoef] = b_denom_data;
  }

  nx = csz_idx_0;
  for (acoef = 0; acoef < nx; acoef++) {
    b_denom_data = t_data[acoef];
    if (t_data[acoef] > 1.0) {
      b_denom_data = 1.0;
    }

    t_data[acoef] = b_denom_data;
  }

  //  difference of coordinates between projected point and base point
  b_loop_ub = csz_idx_0;
  for (k = 0; k < b_loop_ub; k++) {
    b_edge_data[k] = t_data[k] * vl_data[k];
  }

  b_loop_ub = csz_idx_0;
  for (k = 0; k < b_loop_ub; k++) {
    b_edge_data[k + csz_idx_0] = t_data[k] * vl_data[k + loop_ub];
  }

  b_loop_ub = csz_idx_0;
  for (k = 0; k < b_loop_ub; k++) {
    b_edge_data[k + (csz_idx_0 << 1)] = t_data[k] * vl_data[k + (loop_ub << 1)];
  }

  if ((signed char)sck == 1) {
    nx = edge_size[0];
  } else if (edge_size[0] == 1) {
    nx = (signed char)sck;
  } else if (edge_size[0] == (signed char)sck) {
    nx = edge_size[0];
  } else if ((signed char)sck < edge_size[0]) {
    nx = (signed char)sck;
  } else {
    nx = edge_size[0];
  }

  if ((signed char)nx != 0) {
    acoef = (edge_size[0] != 1);
    bcoef = ((signed char)sck != 1);
    for (b_loop_ub = 0; b_loop_ub < 3; b_loop_ub++) {
      for (k = 0; k < (signed char)nx; k++) {
        b_vl_data[k + (signed char)nx * b_loop_ub] = edge_data[acoef * k +
          edge_size[0] * b_loop_ub] + b_edge_data[bcoef * k + csz_idx_0 *
          b_loop_ub];
      }
    }
  }

  b_vl_size[0] = (signed char)nx;
  b_vl_size[1] = 3;
  loop_ub = (signed char)nx * 3;
  if (0 <= loop_ub - 1) {
    memcpy(&vl_data[0], &b_vl_data[0], (unsigned int)(loop_ub * (int)sizeof
            (double)));
  }

  g_bsxfun(point_data, point_size, vl_data, b_vl_size, b_vl_data, vl_size);

  //  compute distance between point and its projection on the edge
  c_vl_size[0] = vl_size[0];
  c_vl_size[1] = 3;
  loop_ub = vl_size[0] * vl_size[1];
  for (k = 0; k < loop_ub; k++) {
    vl_data[k] = b_vl_data[k] * b_vl_data[k];
  }

  b_sum(vl_data, c_vl_size, dist_data, dist_size);
  for (b_loop_ub = 0; b_loop_ub < dist_size[0]; b_loop_ub++) {
    dist_data[b_loop_ub] = sqrt(dist_data[b_loop_ub]);
  }
}

//
// File trailer for distancePointEdge3d.cpp
//
// [EOF]
//
