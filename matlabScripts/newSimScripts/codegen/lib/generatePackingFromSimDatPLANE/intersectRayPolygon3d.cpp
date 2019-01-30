//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: intersectRayPolygon3d.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include <math.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "intersectRayPolygon3d.h"
#include "polygonArea.h"
#include "inpolygon.h"
#include "sum.h"
#include "bsxfun.h"
#include "splitPolygons.h"
#include "dot.h"
#include "vectorNorm3d.h"
#include "projPointOnPlane.h"
#include "crossProduct3d.h"
#include "createPlane.h"

// Function Definitions

//
// INTERSECTRAYPOLYGON3D Intersection point of a 3D ray and a 3D polygon
//
//    INTER = intersectRayPolygon3d(RAY, POLY)
//    Compute coordinates of intersection point between the 3D ray RAY and
//    the 3D polygon POLY. RAY is a 1-by-6 row vector containing origin and
//    direction vector of the ray, POLY is a Np-by-3 array containing
//    coordinates of 3D polygon vertices.
//    INTER is a 1-by-3 row vector containing coordinates of intersection
//    point, or [NaN NaN NaN] if ray and polygon do not intersect.
//
//    INTERS = intersectRayPolygon3d(RAYS, POLY)
//    If RAYS is a N-by-6 array representing several rays, the result
//    INTERS is a N-by-3 array containing coordinates of intersection of each
//    ray with the polygon.
//
//    [INTER INSIDE] = intersectRayPolygon3d(RAY, POLY)
//    Also return a N-by-1 boolean array containing TRUE if both the polygon
//    and the corresponding ray contain the intersection point.
//
//    Example
//      % Compute intersection between a 3D ray and a 3D triangle
//      pts3d = [3 0 0; 0 6 0;0 0 9];
//      ray1 = [0 0 0 3 6 9];
//      inter = intersectRayPolygon3d(ray1, pts3d)
//      inter =
//            1   2   3
//
//      % keep only valid intersections with several rays
//      pts3d = [3 0 0; 0 6 0;0 0 9];
//      rays = [0 0 0 3 6 9;10 0 0 1 2 3;3 6 9 3 6 9];
//      [inter inside] = intersectRayPolygon3d(rays, pts3d);
//      inter(inside, :)
//      ans =
//            1   2   3
//
//    See Also
//    intersectRayPolygon, intersectLinePolygon3d, intersectLineTriangle3d
// Arguments    : const double ray[18]
//                const double poly[12]
//                double inter[9]
//                boolean_T inside[3]
// Return Type  : void
//
void intersectRayPolygon3d(const double ray[18], const double poly[12], double
  inter[9], boolean_T inside[3])
{
  int iy;
  double b_poly[9];
  double plane[9];
  int i2;
  double n[3];
  double dv0[9];
  double denom[3];
  double dv1[9];
  double dv2[3];
  double b_n;
  int N;
  int i1;
  boolean_T par;
  int i;
  boolean_T b_par[3];
  int tmp_data[3];
  double point[12];
  double y;
  double dv3[12];
  double b_plane[12];
  double x[4];
  double polygons_data[4];
  double pts2d[8];
  double c_plane[9];
  double pInt2d[6];
  boolean_T b[8];
  boolean_T exitg1;
  boolean_T inPoly[3];
  emxArray_cell_wrap_3_5 polygons;
  double areas_data[5];
  boolean_T a_data[5];
  boolean_T ccw_data[15];
  int polygons_size[1];
  signed char x_data[15];
  boolean_T in_data[15];
  boolean_T b_y[2];
  static const signed char iv0[4] = { 1, 2, 3, 0 };

  int b_tmp_data[3];

  //  ------
  //  Author: David Legland
  //  e-mail: david.legland@inra.fr
  //  Created: 2011-05-22,    using Matlab 7.9.0.529 (R2009b)
  //  Copyright 2011 INRA - Cepia Software Platform.
  //  supporting plane of polygon vertices
  for (iy = 0; iy < 3; iy++) {
    for (i2 = 0; i2 < 3; i2++) {
      b_poly[i2 + 3 * iy] = poly[i2 + (iy << 2)];
    }
  }

  createPlane(b_poly, plane);

  //  intersection of 3D ray with the plane
  // INTERSECTLINEPLANE Intersection point between a 3D line and a plane
  //
  //    PT = intersectLinePlane(LINE, PLANE)
  //    Returns the intersection point of the given line and the given plane.
  //    LINE:  [x0 y0 z0 dx dy dz]
  //    PLANE: [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2]
  //    PT:    [xi yi zi]
  //    If LINE and PLANE are parallel, return [NaN NaN NaN].
  //    If LINE (or PLANE) is a matrix with 6 (or 9) columns and N rows, result
  //    is an array of points with N rows and 3 columns.
  //
  //    PT = intersectLinePlane(LINE, PLANE, TOL)
  //    Specifies the tolerance factor to test if a line is parallel to a
  //    plane. Default is 1e-14.
  //
  //    Example
  //      % define horizontal plane through origin
  //      plane = [0 0 0   1 0 0   0 1 0];
  //      % intersection with a vertical line
  //      line = [2 3 4  0 0 1];
  //      intersectLinePlane(line, plane)
  //      ans =
  //         2   3   0
  //      % intersection with a line "parallel" to plane
  //      line = [2 3 4  1 2 0];
  //      intersectLinePlane(line, plane)
  //      ans =
  //        NaN  NaN  NaN
  //
  //    See also:
  //    lines3d, planes3d, points3d, clipLine3d
  //
  //    ---------
  //    author : David Legland
  //    INRA - TPV URPOI - BIA IMASTE
  //    created the 17/02/2005.
  //
  //    HISTORY
  //    24/11/2005 add support for multiple input
  //    23/06/2006 correction from Songbai Ji allowing different number of
  //        lines or plane if other input has one row
  //    14/12/2006 correction for parallel lines and plane normals
  //    05/01/2007 fixup for parallel lines and plane normals
  //    24/04/2007 rename as 'intersectLinePlane'
  //    11/19/2010 Added bsxfun functionality for improved speed (Sven Holcombe) 
  //    01/02/2011 code cleanup, add option for tolerance, update doc
  //  extract tolerance if needed
  //  unify sizes of data
  //  N planes and M lines not allowed
  //  plane normal
  crossProduct3d(*(double (*)[3])&plane[3], *(double (*)[3])&plane[6], n);

  //  difference between origins of plane and line
  //  dot product of line direction with plane normal
  b_bsxfun(n, *(double (*)[9])&ray[9], dv0);
  sum(dv0, denom);

  //  relative position of intersection point on line (can be inf in case of a
  //  line parallel to the plane)
  bsxfun(*(double (*)[3])&plane[0], *(double (*)[9])&ray[0], dv0);
  b_bsxfun(n, dv0, dv1);
  sum(dv1, dv2);

  //  compute coord of intersection point
  for (iy = 0; iy < 3; iy++) {
    b_n = dv2[iy] / denom[iy];
    b_poly[iy] = b_n;
    b_poly[3 + iy] = b_n;
    b_poly[6 + iy] = b_n;
  }

  c_bsxfun(b_poly, *(double (*)[9])&ray[9], dv0);
  d_bsxfun(*(double (*)[9])&ray[0], dv0, inter);

  //  set indices of line and plane which are parallel to NaN
  N = 0;
  for (i1 = 0; i1 < 3; i1++) {
    par = (fabs(denom[i1]) < 1.0E-14);
    if (par) {
      N++;
    }

    b_par[i1] = par;
  }

  i2 = 0;
  for (i = 0; i < 3; i++) {
    if (b_par[i]) {
      tmp_data[i2] = i + 1;
      i2++;
    }
  }

  for (iy = 0; iy < 3; iy++) {
    for (i2 = 0; i2 < N; i2++) {
      inter[(tmp_data[i2] + 3 * iy) - 1] = rtNaN;
    }
  }

  //  project all points on reference plane
  memcpy(&point[0], &poly[0], 12U * sizeof(double));
  projPointOnPlane(point, plane);

  // PLANEPOSITION Compute position of a point on a plane
  //
  //    PT2 = planePosition(POINT, PLANE)
  //    POINT has format [X Y Z], and plane has format
  //    [X0 Y0 Z0  DX1 DY1 DZ1  DX2 DY2 DZ2], where :
  //    - (X0, Y0, Z0) is a point belonging to the plane
  //    - (DX1, DY1, DZ1) is a first direction vector
  //    - (DX2, DY2, DZ2) is a second direction vector
  //
  //    Result PT2 has the form [XP YP], with [XP YP] coordinate of the point
  //    in the coordinate system of the plane.
  //
  //
  //    CAUTION:
  //    WORKS ONLY FOR PLANES WITH ORTHOGONAL DIRECTION VECTORS
  //
  //    Example
  //      plane = [10 20 30  1 0 0  0 1 0];
  //      point = [13 24 35];
  //      pos = planePosition(point, plane)
  //      pos =
  //          3   4
  //
  //    See also:
  //    geom3d, planes3d, points3d, planePoint
  //    ---------
  //    author : David Legland
  //    INRA - TPV URPOI - BIA IMASTE
  //    created the 21/02/2005.
  //
  //    HISTORY
  //    24/11/2005 add support for multiple inputs
  //  size of input arguments
  //  check inputs have compatible sizes
  //  origin and direction vectors of the plane
  //  Compute dot products with direction vectors of the plane
  //  we have npl == 1 and npt > 1
  b_n = vectorNorm3d(*(double (*)[3])&plane[3]);
  y = vectorNorm3d(*(double (*)[3])&plane[6]);

  //  % old version:
  //  s = dot(point-p0, d1, 2) ./ vectorNorm3d(d1);
  //  t = dot(point-p0, d2, 2) ./ vectorNorm3d(d2);
  e_bsxfun(point, *(double (*)[3])&plane[0], dv3);
  for (iy = 0; iy < 3; iy++) {
    for (i2 = 0; i2 < 4; i2++) {
      b_plane[i2 + (iy << 2)] = plane[3 + iy] / b_n;
    }
  }

  dot(dv3, b_plane, x);
  e_bsxfun(point, *(double (*)[3])&plane[0], dv3);
  for (iy = 0; iy < 3; iy++) {
    for (i2 = 0; i2 < 4; i2++) {
      b_plane[i2 + (iy << 2)] = plane[6 + iy] / y;
    }
  }

  dot(dv3, b_plane, polygons_data);
  for (iy = 0; iy < 4; iy++) {
    pts2d[iy] = x[iy];
    pts2d[4 + iy] = polygons_data[iy];
  }

  memcpy(&b_poly[0], &inter[0], 9U * sizeof(double));
  b_projPointOnPlane(b_poly, plane);

  // PLANEPOSITION Compute position of a point on a plane
  //
  //    PT2 = planePosition(POINT, PLANE)
  //    POINT has format [X Y Z], and plane has format
  //    [X0 Y0 Z0  DX1 DY1 DZ1  DX2 DY2 DZ2], where :
  //    - (X0, Y0, Z0) is a point belonging to the plane
  //    - (DX1, DY1, DZ1) is a first direction vector
  //    - (DX2, DY2, DZ2) is a second direction vector
  //
  //    Result PT2 has the form [XP YP], with [XP YP] coordinate of the point
  //    in the coordinate system of the plane.
  //
  //
  //    CAUTION:
  //    WORKS ONLY FOR PLANES WITH ORTHOGONAL DIRECTION VECTORS
  //
  //    Example
  //      plane = [10 20 30  1 0 0  0 1 0];
  //      point = [13 24 35];
  //      pos = planePosition(point, plane)
  //      pos =
  //          3   4
  //
  //    See also:
  //    geom3d, planes3d, points3d, planePoint
  //    ---------
  //    author : David Legland
  //    INRA - TPV URPOI - BIA IMASTE
  //    created the 21/02/2005.
  //
  //    HISTORY
  //    24/11/2005 add support for multiple inputs
  //  size of input arguments
  //  check inputs have compatible sizes
  //  origin and direction vectors of the plane
  //  Compute dot products with direction vectors of the plane
  //  we have npl == 1 and npt > 1
  b_n = vectorNorm3d(*(double (*)[3])&plane[3]);
  y = vectorNorm3d(*(double (*)[3])&plane[6]);

  //  % old version:
  //  s = dot(point-p0, d1, 2) ./ vectorNorm3d(d1);
  //  t = dot(point-p0, d2, 2) ./ vectorNorm3d(d2);
  f_bsxfun(b_poly, *(double (*)[3])&plane[0], dv0);
  for (iy = 0; iy < 3; iy++) {
    for (i2 = 0; i2 < 3; i2++) {
      c_plane[i2 + 3 * iy] = plane[3 + iy] / b_n;
    }
  }

  b_dot(dv0, c_plane, dv2);
  f_bsxfun(b_poly, *(double (*)[3])&plane[0], dv0);
  for (iy = 0; iy < 3; iy++) {
    for (i2 = 0; i2 < 3; i2++) {
      c_plane[i2 + 3 * iy] = plane[6 + iy] / y;
    }

    pInt2d[iy] = dv2[iy];
  }

  b_dot(dv0, c_plane, n);
  for (iy = 0; iy < 3; iy++) {
    pInt2d[3 + iy] = n[iy];
  }

  //  need to check polygon orientation
  // ISPOINTINPOLYGON Test if a point is located inside a polygon
  //
  //    B = isPointInPolygon(POINT, POLYGON)
  //    Returns true if the point is located within the given polygon.
  //
  //    This function is simply a wrapper for the function inpolygon, to avoid
  //    decomposition of point and polygon coordinates.
  //
  //    Example
  //      pt1 = [30 20];
  //      pt2 = [30 5];
  //      poly = [10 10;50 10;50 50;10 50];
  //      isPointInPolygon([pt1;pt2], poly)
  //      ans =
  //           1
  //           0
  //
  //      poly = [0 0; 10 0;10 10;0 10;NaN NaN;3 3;3 7;7 7;7 3];
  //      pts = [5 1;5 4];
  //      isPointInPolygon(pts, poly);
  //      ans =
  //           1
  //           0
  //
  //
  //    See also
  //    points2d, polygons2d, inpolygon, isPointInTriangle
  //
  //  ------
  //  Author: David Legland
  //  e-mail: david.legland@inra.fr
  //  Created: 2009-06-19,    using Matlab 7.7.0.471 (R2008b)
  //  Copyright 2009 INRA - Cepia Software Platform.
  //    HISTORY
  //    2013-04-24 add support for multiply connected polygons
  //  In case of a multiple polygon, decompose into a set of contours, and
  //  performs test for each contour
  for (iy = 0; iy < 8; iy++) {
    b[iy] = rtIsNaN(pts2d[iy]);
  }

  par = false;
  i1 = 0;
  exitg1 = false;
  while ((!exitg1) && (i1 < 8)) {
    if (b[i1]) {
      par = true;
      exitg1 = true;
    } else {
      i1++;
    }
  }

  if (par) {
    //  transform as a cell array of simple polygons
    splitPolygons(pts2d, polygons.data, polygons.size);
    N = polygons.size[0] - 1;

    //  compute orientation of polygon, and format to have Np*N matrix
    for (i = 0; i <= N; i++) {
      areas_data[i] = polygonArea(polygons.data[i].f1.data, polygons.data[i].
        f1.size);
    }

    i2 = polygons.size[0];
    for (iy = 0; iy < i2; iy++) {
      a_data[iy] = (areas_data[iy] > 0.0);
    }

    for (i2 = 0; i2 < polygons.size[0]; i2++) {
      iy = i2 * 3;
      for (i1 = 0; i1 < 3; i1++) {
        ccw_data[iy + i1] = a_data[i2];
      }
    }

    //  test if point inside each polygon
    for (i = 0; i <= N; i++) {
      i2 = polygons.data[i].f1.size[0];
      polygons_size[0] = polygons.data[i].f1.size[0];
      for (iy = 0; iy < i2; iy++) {
        x[iy] = polygons.data[i].f1.data[iy];
      }

      i2 = polygons.data[i].f1.size[0];
      for (iy = 0; iy < i2; iy++) {
        polygons_data[iy] = polygons.data[i].f1.data[iy + polygons.data[i].
          f1.size[0]];
      }

      inpolygon(*(double (*)[3])&pInt2d[0], *(double (*)[3])&pInt2d[3], x,
                polygons_size, polygons_data, *(boolean_T (*)[3])&in_data[3 * i]);
    }

    //  count polygons containing point, weighted by polygon orientation
    i2 = 3 * polygons.size[0];
    for (iy = 0; iy < i2; iy++) {
      x_data[iy] = (signed char)(in_data[iy] * ccw_data[iy] - in_data[iy] *
        !ccw_data[iy]);
    }

    for (i2 = 0; i2 < 3; i2++) {
      n[i2] = x_data[i2];
    }

    for (i1 = 2; i1 <= polygons.size[0]; i1++) {
      iy = (i1 - 1) * 3;
      for (i2 = 0; i2 < 3; i2++) {
        n[i2] += (double)x_data[iy + i2];
      }
    }

    for (iy = 0; iy < 3; iy++) {
      inPoly[iy] = (n[iy] > 0.0);
    }
  } else {
    //  standard test for simple polygons
    b_inpolygon(*(double (*)[3])&pInt2d[0], *(double (*)[3])&pInt2d[3], *(double
      (*)[4])&pts2d[0], *(double (*)[4])&pts2d[4], inPoly);
  }

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
  //
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
  //  case of polygons with holes -> computes the sum of areas
  for (iy = 0; iy < 8; iy++) {
    b[iy] = rtIsNaN(pts2d[iy]);
  }

  for (iy = 0; iy < 2; iy++) {
    b_y[iy] = false;
  }

  i2 = 1;
  iy = -1;
  for (i = 0; i < 2; i++) {
    i1 = i2;
    i2 += 4;
    iy++;
    exitg1 = false;
    while ((!exitg1) && (i1 <= i2 - 1)) {
      if (b[i1 - 1]) {
        b_y[iy] = true;
        exitg1 = true;
      } else {
        i1++;
      }
    }
  }

  par = true;
  i1 = 1;
  exitg1 = false;
  while ((!exitg1) && (i1 < 3)) {
    if (!b_y[i1 - 1]) {
      par = false;
      exitg1 = true;
    } else {
      i1++;
    }
  }

  if (par) {
    splitPolygons(pts2d, polygons.data, polygons.size);
    b_n = b_polygonArea(polygons.data, polygons.size);
  } else {
    //  Process single polygons or single rings
    //  extract coordinates
    //  polygon given as N-by-2 array
    //  indices of next vertices
    //  compute area (vectorized version)
    for (iy = 0; iy < 4; iy++) {
      x[iy] = pts2d[iy] * pts2d[4 + iv0[iy]] - pts2d[iv0[iy]] * pts2d[4 + iy];
    }

    b_n = x[0];
    for (i1 = 0; i1 < 3; i1++) {
      b_n += x[i1 + 1];
    }

    b_n /= 2.0;
  }

  par = (b_n < 0.0);

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
  for (i1 = 0; i1 < 3; i1++) {
    for (i2 = 0; i2 < 3; i2++) {
      plane[i2 + 3 * i1] = inter[i2 + 3 * i1] - ray[i2 + 3 * i1];
    }

    inPoly[i1] = ((inPoly[i1] ^ par) != 0);
  }

  //  direction vector of the line
  //  precompute and check validity of denominator
  for (i1 = 0; i1 < 9; i1++) {
    b_poly[i1] = ray[i1 % 3 + 3 * (3 + i1 / 3)] * ray[i1 % 3 + 3 * (3 + i1 / 3)];
  }

  sum(b_poly, denom);

  //  compute position using dot product normalized with norm of line vector.
  c_bsxfun(plane, *(double (*)[9])&ray[9], dv0);
  sum(dv0, n);

  //  position on a degenerated line is set to 0
  //  intersection points outside the polygon are set to NaN
  N = 0;
  for (i1 = 0; i1 < 3; i1++) {
    b_n = denom[i1];
    if (denom[i1] < 2.2204460492503131E-16) {
      b_n = 1.0;
    }

    b_n = n[i1] / b_n;
    if (denom[i1] < 2.2204460492503131E-16) {
      b_n = 0.0;
    }

    par = (b_n >= 0.0);
    inside[i1] = (inPoly[i1] && par);
    if (!(inPoly[i1] && par)) {
      N++;
    }

    b_par[i1] = par;
  }

  i2 = 0;
  for (i = 0; i < 3; i++) {
    if (!(inPoly[i] && b_par[i])) {
      b_tmp_data[i2] = i + 1;
      i2++;
    }
  }

  for (iy = 0; iy < 3; iy++) {
    for (i2 = 0; i2 < N; i2++) {
      inter[(b_tmp_data[i2] + 3 * iy) - 1] = rtNaN;
    }
  }
}

//
// File trailer for intersectRayPolygon3d.cpp
//
// [EOF]
//
