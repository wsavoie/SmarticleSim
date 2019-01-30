//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: createPlane.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "createPlane.h"
#include "crossProduct3d.h"
#include "normalizeVector3d.h"

// Function Definitions

//
// CREATEPLANE Create a plane in parametrized form
//
//    PLANE = createPlane(P1, P2, P3)
//    creates a plane containing the 3 points
//
//    PLANE = createPlane(PTS)
//    The 3 points are packed into a single 3x3 array.
//
//    PLANE = createPlane(P0, N);
//    Creates a plane from a point and from a normal to the plane. The
//    parameter N is given either as a 3D vector (1-by-3 row vector), or as
//    [THETA PHI], where THETA is the colatitute (angle with the vertical
//    axis) and PHI is angle with Ox axis, counted counter-clockwise (both
//    given in radians).
//
//    PLANE = createPlane(P0, Dip, DipDir);
//    Creates a plane from a point and from a dip and dip direction angles
//    of the plane. Parameters Dip and DipDir angles are given as numbers.
//    Dip : maximum inclination to the horizontal.
//    DipDir : direction of the horizontal trace of the line of dip,
//             measured clockwise from north.
//
//    The created plane data has the following format:
//    PLANE = [X0 Y0 Z0  DX1 DY1 DZ1  DX2 DY2 DZ2], with
//    - (X0, Y0, Z0) is a point belonging to the plane
//    - (DX1, DY1, DZ1) is a first direction vector
//    - (DX2, DY2, DZ2) is a second direction vector
//    The 2 direction vectors are normalized and orthogonal.
//
//    See also:
//    planes3d, medianPlane
//
//    ---------
//    author: David Legland
//    INRA - TPV URPOI - BIA IMASTE
//    created the 18/02/2005.
// Arguments    : const double varargin_1[9]
//                double plane[9]
// Return Type  : void
//
void createPlane(const double varargin_1[9], double plane[9])
{
  int k;
  double plane1[9];
  double d1[3];
  double normals[3];
  double dv4[3];
  double dv5[3];
  double d2[3];
  double t;
  double dp[3];
  double y;

  //    HISTORY
  //    24/11/2005 add possibility to pack points for plane creation
  //    21/08/2006 return normalized planes
  //    06/11/2006 update doc for planes created from normal
  //  3 points in a single array
  //  create direction vectors
  //  create plane
  for (k = 0; k < 3; k++) {
    plane1[k] = varargin_1[3 * k];
    plane1[k + 3] = varargin_1[1 + 3 * k] - varargin_1[3 * k];
    plane1[k + 6] = varargin_1[2 + 3 * k] - varargin_1[3 * k];
  }

  // NORMALIZEPLANE Normalize parametric representation of a plane
  //
  //    PLANE2 = normalizePlane(PLANE1);
  //    Transforms the plane PLANE1 in the following format:
  //    [X0 Y0 Z0  DX1 DY1 DZ1  DX2 DY2 DZ2], where:
  //    - (X0, Y0, Z0) is a point belonging to the plane
  //    - (DX1, DY1, DZ1) is a first direction vector
  //    - (DX2, DY2, DZ2) is a second direction vector
  //    into another plane, with the same format, but with:
  //    - (x0 y0 z0) is the closest point of plane to the origin
  //    - (DX1 DY1 DZ1) has norm equal to 1
  //    - (DX2 DY2 DZ2) has norm equal to 1 and is orthogonal to (DX1 DY1 DZ1)
  //
  //    See also:
  //    planes3d, createPlane
  //
  //    ---------
  //    author: David Legland
  //    INRA - TPV URPOI - BIA IMASTE
  //    created the 21/02/2005.
  //
  //    HISTORY
  //    21/08/2009 compute origin after computation of vectors (more precise)
  //        and add support for several planes.
  //  compute first direction vector
  normalizeVector3d(*(double (*)[3])&plane1[3], d1);

  //  compute second direction vector
  // PLANENORMAL Compute the normal to a plane
  //
  //    N = planeNormal(PLANE)
  //    compute the normal of the given plane
  //    PLANE : [x0 y0 z0 dx1 dy1 dz1 dx2 dy2 dz2]
  //    N : [dx dy dz]
  //
  //    See also
  //    geom3d, planes3d, createPlane
  //
  //    ---------
  //    author : David Legland
  //    INRA - TPV URPOI - BIA IMASTE
  //    created the 17/02/2005.
  //
  //    HISTORY
  //    15/04/2013 Extended to N-dim planes by Sven Holcombe
  //  plane normal
  crossProduct3d(*(double (*)[3])&plane1[3], *(double (*)[3])&plane1[6], normals);
  normalizeVector3d(normals, dv4);
  crossProduct3d(d1, dv4, dv5);
  normalizeVector3d(dv5, d2);

  //  compute origin point of the plane
  for (k = 0; k < 3; k++) {
    plane[k] = plane1[k];
    plane[k + 3] = d1[k];
    plane[k + 6] = -d2[k];
    d2[k] = -d2[k];
  }

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
    dp[k] = normals[k] * plane[k];
  }

  t = dp[0];
  for (k = 0; k < 2; k++) {
    t += dp[k + 1];
  }

  for (k = 0; k < 3; k++) {
    dp[k] = normals[k] * normals[k];
  }

  y = dp[0];
  for (k = 0; k < 2; k++) {
    y += dp[k + 1];
  }

  t /= y;

  //  add relative difference to project point back to plane
  for (k = 0; k < 3; k++) {
    dp[k] = t * normals[k];
  }

  //  create the resulting plane
  for (k = 0; k < 3; k++) {
    plane[k] = dp[k];
    plane[k + 3] = d1[k];
    plane[k + 6] = d2[k];
  }
}

//
// File trailer for createPlane.cpp
//
// [EOF]
//
