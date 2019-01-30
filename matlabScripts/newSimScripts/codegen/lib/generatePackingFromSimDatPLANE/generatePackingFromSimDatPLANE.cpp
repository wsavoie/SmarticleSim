//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: generatePackingFromSimDatPLANE.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include <math.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "quat2rotm.h"
#include "distancePointEdge3d.h"
#include "intersectRayPolygon3d.h"
#include "pdist.h"
#include "generatePackingFromSimDatPLANE_emxutil.h"

// Function Definitions

//
// [smartCross]=GENERATEPACKINGFROMSIMDAT(dat,startInd,countDubs)
//
// dat is a 1x1 struct (single version
// countDubs= 1 or 0, flag to count double plane crossings as entanglements
// Arguments    : const struct0_T *dat
//                double countDubs
//                emxArray_real_T *smartCross
// Return Type  : void
//
void generatePackingFromSimDatPLANE(const struct0_T *dat, double countDubs,
  emxArray_real_T *smartCross)
{
  int loop_ub;
  emxArray_real_T *c;
  int b_loop_ub;
  int i0;
  int i1;
  double t2;
  int ii;
  double l;
  double w;
  int i;
  emxArray_real_T *d;
  emxArray_real_T *l2r;
  emxArray_real_T *sp;
  int jj;
  int idx;
  int j;
  double distz_data[3];
  double count;
  double b_l2r[4];
  double L1[9];
  double b_sp[6];
  int k;
  double b_w[3];
  double c_w[3];
  double d_w[3];
  double e_w[3];
  double f_w[3];
  double d0;
  double g_w[3];
  double h_w[3];
  double b_L1[18];
  double iPt[9];
  boolean_T inside[3];
  double L2[9];
  boolean_T y;
  boolean_T exitg1;
  boolean_T guard1 = false;
  int i_data[3];
  int r_data[3];
  int iPt_size[2];
  int L1_size[2];
  double iPt_data[9];
  double L1_data[18];
  int distz_size[1];
  boolean_T x_data[3];

  //  if(PLOTZ)
  //      figure(129);
  //      hold on;
  //  end
  //  startInd=5;
  //  vInd=1;
  // represents ~85 degs
  // %%%
  //  warning('use single struct ind as input param');
  //  c=dat(vInd).smartInfo(:,:,startInd:end); %[x,y,z,e0,e1,e2,e3,ang1,ang2]
  if (1 > dat->smartInfo->size[2]) {
    loop_ub = 0;
  } else {
    loop_ub = dat->smartInfo->size[2];
  }

  emxInit_real_T(&c, 3);
  b_loop_ub = dat->smartInfo->size[0];
  i0 = c->size[0] * c->size[1] * c->size[2];
  c->size[0] = b_loop_ub;
  c->size[1] = 9;
  c->size[2] = loop_ub;
  emxEnsureCapacity_real_T(c, i0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    for (i1 = 0; i1 < 9; i1++) {
      for (ii = 0; ii < b_loop_ub; ii++) {
        c->data[(ii + c->size[0] * i1) + c->size[0] * c->size[1] * i0] =
          dat->smartInfo->data[(ii + dat->smartInfo->size[0] * i1) +
          dat->smartInfo->size[0] * dat->smartInfo->size[1] * i0];
      }
    }
  }

  // [x,y,z,e0,e1,e2,e3,ang1,ang2]
  //  [rho,t1,t2,l,w]=separateVec(dat(1).smartSize,1);
  t2 = dat->smartSize[2];
  l = dat->smartSize[3];
  w = dat->smartSize[4];
  i0 = dat->smartInfo->size[0];
  i1 = smartCross->size[0] * smartCross->size[1];
  smartCross->size[0] = loop_ub;
  smartCross->size[1] = i0;
  emxEnsureCapacity_real_T1(smartCross, i1);
  b_loop_ub = loop_ub * i0;
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    smartCross->data[i0] = 0.0;
  }

  i = 0;
  emxInit_real_T1(&d, 2);
  emxInit_real_T1(&l2r, 2);
  emxInit_real_T(&sp, 3);
  while (i <= loop_ub - 1) {
    b_loop_ub = dat->smartInfo->size[0];
    i0 = d->size[0] * d->size[1];
    d->size[0] = b_loop_ub;
    d->size[1] = 9;
    emxEnsureCapacity_real_T1(d, i0);
    for (i0 = 0; i0 < 9; i0++) {
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        d->data[i1 + d->size[0] * i0] = c->data[(i1 + c->size[0] * i0) + c->
          size[0] * c->size[1] * i];
      }
    }

    // data for current frame
    // final timestep is d
    //      [e0,e1,e2,e3]=separateVec(d(:,4:7),1); %get quaternion data
    b_loop_ub = dat->smartInfo->size[0];
    ii = dat->smartInfo->size[0];
    jj = dat->smartInfo->size[0];
    idx = dat->smartInfo->size[0];
    i0 = l2r->size[0] * l2r->size[1];
    l2r->size[0] = b_loop_ub;
    l2r->size[1] = 4;
    emxEnsureCapacity_real_T1(l2r, i0);
    for (i0 = 0; i0 < b_loop_ub; i0++) {
      l2r->data[i0] = d->data[i0 + d->size[0] * 5];
    }

    for (i0 = 0; i0 < ii; i0++) {
      l2r->data[i0 + l2r->size[0]] = -d->data[i0 + d->size[0] * 6];
    }

    for (i0 = 0; i0 < jj; i0++) {
      l2r->data[i0 + (l2r->size[0] << 1)] = -d->data[i0 + d->size[0] * 3];
    }

    for (i0 = 0; i0 < idx; i0++) {
      l2r->data[i0 + l2r->size[0] * 3] = -d->data[i0 + (d->size[0] << 2)];
    }

    // w,x,y,z->-y,z,w,x
    //      q = [d(:,4:7)];
    //      q=quaternion(l2r); %I need to transpose it before setting to quat for matrix 
    //      angles = EulerAngles(q,'XYZ');
    //      angles2 = quat2eul(l2r);
    //      angles=reshape(angles(3,1,:),[size(angles,3),1]);'
    i0 = dat->smartInfo->size[0];
    i1 = sp->size[0] * sp->size[1] * sp->size[2];
    sp->size[0] = 4;
    sp->size[1] = 3;
    sp->size[2] = i0;
    emxEnsureCapacity_real_T(sp, i1);
    b_loop_ub = 12 * i0;
    for (i0 = 0; i0 < b_loop_ub; i0++) {
      sp->data[i0] = 0.0;
    }

    i0 = dat->smartInfo->size[0];
    for (j = 0; j < i0; j++) {
      for (i1 = 0; i1 < 3; i1++) {
        distz_data[i1] = d->data[j + d->size[0] * i1];
      }

      distz_data[0] = -d->data[j];
      count = l + t2;

      //          RR=RotationMatrix(q(j));
      for (i1 = 0; i1 < 4; i1++) {
        b_l2r[i1] = l2r->data[j + l2r->size[0] * i1];
      }

      quat2rotm(b_l2r, L1);

      // scatter3(smPos(1),smPos(2),smPos(3),'o');
      // v1   v4
      // |    |
      // v2-0-v3
      //          if(PLOTZ)
      //              patch(ptz(:,1),ptz(:,2),ptz(:,3),[1,0,0]);
      //          end
      //          patch(ptz1(:,1),ptz1(:,2),ptz1(:,3),[1,0,0]);
      //          patch(ptz2(:,1),ptz2(:,2),ptz2(:,3),[1,0,0]);
      b_w[0] = -w / 2.0 - count / 2.0 * cos(d->data[j + d->size[0] * 7]);
      b_w[1] = 0.0;
      b_w[2] = -(l - t2) * sin(d->data[j + d->size[0] * 7]);
      c_w[0] = -w / 2.0;
      c_w[1] = 0.0;
      c_w[2] = 0.0;
      d_w[0] = w / 2.0;
      d_w[1] = 0.0;
      d_w[2] = 0.0;
      e_w[0] = w / 2.0 + count / 2.0 * cos(d->data[j + (d->size[0] << 3)]);
      e_w[1] = 0.0;
      e_w[2] = -(l - t2) * sin(d->data[j + (d->size[0] << 3)]);
      for (i1 = 0; i1 < 3; i1++) {
        count = 0.0;
        for (ii = 0; ii < 3; ii++) {
          count += b_w[ii] * L1[ii + 3 * i1];
        }

        d0 = 0.0;
        for (ii = 0; ii < 3; ii++) {
          d0 += c_w[ii] * L1[ii + 3 * i1];
        }

        f_w[i1] = d0 + distz_data[i1];
        d0 = 0.0;
        for (ii = 0; ii < 3; ii++) {
          d0 += d_w[ii] * L1[ii + 3 * i1];
        }

        g_w[i1] = d0 + distz_data[i1];
        d0 = 0.0;
        for (ii = 0; ii < 3; ii++) {
          d0 += e_w[ii] * L1[ii + 3 * i1];
        }

        h_w[i1] = d0 + distz_data[i1];
        sp->data[sp->size[0] * i1 + sp->size[0] * sp->size[1] * j] = count +
          distz_data[i1];
      }

      for (i1 = 0; i1 < 3; i1++) {
        sp->data[(sp->size[0] * i1 + sp->size[0] * sp->size[1] * j) + 1] =
          f_w[i1];
      }

      for (i1 = 0; i1 < 3; i1++) {
        sp->data[(sp->size[0] * i1 + sp->size[0] * sp->size[1] * j) + 2] =
          g_w[i1];
      }

      for (i1 = 0; i1 < 3; i1++) {
        sp->data[(sp->size[0] * i1 + sp->size[0] * sp->size[1] * j) + 3] =
          h_w[i1];
      }

      //
      //          drawEdge3d(createEdge3d(v2,v3));
      //          drawEdge3d(createEdge3d(v1,v2));
      //          drawEdge3d(createEdge3d(v4,v3));
      //          axis equal
    }

    i0 = dat->smartInfo->size[0];
    for (j = 0; j < i0; j++) {
      //          if(PLOTZ)
      //              pH=drawPolygon3d(polyg,'color',[.2,.5,.3],'linewidth',2);
      //          end
      //
      count = 0.0;

      // filter automatically when ang is nearly straight <5deg
      // sin 85deg=.9962
      for (i1 = 0; i1 < 3; i1++) {
        b_sp[i1 << 1] = sp->data[sp->size[0] * i1 + sp->size[0] * sp->size[1] *
          j];
        b_sp[1 + (i1 << 1)] = sp->data[(sp->size[0] * i1 + sp->size[0] *
          sp->size[1] * j) + 3];
      }

      if (!(pdist(b_sp) > 2.0 * l * 0.99 + w)) {
        i1 = dat->smartInfo->size[0];
        for (k = 0; k < i1; k++) {
          if (1 + k != 1 + j) {
            // CREATELINE3D Create a line with various inputs.
            //
            //    Line is represented in a parametric form : [x0 y0 z0 dx dy dz] 
            //        x = x0 + t*dx
            //        y = y0 + t*dy;
            //        z = z0 + t*dz;
            //
            //
            //    L = createLine3d(P1, P2);
            //    Returns the line going through the two given points P1 and P2. 
            //
            //    L = createLine3d(X0, Y0, Z0, DX, DY, DZ);
            //    Returns the line going through the point [x0, y0, z0], and with 
            //    direction vector given by [DX DY DZ].
            //
            //    L = createLine3d(P0, DX, DY, DZ);
            //    Returns the line going through point P0 given by [x0, y0, z0] and with 
            //    direction vector given by [DX DY DZ].
            //
            //    L = createLine3d(THETA, PHI);
            //    Create a line originated at (0,0) and with angles theta and phi. 
            //
            //    L = createLine3d(P0, THETA, PHI);
            //    Create a line with direction given by theta and phi, and which contains 
            //    point P0.
            //
            //
            //    Note : in all cases, parameters can be vertical arrays of the same 
            //    dimension. The result is then an array of lines, of dimensions [N*6]. 
            //
            //    See also:
            //    lines3d
            //
            //    ---------
            //
            //    author : David Legland
            //    INRA - TPV URPOI - BIA IMASTE
            //    created the 17/02/2005.
            //
            //    HISTORY
            //    30/11/2005 add more cases
            //    04/01/2007 remove unused variables
            //    NOTE : A 3d line can also be represented with a 1*7 array :
            //    [x0 y0 z0 dx dy dz t].
            //    whith 't' being one of the following :
            //    - t=0 : line is a singleton (x0,y0)
            //    - t=1 : line is an edge segment, between points (x0,y0) and (x0+dx, 
            //    y0+dy).
            //    - t=Inf : line is a Ray, originated from (x0,y0) and going to infinity 
            //    in the direction(dx,dy).
            //    - t=-Inf : line is a Ray, originated from (x0,y0) and going to infinity 
            //    in the direction(-dx,-dy).
            //    - t=NaN : line is a real straight line, and contains all points 
            //    verifying the above equation.
            //    This seems to be a convenient way to represent uniformly all kind of 
            //    lines (including edges, rays, and even point).
            //
            //    NOTE2 : Another form is the 1*8 array :
            //    [x0 y0 z0 dx dy dz t1 t2].
            //    with t1 and t2 giving first and last position of the edge on the 
            //    supporting straight line, and with t2>t1.
            //  2 input parameters. They can be :
            //  - 2 points, then 2 arrays of 1*2 double.
            //  first input parameter is first point, and second input is the
            //  second point.
            // CREATEEDGE3D Create an edge between two 3D points, or from a 3D line 
            //
            //    E = createEdge3d(P1, P2)
            //    Creates the 3D edge joining the two points P1 and P2.
            //
            //    E = createEdge3d(LIN)
            //    Creates the 3D edge with same origin and same direction vector as the 
            //    3D line LIN.
            //
            //    Example
            //      p1 = [1 1 1];
            //      p2 = [3 4 5];
            //      edge = createEdge3d(p1, p2);
            //      edgeLength3d(edge)
            //      ans =
            //          5.3852
            //
            //    See also
            //      edges3d, drawEdge3d, clipEdge3d, edgelength3d
            //  ------
            //  Author: David Legland
            //  e-mail: david.legland@inra.fr
            //  Created: 2018-08-29,    using Matlab 9.4.0.813654 (R2018a)
            //  Copyright 2018 INRA - Cepia Software Platform.
            //  2 input parameters correspond to two 3D points
            //  extract the two arguments
            //  first input parameter is first point, and second input is the
            //  second point. Allows multiple points.
            //                  l1=createLine3d(armPts(2,:),armPts(1,:)); %left barb 
            //                  l2=createLine3d(armPts(2,:),armPts(3,:)); %spine 
            //                  l3=createLine3d(armPts(3,:),armPts(4,:)); %right barb 
            //
            //                  if(PLOTZ)
            //                      h=drawEdge3d(E,'linewidth',2);
            //                      g=drawLine3d(L,'linewidth',1,'color',[.8,.2,.3]); 
            //                  end
            // intersectionPt
            for (ii = 0; ii < 3; ii++) {
              L1[3 * ii] = sp->data[(sp->size[0] * ii + sp->size[0] * sp->size[1]
                * k) + 1];
              L1[1 + 3 * ii] = sp->data[(sp->size[0] * ii + sp->size[0] *
                sp->size[1] * k) + 1];
              L1[2 + 3 * ii] = sp->data[(sp->size[0] * ii + sp->size[0] *
                sp->size[1] * k) + 2];
              L2[3 * ii] = sp->data[sp->size[0] * ii + sp->size[0] * sp->size[1]
                * k];
              L2[1 + 3 * ii] = sp->data[(sp->size[0] * ii + sp->size[0] *
                sp->size[1] * k) + 2];
              L2[2 + 3 * ii] = sp->data[(sp->size[0] * ii + sp->size[0] *
                sp->size[1] * k) + 3];
              b_sp[ii << 1] = sp->data[sp->size[0] * ii + sp->size[0] * sp->
                size[1] * k];
              b_sp[1 + (ii << 1)] = sp->data[(sp->size[0] * ii + sp->size[0] *
                sp->size[1] * k) + 3];
            }

            if (!(pdist(b_sp) > 2.0 * l * 0.99 + w)) {
              for (ii = 0; ii < 3; ii++) {
                b_L1[ii] = L1[ii];
                b_L1[3 + ii] = L1[3 + ii];
                b_L1[6 + ii] = L1[6 + ii];
                b_L1[9 + ii] = L2[ii] - L1[ii];
                b_L1[12 + ii] = L2[3 + ii] - L1[3 + ii];
                b_L1[15 + ii] = L2[6 + ii] - L1[6 + ii];
              }

              intersectRayPolygon3d(b_L1, *(double (*)[12])&sp->data[sp->size[0]
                                    * sp->size[1] * j], iPt, inside);

              // straight shapes are not accurate
              y = false;
              ii = 0;
              exitg1 = false;
              while ((!exitg1) && (ii < 3)) {
                if (inside[ii]) {
                  y = true;
                  exitg1 = true;
                } else {
                  ii++;
                }
              }

              if (y) {
                //                      if(sum(inside)>1)
                //                          warning('multiple inside')
                //                          pH=drawPolygon3d(polyg,'color',[.2,.5,.3],'linewidth',2); 
                //                          h=drawEdge3d(E,'linewidth',2);
                //                          g=drawLine3d(L,'linewidth',1,'color',[.8,.2,.3]); 
                //                      end
                idx = 0;
                ii = 1;
                jj = 1;
                exitg1 = false;
                while ((!exitg1) && (jj <= 1)) {
                  guard1 = false;
                  if (inside[ii - 1]) {
                    idx++;
                    i_data[idx - 1] = ii;
                    if (idx >= 3) {
                      exitg1 = true;
                    } else {
                      guard1 = true;
                    }
                  } else {
                    guard1 = true;
                  }

                  if (guard1) {
                    ii++;
                    if (ii > 3) {
                      ii = 1;
                      jj = 2;
                    }
                  }
                }

                if (1 > idx) {
                  jj = 0;
                  b_loop_ub = 0;
                } else {
                  jj = idx;
                  b_loop_ub = idx;
                }

                if (0 <= b_loop_ub - 1) {
                  memcpy(&r_data[0], &i_data[0], (unsigned int)(b_loop_ub * (int)
                          sizeof(int)));
                }

                // point may on actual edge only on line
                // so check that it exists on the actual edge
                if (((jj > 1) && (countDubs != 0.0)) || (jj == 1)) {
                  iPt_size[0] = jj;
                  iPt_size[1] = 3;
                  for (ii = 0; ii < 3; ii++) {
                    for (idx = 0; idx < jj; idx++) {
                      iPt_data[idx + jj * ii] = iPt[(r_data[idx] + 3 * ii) - 1];
                    }

                    for (idx = 0; idx < 3; idx++) {
                      b_L1[idx + 3 * ii] = L1[idx + 3 * ii];
                      b_L1[idx + 3 * (ii + 3)] = L2[idx + 3 * ii];
                    }
                  }

                  L1_size[0] = jj;
                  L1_size[1] = 6;
                  for (ii = 0; ii < 6; ii++) {
                    for (idx = 0; idx < jj; idx++) {
                      L1_data[idx + jj * ii] = b_L1[(r_data[idx] + 3 * ii) - 1];
                    }
                  }

                  distancePointEdge3d(iPt_data, iPt_size, L1_data, L1_size,
                                      distz_data, distz_size);
                  b_loop_ub = distz_size[0];
                  for (ii = 0; ii < b_loop_ub; ii++) {
                    x_data[ii] = (distz_data[ii] < 1.0E-8);
                  }

                  y = !(distz_size[0] == 0);
                  if (y) {
                    ii = 1;
                    exitg1 = false;
                    while ((!exitg1) && (ii <= distz_size[0])) {
                      if (!x_data[ii - 1]) {
                        y = false;
                        exitg1 = true;
                      } else {
                        ii++;
                      }
                    }
                  }

                  if (y) {
                    //                              if(PLOTZ)
                    //                                  scatter3(iPt(r,1),iPt(r,2),iPt(r,3),60*ones(size(ptArr,1),1),'filled') 
                    //                              end
                    count += (double)jj;

                    //                                      count=count+1;
                  }
                }
              }
            }

            //                  if(PLOTZ)
            //                      delete(h);
            //                      delete(g);
            //                  end
          }
        }
      }

      smartCross->data[i + smartCross->size[0] * j] = count;

      //          if(PLOTZ)
      //              delete(pH);
      //          end
    }

    i++;
  }

  emxFree_real_T(&sp);
  emxFree_real_T(&l2r);
  emxFree_real_T(&d);
  emxFree_real_T(&c);

  //  beep;
}

//
// File trailer for generatePackingFromSimDatPLANE.cpp
//
// [EOF]
//
