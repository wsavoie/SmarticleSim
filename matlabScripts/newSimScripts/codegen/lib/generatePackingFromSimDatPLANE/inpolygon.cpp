//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: inpolygon.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include <math.h>
#include "rt_nonfinite.h"
#include <string.h>
#include "generatePackingFromSimDatPLANE.h"
#include "inpolygon.h"

// Function Declarations
static void b_computeRange(const double xv[4], int nloops, const int first[4],
  const int last[4], double *minxv, double *maxxv);
static void computeRange(const double xv_data[], int nloops, const int
  first_data[], const int last_data[], double *minxv, double *maxxv);
static void contrib(double x1, double b_y1, double x2, double y2, signed char
                    quad1, signed char quad2, double scale, signed char
                    *diffQuad, boolean_T *onj);

// Function Definitions

//
// Arguments    : const double xv[4]
//                int nloops
//                const int first[4]
//                const int last[4]
//                double *minxv
//                double *maxxv
// Return Type  : void
//
static void b_computeRange(const double xv[4], int nloops, const int first[4],
  const int last[4], double *minxv, double *maxxv)
{
  int k;
  int j;
  *minxv = xv[first[0] - 1];
  *maxxv = xv[first[0] - 1];
  for (k = 0; k < nloops; k++) {
    for (j = first[k] - 1; j < last[k]; j++) {
      if (xv[j] < *minxv) {
        *minxv = xv[j];
      } else {
        if (xv[j] > *maxxv) {
          *maxxv = xv[j];
        }
      }
    }
  }
}

//
// Arguments    : const double xv_data[]
//                int nloops
//                const int first_data[]
//                const int last_data[]
//                double *minxv
//                double *maxxv
// Return Type  : void
//
static void computeRange(const double xv_data[], int nloops, const int
  first_data[], const int last_data[], double *minxv, double *maxxv)
{
  int k;
  int j;
  *minxv = xv_data[first_data[0] - 1];
  *maxxv = xv_data[first_data[0] - 1];
  for (k = 0; k < nloops; k++) {
    for (j = first_data[k] - 1; j < last_data[k]; j++) {
      if (xv_data[j] < *minxv) {
        *minxv = xv_data[j];
      } else {
        if (xv_data[j] > *maxxv) {
          *maxxv = xv_data[j];
        }
      }
    }
  }
}

//
// Arguments    : double x1
//                double b_y1
//                double x2
//                double y2
//                signed char quad1
//                signed char quad2
//                double scale
//                signed char *diffQuad
//                boolean_T *onj
// Return Type  : void
//
static void contrib(double x1, double b_y1, double x2, double y2, signed char
                    quad1, signed char quad2, double scale, signed char
                    *diffQuad, boolean_T *onj)
{
  double cp;
  *onj = false;
  *diffQuad = (signed char)(quad2 - quad1);
  cp = x1 * y2 - x2 * b_y1;
  if (fabs(cp) < scale) {
    *onj = (x1 * x2 + b_y1 * y2 <= 0.0);
    if ((*diffQuad == 2) || (*diffQuad == -2)) {
      *diffQuad = 0;
    } else if (*diffQuad == -3) {
      *diffQuad = 1;
    } else {
      if (*diffQuad == 3) {
        *diffQuad = -1;
      }
    }
  } else if (cp < 0.0) {
    if (*diffQuad == 2) {
      *diffQuad = -2;
    } else if (*diffQuad == -3) {
      *diffQuad = 1;
    } else {
      if (*diffQuad == 3) {
        *diffQuad = -1;
      }
    }
  } else if (*diffQuad == -2) {
    *diffQuad = 2;
  } else if (*diffQuad == -3) {
    *diffQuad = 1;
  } else {
    if (*diffQuad == 3) {
      *diffQuad = -1;
    }
  }
}

//
// Arguments    : const double x[3]
//                const double y[3]
//                const double xv[4]
//                const double yv[4]
//                boolean_T in[3]
// Return Type  : void
//
void b_inpolygon(const double x[3], const double y[3], const double xv[4], const
                 double yv[4], boolean_T in[3])
{
  int i;
  int nloops;
  int k;
  int first[4];
  int last[4];
  double minxv;
  double maxxv;
  double minyv;
  double maxyv;
  boolean_T exitg1;
  int j;
  double scale[4];
  double avxi;
  double avyi;
  boolean_T inj;
  signed char sdq;
  int exitg3;
  double xvFirst;
  double yvFirst;
  signed char quadFirst;
  double xv2;
  double yv2;
  signed char quad2;
  int exitg2;
  signed char dquad;
  double xv1;
  double yv1;
  signed char quad1;
  for (i = 0; i < 3; i++) {
    in[i] = false;
  }

  nloops = 0;
  for (i = 0; i < 4; i++) {
    first[i] = 0;
    last[i] = 0;
  }

  k = 0;
  while ((k + 1 <= 4) && rtIsNaN(xv[k])) {
    k++;
  }

  while (k + 1 <= 4) {
    nloops++;
    i = k;
    first[nloops - 1] = k + 1;
    exitg1 = false;
    while ((!exitg1) && (k + 1 < 4)) {
      k++;
      if (rtIsNaN(xv[k]) || rtIsNaN(yv[k])) {
        k--;
        exitg1 = true;
      }
    }

    if ((xv[k] == xv[i]) && (yv[k] == yv[i])) {
      last[nloops - 1] = k;
    } else {
      last[nloops - 1] = k + 1;
    }

    k += 2;
    while ((k + 1 <= 4) && rtIsNaN(xv[k])) {
      k++;
    }
  }

  if (nloops != 0) {
    b_computeRange(xv, nloops, first, last, &minxv, &maxxv);
    b_computeRange(yv, nloops, first, last, &minyv, &maxyv);
    for (i = 0; i < 4; i++) {
      scale[i] = 0.0;
    }

    for (j = 0; j < nloops; j++) {
      for (i = first[j]; i < last[j]; i++) {
        avxi = fabs(0.5 * (xv[i - 1] + xv[i]));
        avyi = fabs(0.5 * (yv[i - 1] + yv[i]));
        if ((avxi > 1.0) && (avyi > 1.0)) {
          avxi *= avyi;
        } else {
          if ((avyi > avxi) || rtIsNaN(avxi)) {
            avxi = avyi;
          }
        }

        scale[i - 1] = avxi * 6.6613381477509392E-16;
      }

      avxi = fabs(0.5 * (xv[last[j] - 1] + xv[first[j] - 1]));
      avyi = fabs(0.5 * (yv[last[j] - 1] + yv[first[j] - 1]));
      if ((avxi > 1.0) && (avyi > 1.0)) {
        avxi *= avyi;
      } else {
        if ((avyi > avxi) || rtIsNaN(avxi)) {
          avxi = avyi;
        }
      }

      scale[last[j] - 1] = avxi * 6.6613381477509392E-16;
    }

    for (j = 0; j < 3; j++) {
      avxi = x[j];
      avyi = y[j];
      inj = false;
      if ((avxi >= minxv) && (avxi <= maxxv) && (avyi >= minyv) && (avyi <=
           maxyv)) {
        sdq = 0;
        k = 0;
        do {
          exitg3 = 0;
          if (k + 1 <= nloops) {
            xvFirst = xv[first[k] - 1] - avxi;
            yvFirst = yv[first[k] - 1] - avyi;
            if (xvFirst > 0.0) {
              if (yvFirst > 0.0) {
                quadFirst = 0;
              } else {
                quadFirst = 3;
              }
            } else if (yvFirst > 0.0) {
              quadFirst = 1;
            } else {
              quadFirst = 2;
            }

            xv2 = xvFirst;
            yv2 = yvFirst;
            quad2 = quadFirst;
            i = first[k];
            do {
              exitg2 = 0;
              if (i <= last[k] - 1) {
                xv1 = xv2;
                yv1 = yv2;
                xv2 = xv[i] - avxi;
                yv2 = yv[i] - avyi;
                quad1 = quad2;
                if (xv2 > 0.0) {
                  if (yv2 > 0.0) {
                    quad2 = 0;
                  } else {
                    quad2 = 3;
                  }
                } else if (yv2 > 0.0) {
                  quad2 = 1;
                } else {
                  quad2 = 2;
                }

                contrib(xv1, yv1, xv2, yv2, quad1, quad2, scale[i - 1], &dquad,
                        &inj);
                if (inj) {
                  inj = true;
                  exitg2 = 1;
                } else {
                  sdq += dquad;
                  i++;
                }
              } else {
                contrib(xv2, yv2, xvFirst, yvFirst, quad2, quadFirst,
                        scale[last[k] - 1], &dquad, &inj);
                exitg2 = 2;
              }
            } while (exitg2 == 0);

            if (exitg2 == 1) {
              exitg3 = 1;
            } else if (inj) {
              inj = true;
              exitg3 = 1;
            } else {
              sdq += dquad;
              k++;
            }
          } else {
            inj = (sdq != 0);
            exitg3 = 1;
          }
        } while (exitg3 == 0);
      }

      in[j] = inj;
    }
  }
}

//
// Arguments    : const double x[3]
//                const double y[3]
//                const double xv_data[]
//                const int xv_size[1]
//                const double yv_data[]
//                boolean_T in[3]
// Return Type  : void
//
void inpolygon(const double x[3], const double y[3], const double xv_data[],
               const int xv_size[1], const double yv_data[], boolean_T in[3])
{
  int i;
  int n;
  int nloops;
  int first_data[4];
  int last_data[4];
  int k;
  double minxv;
  double maxxv;
  double minyv;
  double maxyv;
  boolean_T exitg1;
  double scale_data[4];
  double avxi;
  double avyi;
  boolean_T inj;
  long long sdq;
  int exitg3;
  double xvFirst;
  double yvFirst;
  signed char quadFirst;
  double xv2;
  double yv2;
  signed char quad2;
  int exitg2;
  signed char dquad;
  double xv1;
  double yv1;
  signed char quad1;
  for (i = 0; i < 3; i++) {
    in[i] = false;
  }

  if (xv_size[0] != 0) {
    n = xv_size[0];
    nloops = 0;
    i = xv_size[0];
    if (0 <= i - 1) {
      memset(&first_data[0], 0, (unsigned int)(i * (int)sizeof(int)));
    }

    i = xv_size[0];
    if (0 <= i - 1) {
      memset(&last_data[0], 0, (unsigned int)(i * (int)sizeof(int)));
    }

    k = 0;
    while ((k + 1 <= n) && rtIsNaN(xv_data[k])) {
      k++;
    }

    while (k + 1 <= n) {
      nloops++;
      i = k;
      first_data[nloops - 1] = k + 1;
      exitg1 = false;
      while ((!exitg1) && (k + 1 < n)) {
        k++;
        if (rtIsNaN(xv_data[k]) || rtIsNaN(yv_data[k])) {
          k--;
          exitg1 = true;
        }
      }

      if ((xv_data[k] == xv_data[i]) && (yv_data[k] == yv_data[i])) {
        last_data[nloops - 1] = k;
      } else {
        last_data[nloops - 1] = k + 1;
      }

      k += 2;
      while ((k + 1 <= n) && rtIsNaN(xv_data[k])) {
        k++;
      }
    }

    if (nloops != 0) {
      computeRange(xv_data, nloops, first_data, last_data, &minxv, &maxxv);
      computeRange(yv_data, nloops, first_data, last_data, &minyv, &maxyv);
      i = xv_size[0];
      if (0 <= i - 1) {
        memset(&scale_data[0], 0, (unsigned int)(i * (int)sizeof(double)));
      }

      for (n = 0; n < nloops; n++) {
        for (i = first_data[n]; i < last_data[n]; i++) {
          avxi = fabs(0.5 * (xv_data[i - 1] + xv_data[i]));
          avyi = fabs(0.5 * (yv_data[i - 1] + yv_data[i]));
          if ((avxi > 1.0) && (avyi > 1.0)) {
            avxi *= avyi;
          } else {
            if ((avyi > avxi) || rtIsNaN(avxi)) {
              avxi = avyi;
            }
          }

          scale_data[i - 1] = avxi * 6.6613381477509392E-16;
        }

        avxi = fabs(0.5 * (xv_data[last_data[n] - 1] + xv_data[first_data[n] - 1]));
        avyi = fabs(0.5 * (yv_data[last_data[n] - 1] + yv_data[first_data[n] - 1]));
        if ((avxi > 1.0) && (avyi > 1.0)) {
          avxi *= avyi;
        } else {
          if ((avyi > avxi) || rtIsNaN(avxi)) {
            avxi = avyi;
          }
        }

        scale_data[last_data[n] - 1] = avxi * 6.6613381477509392E-16;
      }

      for (n = 0; n < 3; n++) {
        avxi = x[n];
        avyi = y[n];
        inj = false;
        if ((avxi >= minxv) && (avxi <= maxxv) && (avyi >= minyv) && (avyi <=
             maxyv)) {
          sdq = 0LL;
          k = 0;
          do {
            exitg3 = 0;
            if (k + 1 <= nloops) {
              xvFirst = xv_data[first_data[k] - 1] - avxi;
              yvFirst = yv_data[first_data[k] - 1] - avyi;
              if (xvFirst > 0.0) {
                if (yvFirst > 0.0) {
                  quadFirst = 0;
                } else {
                  quadFirst = 3;
                }
              } else if (yvFirst > 0.0) {
                quadFirst = 1;
              } else {
                quadFirst = 2;
              }

              xv2 = xvFirst;
              yv2 = yvFirst;
              quad2 = quadFirst;
              i = first_data[k];
              do {
                exitg2 = 0;
                if (i <= last_data[k] - 1) {
                  xv1 = xv2;
                  yv1 = yv2;
                  xv2 = xv_data[i] - avxi;
                  yv2 = yv_data[i] - avyi;
                  quad1 = quad2;
                  if (xv2 > 0.0) {
                    if (yv2 > 0.0) {
                      quad2 = 0;
                    } else {
                      quad2 = 3;
                    }
                  } else if (yv2 > 0.0) {
                    quad2 = 1;
                  } else {
                    quad2 = 2;
                  }

                  contrib(xv1, yv1, xv2, yv2, quad1, quad2, scale_data[i - 1],
                          &dquad, &inj);
                  if (inj) {
                    inj = true;
                    exitg2 = 1;
                  } else {
                    sdq += dquad;
                    i++;
                  }
                } else {
                  contrib(xv2, yv2, xvFirst, yvFirst, quad2, quadFirst,
                          scale_data[last_data[k] - 1], &dquad, &inj);
                  exitg2 = 2;
                }
              } while (exitg2 == 0);

              if (exitg2 == 1) {
                exitg3 = 1;
              } else if (inj) {
                inj = true;
                exitg3 = 1;
              } else {
                sdq += dquad;
                k++;
              }
            } else {
              inj = (sdq != 0LL);
              exitg3 = 1;
            }
          } while (exitg3 == 0);
        }

        in[n] = inj;
      }
    }
  }
}

//
// File trailer for inpolygon.cpp
//
// [EOF]
//
