/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_generatePackingFromSimDatPLANE_api.h
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 20-Dec-2018 12:28:38
 */

#ifndef _CODER_GENERATEPACKINGFROMSIMDATPLANE_API_H
#define _CODER_GENERATEPACKINGFROMSIMDATPLANE_API_H

/* Include Files */
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include <stddef.h>
#include <stdlib.h>
#include "_coder_generatePackingFromSimDatPLANE_api.h"

/* Type Definitions */
#ifndef struct_emxArray_char_T
#define struct_emxArray_char_T

struct emxArray_char_T
{
  char_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_char_T*/

#ifndef typedef_emxArray_char_T
#define typedef_emxArray_char_T

typedef struct emxArray_char_T emxArray_char_T;

#endif                                 /*typedef_emxArray_char_T*/

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  real_T *data;
  int32_T *size;
  int32_T allocatedSize;
  int32_T numDimensions;
  boolean_T canFreeData;
};

#endif                                 /*struct_emxArray_real_T*/

#ifndef typedef_emxArray_real_T
#define typedef_emxArray_real_T

typedef struct emxArray_real_T emxArray_real_T;

#endif                                 /*typedef_emxArray_real_T*/

#ifndef typedef_struct0_T
#define typedef_struct0_T

typedef struct {
  real_T lw;
  real_T nl;
  real_T npl;
  real_T vib;
  real_T N;
  real_T v;
  emxArray_char_T *fold;
  real_T pars[6];
  emxArray_real_T *ballOut;
  emxArray_real_T *smartInfo;
  emxArray_real_T *t;
  emxArray_real_T *gui;
  emxArray_real_T *totTorque;
  real_T smartVol;
  real_T smartSize[5];
} struct0_T;

#endif                                 /*typedef_struct0_T*/

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern emlrtContext emlrtContextGlobal;

/* Function Declarations */
extern void generatePackingFromSimDatPLANE(struct0_T *dat, real_T countDubs,
  emxArray_real_T *smartCross);
extern void generatePackingFromSimDatPLANE_api(const mxArray * const prhs[2],
  int32_T nlhs, const mxArray *plhs[1]);
extern void generatePackingFromSimDatPLANE_atexit(void);
extern void generatePackingFromSimDatPLANE_initialize(void);
extern void generatePackingFromSimDatPLANE_terminate(void);
extern void generatePackingFromSimDatPLANE_xil_terminate(void);

#endif

/*
 * File trailer for _coder_generatePackingFromSimDatPLANE_api.h
 *
 * [EOF]
 */
