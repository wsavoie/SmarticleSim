//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: generatePackingFromSimDatPLANE_types.h
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//
#ifndef GENERATEPACKINGFROMSIMDATPLANE_TYPES_H
#define GENERATEPACKINGFROMSIMDATPLANE_TYPES_H

// Include Files
#include "rtwtypes.h"

// Type Definitions
#ifndef struct_emxArray_real_T_4x2
#define struct_emxArray_real_T_4x2

struct emxArray_real_T_4x2
{
  double data[8];
  int size[2];
};

#endif                                 //struct_emxArray_real_T_4x2

#ifndef struct_sN4Iv1dHUk7K5DirCqGMRGD_tag
#define struct_sN4Iv1dHUk7K5DirCqGMRGD_tag

struct sN4Iv1dHUk7K5DirCqGMRGD_tag
{
  emxArray_real_T_4x2 f1;
};

#endif                                 //struct_sN4Iv1dHUk7K5DirCqGMRGD_tag

typedef sN4Iv1dHUk7K5DirCqGMRGD_tag cell_wrap_3;

#ifndef struct_c_emxArray_sN4Iv1dHUk7K5DirCqGM
#define struct_c_emxArray_sN4Iv1dHUk7K5DirCqGM

struct c_emxArray_sN4Iv1dHUk7K5DirCqGM
{
  cell_wrap_3 data[5];
  int size[1];
};

#endif                                 //struct_c_emxArray_sN4Iv1dHUk7K5DirCqGM

typedef c_emxArray_sN4Iv1dHUk7K5DirCqGM emxArray_cell_wrap_3_5;

#ifndef struct_emxArray_char_T
#define struct_emxArray_char_T

struct emxArray_char_T
{
  char *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_char_T

#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T

struct emxArray_real_T
{
  double *data;
  int *size;
  int allocatedSize;
  int numDimensions;
  boolean_T canFreeData;
};

#endif                                 //struct_emxArray_real_T

typedef struct {
  double lw;
  double nl;
  double npl;
  double vib;
  double N;
  double v;
  emxArray_char_T *fold;
  double pars[6];
  emxArray_real_T *ballOut;
  emxArray_real_T *smartInfo;
  emxArray_real_T *t;
  emxArray_real_T *gui;
  emxArray_real_T *totTorque;
  double smartVol;
  double smartSize[5];
} struct0_T;

#endif

//
// File trailer for generatePackingFromSimDatPLANE_types.h
//
// [EOF]
//
