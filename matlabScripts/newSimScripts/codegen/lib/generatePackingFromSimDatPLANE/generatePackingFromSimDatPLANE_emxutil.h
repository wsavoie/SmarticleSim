//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: generatePackingFromSimDatPLANE_emxutil.h
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//
#ifndef GENERATEPACKINGFROMSIMDATPLANE_EMXUTIL_H
#define GENERATEPACKINGFROMSIMDATPLANE_EMXUTIL_H

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "omp.h"
#include "generatePackingFromSimDatPLANE_types.h"

// Function Declarations
extern void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel);
extern void emxEnsureCapacity_real_T1(emxArray_real_T *emxArray, int oldNumel);
extern void emxFreeStruct_struct0_T(struct0_T *pStruct);
extern void emxFree_char_T(emxArray_char_T **pEmxArray);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInitStruct_struct0_T(struct0_T *pStruct);
extern void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions);
extern void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions);

#endif

//
// File trailer for generatePackingFromSimDatPLANE_emxutil.h
//
// [EOF]
//
