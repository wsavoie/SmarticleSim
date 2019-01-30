//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: generatePackingFromSimDatPLANE_emxutil.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "generatePackingFromSimDatPLANE_emxutil.h"

// Function Definitions

//
// Arguments    : emxArray_real_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_real_T(emxArray_real_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i <<= 1;
      }
    }

    newData = calloc((unsigned int)i, sizeof(double));
    if (emxArray->data != NULL) {
      memcpy(newData, (void *)emxArray->data, sizeof(double) * oldNumel);
      if (emxArray->canFreeData) {
        free((void *)emxArray->data);
      }
    }

    emxArray->data = (double *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

//
// Arguments    : emxArray_real_T *emxArray
//                int oldNumel
// Return Type  : void
//
void emxEnsureCapacity_real_T1(emxArray_real_T *emxArray, int oldNumel)
{
  int newNumel;
  int i;
  void *newData;
  if (oldNumel < 0) {
    oldNumel = 0;
  }

  newNumel = 1;
  for (i = 0; i < emxArray->numDimensions; i++) {
    newNumel *= emxArray->size[i];
  }

  if (newNumel > emxArray->allocatedSize) {
    i = emxArray->allocatedSize;
    if (i < 16) {
      i = 16;
    }

    while (i < newNumel) {
      if (i > 1073741823) {
        i = MAX_int32_T;
      } else {
        i <<= 1;
      }
    }

    newData = calloc((unsigned int)i, sizeof(double));
    if (emxArray->data != NULL) {
      memcpy(newData, (void *)emxArray->data, sizeof(double) * oldNumel);
      if (emxArray->canFreeData) {
        free((void *)emxArray->data);
      }
    }

    emxArray->data = (double *)newData;
    emxArray->allocatedSize = i;
    emxArray->canFreeData = true;
  }
}

//
// Arguments    : struct0_T *pStruct
// Return Type  : void
//
void emxFreeStruct_struct0_T(struct0_T *pStruct)
{
  emxFree_char_T(&pStruct->fold);
  emxFree_real_T(&pStruct->ballOut);
  emxFree_real_T(&pStruct->smartInfo);
  emxFree_real_T(&pStruct->t);
  emxFree_real_T(&pStruct->gui);
  emxFree_real_T(&pStruct->totTorque);
}

//
// Arguments    : emxArray_char_T **pEmxArray
// Return Type  : void
//
void emxFree_char_T(emxArray_char_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_char_T *)NULL) {
    if (((*pEmxArray)->data != (char *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_char_T *)NULL;
  }
}

//
// Arguments    : emxArray_real_T **pEmxArray
// Return Type  : void
//
void emxFree_real_T(emxArray_real_T **pEmxArray)
{
  if (*pEmxArray != (emxArray_real_T *)NULL) {
    if (((*pEmxArray)->data != (double *)NULL) && (*pEmxArray)->canFreeData) {
      free((void *)(*pEmxArray)->data);
    }

    free((void *)(*pEmxArray)->size);
    free((void *)*pEmxArray);
    *pEmxArray = (emxArray_real_T *)NULL;
  }
}

//
// Arguments    : struct0_T *pStruct
// Return Type  : void
//
void emxInitStruct_struct0_T(struct0_T *pStruct)
{
  emxInit_char_T(&pStruct->fold, 2);
  emxInit_real_T(&pStruct->ballOut, 3);
  emxInit_real_T(&pStruct->smartInfo, 3);
  emxInit_real_T1(&pStruct->t, 2);
  emxInit_real_T1(&pStruct->gui, 2);
  emxInit_real_T1(&pStruct->totTorque, 2);
}

//
// Arguments    : emxArray_char_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_char_T(emxArray_char_T **pEmxArray, int numDimensions)
{
  emxArray_char_T *emxArray;
  int i;
  *pEmxArray = (emxArray_char_T *)malloc(sizeof(emxArray_char_T));
  emxArray = *pEmxArray;
  emxArray->data = (char *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_real_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_real_T(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// Arguments    : emxArray_real_T **pEmxArray
//                int numDimensions
// Return Type  : void
//
void emxInit_real_T1(emxArray_real_T **pEmxArray, int numDimensions)
{
  emxArray_real_T *emxArray;
  int i;
  *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
  emxArray = *pEmxArray;
  emxArray->data = (double *)NULL;
  emxArray->numDimensions = numDimensions;
  emxArray->size = (int *)malloc(sizeof(int) * numDimensions);
  emxArray->allocatedSize = 0;
  emxArray->canFreeData = true;
  for (i = 0; i < numDimensions; i++) {
    emxArray->size[i] = 0;
  }
}

//
// File trailer for generatePackingFromSimDatPLANE_emxutil.cpp
//
// [EOF]
//
