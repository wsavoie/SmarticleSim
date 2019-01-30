//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: generatePackingFromSimDatPLANE_initialize.cpp
//
// MATLAB Coder version            : 4.0
// C/C++ source code generated on  : 20-Dec-2018 12:28:38
//

// Include Files
#include "rt_nonfinite.h"
#include "generatePackingFromSimDatPLANE.h"
#include "generatePackingFromSimDatPLANE_initialize.h"
#include "generatePackingFromSimDatPLANE_data.h"

// Function Definitions

//
// Arguments    : void
// Return Type  : void
//
void generatePackingFromSimDatPLANE_initialize()
{
  rt_InitInfAndNaN(8U);
  omp_init_nest_lock(&emlrtNestLockGlobal);
}

//
// File trailer for generatePackingFromSimDatPLANE_initialize.cpp
//
// [EOF]
//
