/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: _coder_generatePackingFromSimDatPLANE_mex.cpp
 *
 * MATLAB Coder version            : 4.0
 * C/C++ source code generated on  : 20-Dec-2018 12:28:38
 */

/* Include Files */
#include "_coder_generatePackingFromSimDatPLANE_api.h"
#include "_coder_generatePackingFromSimDatPLANE_mex.h"

/* Function Declarations */
static void c_generatePackingFromSimDatPLAN(int32_T nlhs, mxArray *plhs[1],
  int32_T nrhs, const mxArray *prhs[2]);

/* Function Definitions */

/*
 * Arguments    : int32_T nlhs
 *                mxArray *plhs[1]
 *                int32_T nrhs
 *                const mxArray *prhs[2]
 * Return Type  : void
 */
static void c_generatePackingFromSimDatPLAN(int32_T nlhs, mxArray *plhs[1],
  int32_T nrhs, const mxArray *prhs[2])
{
  const mxArray *outputs[1];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 2, 4,
                        30, "generatePackingFromSimDatPLANE");
  }

  if (nlhs > 1) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 30,
                        "generatePackingFromSimDatPLANE");
  }

  /* Call the function. */
  generatePackingFromSimDatPLANE_api(prhs, nlhs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  generatePackingFromSimDatPLANE_terminate();
}

/*
 * Arguments    : int32_T nlhs
 *                mxArray * const plhs[]
 *                int32_T nrhs
 *                const mxArray * const prhs[]
 * Return Type  : void
 */
void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(generatePackingFromSimDatPLANE_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  generatePackingFromSimDatPLANE_initialize();

  /* Dispatch the entry-point. */
  c_generatePackingFromSimDatPLAN(nlhs, plhs, nrhs, prhs);
}

/*
 * Arguments    : void
 * Return Type  : emlrtCTX
 */
emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/*
 * File trailer for _coder_generatePackingFromSimDatPLANE_mex.cpp
 *
 * [EOF]
 */
