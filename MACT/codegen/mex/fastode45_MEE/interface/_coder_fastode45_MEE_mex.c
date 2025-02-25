/*
 * _coder_fastode45_MEE_mex.c
 *
 * Code generation for function '_coder_fastode45_MEE_mex'
 *
 */

/* Include files */
#include "fastode45_MEE.h"
#include "_coder_fastode45_MEE_mex.h"
#include "fastode45_MEE_terminate.h"
#include "_coder_fastode45_MEE_api.h"
#include "fastode45_MEE_initialize.h"
#include "fastode45_MEE_data.h"

/* Function Declarations */
static void fastode45_MEE_mexFunction(int32_T nlhs, mxArray *plhs[2], int32_T
  nrhs, const mxArray *prhs[6]);

/* Function Definitions */
static void fastode45_MEE_mexFunction(int32_T nlhs, mxArray *plhs[2], int32_T
  nrhs, const mxArray *prhs[6])
{
  int32_T n;
  const mxArray *inputs[6];
  const mxArray *outputs[2];
  int32_T b_nlhs;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;

  /* Check for proper number of arguments. */
  if (nrhs != 6) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:WrongNumberOfInputs", 5, 12, 6, 4,
                        13, "fastode45_MEE");
  }

  if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(&st, "EMLRT:runTime:TooManyOutputArguments", 3, 4, 13,
                        "fastode45_MEE");
  }

  /* Temporary copy for mex inputs. */
  for (n = 0; n < nrhs; n++) {
    inputs[n] = prhs[n];
    if (*emlrtBreakCheckR2012bFlagVar != 0) {
      emlrtBreakCheckR2012b(&st);
    }
  }

  /* Call the function. */
  fastode45_MEE_api(inputs, outputs);

  /* Copy over outputs to the caller. */
  if (nlhs < 1) {
    b_nlhs = 1;
  } else {
    b_nlhs = nlhs;
  }

  emlrtReturnArrays(b_nlhs, plhs, outputs);

  /* Module termination. */
  fastode45_MEE_terminate();
}

void mexFunction(int32_T nlhs, mxArray *plhs[], int32_T nrhs, const mxArray
                 *prhs[])
{
  mexAtExit(fastode45_MEE_atexit);

  /* Initialize the memory manager. */
  /* Module initialization. */
  fastode45_MEE_initialize();

  /* Dispatch the entry-point. */
  fastode45_MEE_mexFunction(nlhs, plhs, nrhs, prhs);
}

emlrtCTX mexFunctionCreateRootTLS(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  return emlrtRootTLSGlobal;
}

/* End of code generation (_coder_fastode45_MEE_mex.c) */
