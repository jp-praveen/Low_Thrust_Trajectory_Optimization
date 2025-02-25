/*
 * fastode45_MEE_terminate.c
 *
 * Code generation for function 'fastode45_MEE_terminate'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fastode45_MEE.h"
#include "fastode45_MEE_terminate.h"
#include "_coder_fastode45_MEE_mex.h"
#include "fastode45_MEE_data.h"

/* Function Definitions */
void fastode45_MEE_atexit(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  st.tls = emlrtRootTLSGlobal;
  emlrtEnterRtStackR2012b(&st);
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void fastode45_MEE_terminate(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  emlrtLeaveRtStackR2012b(&st);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (fastode45_MEE_terminate.c) */
