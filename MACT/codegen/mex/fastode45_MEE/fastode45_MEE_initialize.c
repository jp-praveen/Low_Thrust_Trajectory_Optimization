/*
 * fastode45_MEE_initialize.c
 *
 * Code generation for function 'fastode45_MEE_initialize'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fastode45_MEE.h"
#include "fastode45_MEE_initialize.h"
#include "_coder_fastode45_MEE_mex.h"
#include "fastode45_MEE_data.h"

/* Named Constants */
#define c_atol                         (1.0E-10)
#define b_rtol                         (1.0E-10)
#define c_pow                          (0.2)

/* Function Declarations */
static void fastode45_MEE_once(void);

/* Function Definitions */
static void fastode45_MEE_once(void)
{
  int32_T i;
  static const real_T dv0[42] = { 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.075,
    0.225, 0.0, 0.0, 0.0, 0.0, 0.0, 0.97777777777777775, -3.7333333333333334,
    3.5555555555555554, 0.0, 0.0, 0.0, 0.0, 2.9525986892242035,
    -11.595793324188385, 9.8228928516994358, -0.29080932784636487, 0.0, 0.0, 0.0,
    2.8462752525252526, -10.757575757575758, 8.9064227177434727,
    0.27840909090909088, -0.2735313036020583, 0.0, 0.0, 0.091145833333333329,
    0.0, 0.44923629829290207, 0.65104166666666663, -0.322376179245283,
    0.13095238095238096, 0.0 };

  static const real_T dv1[7] = { 0.0012326388888888888, 0.0,
    -0.0042527702905061394, 0.036979166666666667, -0.05086379716981132,
    0.0419047619047619, -0.025 };

  static const real_T dv2[6] = { 0.2, 0.3, 0.8, 0.88888888888888884, 1.0, 1.0 };

  for (i = 0; i < 7; i++) {
    E[i] = dv1[i];
  }

  memcpy(&B[0], &dv0[0], 42U * sizeof(real_T));
  for (i = 0; i < 6; i++) {
    A[i] = dv2[i];
  }

  E_dirty = 0U;
  B_dirty = 0U;
  A_dirty = 0U;
  pow_dirty = 0U;
  rtol_dirty = 0U;
  atol_dirty = 0U;
  b_pow = c_pow;
  rtol = b_rtol;
  b_atol = c_atol;
}

void fastode45_MEE_initialize(void)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  mexFunctionCreateRootTLS();
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  st.tls = emlrtRootTLSGlobal;
  emlrtClearAllocCountR2012b(&st, false, 0U, 0);
  emlrtEnterRtStackR2012b(&st);
  if (emlrtFirstTimeR2012b(emlrtRootTLSGlobal)) {
    fastode45_MEE_once();
  }
}

/* End of code generation (fastode45_MEE_initialize.c) */
