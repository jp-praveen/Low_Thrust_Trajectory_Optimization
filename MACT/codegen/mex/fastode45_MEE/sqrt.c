/*
 * sqrt.c
 *
 * Code generation for function 'sqrt'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fastode45_MEE.h"
#include "sqrt.h"
#include "error.h"
#include "fastode45_MEE_data.h"

/* Function Definitions */
void b_sqrt(const emlrtStack *sp, real_T *x)
{
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  if (*x < 0.0) {
    st.site = &fb_emlrtRSI;
    error(&st);
  }

  *x = muDoubleScalarSqrt(*x);
}

/* End of code generation (sqrt.c) */
