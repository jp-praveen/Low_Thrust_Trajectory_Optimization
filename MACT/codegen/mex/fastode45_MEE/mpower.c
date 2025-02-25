/*
 * mpower.c
 *
 * Code generation for function 'mpower'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fastode45_MEE.h"
#include "mpower.h"
#include "error.h"
#include "fastode45_MEE_data.h"

/* Variable Definitions */
static emlrtRSInfo ib_emlrtRSI = { 37, /* lineNo */
  "mpower",                            /* fcnName */
  "C:\\Program Files\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\ops\\mpower.m"/* pathName */
};

/* Function Definitions */
real_T mpower(const emlrtStack *sp, real_T a, real_T b)
{
  real_T c;
  boolean_T p;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &ib_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  b_st.site = &gb_emlrtRSI;
  c = muDoubleScalarPower(a, b);
  p = false;
  if (a < 0.0) {
    if (muDoubleScalarIsNaN(b) || (muDoubleScalarFloor(b) == b)) {
      p = true;
    } else {
      p = false;
    }

    p = !p;
  }

  if (p) {
    c_st.site = &hb_emlrtRSI;
    b_error(&c_st);
  }

  return c;
}

/* End of code generation (mpower.c) */
