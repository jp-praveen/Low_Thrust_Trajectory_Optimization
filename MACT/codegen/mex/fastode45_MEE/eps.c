/*
 * eps.c
 *
 * Code generation for function 'eps'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fastode45_MEE.h"
#include "eps.h"

/* Function Definitions */
real_T eps(real_T x)
{
  real_T r;
  real_T absxk;
  int32_T exponent;
  absxk = muDoubleScalarAbs(x);
  if ((!muDoubleScalarIsInf(absxk)) && (!muDoubleScalarIsNaN(absxk))) {
    if (absxk <= 2.2250738585072014E-308) {
      r = 4.94065645841247E-324;
    } else {
      frexp(absxk, &exponent);
      r = ldexp(1.0, exponent - 53);
    }
  } else {
    r = rtNaN;
  }

  return r;
}

/* End of code generation (eps.c) */
