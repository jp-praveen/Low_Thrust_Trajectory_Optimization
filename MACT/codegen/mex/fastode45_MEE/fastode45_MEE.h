/*
 * fastode45_MEE.h
 *
 * Code generation for function 'fastode45_MEE'
 *
 */

#ifndef FASTODE45_MEE_H
#define FASTODE45_MEE_H

/* Include files */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"
#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "rtwtypes.h"
#include "fastode45_MEE_types.h"

/* Function Declarations */
extern void fastode45_MEE(const emlrtStack *sp, const real_T tspan[2], const
  real_T b_y0[14], real_T c, real_T Thr, real_T rho, real_T SI2CAN, real_T
  *t_out, real_T y_out[14]);

#endif

/* End of code generation (fastode45_MEE.h) */
