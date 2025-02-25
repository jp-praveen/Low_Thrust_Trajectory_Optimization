/*
 * fastode45_MEE_data.h
 *
 * Code generation for function 'fastode45_MEE_data'
 *
 */

#ifndef FASTODE45_MEE_DATA_H
#define FASTODE45_MEE_DATA_H

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

/* Variable Declarations */
extern emlrtCTX emlrtRootTLSGlobal;
extern const volatile char_T *emlrtBreakCheckR2012bFlagVar;
extern real_T b_atol;
extern real_T rtol;
extern real_T b_pow;
extern real_T A[6];
extern real_T B[42];
extern real_T E[7];
extern uint32_T atol_dirty;
extern uint32_T rtol_dirty;
extern uint32_T pow_dirty;
extern uint32_T A_dirty;
extern uint32_T B_dirty;
extern uint32_T E_dirty;
extern emlrtContext emlrtContextGlobal;
extern emlrtRSInfo fb_emlrtRSI;
extern emlrtRSInfo gb_emlrtRSI;
extern emlrtRSInfo hb_emlrtRSI;

#endif

/* End of code generation (fastode45_MEE_data.h) */
