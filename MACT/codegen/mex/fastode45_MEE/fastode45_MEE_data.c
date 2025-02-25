/*
 * fastode45_MEE_data.c
 *
 * Code generation for function 'fastode45_MEE_data'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fastode45_MEE.h"
#include "fastode45_MEE_data.h"

/* Variable Definitions */
emlrtCTX emlrtRootTLSGlobal = NULL;
const volatile char_T *emlrtBreakCheckR2012bFlagVar = NULL;
real_T b_atol;
real_T rtol;
real_T b_pow;
real_T A[6];
real_T B[42];
real_T E[7];
uint32_T atol_dirty;
uint32_T rtol_dirty;
uint32_T pow_dirty;
uint32_T A_dirty;
uint32_T B_dirty;
uint32_T E_dirty;
emlrtContext emlrtContextGlobal = { true,/* bFirstTime */
  false,                               /* bInitialized */
  131450U,                             /* fVersionInfo */
  NULL,                                /* fErrorFunction */
  "fastode45_MEE",                     /* fFunctionName */
  NULL,                                /* fRTCallStack */
  false,                               /* bDebugMode */
  { 2045744189U, 2170104910U, 2743257031U, 4284093946U },/* fSigWrd */
  NULL                                 /* fSigMem */
};

emlrtRSInfo fb_emlrtRSI = { 12,        /* lineNo */
  "sqrt",                              /* fcnName */
  "C:\\Program Files\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\elfun\\sqrt.m"/* pathName */
};

emlrtRSInfo gb_emlrtRSI = { 49,        /* lineNo */
  "power",                             /* fcnName */
  "C:\\Program Files\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\ops\\power.m"/* pathName */
};

emlrtRSInfo hb_emlrtRSI = { 61,        /* lineNo */
  "power",                             /* fcnName */
  "C:\\Program Files\\MATLAB\\R2017a\\toolbox\\eml\\lib\\matlab\\ops\\power.m"/* pathName */
};

/* End of code generation (fastode45_MEE_data.c) */
