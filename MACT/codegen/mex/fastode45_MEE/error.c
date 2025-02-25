/*
 * error.c
 *
 * Code generation for function 'error'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fastode45_MEE.h"
#include "error.h"

/* Variable Definitions */
static emlrtRTEInfo emlrtRTEI = { 17,  /* lineNo */
  9,                                   /* colNo */
  "error",                             /* fName */
  "C:\\Program Files\\MATLAB\\R2017a\\toolbox\\eml\\eml\\+coder\\+internal\\error.m"/* pName */
};

/* Function Definitions */
void b_error(const emlrtStack *sp)
{
  emlrtErrorWithMessageIdR2012b(sp, &emlrtRTEI,
    "Coder:toolbox:power_domainError", 0);
}

void error(const emlrtStack *sp)
{
  emlrtErrorWithMessageIdR2012b(sp, &emlrtRTEI, "Coder:toolbox:ElFunDomainError",
    3, 4, 4, "sqrt");
}

/* End of code generation (error.c) */
