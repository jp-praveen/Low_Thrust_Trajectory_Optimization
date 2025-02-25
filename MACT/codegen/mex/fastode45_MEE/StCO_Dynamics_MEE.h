/*
 * StCO_Dynamics_MEE.h
 *
 * Code generation for function 'StCO_Dynamics_MEE'
 *
 */

#ifndef STCO_DYNAMICS_MEE_H
#define STCO_DYNAMICS_MEE_H

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
extern void StCO_Dynamics_MEE(const emlrtStack *sp, const real_T in2[14], real_T
  c, real_T Thr, real_T rho, real_T SI2CAN, real_T Fdot[14]);

#endif

/* End of code generation (StCO_Dynamics_MEE.h) */
