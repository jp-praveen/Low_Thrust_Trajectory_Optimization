/*
 * eps.h
 *
 * Code generation for function 'eps'
 *
 */

#ifndef EPS_H
#define EPS_H

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
extern real_T eps(real_T x);

#ifdef __WATCOMC__

#pragma aux eps value [8087];

#endif
#endif

/* End of code generation (eps.h) */
