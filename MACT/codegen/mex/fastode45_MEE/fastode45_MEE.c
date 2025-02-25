/*
 * fastode45_MEE.c
 *
 * Code generation for function 'fastode45_MEE'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fastode45_MEE.h"
#include "mpower.h"
#include "StCO_Dynamics_MEE.h"
#include "eps.h"
#include "fastode45_MEE_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 6,     /* lineNo */
  "fastode45_MEE",                     /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\fastode45_MEE.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 22,  /* lineNo */
  "fastode45_MEE",                     /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\fastode45_MEE.m"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 54,  /* lineNo */
  "fastode45_MEE",                     /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\fastode45_MEE.m"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 55,  /* lineNo */
  "fastode45_MEE",                     /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\fastode45_MEE.m"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 56,  /* lineNo */
  "fastode45_MEE",                     /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\fastode45_MEE.m"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 57,  /* lineNo */
  "fastode45_MEE",                     /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\fastode45_MEE.m"/* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 58,  /* lineNo */
  "fastode45_MEE",                     /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\fastode45_MEE.m"/* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 67,  /* lineNo */
  "fastode45_MEE",                     /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\fastode45_MEE.m"/* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 87,  /* lineNo */
  "fastode45_MEE",                     /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\fastode45_MEE.m"/* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 108, /* lineNo */
  "fastode45_MEE",                     /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\fastode45_MEE.m"/* pathName */
};

/* Function Definitions */
void fastode45_MEE(const emlrtStack *sp, const real_T tspan[2], const real_T
                   b_y0[14], real_T c, real_T Thr, real_T rho, real_T SI2CAN,
                   real_T *t_out, real_T y_out[14])
{
  real_T tfinal;
  real_T f0[14];
  real_T threshold;
  real_T hmax;
  real_T t;
  real_T y[14];
  real_T f[98];
  real_T absh;
  int32_T i;
  real_T b_y;
  real_T maxval[14];
  boolean_T exitg1;
  real_T absx;
  boolean_T done;
  int32_T exitg3;
  real_T hmin;
  boolean_T nofailed;
  int32_T exitg2;
  int32_T i0;
  real_T hA[6];
  real_T hB[42];
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  tfinal = tspan[1];
  st.site = &emlrtRSI;
  StCO_Dynamics_MEE(&st, b_y0, c, Thr, rho, SI2CAN, f0);
  threshold = b_atol / rtol;

  /*  By default, hmax is 1/10 of the interval. */
  hmax = 0.1 * (tspan[1] - tspan[0]);
  t = tspan[0];
  memcpy(&y[0], &b_y0[0], 14U * sizeof(real_T));

  /* <---------------- */
  memset(&f[0], 0, 98U * sizeof(real_T));

  /*  the first number is number of equations */
  /*  Compute an initial step size h using y'(t). */
  absh = muDoubleScalarMin(hmax, tspan[1] - tspan[0]);
  for (i = 0; i < 14; i++) {
    maxval[i] = f0[i] / muDoubleScalarMax(muDoubleScalarAbs(b_y0[i]), threshold);
  }

  b_y = 0.0;
  i = 0;
  exitg1 = false;
  while ((!exitg1) && (i < 14)) {
    absx = muDoubleScalarAbs(maxval[i]);
    if (muDoubleScalarIsNaN(absx)) {
      b_y = rtNaN;
      exitg1 = true;
    } else {
      if (absx > b_y) {
        b_y = absx;
      }

      i++;
    }
  }

  st.site = &b_emlrtRSI;
  absx = 0.8 * mpower(&st, rtol, b_pow);
  absx = b_y / absx;
  if (absh * absx > 1.0) {
    absh = 1.0 / absx;
  }

  absh = muDoubleScalarMax(absh, 16.0 * eps(tspan[0]));
  memcpy(&f[0], &f0[0], 14U * sizeof(real_T));

  /*  THE MAIN LOOP */
  done = false;
  do {
    exitg3 = 0;

    /*  By default, hmin is a small number such that t+hmin is only slightly */
    /*  different than t.  It might be 0 if t is 0. */
    hmin = 16.0 * eps(t);
    absh = muDoubleScalarMin(hmax, muDoubleScalarMax(hmin, absh));

    /*  couldn't limit absh until new hmin */
    absx = absh;

    /*  Stretch the step if within 10% of tfinal-t. */
    if (1.1 * absh >= muDoubleScalarAbs(tfinal - t)) {
      absx = tfinal - t;
      absh = muDoubleScalarAbs(absx);
      done = true;
    }

    /*  LOOP FOR ADVANCING ONE STEP. */
    nofailed = true;

    /*  no failed attempts */
    do {
      exitg2 = 0;
      for (i0 = 0; i0 < 6; i0++) {
        hA[i0] = absx * A[i0];
      }

      for (i0 = 0; i0 < 42; i0++) {
        hB[i0] = absx * B[i0];
      }

      for (i0 = 0; i0 < 14; i0++) {
        absx = 0.0;
        for (i = 0; i < 7; i++) {
          absx += f[i0 + 14 * i] * hB[i];
        }

        f0[i0] = y[i0] + absx;
      }

      st.site = &c_emlrtRSI;
      StCO_Dynamics_MEE(&st, f0, c, Thr, rho, SI2CAN, *(real_T (*)[14])&f[14]);
      for (i0 = 0; i0 < 14; i0++) {
        absx = 0.0;
        for (i = 0; i < 7; i++) {
          absx += f[i0 + 14 * i] * hB[7 + i];
        }

        f0[i0] = y[i0] + absx;
      }

      st.site = &d_emlrtRSI;
      StCO_Dynamics_MEE(&st, f0, c, Thr, rho, SI2CAN, *(real_T (*)[14])&f[28]);
      for (i0 = 0; i0 < 14; i0++) {
        absx = 0.0;
        for (i = 0; i < 7; i++) {
          absx += f[i0 + 14 * i] * hB[14 + i];
        }

        f0[i0] = y[i0] + absx;
      }

      st.site = &e_emlrtRSI;
      StCO_Dynamics_MEE(&st, f0, c, Thr, rho, SI2CAN, *(real_T (*)[14])&f[42]);
      for (i0 = 0; i0 < 14; i0++) {
        absx = 0.0;
        for (i = 0; i < 7; i++) {
          absx += f[i0 + 14 * i] * hB[21 + i];
        }

        f0[i0] = y[i0] + absx;
      }

      st.site = &f_emlrtRSI;
      StCO_Dynamics_MEE(&st, f0, c, Thr, rho, SI2CAN, *(real_T (*)[14])&f[56]);
      for (i0 = 0; i0 < 14; i0++) {
        absx = 0.0;
        for (i = 0; i < 7; i++) {
          absx += f[i0 + 14 * i] * hB[28 + i];
        }

        f0[i0] = y[i0] + absx;
      }

      st.site = &g_emlrtRSI;
      StCO_Dynamics_MEE(&st, f0, c, Thr, rho, SI2CAN, *(real_T (*)[14])&f[70]);
      *t_out = t + hA[5];
      if (done) {
        *t_out = tfinal;

        /*  Hit end point exactly. */
      }

      /*        h = tnew - t;      % Purify h. */
      for (i0 = 0; i0 < 14; i0++) {
        absx = 0.0;
        for (i = 0; i < 7; i++) {
          absx += f[i0 + 14 * i] * hB[35 + i];
        }

        y_out[i0] = y[i0] + absx;
      }

      st.site = &h_emlrtRSI;
      StCO_Dynamics_MEE(&st, y_out, c, Thr, rho, SI2CAN, *(real_T (*)[14])&f[84]);

      /*  Estimate the error. */
      for (i = 0; i < 14; i++) {
        maxval[i] = muDoubleScalarMax(muDoubleScalarAbs(y[i]), muDoubleScalarAbs
          (y_out[i]));
      }

      for (i = 0; i < 14; i++) {
        absx = 0.0;
        for (i0 = 0; i0 < 7; i0++) {
          absx += f[i + 14 * i0] * E[i0];
        }

        f0[i] = absx / muDoubleScalarMax(maxval[i], threshold);
      }

      b_y = 0.0;
      i = 0;
      exitg1 = false;
      while ((!exitg1) && (i < 14)) {
        absx = muDoubleScalarAbs(f0[i]);
        if (muDoubleScalarIsNaN(absx)) {
          b_y = rtNaN;
          exitg1 = true;
        } else {
          if (absx > b_y) {
            b_y = absx;
          }

          i++;
        }
      }

      absx = absh * b_y;

      /*  Accept the solution only if the weighted error is no more than the */
      /*  tolerance rtol.  Estimate an h that will yield an error of rtol on */
      /*  the next step or the next try at taking this step, as the case may be, */
      /*  and use 0.8 of this value to avoid failures. */
      if (absx > rtol) {
        /*  Failed step */
        if (absh <= hmin) {
          *t_out = 0.0;
          memset(&y_out[0], 0, 14U * sizeof(real_T));

          /*              sprintf('tolerance not satisfied'); */
          exitg2 = 1;
        } else {
          if (nofailed) {
            nofailed = false;
            absx = rtol / absx;
            st.site = &i_emlrtRSI;
            absx = 0.8 * mpower(&st, absx, b_pow);
            absh = muDoubleScalarMax(hmin, absh * muDoubleScalarMax(0.1, absx));
          } else {
            absh = muDoubleScalarMax(hmin, 0.5 * absh);
          }

          absx = absh;
          done = false;
          if (*emlrtBreakCheckR2012bFlagVar != 0) {
            emlrtBreakCheckR2012b(sp);
          }
        }
      } else {
        exitg2 = 2;
      }
    } while (exitg2 == 0);

    if (exitg2 == 1) {
      exitg3 = 1;
    } else if (done) {
      exitg3 = 1;
    } else {
      /*  If there were no failures compute a new h. */
      if (nofailed) {
        /*  Note that absh may shrink by 0.8, and that err may be 0. */
        absx /= rtol;
        st.site = &j_emlrtRSI;
        absx = 1.25 * mpower(&st, absx, b_pow);
        if (absx > 0.2) {
          absh /= absx;
        } else {
          absh *= 5.0;
        }
      }

      /*  Advance the integration one step. */
      t = *t_out;
      for (i = 0; i < 14; i++) {
        y[i] = y_out[i];
        f[i] = f[84 + i];
      }

      /*  Already have f(tnew,ynew) */
      if (*emlrtBreakCheckR2012bFlagVar != 0) {
        emlrtBreakCheckR2012b(sp);
      }
    }
  } while (exitg3 == 0);
}

/* End of code generation (fastode45_MEE.c) */
