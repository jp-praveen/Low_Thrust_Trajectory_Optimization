/*
 * StCO_Dynamics_MEE.c
 *
 * Code generation for function 'StCO_Dynamics_MEE'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fastode45_MEE.h"
#include "StCO_Dynamics_MEE.h"
#include "error.h"
#include "sqrt.h"
#include "fastode45_MEE_data.h"

/* Variable Definitions */
static emlrtRSInfo k_emlrtRSI = { 25,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 30,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 31,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 47,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 61,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 63,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 64,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo r_emlrtRSI = { 66,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo s_emlrtRSI = { 72,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo t_emlrtRSI = { 74,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo u_emlrtRSI = { 75,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 82,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo w_emlrtRSI = { 83,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo x_emlrtRSI = { 84,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo y_emlrtRSI = { 85,  /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo ab_emlrtRSI = { 86, /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo bb_emlrtRSI = { 87, /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo cb_emlrtRSI = { 89, /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo db_emlrtRSI = { 99, /* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

static emlrtRSInfo eb_emlrtRSI = { 102,/* lineNo */
  "StCO_Dynamics_MEE",                 /* fcnName */
  "C:\\Praveen\\Matlab_2017\\graduate\\Assignment_7\\StCO_Dynamics_MEE.m"/* pathName */
};

/* Function Definitions */
void StCO_Dynamics_MEE(const emlrtStack *sp, const real_T in2[14], real_T c,
  real_T Thr, real_T rho, real_T SI2CAN, real_T Fdot[14])
{
  real_T t2;
  real_T t3;
  real_T t4;
  real_T t5;
  real_T t8;
  real_T t9;
  real_T t12;
  real_T t13;
  real_T t14;
  real_T t15;
  real_T t16;
  real_T t17;
  real_T t31;
  real_T t18;
  real_T t32;
  real_T t33;
  real_T t34;
  real_T t19;
  real_T t36;
  real_T t37;
  real_T t20;
  real_T t21;
  boolean_T p;
  real_T t25;
  real_T t26;
  real_T t28;
  real_T t29;
  real_T t23;
  real_T t38;
  real_T t39;
  real_T t40;
  real_T t45;
  real_T t46;
  real_T t47;
  real_T t48;
  real_T t54;
  real_T t60;
  real_T t62;
  real_T t63;
  real_T t64;
  real_T t65;
  real_T t66;
  real_T t79;
  real_T t81;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;

  /* STCO_DYNAMICS_MEE */
  /*     FDOT = STCO_DYNAMICS_MEE(T,IN2,C,THR,RHO,SI2CAN) */
  /*     This function was generated by the Symbolic Math Toolbox version 7.2. */
  /*     08-Oct-2020 19:30:47 */
  t2 = 1.0 / in2[6];
  t3 = muDoubleScalarCos(in2[5]);
  t4 = muDoubleScalarSin(in2[5]);
  t5 = in2[0];
  st.site = &k_emlrtRSI;
  b_sqrt(&st, &t5);
  t8 = (in2[1] * t3 * 2.0 + in2[2] * t4 * 2.0) + 2.0;
  t9 = 1.0 / t8;
  st.site = &l_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  st.site = &m_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  t12 = (in2[3] * in2[3] + in2[4] * in2[4]) + 1.0;
  t13 = in2[3] * t4;
  t14 = in2[1] * t3;
  t15 = in2[2] * t4;
  t16 = (t14 + t15) + 1.0;
  t17 = 1.0 / t16;
  t31 = in2[4] * t3;
  t18 = t13 - t31;
  t32 = in2[10] * t3 * t5 * t9 * t12;
  t33 = in2[11] * t4 * t5 * t9 * t12;
  t34 = in2[8] * in2[2] * t5 * t17 * t18;
  t19 = (((t32 + t33) - t34) + in2[12] * t5 * t17 * t18) + in2[9] * in2[1] * t5 *
    t17 * t18;
  t36 = in2[9] * t3 * t5;
  t37 = in2[8] * t4 * t5;
  t20 = -t36 + t37;
  st.site = &n_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  t21 = muDoubleScalarPower(in2[0], 1.5);
  p = false;
  if (in2[0] < 0.0) {
    p = true;
  }

  if (p) {
    c_st.site = &hb_emlrtRSI;
    b_error(&c_st);
  }

  t18 = (t14 + t15) + 2.0;
  t25 = t3 * t18;
  t26 = in2[1] + t25;
  t28 = t4 * t18;
  t29 = in2[2] + t28;
  t23 = (in2[7] * t17 * t21 * 2.0 + in2[8] * t5 * t17 * t26) + in2[9] * t5 * t17
    * t29;
  t18 = (((t32 + t33) - t34) + in2[12] * t5 * t17 * (t13 - t31)) + in2[9] * in2
    [1] * t5 * t17 * (t13 - t31);
  t38 = t36 - t37;
  st.site = &o_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  t39 = t23 * t23;
  t40 = 1.0 / rho;
  st.site = &p_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  st.site = &q_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  t14 = t38 * t38;
  t45 = (t39 + t18 * t18) + t14;
  t15 = t45;
  st.site = &r_emlrtRSI;
  b_sqrt(&st, &t15);
  t46 = 1.0 / t15;
  t47 = -t13 + t31;
  t48 = (((t32 + t33) - in2[12] * t5 * t17 * t47) - in2[9] * in2[1] * t5 * t17 *
         t47) + in2[8] * in2[2] * t5 * t17 * t47;
  st.site = &s_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  t18 = (t39 + t14) + t48 * t48;
  t15 = t18;
  st.site = &t_emlrtRSI;
  b_sqrt(&st, &t15);
  t54 = 1.0 / t15;
  st.site = &u_emlrtRSI;
  b_sqrt(&st, &t18);
  t60 = muDoubleScalarTanh(t40 * (in2[13] + SI2CAN * c * t2 * t18)) * 0.5;
  st.site = &v_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  t62 = t16 * t16;
  t15 = in2[0];
  st.site = &w_emlrtRSI;
  b_sqrt(&st, &t15);
  t63 = 1.0 / t15;
  st.site = &x_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  p = false;
  if (in2[0] < 0.0) {
    p = true;
  }

  if (p) {
    c_st.site = &hb_emlrtRSI;
    b_error(&c_st);
  }

  t64 = 1.0 / muDoubleScalarPower(in2[0], 1.5);
  st.site = &y_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  t65 = 1.0 / (t16 * t16);
  st.site = &ab_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  t66 = t3 * t3;
  st.site = &bb_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  t14 = 1.0 / (t8 * t8);
  t15 = t5 * t17 * t47 * t48 * t54;
  st.site = &cb_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  t34 = t4 * t4;
  t36 = in2[1] * t4;
  t37 = in2[2] * t3;
  t13 = t36 - t37;
  t31 = in2[3] * t3 + in2[4] * t4;
  t32 = t3 * t5 * t38 * t54;
  t33 = in2[1] * t4 * 2.0;
  t8 = in2[2] * t5 * t17 * t47 * t48 * t54;
  st.site = &db_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  t79 = 1.0 / (in2[6] * in2[6]);
  t81 = (t32 + t5 * t17 * t23 * t29 * t54) - in2[1] * t5 * t17 * t47 * t48 * t54;
  st.site = &eb_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  st.site = &eb_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  st.site = &eb_emlrtRSI;
  t18 = (t39 + t19 * t19) + t20 * t20;
  if (t18 < 0.0) {
    b_st.site = &fb_emlrtRSI;
    error(&b_st);
  }

  st.site = &eb_emlrtRSI;
  if (t45 < 0.0) {
    b_st.site = &fb_emlrtRSI;
    error(&b_st);
  }

  st.site = &eb_emlrtRSI;
  b_st.site = &gb_emlrtRSI;
  p = false;
  if (in2[0] < 0.0) {
    p = true;
  }

  if (p) {
    c_st.site = &hb_emlrtRSI;
    b_error(&c_st);
  }

  Fdot[0] = SI2CAN * Thr * t2 * t17 * t21 * t23 * t46 * (muDoubleScalarTanh(t40 *
    (in2[13] + SI2CAN * c * t2 * muDoubleScalarSqrt(t18))) * 0.5 + 0.5) * -2.0;
  Fdot[1] = -SI2CAN * Thr * t2 * (muDoubleScalarTanh(t40 * (in2[13] + SI2CAN * c
    * t2 * muDoubleScalarSqrt(t45))) * 0.5 + 0.5) * ((t8 - t4 * t5 * t38 * t46)
    + t5 * t17 * t23 * t26 * t46);
  Fdot[2] = -SI2CAN * Thr * t2 * (t60 + 0.5) * t81;
  Fdot[3] = -SI2CAN * Thr * t2 * t3 * t5 * t9 * t12 * t48 * t54 * (t60 + 0.5);
  Fdot[4] = -SI2CAN * Thr * t2 * t4 * t5 * t9 * t12 * t48 * t54 * (t60 + 0.5);
  Fdot[5] = t62 * t64 + SI2CAN * Thr * t2 * t5 * t17 * t47 * t48 * t54 * (t60 +
    0.5);
  Fdot[6] = -(Thr * (t60 + 0.5)) / c;
  Fdot[7] = ((((in2[12] * (1.0 / muDoubleScalarPower(in2[0], 2.5) * t62 * 1.5 -
    SI2CAN * Thr * t2 * t17 * t47 * t48 * t54 * (t60 + 0.5) * t63 * 0.5) + in2[9]
                * SI2CAN * Thr * t2 * (t60 + 0.5) * ((t3 * t38 * t54 * t63 * 0.5
    + t17 * t23 * t29 * t54 * t63 * 0.5) - in2[1] * t17 * t47 * t48 * t54 * t63 *
    0.5)) + in2[8] * SI2CAN * Thr * t2 * (t60 + 0.5) * ((t4 * t38 * t54 * t63 *
    -0.5 + t17 * t23 * t26 * t54 * t63 * 0.5) + in2[2] * t17 * t47 * t48 * t54 *
    t63 * 0.5)) + in2[7] * SI2CAN * Thr * t2 * t5 * t17 * t23 * t54 * (t60 + 0.5)
              * 3.0) + in2[10] * SI2CAN * Thr * t2 * t3 * t9 * t12 * t48 * t54 *
             (t60 + 0.5) * t63 * 0.5) + in2[11] * SI2CAN * Thr * t2 * t4 * t9 *
    t12 * t48 * t54 * (t60 + 0.5) * t63 * 0.5;
  Fdot[8] = ((((-in2[12] * (t3 * t16 * t64 * 2.0 - SI2CAN * Thr * t2 * t3 * t5 *
    t47 * t48 * t54 * (t60 + 0.5) * t65) - in2[8] * SI2CAN * Thr * t2 * (t60 +
    0.5) * ((-t5 * t17 * t23 * t54 * (t66 + 1.0) + t3 * t5 * t23 * t26 * t54 *
             t65) + in2[2] * t3 * t5 * t47 * t48 * t54 * t65)) - in2[9] * SI2CAN
               * Thr * t2 * (t60 + 0.5) * (((t15 - t3 * t4 * t5 * t17 * t23 *
    t54) + t3 * t5 * t23 * t29 * t54 * t65) - in2[1] * t3 * t5 * t47 * t48 * t54
    * t65)) - in2[7] * SI2CAN * Thr * t2 * t3 * t21 * t23 * t54 * (t60 + 0.5) *
              t65 * 2.0) - in2[10] * SI2CAN * Thr * t2 * t5 * t12 * t48 * t54 *
             (t60 + 0.5) * t66 * t14 * 2.0) - in2[11] * SI2CAN * Thr * t2 * t3 *
    t4 * t5 * t12 * t48 * t54 * (t60 + 0.5) * t14 * 2.0;
  Fdot[9] = ((((-in2[12] * (t4 * t16 * t64 * 2.0 - SI2CAN * Thr * t2 * t4 * t5 *
    t47 * t48 * t54 * (t60 + 0.5) * t65) + in2[9] * SI2CAN * Thr * t2 * (t60 +
    0.5) * ((t5 * t17 * t23 * t54 * (t34 + 1.0) - t4 * t5 * t23 * t29 * t54 *
             t65) + in2[1] * t4 * t5 * t47 * t48 * t54 * t65)) + in2[8] * SI2CAN
               * Thr * t2 * (t60 + 0.5) * (((t15 + t3 * t4 * t5 * t17 * t23 *
    t54) - t4 * t5 * t23 * t26 * t54 * t65) - in2[2] * t4 * t5 * t47 * t48 * t54
    * t65)) - in2[7] * SI2CAN * Thr * t2 * t4 * t21 * t23 * t54 * (t60 + 0.5) *
              t65 * 2.0) - in2[11] * SI2CAN * Thr * t2 * t5 * t12 * t48 * t54 *
             (t60 + 0.5) * t14 * t34 * 2.0) - in2[10] * SI2CAN * Thr * t2 * t3 *
    t4 * t5 * t12 * t48 * t54 * (t60 + 0.5) * t14 * 2.0;
  Fdot[10] = (((in2[12] * SI2CAN * Thr * t2 * t4 * t5 * t17 * t48 * t54 * (t60 +
    0.5) + in2[9] * SI2CAN * Thr * in2[1] * t2 * t4 * t5 * t17 * t48 * t54 *
                (t60 + 0.5)) - in2[8] * SI2CAN * Thr * in2[2] * t2 * t4 * t5 *
               t17 * t48 * t54 * (t60 + 0.5)) + in2[10] * SI2CAN * Thr * in2[3] *
              t2 * t3 * t5 * t9 * t48 * t54 * (t60 + 0.5) * 2.0) + in2[11] *
    SI2CAN * Thr * in2[3] * t2 * t4 * t5 * t9 * t48 * t54 * (t60 + 0.5) * 2.0;
  Fdot[11] = (((-in2[12] * SI2CAN * Thr * t2 * t3 * t5 * t17 * t48 * t54 * (t60
    + 0.5) - in2[9] * SI2CAN * Thr * in2[1] * t2 * t3 * t5 * t17 * t48 * t54 *
                (t60 + 0.5)) + in2[8] * SI2CAN * Thr * in2[2] * t2 * t3 * t5 *
               t17 * t48 * t54 * (t60 + 0.5)) + in2[10] * SI2CAN * Thr * in2[4] *
              t2 * t3 * t5 * t9 * t48 * t54 * (t60 + 0.5) * 2.0) + in2[11] *
    SI2CAN * Thr * in2[4] * t2 * t4 * t5 * t9 * t48 * t54 * (t60 + 0.5) * 2.0;
  Fdot[12] = ((((((in2[12] * ((t16 * t64 * t13 * 2.0 + SI2CAN * Thr * t2 * t5 *
    t17 * t48 * t54 * (t60 + 0.5) * t31) - SI2CAN * Thr * t2 * t5 * t47 * t48 *
    t54 * (t60 + 0.5) * t65 * t13) + in2[9] * SI2CAN * Thr * t2 * (t60 + 0.5) *
                   ((((-t4 * t5 * t38 * t54 + t5 * t17 * t23 * t54 * (t25 - t4 *
    t13)) + in2[1] * t5 * t17 * t48 * t54 * t31) + t5 * t23 * t29 * t54 * t65 *
                     (t36 - t37)) - in2[1] * t5 * t47 * t48 * t54 * t65 * t13))
                  - in2[8] * SI2CAN * Thr * t2 * (t60 + 0.5) * ((((t32 + t5 *
    t17 * t23 * t54 * (t28 + t3 * t13)) + in2[2] * t5 * t17 * t48 * t54 * t31) -
    t5 * t23 * t26 * t54 * t65 * t13) - in2[2] * t5 * t47 * t48 * t54 * t65 *
    t13)) + in2[7] * SI2CAN * Thr * t2 * t21 * t23 * t54 * (t60 + 0.5) * t65 *
                 (t36 - t37) * 2.0) - in2[10] * SI2CAN * Thr * t2 * t4 * t5 * t9
                * t12 * t48 * t54 * (t60 + 0.5)) + in2[11] * SI2CAN * Thr * t2 *
               t3 * t5 * t9 * t12 * t48 * t54 * (t60 + 0.5)) + in2[10] * SI2CAN *
              Thr * t2 * t3 * t5 * t12 * t48 * t54 * (t60 + 0.5) * t14 * (t33 -
    in2[2] * t3 * 2.0)) + in2[11] * SI2CAN * Thr * t2 * t4 * t5 * t12 * t48 *
    t54 * (t60 + 0.5) * t14 * (t33 - in2[2] * t3 * 2.0);
  Fdot[13] = ((((-in2[9] * SI2CAN * Thr * (t60 + 0.5) * t79 * t81 - in2[8] *
                 SI2CAN * Thr * (t60 + 0.5) * t79 * ((t8 - t4 * t5 * t38 * t54)
    + t5 * t17 * t23 * t26 * t54)) - in2[7] * SI2CAN * Thr * t17 * t21 * t23 *
                t54 * (t60 + 0.5) * t79 * 2.0) + in2[12] * SI2CAN * Thr * t5 *
               t17 * t47 * t48 * t54 * (t60 + 0.5) * t79) - in2[10] * SI2CAN *
              Thr * t3 * t5 * t9 * t12 * t48 * t54 * (t60 + 0.5) * t79) - in2[11]
    * SI2CAN * Thr * t4 * t5 * t9 * t12 * t48 * t54 * (t60 + 0.5) * t79;
}

/* End of code generation (StCO_Dynamics_MEE.c) */
