/*
 * _coder_fastode45_MEE_api.c
 *
 * Code generation for function '_coder_fastode45_MEE_api'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fastode45_MEE.h"
#include "_coder_fastode45_MEE_api.h"
#include "fastode45_MEE_data.h"

/* Function Declarations */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[2];
static const mxArray *b_emlrt_marshallOut(const real_T u[6]);
static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_y0,
  const char_T *identifier))[14];
static const mxArray *c_emlrt_marshallOut(const real_T u[42]);
static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[14];
static const mxArray *d_emlrt_marshallOut(const real_T u[7]);
static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *c, const
  char_T *identifier);
static const mxArray *e_emlrt_marshallOut(const real_T u[14]);
static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *tspan,
  const char_T *identifier))[2];
static const mxArray *emlrt_marshallOut(const real_T u);
static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_A, const
  char_T *identifier, real_T y[6]);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[6]);
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_B, const
  char_T *identifier, real_T y[42]);
static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[42]);
static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_E, const
  char_T *identifier, real_T y[7]);
static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[7]);
static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[2];
static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[14];
static real_T o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static void p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[6]);
static void q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[42]);
static void r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[7]);

/* Function Definitions */
static real_T (*b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[2]
{
  real_T (*y)[2];
  y = m_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static const mxArray *b_emlrt_marshallOut(const real_T u[6])
{
  const mxArray *y;
  const mxArray *m1;
  static const int32_T iv0[2] = { 1, 6 };

  real_T *pData;
  int32_T i;
  y = NULL;
  m1 = emlrtCreateNumericArray(2, iv0, mxDOUBLE_CLASS, mxREAL);
  pData = (real_T *)mxGetPr(m1);
  for (i = 0; i < 6; i++) {
    pData[i] = u[i];
  }

  emlrtAssign(&y, m1);
  return y;
}

static real_T (*c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_y0,
  const char_T *identifier))[14]
{
  real_T (*y)[14];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = d_emlrt_marshallIn(sp, emlrtAlias(b_y0), &thisId);
  emlrtDestroyArray(&b_y0);
  return y;
}
  static const mxArray *c_emlrt_marshallOut(const real_T u[42])
{
  const mxArray *y;
  const mxArray *m2;
  static const int32_T iv1[2] = { 7, 6 };

  real_T *pData;
  int32_T i;
  y = NULL;
  m2 = emlrtCreateNumericArray(2, iv1, mxDOUBLE_CLASS, mxREAL);
  pData = (real_T *)mxGetPr(m2);
  for (i = 0; i < 42; i++) {
    pData[i] = u[i];
  }

  emlrtAssign(&y, m2);
  return y;
}

static real_T (*d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId))[14]
{
  real_T (*y)[14];
  y = n_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static const mxArray *d_emlrt_marshallOut(const real_T u[7])
{
  const mxArray *y;
  const mxArray *m3;
  static const int32_T iv2[1] = { 7 };

  real_T *pData;
  int32_T i;
  y = NULL;
  m3 = emlrtCreateNumericArray(1, iv2, mxDOUBLE_CLASS, mxREAL);
  pData = (real_T *)mxGetPr(m3);
  for (i = 0; i < 7; i++) {
    pData[i] = u[i];
  }

  emlrtAssign(&y, m3);
  return y;
}

static real_T e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *c, const
  char_T *identifier)
{
  real_T y;
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = f_emlrt_marshallIn(sp, emlrtAlias(c), &thisId);
  emlrtDestroyArray(&c);
  return y;
}

static const mxArray *e_emlrt_marshallOut(const real_T u[14])
{
  const mxArray *y;
  const mxArray *m4;
  static const int32_T iv3[1] = { 0 };

  static const int32_T iv4[1] = { 14 };

  y = NULL;
  m4 = emlrtCreateNumericArray(1, iv3, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m4, (void *)&u[0]);
  emlrtSetDimensions((mxArray *)m4, iv4, 1);
  emlrtAssign(&y, m4);
  return y;
}

static real_T (*emlrt_marshallIn(const emlrtStack *sp, const mxArray *tspan,
  const char_T *identifier))[2]
{
  real_T (*y)[2];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = b_emlrt_marshallIn(sp, emlrtAlias(tspan), &thisId);
  emlrtDestroyArray(&tspan);
  return y;
}
  static const mxArray *emlrt_marshallOut(const real_T u)
{
  const mxArray *y;
  const mxArray *m0;
  y = NULL;
  m0 = emlrtCreateDoubleScalar(u);
  emlrtAssign(&y, m0);
  return y;
}

static real_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = o_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_A, const
  char_T *identifier, real_T y[6])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  h_emlrt_marshallIn(sp, emlrtAlias(b_A), &thisId, y);
  emlrtDestroyArray(&b_A);
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[6])
{
  p_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_B, const
  char_T *identifier, real_T y[42])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  j_emlrt_marshallIn(sp, emlrtAlias(b_B), &thisId, y);
  emlrtDestroyArray(&b_B);
}

static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[42])
{
  q_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *b_E, const
  char_T *identifier, real_T y[7])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  l_emlrt_marshallIn(sp, emlrtAlias(b_E), &thisId, y);
  emlrtDestroyArray(&b_E);
}

static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T y[7])
{
  r_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T (*m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[2]
{
  real_T (*ret)[2];
  static const int32_T dims[2] = { 1, 2 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  ret = (real_T (*)[2])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  static real_T (*n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId))[14]
{
  real_T (*ret)[14];
  static const int32_T dims[1] = { 14 };

  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  ret = (real_T (*)[14])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static real_T o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  real_T ret;
  static const int32_T dims = 0;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 0U, &dims);
  ret = *(real_T *)mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[6])
{
  static const int32_T dims[2] = { 1, 6 };

  int32_T i1;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  for (i1 = 0; i1 < 6; i1++) {
    ret[i1] = (*(real_T (*)[6])mxGetData(src))[i1];
  }

  emlrtDestroyArray(&src);
}

static void q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[42])
{
  static const int32_T dims[2] = { 7, 6 };

  int32_T i2;
  int32_T i3;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 2U, dims);
  for (i2 = 0; i2 < 6; i2++) {
    for (i3 = 0; i3 < 7; i3++) {
      ret[i3 + 7 * i2] = (*(real_T (*)[42])mxGetData(src))[i3 + 7 * i2];
    }
  }

  emlrtDestroyArray(&src);
}

static void r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T ret[7])
{
  static const int32_T dims[1] = { 7 };

  int32_T i4;
  emlrtCheckBuiltInR2012b(sp, msgId, src, "double", false, 1U, dims);
  for (i4 = 0; i4 < 7; i4++) {
    ret[i4] = (*(real_T (*)[7])mxGetData(src))[i4];
  }

  emlrtDestroyArray(&src);
}

void fastode45_MEE_api(const mxArray * const prhs[6], const mxArray *plhs[2])
{
  real_T (*y_out)[14];
  real_T (*tspan)[2];
  real_T (*b_y0)[14];
  real_T c;
  real_T Thr;
  real_T rho;
  real_T SI2CAN;
  const mxArray *tmp;
  real_T t_out;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  st.tls = emlrtRootTLSGlobal;
  y_out = (real_T (*)[14])mxMalloc(sizeof(real_T [14]));

  /* Marshall function inputs */
  tspan = emlrt_marshallIn(&st, emlrtAlias(prhs[0]), "tspan");
  b_y0 = c_emlrt_marshallIn(&st, emlrtAlias(prhs[1]), "y0");
  c = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[2]), "c");
  Thr = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[3]), "Thr");
  rho = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[4]), "rho");
  SI2CAN = e_emlrt_marshallIn(&st, emlrtAliasP(prhs[5]), "SI2CAN");

  /* Marshall in global variables */
  tmp = emlrtGetGlobalVariable("atol");
  if (tmp != NULL) {
    b_atol = e_emlrt_marshallIn(&st, tmp, "atol");
    atol_dirty = 0U;
  }

  tmp = emlrtGetGlobalVariable("rtol");
  if (tmp != NULL) {
    rtol = e_emlrt_marshallIn(&st, tmp, "rtol");
    rtol_dirty = 0U;
  }

  tmp = emlrtGetGlobalVariable("pow");
  if (tmp != NULL) {
    b_pow = e_emlrt_marshallIn(&st, tmp, "pow");
    pow_dirty = 0U;
  }

  tmp = emlrtGetGlobalVariable("A");
  if (tmp != NULL) {
    g_emlrt_marshallIn(&st, tmp, "A", A);
    A_dirty = 0U;
  }

  tmp = emlrtGetGlobalVariable("B");
  if (tmp != NULL) {
    i_emlrt_marshallIn(&st, tmp, "B", B);
    B_dirty = 0U;
  }

  tmp = emlrtGetGlobalVariable("E");
  if (tmp != NULL) {
    k_emlrt_marshallIn(&st, tmp, "E", E);
    E_dirty = 0U;
  }

  /* Invoke the target function */
  fastode45_MEE(&st, *tspan, *b_y0, c, Thr, rho, SI2CAN, &t_out, *y_out);

  /* Marshall out global variables */
  emlrtPutGlobalVariable("atol", emlrt_marshallOut(b_atol));
  emlrtPutGlobalVariable("rtol", emlrt_marshallOut(rtol));
  emlrtPutGlobalVariable("pow", emlrt_marshallOut(b_pow));
  emlrtPutGlobalVariable("A", b_emlrt_marshallOut(A));
  emlrtPutGlobalVariable("B", c_emlrt_marshallOut(B));
  emlrtPutGlobalVariable("E", d_emlrt_marshallOut(E));

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(t_out);
  plhs[1] = e_emlrt_marshallOut(*y_out);
}

/* End of code generation (_coder_fastode45_MEE_api.c) */
