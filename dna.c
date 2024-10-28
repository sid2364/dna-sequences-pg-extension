/*
 * dna.C
 *
 * PostgreSQL DNA Sequences Type:
 *
 * dna '(ACGT)'
 *
 */

#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "postgres.h"
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"

PG_MODULE_MAGIC;

/**
* TODO: This entire file!
*/

#define EPSILON         1.0E-06

#define FPzero(A)       (fabs(A) <= EPSILON)
#define FPeq(A,B)       (fabs((A) - (B)) <= EPSILON)
#define FPne(A,B)       (fabs((A) - (B)) > EPSILON)
#define FPlt(A,B)       ((B) - (A) > EPSILON)
#define FPle(A,B)       ((A) - (B) <= EPSILON)
#define FPgt(A,B)       ((A) - (B) > EPSILON)
#define FPge(A,B)       ((B) - (A) <= EPSILON)

/*****************************************************************************/

/* Structure to represent complex numbers */
typedef struct
{
  double    a,
            b;
} Complex;

/* fmgr macros complex type */

#define DatumGetComplexP(X)  ((Complex *) DatumGetPointer(X))
#define ComplexPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_COMPLEX_P(n) DatumGetComplexP(PG_GETARG_DATUM(n))
#define PG_RETURN_COMPLEX_P(x) return ComplexPGetDatum(x)

/*****************************************************************************/

static Complex *
complex_make(double a, double b)
{
  Complex *c = palloc0(sizeof(Complex));
  c->a = a;
  c->b = b;
  /* this is not superflous code, we need
  to account for negative zeroes here */
  if (c->a == 0)
    c->a = 0;
  if (c->b == 0)
    c->b = 0;
  return c;
}

/*****************************************************************************/

static void
p_whitespace(char **str)
{
  while (**str == ' ' || **str == '\n' || **str == '\r' || **str == '\t')
    *str += 1;
}

static void
ensure_end_input(char **str, bool end)
{
  if (end)
  {
    p_whitespace(str);
    if (**str != 0)
      ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
        errmsg("Could not parse temporal value")));
  }
}

static bool
p_oparen(char **str)
{
  p_whitespace(str);
  if (**str == '(')
  {
    *str += 1;
    return true;
  }
  return false;
}

static bool
p_cparen(char **str)
{
  p_whitespace(str);
  if (**str == ')')
  {
    *str += 1;
    return true;
  }
  return false;
}

static bool
p_comma(char **str)
{
  p_whitespace(str);
  if (**str == ',')
  {
    *str += 1;
    return true;
  }
  return false;
}

static double
double_parse(char **str)
{
  char *nextstr = *str;
  double result = strtod(*str, &nextstr);
  if (*str == nextstr)
    ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
      errmsg("Invalid input syntax for type double")));
  *str = nextstr;
  return result;
}

static Complex *
complex_parse(char **str)
{
  double a, b;
  if (!p_oparen(str))
    ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
      errmsg("Invalid input syntax for type complex")));
  a = double_parse(str);
  if (!p_comma(str))
    ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
      errmsg("Invalid input syntax for type complex")));
  b = double_parse(str);
  if (!p_cparen(str))
    ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
      errmsg("Invalid input syntax for type complex")));
  ensure_end_input(str, true);
  return complex_make(a, b);
}

static char *
complex_to_str(const Complex *c)
{
  char *result = psprintf("(%.*g, %.*g)",
    DBL_DIG, c->a,
    DBL_DIG, c->b);
  return result;
}

/*****************************************************************************/

PG_FUNCTION_INFO_V1(complex_in);
Datum
complex_in(PG_FUNCTION_ARGS)
{
  char *str = PG_GETARG_CSTRING(0);
  PG_RETURN_COMPLEX_P(complex_parse(&str));
}

PG_FUNCTION_INFO_V1(complex_out);
Datum
complex_out(PG_FUNCTION_ARGS)
{
  Complex *c = PG_GETARG_COMPLEX_P(0);
  char *result = complex_to_str(c);
  PG_FREE_IF_COPY(c, 0);
  PG_RETURN_CSTRING(result);
}

PG_FUNCTION_INFO_V1(complex_recv);
Datum
complex_recv(PG_FUNCTION_ARGS)
{
  StringInfo  buf = (StringInfo) PG_GETARG_POINTER(0);
  Complex *c = (Complex *) palloc(sizeof(Complex));
  c->a = pq_getmsgfloat8(buf);
  c->b = pq_getmsgfloat8(buf);
  PG_RETURN_COMPLEX_P(c);
}

PG_FUNCTION_INFO_V1(complex_send);
Datum
complex_send(PG_FUNCTION_ARGS)
{
  Complex *c = PG_GETARG_COMPLEX_P(0);
  StringInfoData buf;
  pq_begintypsend(&buf);
  pq_sendfloat8(&buf, c->a);
  pq_sendfloat8(&buf, c->b);
  PG_FREE_IF_COPY(c, 0);
  PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

PG_FUNCTION_INFO_V1(complex_cast_from_text);
Datum
complex_cast_from_text(PG_FUNCTION_ARGS)
{
  text *txt = PG_GETARG_TEXT_P(0);
  char *str = DatumGetCString(DirectFunctionCall1(textout,
               PointerGetDatum(txt)));
  PG_RETURN_COMPLEX_P(complex_parse(&str));
}

PG_FUNCTION_INFO_V1(complex_cast_to_text);
Datum
complex_cast_to_text(PG_FUNCTION_ARGS)
{
  Complex *c  = PG_GETARG_COMPLEX_P(0);
  text *out = (text *)DirectFunctionCall1(textin,
            PointerGetDatum(complex_to_str(c)));
  PG_FREE_IF_COPY(c, 0);
  PG_RETURN_TEXT_P(out);
}

/*****************************************************************************/

PG_FUNCTION_INFO_V1(complex_constructor);
Datum
complex_constructor(PG_FUNCTION_ARGS)
{
  double a = PG_GETARG_FLOAT8(0);
  double b = PG_GETARG_FLOAT8(1);
  PG_RETURN_COMPLEX_P(complex_make(a, b));
}

/*****************************************************************************/

PG_FUNCTION_INFO_V1(complex_re);
Datum
complex_re(PG_FUNCTION_ARGS)
{
  Complex *c = PG_GETARG_COMPLEX_P(0);
  double result = c->a;
  PG_FREE_IF_COPY(c, 0);
  PG_RETURN_FLOAT8(result);
}

PG_FUNCTION_INFO_V1(complex_im);
Datum
complex_im(PG_FUNCTION_ARGS)
{
  Complex *c = PG_GETARG_COMPLEX_P(0);
  double result = c->b;
  PG_FREE_IF_COPY(c, 0);
  PG_RETURN_FLOAT8(result);
}

/* Complete the following code by filling in the dots (...) 
 * This function returns the conjugate of a complex number
 */
PG_FUNCTION_INFO_V1(complex_conj);
Datum
complex_conj(PG_FUNCTION_ARGS)
{
  Complex *c = PG_GETARG_COMPLEX_P(0);
  double negative_b = c->b;
  PG_RETURN_COMPLEX_P(complex_make(c->a, -negative_b));
}

/*****************************************************************************/

static bool
complex_eq_internal(Complex *c, Complex *d)
{
  return FPeq(c->a, d->a) && FPeq(c->b, d->b);
}

PG_FUNCTION_INFO_V1(complex_eq);
Datum
complex_eq(PG_FUNCTION_ARGS)
{
  Complex *c = PG_GETARG_COMPLEX_P(0);
  Complex *d = PG_GETARG_COMPLEX_P(1);
  bool result = complex_eq_internal(c, d);
  PG_FREE_IF_COPY(c, 0);
  PG_FREE_IF_COPY(d, 1);
  PG_RETURN_BOOL(result);
}

PG_FUNCTION_INFO_V1(complex_ne);
Datum
complex_ne(PG_FUNCTION_ARGS)
{
  Complex *c = PG_GETARG_COMPLEX_P(0);
  Complex *d = PG_GETARG_COMPLEX_P(1);
  bool result = !complex_eq_internal(c, d);
  PG_FREE_IF_COPY(c, 0);
  PG_FREE_IF_COPY(d, 1);
  PG_RETURN_BOOL(result);
}

PG_FUNCTION_INFO_V1(complex_left);
Datum
complex_left(PG_FUNCTION_ARGS)
{
  Complex *c = PG_GETARG_COMPLEX_P(0);
  Complex *d = PG_GETARG_COMPLEX_P(1);
  bool result = FPlt(c->a, d->a);
  PG_FREE_IF_COPY(c, 0);
  PG_FREE_IF_COPY(d, 1);
  PG_RETURN_BOOL(result);
}

/* Complete the following code by filling in the dots (...) 
 * This is the function for the strictly right '>>' operator
 */
/*
PG_FUNCTION_INFO_V1(complex_right);
Datum
complex_right(PG_FUNCTION_ARGS)
{
  Complex *c = PG_GETARG_COMPLEX_P(0);
  Complex *d = PG_GETARG_COMPLEX_P(1);
  ...
  PG_FREE_IF_COPY(c, 0);
  PG_FREE_IF_COPY(d, 1);
  PG_RETURN_BOOL(...);
}
*/

PG_FUNCTION_INFO_V1(complex_below);
Datum
complex_below(PG_FUNCTION_ARGS)
{
  Complex *c = PG_GETARG_COMPLEX_P(0);
  Complex *d = PG_GETARG_COMPLEX_P(1);
  bool result = FPlt(c->b, d->b);
  PG_FREE_IF_COPY(c, 0);
  PG_FREE_IF_COPY(d, 1);
  PG_RETURN_BOOL(result);
}

/* Write the code for the following function 
 * This is the function for the strictly above '|>>' operator
 */
/*
PG_FUNCTION_INFO_V1(complex_above);
Datum
complex_above(PG_FUNCTION_ARGS)
{
  ...
}
*/

/*****************************************************************************/

PG_FUNCTION_INFO_V1(complex_add);
Datum
complex_add(PG_FUNCTION_ARGS)
{
  Complex *c = PG_GETARG_COMPLEX_P(0);
  Complex *d = PG_GETARG_COMPLEX_P(1);
  Complex *result = complex_make(c->a + d->a, c->b + d->b);
  PG_FREE_IF_COPY(c, 0);
  PG_FREE_IF_COPY(d, 1);
  PG_RETURN_COMPLEX_P(result);
}

/* Complete the following code by filling in the dots (...) 
 * This is the function for the subtract '-' operator
 */
/*
PG_FUNCTION_INFO_V1(complex_sub);
Datum
complex_sub(PG_FUNCTION_ARGS)
{
  ...
}
*/

/* Write the code for the following function 
 * This is the function for the multiplication '*' operator
 */
/*
...
*/

PG_FUNCTION_INFO_V1(complex_div);
Datum
complex_div(PG_FUNCTION_ARGS)
{
  Complex *c = PG_GETARG_COMPLEX_P(0);
  Complex *d = PG_GETARG_COMPLEX_P(1);
  double div;
  Complex *result;
  if (FPzero(d->a) && FPzero(d->b))
  {
    ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE),
      errmsg("Can only divide by a non-zero complex number")));
  }
  div = (d->a * d->a) + (d->b * d->b);
  result = complex_make(
    (c->a*d->a + c->b*d->b) / div, 
    (c->b*d->a - c->a*d->b) / div);
  PG_FREE_IF_COPY(c, 0);
  PG_FREE_IF_COPY(d, 1);
  PG_RETURN_COMPLEX_P(result);
}

/*****************************************************************************/

static double
complex_dist_internal(Complex *c, Complex *d)
{
  double x = fabs(c->a - d->a);
  double y = fabs(c->b - d->b);
  double yx, result;
  if (x < y)
  {
    double temp = x;
    x = y;
    y = temp;
  }
  if (FPzero(y))
    return x;
  yx = y / x;
  result = sqrt(1.0 + (yx * yx));
  return result;
}

PG_FUNCTION_INFO_V1(complex_dist);
Datum
complex_dist(PG_FUNCTION_ARGS)
{
  Complex *c = PG_GETARG_COMPLEX_P(0);
  Complex *d = PG_GETARG_COMPLEX_P(1);
  double result = complex_dist_internal(c, d);
  PG_FREE_IF_COPY(c, 0);
  PG_FREE_IF_COPY(d, 1);
  PG_RETURN_FLOAT8(result);
}

/*****************************************************************************/