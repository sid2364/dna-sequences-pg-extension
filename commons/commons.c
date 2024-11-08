#include "postgres.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"

PG_MODULE_MAGIC;

static void p_whitespace(char **str)
{
  while (**str == ' ' || **str == '\n' || **str == '\r' || **str == '\t')
    *str += 1;
}

static void ensure_end_input(char **str, bool end)
{
  if (end)
  {
    p_whitespace(str);
    if (**str != 0)
      ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
        errmsg("Could not parse temporal value")));
  }
}

static double double_parse(char **str)
{
  char *nextstr = *str;
  double result = strtod(*str, &nextstr);
  if (*str == nextstr)
    ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
      errmsg("Invalid input syntax for type double")));
  *str = nextstr;
  return result;
}

