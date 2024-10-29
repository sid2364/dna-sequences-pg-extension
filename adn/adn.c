#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "postgres.h"
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"

PG_MODULE_MAGIC;

typedef struct
{
  char* sequence;
} Adn;

// In simple words, datum is like void * with additional size header and here we define macros.
#define DatumGetAdnP(X)  ((Adn *) DatumGetPointer(X)) // We convert the datum pointer into a adn pointer
#define AdnPGetDatum(X)  PointerGetDatum(X) // We covert the adn pointer into a Datum pointer
#define PG_GETARG_ADN_P(n) DatumGetAdnP(PG_GETARG_DATUM(n)) // We get the nth argument given to a function
#define PG_RETURN_ADN_P(x) return AdnPGetDatum(x) // ¯\_(ツ)_/¯

static Adn * adn_make(const char *sequence)
{
  Adn *adn = palloc0(sizeof(Adn));        
  if (sequence != NULL){
    adn->sequence = pstrdup(sequence);     
  }else{
    adn->sequence = pstrdup("");           
  
} return adn;
}


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

static Adn * adn_parse(char **str) {
    p_whitespace(str); 
    char *sequence = *str; 
    char *end = sequence;
    while (*end && *end != ' ' && *end != '\n' && *end != '\r' && *end != '\t') {
        end++;
    }

    *end = '\0'; 
    Adn *result = adn_make(sequence);
    *str = end;
    ensure_end_input(str, true); 
    return result;
}

static char * adn_to_str(const Adn *adn)
{
  char *result = psprintf("%s", adn->sequence);
  return result;
}


PG_FUNCTION_INFO_V1(adn_in);
Datum
adn_in(PG_FUNCTION_ARGS)
{
  char *str = PG_GETARG_CSTRING(0);
  PG_RETURN_ADN_P(adn_parse(&str));
}

PG_FUNCTION_INFO_V1(adn_out);
Datum
adn_out(PG_FUNCTION_ARGS)
{
  Adn *adn = PG_GETARG_ADN_P(0);
  char *result = adn_to_str(adn);
  PG_FREE_IF_COPY(adn, 0);
  PG_RETURN_CSTRING(result);
}

PG_FUNCTION_INFO_V1(adn_recv);
Datum
adn_recv(PG_FUNCTION_ARGS)
{
  StringInfo  buf = (StringInfo) PG_GETARG_POINTER(0);
  Adn *adn = (Adn *) palloc(sizeof(Adn));
  adn->sequence = pq_getmsgstring(buf);
  PG_RETURN_ADN_P(adn);
}

PG_FUNCTION_INFO_V1(adn_send);
Datum
adn_send(PG_FUNCTION_ARGS)
{
  Adn *adn = PG_GETARG_ADN_P(0);
  StringInfoData buf;
  pq_begintypsend(&buf);
  pq_sendstring(&buf, adn->sequence);
  PG_FREE_IF_COPY(adn, 0);
  PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

PG_FUNCTION_INFO_V1(adn_cast_from_text);
Datum
adn_cast_from_text(PG_FUNCTION_ARGS)
{
    text *txt = PG_GETARG_TEXT_P(0);
    char *str = DatumGetCString(DirectFunctionCall1(textout, PointerGetDatum(txt)));
    Adn *adn = adn_parse(&str);
    PG_RETURN_ADN_P(adn);
}

PG_FUNCTION_INFO_V1(adn_cast_to_text);
Datum
adn_cast_to_text(PG_FUNCTION_ARGS)
{
  Adn *adn  = PG_GETARG_ADN_P(0);
  text *out = (text *)DirectFunctionCall1(textin,
            PointerGetDatum(adn_to_str(adn)));
  PG_FREE_IF_COPY(adn, 0);
  PG_RETURN_TEXT_P(out);
}


/*****************************************************************************/
PG_FUNCTION_INFO_V1(adn_constructor);
Datum
adn_constructor(PG_FUNCTION_ARGS){
  char *sequence = PG_GETARG_CSTRING(0); 
  PG_RETURN_ADN_P(adn_make(sequence));   
}
/*****************************************************************************/

PG_FUNCTION_INFO_V1(adn_to_string);
Datum
adn_to_string(PG_FUNCTION_ARGS)
{
  Adn *adn = PG_GETARG_ADN_P(0);        
  char *result = psprintf("%s", adn->sequence); 
  PG_FREE_IF_COPY(adn, 0);             
  PG_RETURN_CSTRING(result);           
}

static bool
adn_eq_internal(Adn *adn1, Adn *adn2){
  return strcmp(adn1->sequence, adn2->sequence) == 0;
}

PG_FUNCTION_INFO_V1(adn_eq);
Datum
adn_eq(PG_FUNCTION_ARGS)
{
  Adn *adn1 = PG_GETARG_ADN_P(0);
  Adn *adn2 = PG_GETARG_ADN_P(1);
  bool result = adn_eq_internal(adn1, adn2);
  PG_FREE_IF_COPY(adn1, 0);
  PG_FREE_IF_COPY(adn2, 1);
  PG_RETURN_BOOL(result);
}


PG_FUNCTION_INFO_V1(adn_ne);
Datum
adn_ne(PG_FUNCTION_ARGS)
{
  Adn *adn1 = PG_GETARG_ADN_P(0);
  Adn *adn2 = PG_GETARG_ADN_P(1);
  bool result = !adn_eq_internal(adn1, adn2);
  PG_FREE_IF_COPY(adn1, 0);
  PG_FREE_IF_COPY(adn2, 1);
  PG_RETURN_BOOL(result);
}


static double
adn_dist_internal(Adn *adn1, Adn *adn2)
{
  double result = 1;
  return result;
}

PG_FUNCTION_INFO_V1(adn_dist);
Datum
adn_dist(PG_FUNCTION_ARGS)
{
  Adn *adn1 = PG_GETARG_ADN_P(0);
  Adn *adn2 = PG_GETARG_ADN_P(1);
  double result = adn_dist_internal(adn1, adn2);
  PG_FREE_IF_COPY(adn1, 0);
  PG_FREE_IF_COPY(adn2, 1);
  PG_RETURN_FLOAT8(result);
}
