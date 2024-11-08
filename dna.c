#include "commons/commons.c"
#include "kmer/kmer.c"

typedef struct Dna
{
  char* sequence;
} Dna;

// In simple words, datum is like void * with additional size header and here we define macros.
#define DatumGetDnaP(X)  ((Dna *) DatumGetPointer(X)) // We convert the datum pointer into a dna pointer
#define DnaPGetDatum(X)  PointerGetDatum(X) // We covert the dna pointer into a Datum pointer
#define PG_GETARG_DNA_P(n) DatumGetDnaP(PG_GETARG_DATUM(n)) // We get the nth argument given to a function
#define PG_RETURN_DNA_P(x) return DnaPGetDatum(x) // ¯\_(ツ)_/¯

static bool validate_dna_sequence(const char *sequence) {
    if (sequence == NULL || *sequence == '\0') {
        ereport(ERROR, (errmsg("DNA sequence cannot be empty")));
        return false;
    }
    for (const char *p = sequence; *p; p++) {
        if (*p != 'A' && *p != 'T' && *p != 'C' && *p != 'G') {
            ereport(ERROR, (errmsg("Invalid character in DNA sequence: %c", *p)));
            return false;
        }
    }
    return true;
}

static Dna * dna_make(const char *sequence)
{
  Dna *dna = palloc0(sizeof(Dna));
  if (sequence != NULL) {
    if (!validate_dna_sequence(sequence)) {
         return NULL;  // Returns NULL if the sequence contains invalid characters
     }
    dna->sequence = pstrdup(sequence);
  } else {
    dna->sequence = pstrdup("");
  }
  return dna;
}


static Dna * dna_parse(char **str) {
    p_whitespace(str); 
    char *sequence = *str; 
    char *end = sequence;
    while (*end && *end != ' ' && *end != '\n' && *end != '\r' && *end != '\t') {
        end++;
    }

    *end = '\0'; 
    Dna *result = dna_make(sequence);
    *str = end;
    ensure_end_input(str, true); 
    return result;
}

static char * dna_to_str(const Dna *dna)
{
  char *result = psprintf("%s", dna->sequence);
  return result;
}

PG_FUNCTION_INFO_V1(dna_in);
Datum
dna_in(PG_FUNCTION_ARGS)
{
  char *str = PG_GETARG_CSTRING(0);
  PG_RETURN_DNA_P(dna_parse(&str));
}

PG_FUNCTION_INFO_V1(dna_out);
Datum
dna_out(PG_FUNCTION_ARGS)
{
  Dna *dna = PG_GETARG_DNA_P(0);
  char *result = dna_to_str(dna);
  PG_FREE_IF_COPY(dna, 0);
  PG_RETURN_CSTRING(result);
}

PG_FUNCTION_INFO_V1(dna_recv);
Datum
dna_recv(PG_FUNCTION_ARGS)
{
  StringInfo  buf = (StringInfo) PG_GETARG_POINTER(0);
  Dna *dna = (Dna *) palloc(sizeof(Dna));
  dna->sequence = pq_getmsgstring(buf);
  PG_RETURN_DNA_P(dna);
}

PG_FUNCTION_INFO_V1(dna_send);
Datum
dna_send(PG_FUNCTION_ARGS)
{
  Dna *dna = PG_GETARG_DNA_P(0);
  StringInfoData buf;
  pq_begintypsend(&buf);
  pq_sendstring(&buf, dna->sequence);
  PG_FREE_IF_COPY(dna, 0);
  PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

PG_FUNCTION_INFO_V1(dna_cast_from_text);
Datum
dna_cast_from_text(PG_FUNCTION_ARGS)
{
    text *txt = PG_GETARG_TEXT_P(0);
    char *str = DatumGetCString(DirectFunctionCall1(textout, PointerGetDatum(txt)));
    Dna *dna = dna_parse(&str);
    PG_RETURN_DNA_P(dna);
}

PG_FUNCTION_INFO_V1(dna_cast_to_text);
Datum
dna_cast_to_text(PG_FUNCTION_ARGS)
{
  Dna *dna  = PG_GETARG_DNA_P(0);
  text *out = (text *)DirectFunctionCall1(textin,
            PointerGetDatum(dna_to_str(dna)));
  PG_FREE_IF_COPY(dna, 0);
  PG_RETURN_TEXT_P(out);
}


/*****************************************************************************/
PG_FUNCTION_INFO_V1(dna_constructor);
Datum
dna_constructor(PG_FUNCTION_ARGS){
  char *sequence = PG_GETARG_CSTRING(0); 
  PG_RETURN_DNA_P(dna_make(sequence));
}
/*****************************************************************************/

PG_FUNCTION_INFO_V1(dna_to_string);
Datum
dna_to_string(PG_FUNCTION_ARGS)
{
  Dna *dna = PG_GETARG_DNA_P(0);
  char *result = psprintf("%s", dna->sequence);
  PG_FREE_IF_COPY(dna, 0);
  PG_RETURN_CSTRING(result);           
}

static bool
dna_eq_internal(Dna *dna1, Dna *dna2){
  return strcmp(dna1->sequence, dna2->sequence) == 0;
}

PG_FUNCTION_INFO_V1(equals);
Datum
dna_eq(PG_FUNCTION_ARGS)
{
  Dna *dna1 = PG_GETARG_DNA_P(0);
  Dna *dna2 = PG_GETARG_DNA_P(1);
  bool result = dna_eq_internal(dna1, dna2);
  PG_FREE_IF_COPY(dna1, 0);
  PG_FREE_IF_COPY(dna2, 1);
  PG_RETURN_BOOL(result);
}



static int
dna_length_internal(Dna *dna)
{
  int length; 
  if (dna->sequence == NULL){
    length = 0;
  }
  else{
  length = strlen(dna->sequence);
  }
  return length;
}

PG_FUNCTION_INFO_V1(dna_length);
Datum
dna_length(PG_FUNCTION_ARGS)
{
  Dna *dna = PG_GETARG_DNA_P(0);
  int length = dna_length_internal(dna);
  PG_FREE_IF_COPY(dna, 0);
  PG_RETURN_INT32(length); 
}


PG_FUNCTION_INFO_V1(dna_ne);
Datum
dna_ne(PG_FUNCTION_ARGS)
{
  Dna *dna1 = PG_GETARG_DNA_P(0);
  Dna *dna2 = PG_GETARG_DNA_P(1);
  bool result = !dna_eq_internal(dna1, dna2);
  PG_FREE_IF_COPY(dna1, 0);
  PG_FREE_IF_COPY(dna2, 1);
  PG_RETURN_BOOL(result);
}

// We'll have to remove this function
static double
dna_dist_internal(Dna *dna1, Dna *dna2)
{
  double result = 1;
  return result;
}




PG_FUNCTION_INFO_V1(dna_dist);
Datum
dna_dist(PG_FUNCTION_ARGS)
{
  Dna *dna1 = PG_GETARG_DNA_P(0);
  Dna *dna2 = PG_GETARG_DNA_P(1);
  double result = dna_dist_internal(dna1, dna2);
  PG_FREE_IF_COPY(dna1, 0);
  PG_FREE_IF_COPY(dna2, 1);
  PG_RETURN_FLOAT8(result);
}
