typedef struct Kmer
{
  char* sequence;
} Kmer;

// In simple words, datum is like void * with additional size header and here we define macros.
#define DatumGetKmerP(X)  ((Kmer *) DatumGetPointer(X)) // We convert the datum pointer into a kmer pointer
#define KmerPGetDatum(X)  PointerGetDatum(X) // We covert the kmer pointer into a Datum pointer
#define PG_GETARG_KMER_P(n) DatumGetKmerP(PG_GETARG_DATUM(n)) // We get the nth argument given to a function
#define PG_RETURN_KMER_P(x) return KmerPGetDatum(x) // ¯\_(ツ)_/¯

static bool validate_kmer_sequence(const char *sequence) {
    if (sequence == NULL || *sequence == '\0') {
        ereport(ERROR, (errmsg("KMER sequence cannot be empty")));
        return false;
    }
	const char *length_bound = sequence + 32;
    for (const char *p = sequence; *p; p++) {
		if (p == length_bound) {
			ereport(ERROR, (errmsg("Kmer sequence length must not exceed 32", *p)));
			return false;
		}
        if (*p != 'A' && *p != 'T' && *p != 'C' && *p != 'G') {
            ereport(ERROR, (errmsg("Invalid character in kmer sequence: %c", *p)));
            return false;
        }
    }
    return true;
}

static Kmer * kmer_make(const char *sequence)
{
  Kmer *kmer = palloc0(sizeof(Kmer));
  if (sequence != NULL) {
    if (!validate_kmer_sequence(sequence)) {
         return NULL;  // Returns NULL if the sequence contains invalid characters
     }
    kmer->sequence = pstrdup(sequence);
  } else {
    kmer->sequence = pstrdup("");
  }
  return kmer;
}

static Kmer * kmer_parse(char **str) {
    p_whitespace(str); 
    char *sequence = *str; 
    char *end = sequence;
    while (*end && *end != ' ' && *end != '\n' && *end != '\r' && *end != '\t') {
        end++;
    }

    *end = '\0'; 
    Kmer *result = kmer_make(sequence);
    *str = end;
    ensure_end_input(str, true); 
    return result;
}

static char * kmer_to_str(const Kmer *kmer)
{
  char *result = psprintf("%s", kmer->sequence);
  return result;
}

PG_FUNCTION_INFO_V1(kmer_in);
Datum
kmer_in(PG_FUNCTION_ARGS)
{
  char *str = PG_GETARG_CSTRING(0);
  PG_RETURN_KMER_P(kmer_parse(&str));
}

PG_FUNCTION_INFO_V1(kmer_out);
Datum
kmer_out(PG_FUNCTION_ARGS)
{
  Kmer *kmer = PG_GETARG_KMER_P(0);
  char *result = kmer_to_str(kmer);
  PG_FREE_IF_COPY(kmer, 0);
  PG_RETURN_CSTRING(result);
}

PG_FUNCTION_INFO_V1(kmer_recv);
Datum
kmer_recv(PG_FUNCTION_ARGS)
{
  StringInfo  buf = (StringInfo) PG_GETARG_POINTER(0);
  Kmer *kmer = (Kmer *) palloc(sizeof(Kmer));
  kmer->sequence = pq_getmsgstring(buf);
  PG_RETURN_KMER_P(kmer);
}

PG_FUNCTION_INFO_V1(kmer_send);
Datum
kmer_send(PG_FUNCTION_ARGS)
{
  Kmer *kmer = PG_GETARG_KMER_P(0);
  StringInfoData buf;
  pq_begintypsend(&buf);
  pq_sendstring(&buf, kmer->sequence);
  PG_FREE_IF_COPY(kmer, 0);
  PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

PG_FUNCTION_INFO_V1(kmer_cast_from_text);
Datum
kmer_cast_from_text(PG_FUNCTION_ARGS)
{
    text *txt = PG_GETARG_TEXT_P(0);
    char *str = DatumGetCString(DirectFunctionCall1(textout, PointerGetDatum(txt)));
    Kmer *kmer = kmer_parse(&str);
    PG_RETURN_KMER_P(kmer);
}

PG_FUNCTION_INFO_V1(kmer_cast_to_text);
Datum
kmer_cast_to_text(PG_FUNCTION_ARGS)
{
  Kmer *kmer  = PG_GETARG_KMER_P(0);
  text *out = (text *)DirectFunctionCall1(textin,
            PointerGetDatum(kmer_to_str(kmer)));
  PG_FREE_IF_COPY(kmer, 0);
  PG_RETURN_TEXT_P(out);
}


/*****************************************************************************/
PG_FUNCTION_INFO_V1(kmer_constructor);
Datum
kmer_constructor(PG_FUNCTION_ARGS){
  char *sequence = PG_GETARG_CSTRING(0); 
  PG_RETURN_KMER_P(kmer_make(sequence));
}
/*****************************************************************************/

PG_FUNCTION_INFO_V1(kmer_to_string);
Datum
kmer_to_string(PG_FUNCTION_ARGS)
{
  Kmer *kmer = PG_GETARG_KMER_P(0);
  char *result = psprintf("%s", kmer->sequence);
  PG_FREE_IF_COPY(kmer, 0);
  PG_RETURN_CSTRING(result);           
}

static bool
kmer_eq_internal(Kmer *kmer1, Kmer *kmer2){
  return strcmp(kmer1->sequence, kmer2->sequence) == 0;
}

PG_FUNCTION_INFO_V1(kmer_eq);
Datum
kmer_eq(PG_FUNCTION_ARGS)
{
  Kmer *kmer1 = PG_GETARG_KMER_P(0);
  Kmer *kmer2 = PG_GETARG_KMER_P(1);
  bool result = kmer_eq_internal(kmer1, kmer2);
  PG_FREE_IF_COPY(kmer1, 0);
  PG_FREE_IF_COPY(kmer2, 1);
  PG_RETURN_BOOL(result);
}



static int
kmer_length_internal(Kmer *kmer)
{
  int length; 
  if (kmer->sequence == NULL){
    length = 0;
  }
  else{
  length = strlen(kmer->sequence);
  }
  return length;
}

PG_FUNCTION_INFO_V1(kmer_length);
Datum
kmer_length(PG_FUNCTION_ARGS)
{
  Kmer *kmer = PG_GETARG_KMER_P(0);
  int length = kmer_length_internal(kmer);
  PG_FREE_IF_COPY(kmer, 0);
  PG_RETURN_INT32(length); 
}


PG_FUNCTION_INFO_V1(kmer_ne);
Datum
kmer_ne(PG_FUNCTION_ARGS)
{
  Kmer *kmer1 = PG_GETARG_KMER_P(0);
  Kmer *kmer2 = PG_GETARG_KMER_P(1);
  bool result = !kmer_eq_internal(kmer1, kmer2);
  PG_FREE_IF_COPY(kmer1, 0);
  PG_FREE_IF_COPY(kmer2, 1);
  PG_RETURN_BOOL(result);
}

// We'll have to remove this function
static double
kmer_dist_internal(Kmer *kmer1, Kmer *kmer2)
{
  double result = 1;
  return result;
}




PG_FUNCTION_INFO_V1(kmer_dist);
Datum
kmer_dist(PG_FUNCTION_ARGS)
{
  Kmer *kmer1 = PG_GETARG_KMER_P(0);
  Kmer *kmer2 = PG_GETARG_KMER_P(1);
  double result = kmer_dist_internal(kmer1, kmer2);
  PG_FREE_IF_COPY(kmer1, 0);
  PG_FREE_IF_COPY(kmer2, 1);
  PG_RETURN_FLOAT8(result);
}
