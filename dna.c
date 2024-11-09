#include "postgres.h"
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdint.h>

#include "utils/varlena.h"
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"

PG_MODULE_MAGIC;

/*
Great reference for pg functions: https://doxygen.postgresql.org/pqformat_8c.html
*/

typedef struct Dna
{
    int32 vl_len_;                    // Required header for PostgreSQL variable-length types
    int32 length;                     // Length of the DNA sequence in nucleotides
    uint64_t bit_sequence[FLEXIBLE_ARRAY_MEMBER];  // Array to store packed bits
} Dna;

// In simple words, datum is like void * with additional size header and here we define macros.
#define DatumGetDnaP(X)  ((Dna *) DatumGetPointer(X)) // We convert the datum pointer into a dna pointer
#define DnaPGetDatum(X)  PointerGetDatum(X) // We covert the dna pointer into a Datum pointer
#define PG_GETARG_DNA_P(n) DatumGetDnaP(PG_GETARG_DATUM(n)) // We get the nth argument given to a function
#define PG_RETURN_DNA_P(x) return DnaPGetDatum(x) // ¯\_(ツ)_/¯

/**
Encoding function

Calculate number of 64-bit chunks we need, since rest of the int will be padded with zeros,
 we also store the length so that we can decode it later properly without decoding extra "00"s as "A"s
*/
static void encode_dna(const char *sequence, uint64_t *bit_sequence, int length) {
    int bit_length = (length + 31) / 32;  // Number of 64-bit chunks we need

    for (int i = 0; i < length; i++) {
        int offset = (i * 2) % 64;
        int index = i / 32;

        switch (sequence[i]) {
            case 'A': /* 00 */ break;
            case 'T': bit_sequence[index] |= ((uint64_t)0x1 << offset); break; // 01
            case 'C': bit_sequence[index] |= ((uint64_t)0x2 << offset); break; // 10
            case 'G': bit_sequence[index] |= ((uint64_t)0x3 << offset); break; // 11
            default:
                ereport(ERROR, (errmsg("Invalid character in DNA sequence: %c", sequence[i])));
        }
    }
}


/**
Decoding function

We decode the sequence by shifting the bits to the right and then reading the last 2 bits to get the nucleotide
*/
static char* decode_dna(const uint64_t *bit_sequence, int length) {
    char *sequence = palloc0(length + 1);  // +1 for the null terminator which is added automatically by palloc0

    for (int i = 0; i < length; i++) {
        int offset = (i * 2) % 64; // Offset of the current nucleotide since it's 2 bits long
        int index = i / 32;
        uint64_t bits = (bit_sequence[index] >> offset) & 0x3;

        switch (bits) {
            case 0x0: sequence[i] = 'A'; break;
            case 0x1: sequence[i] = 'T'; break;
            case 0x2: sequence[i] = 'C'; break;
            case 0x3: sequence[i] = 'G'; break;
        }
    }

    return sequence;
}

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

/*
Creates and returns a new Dna struct by encoding the provided DNA sequence string "ATCG" into binary format (2 bits per nucleotide)

It checks the input, calculates required memory, and calls encode_dna which stores the enocded sequence in bit_sequence
 */
static Dna * dna_make(const char *sequence)
{
    int length = strlen(sequence);

    if (sequence != NULL) {
        if (!validate_dna_sequence(sequence)) {
            ereport(ERROR, (errmsg("Invalid DNA sequence: must contain only A, T, C, G")));
            return NULL;
        }
    }

    int bit_length = (length + 31) / 32;  // Number of 64-bit chunks we need, rest will be padded with zeros
    int size = VARHDRSZ + sizeof(int32) + bit_length * sizeof(uint64_t);  // Total size with varlena header

    // Allocate memory for Dna struct and bit_sequence
    Dna *dna = (Dna *) palloc0(size);
    dna->vl_len_ = size;  // Manually set the size in the varlena header, SET_VARSIZE is somehow not working even with the import
    dna->length = length;

    // Encode the DNA sequence directly into bit_sequence, pointer magic
    encode_dna(sequence, dna->bit_sequence, length);

    return dna;
}

// TODO: Check if we really need all this
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
    if (dna->bit_sequence == NULL || dna->length == 0) {
        return pstrdup("");
    }
    return decode_dna(dna->bit_sequence, dna->length);
}

PG_FUNCTION_INFO_V1(dna_in);
Datum
dna_in(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);
    Dna *dna = dna_make(str);  // Use dna_make to create the encoded binary sequence

    PG_RETURN_DNA_P(dna);
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

/*
This function is supposed to take in an existing DNA sequence and return a new DNA sequence with the same values!
An existing sequence means it's a binary encoded sequence in String format; just Postgres things!
*/
PG_FUNCTION_INFO_V1(dna_recv);
Datum
dna_recv(PG_FUNCTION_ARGS)
{
    StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);

    int length = pq_getmsgint(buf, sizeof(int));
    int bit_length = (length + 31) / 32;
    int size = sizeof(Dna) + bit_length * sizeof(uint64_t);

    Dna *dna = (Dna *) palloc0(size);
    dna->vl_len_ = size;  // Similar to dna_make
    dna->length = length;

    for (int i = 0; i < bit_length; i++) {
        dna->bit_sequence[i] = pq_getmsgint64(buf);
    }

    PG_RETURN_POINTER(dna);
}

/*
Does the same but in reverse, takes a DNA sequence and encodes it into a binary format

We can't use dna_to_str here because
dna_to_str converts the binary data back to a string format which is not what we want here!
*/
PG_FUNCTION_INFO_V1(dna_send);
Datum
dna_send(PG_FUNCTION_ARGS)
{
    Dna *dna = (Dna *) PG_GETARG_POINTER(0);
    StringInfoData buf;
    pq_begintypsend(&buf);

    pq_sendint(&buf, dna->length, sizeof(int));

    int bit_length = (dna->length + 31) / 32;
    for (int i = 0; i < bit_length; i++) {
        pq_sendint64(&buf, dna->bit_sequence[i]);
    }

    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}


PG_FUNCTION_INFO_V1(dna_cast_from_text);
Datum
dna_cast_from_text(PG_FUNCTION_ARGS)
{
    text *txt = PG_GETARG_TEXT_P(0);
    char *str = DatumGetCString(DirectFunctionCall1(textout, PointerGetDatum(txt)));
    Dna *dna = dna_make(str);  // Use dna_make for binary encoding
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
    char *result = decode_dna(dna->bit_sequence, dna->length);  // Decode bit_sequence to a readable string
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_CSTRING(result);
}

/**
Magically faster than strcmp!
*/
static bool dna_eq_internal(Dna *dna1, Dna *dna2)
{
    if (dna1->length != dna2->length) {
        return false;  // Different lengths mean they can't be equal
    }

    int bit_length = (dna1->length + 31) / 32;  // Number of 64-bit chunks

    // Compare each 64-bit chunk in the bit_sequence array
    for (int i = 0; i < bit_length; i++) {
        if (dna1->bit_sequence[i] != dna2->bit_sequence[i]) {
            return false;  // If any chunk differs, the sequences are not equal
        }
    }

    return true;  // All chunks are equal, phew!
}


PG_FUNCTION_INFO_V1(equals);
Datum
equals(PG_FUNCTION_ARGS)
{
  Dna *dna1 = PG_GETARG_DNA_P(0);
  Dna *dna2 = PG_GETARG_DNA_P(1);
  bool result = dna_eq_internal(dna1, dna2);
  PG_FREE_IF_COPY(dna1, 0);
  PG_FREE_IF_COPY(dna2, 1);
  PG_RETURN_BOOL(result);
}

static int dna_length_internal(Dna *dna)
{
    return dna->length;  // We store the length already
}

PG_FUNCTION_INFO_V1(length);
Datum
length(PG_FUNCTION_ARGS)
{
    Dna *dna = PG_GETARG_DNA_P(0);
    int length = dna->length;  // Directly get the length field
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_INT32(length);
}

PG_FUNCTION_INFO_V1(dna_ne);
Datum
dna_ne(PG_FUNCTION_ARGS)
{
    Dna *dna1 = PG_GETARG_DNA_P(0);
    Dna *dna2 = PG_GETARG_DNA_P(1);
    bool result = !dna_eq_internal(dna1, dna2);  // ~ the result of dna_eq_internal, viola!
    PG_FREE_IF_COPY(dna1, 0);
    PG_FREE_IF_COPY(dna2, 1);
    PG_RETURN_BOOL(result);
}
/*
TODO: Figure out how to implement distance functions (or if they are even needed)
*/
static double dna_dist_internal(Dna *dna1, Dna *dna2)
{
    return 1.0;
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