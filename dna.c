#include "postgres.h"
#include "utils/varlena.h"
#include "varatt.h" // For SET_VARSIZE!
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"
#include <stddef.h>  // Include for offsetof
#include "funcapi.h"
#include "utils/builtins.h" // For cstring_to_text

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdint.h>

PG_MODULE_MAGIC;

/*
Great reference for pg functions:
https://doxygen.postgresql.org/pqformat_8c.html
https://www.postgresql.org/docs/9.5/xfunc-c.html
*/

typedef struct Dna
{
    char vl_len_[4];                    // Required header for PostgreSQL variable-length types
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
    for (int i = 0; i < length; i++) {
        int offset = (i * 2) % 64;      // Offset for 2 bits per base
        int index = i / 32;             // Each uint64_t holds 32 bases (64 bits)

        switch (sequence[i]) {
            case 'A': /* 00 */ break;
            case 'T': bit_sequence[index] |= ((uint64_t)0x1 << offset); break; // 01 for T
            case 'C': bit_sequence[index] |= ((uint64_t)0x2 << offset); break; // 10 for C
            case 'G': bit_sequence[index] |= ((uint64_t)0x3 << offset); break; // 11 for G
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

    int num_bits = length * 2;  // 2 bits per nucleotide
    int bit_length = (num_bits + 63) / 64;  // Number of 64-bit chunks we need, rest will be padded with zeros
    Size dna_size = offsetof(Dna, bit_sequence) + bit_length * sizeof(uint64_t);

    // Allocate memory for Dna struct and bit_sequence, set all bits to 0 (which is why we use palloc0 and not palloc)
    Dna *dna = (Dna *) palloc0(dna_size);

    if (sequence != NULL) {
        if (!validate_dna_sequence(sequence)) {
            ereport(ERROR, (errmsg("Invalid DNA sequence: must contain only A, T, C, G")));
            return NULL;
        }
    }

    SET_VARSIZE(dna, dna_size); // No need to add VARHDRSZ since the library does it for us!
    //dna->vl_len_ = VARHDRSZ + dna_size;
    dna->length = length;

    // Encode the DNA sequence directly into bit_sequence, pointer magic
    encode_dna(sequence, dna->bit_sequence, length);

    return dna;
}

static char * dna_to_str(const Dna *dna)
{
    if (dna->length == 0) {
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
    int num_bits = length * 2;  // Each DNA nucleotide is 2 bits
    int bit_length = (num_bits + 63) / 64;  // Number of 64-bit chunks needed, rounding up

    Size dna_size = offsetof(Dna, bit_sequence) + bit_length * sizeof(uint64_t);

    Dna *dna = (Dna *) palloc0(dna_size);
    SET_VARSIZE(dna, dna_size);
    //dna->vl_len_ = VARHDRSZ + dna_size;
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

    int bit_length = (dna->length * 2 + 63) / 64;

    pq_begintypsend(&buf);
    pq_sendint(&buf, dna->length, sizeof(int));

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
    int bit_length = (dna1->length * 2 + 63) / 64;  // Number of 64-bit chunks needed

    if (dna1->length != dna2->length) {
        return false;  // Different lengths mean they can't be equal
    }

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


/*****************************************************************************
* K-mer type and functions
*****************************************************************************/

/*
This is a set returning function that generates all possible k-mers from a given DNA sequence

We don't return all kmers at once, we return them one by one, this is why we use SRF_RETURN_NEXT
*/
PG_FUNCTION_INFO_V1(generate_kmers);
Datum
generate_kmers(PG_FUNCTION_ARGS)
{
    FuncCallContext *funcctx;
    MemoryContext oldcontext;

    // Declare state
    struct {
        Dna *dna;
        int k;
    } *state;
    Dna *dna;
    int k;
    int current_index;
    char *kmer;

    // First call initialization
    if (SRF_IS_FIRSTCALL())
    {
        funcctx = SRF_FIRSTCALL_INIT();
        oldcontext = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

        // Extract arguments
        dna = PG_GETARG_DNA_P(0); // We know the first argument is a DNA sequence (and not just text)
        k = PG_GETARG_INT32(1);

        // Validate k
        if (k <= 0 || k > 32)  // k should not exceed 32, that's usually the limit of the usefulness of k in dna sequences in practice
            ereport(ERROR, (errmsg("Invalid k value: must be between 1 and 32")));

        // Save DNA and k in funcctx
        state = palloc(sizeof(*state));
        state->dna = dna;
        state->k = k;

        funcctx->user_fctx = state;
        funcctx->max_calls = dna->length - k + 1;  // Total number of kmers to generate

        MemoryContextSwitchTo(oldcontext);
    }

    // Per-call processing
    funcctx = SRF_PERCALL_SETUP();

    // Get state
    state = funcctx->user_fctx;
    // Current kmer index, we maintain "state" this way - the number of times we/they have called the function
    current_index = funcctx->call_cntr;

    if (current_index < funcctx->max_calls)
    {
        dna = state->dna;
        k = state->k;

        // Allocate memory for the kmer we will return
        kmer = palloc(k + 1);
        kmer[k] = '\0';  // Null-terminate the kmer, since it's a string

        // Decode kmer directly from the bit_sequence
        for (int i = 0; i < k; i++)
        {
            int nucleotide_index = current_index + i;  // Nucleotide/character index in the DNA sequence
            int bit_offset = (nucleotide_index * 2) % 64;  // Offset within the uint64_t
            int chunk_index = nucleotide_index / 32;      // Where is it in the bit_sequence array
            uint64_t bits = (dna->bit_sequence[chunk_index] >> bit_offset) & 0x3;  // AND with 0x3 to filter out just the last 2 bits

            // Decode bits to string nucleotide for n00b users who can't read binary B)
            switch (bits)
            {
                case 0x0: kmer[i] = 'A'; break;
                case 0x1: kmer[i] = 'T'; break;
                case 0x2: kmer[i] = 'C'; break;
                case 0x3: kmer[i] = 'G'; break;
                default:
                    ereport(ERROR, (errmsg("Unexpected bit pattern in DNA sequence - this shouldn't ever happen!")));
            }
        }

        // Return just this kmer
        SRF_RETURN_NEXT(funcctx, PointerGetDatum(cstring_to_text(kmer)));

        // TODO: Check whether we need to store kmers in a hash table (or something) and only return unique
    }
    else
    {
        // Finish and cleanup
        pfree(state);
        SRF_RETURN_DONE(funcctx);
    }
}