#include "postgres.h"
#include "utils/varlena.h"
#include "varatt.h"
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"
#include <stddef.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdint.h>


#define MAX_KMER_LENGTH 32

typedef struct Kmer
{
    char vl_len_[4];                    // Required header for PostgreSQL variable-length types
    int32 length;                       // Length of the kmer in nucleotides
    uint64_t bit_sequence[(MAX_KMER_LENGTH * 2 + 63) / 64];  // Array to store packed bits
} Kmer;

#define DatumGetKmerP(X)  ((Kmer *) DatumGetPointer(X))
#define KmerPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_KMER_P(n) DatumGetKmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_KMER_P(x) return KmerPGetDatum(x)

static void encode_kmer(const char *sequence, uint64_t *bit_sequence, int length) {
    for (int i = 0; i < length; i++) {
        int offset = (i * 2) % 64;      // Offset for 2 bits per base
        int index = i / 32;             // Each uint64_t holds 32 bases (64 bits)

        switch (sequence[i]) {
            case 'A': /* 00 */ break;
            case 'T': bit_sequence[index] |= ((uint64_t)0x1 << offset); break; // 01 for T
            case 'C': bit_sequence[index] |= ((uint64_t)0x2 << offset); break; // 10 for C
            case 'G': bit_sequence[index] |= ((uint64_t)0x3 << offset); break; // 11 for G
            default:
                ereport(ERROR, (errmsg("Invalid character in kmer sequence: %c", sequence[i])));
        }
    }
}

static char* decode_kmer(const uint64_t *bit_sequence, int length) {
    char *sequence = palloc0(length + 1);  // +1 for the null terminator

    for (int i = 0; i < length; i++) {
        int offset = (i * 2) % 64; // Offset of the current nucleotide
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

static bool validate_kmer_sequence(const char *sequence) {
    if (sequence == NULL || *sequence == '\0') {
        ereport(ERROR, (errmsg("Kmer sequence cannot be empty")));
        return false;
    }
    for (const char *p = sequence; *p; p++) {
        if (*p != 'A' && *p != 'T' && *p != 'C' && *p != 'G') {
            ereport(ERROR, (errmsg("Invalid character in kmer sequence: %c", *p)));
            return false;
        }
    }
    return true;
}

static Kmer * kmer_make(const char *sequence)
{
    int length = strlen(sequence);

    if (length > MAX_KMER_LENGTH) {
        ereport(ERROR, (errmsg("Kmer sequence exceeds maximum length of %d", MAX_KMER_LENGTH)));
        return NULL;
    }

    int num_bits = length * 2;  // 2 bits per nucleotide
    int bit_length = (num_bits + 63) / 64;  // Number of 64-bit chunks needed
    Size kmer_size = offsetof(Kmer, bit_sequence) + bit_length * sizeof(uint64_t);

    Kmer *kmer = (Kmer *) palloc0(kmer_size);

    if (!validate_kmer_sequence(sequence)) {
        return NULL;
    }

    SET_VARSIZE(kmer, kmer_size);
    kmer->length = length;

    encode_kmer(sequence, kmer->bit_sequence, length);

    return kmer;
}

static char * kmer_to_str(const Kmer *kmer)
{
    if (kmer->length == 0) {
        return pstrdup("");
    }
    return decode_kmer(kmer->bit_sequence, kmer->length);
}

PG_FUNCTION_INFO_V1(kmer_in);
Datum
kmer_in(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);
    Kmer *kmer = kmer_make(str);

    PG_RETURN_KMER_P(kmer);
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

PG_FUNCTION_INFO_V1(kmer_constructor);
Datum
kmer_constructor(PG_FUNCTION_ARGS)
{
    char *sequence = PG_GETARG_CSTRING(0);
    PG_RETURN_KMER_P(kmer_make(sequence));
}

PG_FUNCTION_INFO_V1(kmer_to_string);
Datum
kmer_to_string(PG_FUNCTION_ARGS)
{
    Kmer *kmer = PG_GETARG_KMER_P(0);
    char *result = kmer_to_str(kmer);
    PG_FREE_IF_COPY(kmer, 0);
    PG_RETURN_CSTRING(result);
}
