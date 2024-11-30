#include "postgres.h"
#include "utils/varlena.h"
#include "varatt.h" // For SET_VARSIZE!
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"
#include <stddef.h>  // Include for offsetof
#include "funcapi.h"
#include "utils/builtins.h" // For cstring_to_text
#include "common/hashfn.h" // For hash_any() in kmer_hash
#include <inttypes.h>
#include "access/spgist.h"
#include "parser/parse_type.h"
#include "utils/lsyscache.h"
#include "nodes/nodes.h"
#include "nodes/makefuncs.h"
#include "nodes/parsenodes.h"
#include "common/int.h" // For pg_cmp_s16()

#include "utils/pg_locale.h"
#include "utils/datum.h"
#include "utils/fmgrprotos.h"
#include "mb/pg_wchar.h"
#include "utils/sortsupport.h"

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdint.h>

PG_MODULE_MAGIC;

/**
 * DNA structure
 *
 * We store the DNA sequence as a bit sequence, each nucleotide takes 2 bits
 * This is done to save space since DNA sequences can be very long (~3 billion nucleotides in the human genome)
 *
 * Also, we store the length of the DNA sequence in nucleotides so that we can decode it properly
 * The vl_len_ header is required for PostgreSQL variable-length types
 */
typedef struct Dna
{
    char vl_len_[4];                    // Required header for PostgreSQL variable-length types
    uint64_t length;                     // Length of the DNA sequence in nucleotides
    uint64_t bit_sequence[FLEXIBLE_ARRAY_MEMBER];  // Array to store packed bits
} Dna;

// In simple words, datum is like void * with additional size header and here we define macros.
#define DatumGetDnaP(X)  ((Dna *) DatumGetPointer(X)) // We convert the datum pointer into a dna pointer
#define DnaPGetDatum(X)  PointerGetDatum(X) // We covert the dna pointer into a Datum pointer
#define PG_GETARG_DNA_P(n) DatumGetDnaP(PG_GETARG_DATUM(n)) // We get the nth argument given to a function
#define PG_RETURN_DNA_P(x) return DnaPGetDatum(x) // ¯\_(ツ)_/¯

/**
 * K-mer structure
 *
 * We store the K-mer as a bit sequence, each nucleotide takes 2 bits
 * Maximum length of a K-mer is 32 nucleotides, which is 64 bits, so one single 64-bit integer is enough
 */
typedef struct Kmer {
    int32 length;     // Length of the K-mer in nucleotides (each is 2 bits)
    uint64_t bit_sequence;  // Encoded bit sequence K-mer of max length 64 bits
    // Here, we don't need vl_len_ since the max length of a kmer is 32 nucleotides, which is 64 bits!
} Kmer;

#define DatumGetKmerP(X)  ((Kmer *) DatumGetPointer(X)) // We convert the datum pointer into a dna pointer
#define KmerPGetDatum(X)  PointerGetDatum(X) // We covert the dna pointer into a Datum pointer
#define PG_GETARG_KMER_P(n) DatumGetKmerP(PG_GETARG_DATUM(n)) // We get the nth argument given to a function
#define PG_RETURN_KMER_P(x) return KmerPGetDatum(x) // ¯\_(ツ)_/¯

/**
 * Qkmer structure
 *
 * We store the Qkmer as a string, each character is a nucleotide
 * We later enforce that the Qkmer must contain valid IUPAC nucleotide codes, and is at most 32 characters long
 * (since the maximum length of a K-mer is 32 nucleotides)
 *
 * The vl_len_ header is required for PostgreSQL variable-length types
 */
typedef struct Qkmer {
    char vl_len_[4];    // Again, required header for PostgreSQL variable-length types
    char sequence[FLEXIBLE_ARRAY_MEMBER]; // Flexible array for storing the sequence
} Qkmer;

// Macros for Qkmer
#define DatumGetQkmerP(X) ((Qkmer *) DatumGetPointer(X))
#define QkmerPGetDatum(X) PointerGetDatum(X)
#define PG_GETARG_QKMER_P(n) DatumGetQkmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_QKMER_P(x) return QkmerPGetDatum(x)


/**
Great reference for pg functions:
https://doxygen.postgresql.org/varatt_8h.html (For SET_VARSIZE)
https://www.postgresql.org/docs/16/index.html

SPGiST refs:
https://www.postgresql.org/docs/16/indexes-types.html#INDEXES-TYPE-SPGIST
https://www.postgresql.org/docs/current/spgist.html
https://doxygen.postgresql.org/spgtextproc_8c_source.html
*/

/********************************************************************************************
* DNA functions
********************************************************************************************/

/**
 * Encoding function
 *
 * Calculate number of 64-bit chunks we need, since rest of the int will be padded with zeros,
 * we also store the length so that we can decode it later properly without decoding extra "00"s as "A"s
 */
static void encode_dna(const char *sequence, uint64_t *bit_sequence, uint64_t length) {
    for (uint64_t i = 0; i < length; i++) {
        uint64_t offset = (i * 2) % 64;      // Offset for 2 bits per base
        uint64_t index = i / 32;             // Each uint64_t holds 32 bases (64 bits)

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
 * Decoding function
 *
 * We decode the sequence by shifting the bits to the right and then reading the last 2 bits to get the nucleotide
 */
static char* decode_dna(const uint64_t *bit_sequence, uint64_t length) {
    char *sequence = palloc0(length + 1);  // +1 for the null terminator which is added automatically by palloc0

    for (uint64_t i = 0; i < length; i++) {
        uint64_t offset = (i * 2) % 64; // Offset of the current nucleotide since it's 2 bits long
        uint64_t index = i / 32;
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

/**
 * Validates the DNA sequence
 *
 * Checks if the sequence is not empty and contains only A, T, C, G
 */
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

/**
 * Creates and returns a new Dna struct by encoding the provided DNA sequence string "ATCG" into binary format (2 bits per nucleotide)
 *
 * It checks the input, calculates required memory, and calls encode_dna which stores the enocded sequence in bit_sequence
 */
static Dna * dna_make(const char *sequence)
{
    uint64_t length = (uint64_t) strlen(sequence);
    uint64_t num_bits = length * 2;  // 2 bits per nucleotide
    uint64_t bit_length = (num_bits + 63) / 64;  // Number of 64-bit chunks we need, rest will be padded with zeros
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

/**
 * Converts a Dna struct to a string
 *
 * We decode the bit_sequence to a string and return it
 */
static char * dna_to_str(const Dna *dna)
{
    if (dna->length == 0) {
        return pstrdup("");
    }
    return decode_dna(dna->bit_sequence, dna->length);
}

/**
 * General functions for DNA
 */
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
  Dna *dna = (Dna *) PG_GETARG_VARLENA_P(0);
  char *result = dna_to_str(dna);
  PG_FREE_IF_COPY(dna, 0);
  PG_RETURN_CSTRING(result);
}

/*
 * This function is supposed to take in an existing DNA sequence and return a new DNA sequence with the same values!
 * An existing sequence means it's a binary encoded sequence in String format; just Postgres things!
 */
PG_FUNCTION_INFO_V1(dna_recv);
Datum
dna_recv(PG_FUNCTION_ARGS)
{
    StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);

    uint64_t length = pq_getmsgint(buf, sizeof(uint64_t));
    uint64_t num_bits = length * 2;  // Each DNA nucleotide is 2 bits
    uint64_t bit_length = (num_bits + 63) / 64;  // Number of 64-bit chunks needed, rounding up

    Size dna_size = offsetof(Dna, bit_sequence) + bit_length * sizeof(uint64_t);

    Dna *dna = (Dna *) palloc0(dna_size);
    SET_VARSIZE(dna, dna_size);
    //dna->vl_len_ = VARHDRSZ + dna_size;
	
    dna->length = length;

    for (uint64_t i = 0; i < bit_length; i++) {
        dna->bit_sequence[i] = pq_getmsgint64(buf);
    }

    PG_RETURN_POINTER(dna);
}

/*
 * Does the same but in reverse, takes a DNA sequence and encodes it into a binary format
 *
 * We can't use dna_to_str here because
 * dna_to_str converts the binary data back to a string format which is not what we want here!
 */
PG_FUNCTION_INFO_V1(dna_send);
Datum
dna_send(PG_FUNCTION_ARGS)
{
    Dna *dna = (Dna *) PG_GETARG_POINTER(0);
    StringInfoData buf;

    uint64_t bit_length = (dna->length * 2 + 63) / 64;

    pq_begintypsend(&buf);
    pq_sendint(&buf, dna->length, sizeof(uint64_t));

    for (uint64_t i = 0; i < bit_length; i++) {
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
  Dna *dna  = (Dna *) PG_GETARG_VARLENA_P(0);
  text *out = (text *)DirectFunctionCall1(textin,
            PointerGetDatum(dna_to_str(dna)));
  PG_FREE_IF_COPY(dna, 0);
  PG_RETURN_TEXT_P(out);
}

PG_FUNCTION_INFO_V1(dna_constructor);
Datum
dna_constructor(PG_FUNCTION_ARGS){
  char *sequence = PG_GETARG_CSTRING(0); 
  PG_RETURN_DNA_P(dna_make(sequence));
}

PG_FUNCTION_INFO_V1(dna_to_string);
Datum
dna_to_string(PG_FUNCTION_ARGS)
{
    Dna *dna = (Dna *) PG_GETARG_VARLENA_P(0);
    char *result = decode_dna(dna->bit_sequence, dna->length);  // Decode bit_sequence to a readable string
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_CSTRING(result);
}

/**
* Magically faster than strcmp!
*/
static bool dna_eq_internal(Dna *dna1, Dna *dna2)
{
    uint64_t bit_length = (dna1->length * 2 + 63) / 64;  // Number of 64-bit chunks needed

    if (dna1->length != dna2->length) {
        return false;  // Different lengths mean they can't be equal
    }

    // Compare each 64-bit chunk in the bit_sequence array
    for (uint64_t i = 0; i < bit_length; i++) {
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
  Dna *dna1 = (Dna *) PG_GETARG_VARLENA_P(0);
  Dna *dna2 = (Dna *) PG_GETARG_VARLENA_P(1);
  bool result = dna_eq_internal(dna1, dna2);
  PG_FREE_IF_COPY(dna1, 0);
  PG_FREE_IF_COPY(dna2, 1);
  PG_RETURN_BOOL(result);
}

PG_FUNCTION_INFO_V1(length);
Datum
length(PG_FUNCTION_ARGS)
{
    Dna *dna = (Dna *) PG_GETARG_VARLENA_P(0);
    uint64_t length = dna->length;  // Directly get the length field
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_INT32(length);
}

PG_FUNCTION_INFO_V1(dna_ne);
Datum
dna_ne(PG_FUNCTION_ARGS)
{
    Dna *dna1 = (Dna *) PG_GETARG_VARLENA_P(0);
    Dna *dna2 = (Dna *) PG_GETARG_VARLENA_P(1);
    bool result = !dna_eq_internal(dna1, dna2);  // ~ the result of dna_eq_internal, viola!
    PG_FREE_IF_COPY(dna1, 0);
    PG_FREE_IF_COPY(dna2, 1);
    PG_RETURN_BOOL(result);
}


/********************************************************************************************
* Kmer functions
********************************************************************************************/

/**
 * Encoding function for K-mers
 *
 * Store the k-mer in a single 64-bit variable, each nucleotide takes 2 bits
 * Encoding 'A' as 00, 'T' as 01, 'C' as 10, and 'G' as 11
*/
static uint64_t encode_kmer(const char *sequence, int length) {
    uint64_t bit_sequence = 0;  // Single 64-bit variable to store the k-mer!

    if (length <= -1 || length > 32) { // Literally cannot store a k-mer longer than 32 nucleotides
        ereport(ERROR, (errmsg("K-mer length must be between 1 and 32 nucleotides")));
    }

    for (int i = 0; i < length; i++) {
        int offset = i * 2;  // Each nucleotide takes 2 bits

        switch (sequence[i]) {
            case 'A': /* 00 for A */ break;  // No operation needed for 'A'
            case 'T': bit_sequence |= ((uint64_t)0x1 << offset); break; // 01 for T
            case 'C': bit_sequence |= ((uint64_t)0x2 << offset); break; // 10 for C
            case 'G': bit_sequence |= ((uint64_t)0x3 << offset); break; // 11 for G
            default:
                ereport(ERROR, (errmsg("Invalid character in K-mer: %c", sequence[i])));
        }
    }

    return bit_sequence;  // Return the encoded k-mer
}


/*
 * Decoding function for K-mers
 *
 * We decode the sequence by shifting the bits to the right and then reading the last 2 bits to get the nucleotide
 */
static char* decode_kmer(uint64_t bit_sequence, int length) {
    // Allocate memory for the k-mer string (+1 for null terminator)
    char *sequence = palloc(length + 1);
    sequence[length] = '\0';

    if (length <= 0 || length > 32) { // Ideally, nobody should pass invalid lengths, but just in case
        ereport(ERROR, (errmsg("K-mer length must be between 1 and 32 nucleotides")));
    }

    for (int i = 0; i < length; i++) {
        int offset = i * 2;  // Each nucleotide is 2 bits
        uint64_t bits = (bit_sequence >> offset) & 0x3;  // Get the 2 bits we care about in this iteration

        switch (bits) {
            case 0x0: sequence[i] = 'A'; break;  // 00 -> A
            case 0x1: sequence[i] = 'T'; break;  // 01 -> T
            case 0x2: sequence[i] = 'C'; break;  // 10 -> C
            case 0x3: sequence[i] = 'G'; break;  // 11 -> G
            default:
                ereport(ERROR, (errmsg("Unexpected bit pattern in K-mer: %lu", bits)));
        }
    }

    return sequence;
}

/*
 * Only difference here (from DNA) is that we also check the length of the k-mer
 */
static bool validate_kmer_sequence(const char *sequence) {
    int length = strlen(sequence);

    if (sequence == NULL /*|| *sequence == '\0'*/) {
        ereport(ERROR, (errmsg("K-mer sequence cannot be empty")));
        return false;
    }

    if (length > 32) {
        ereport(ERROR, (errmsg("K-mer length cannot exceed 32 nucleotides")));
        return false;
    }

    for (const char *p = sequence; *p; p++) {
        if (*p != 'A' && *p != 'T' && *p != 'C' && *p != 'G') {
            ereport(ERROR, (errmsg("Invalid character in K-mer sequence: %c", *p)));
            return false;
        }
    }

    return true;
}


/*
 * Creates and returns a new Kmer struct by encoding the provided DNA/K-mer sequence string "ATCG" into binary format (2 bits per nucleotide)
 *
 * It checks the input, calculates required memory, and calls encode_kmer which stores the encoded sequence in bit_sequence
 */
static Kmer *kmer_make(const char *sequence)
{
   int length = strlen(sequence);

   // Allocate memory for Kmer struct
   Kmer *kmer = (Kmer *) palloc0(sizeof(Kmer));
   kmer->length = length;

   // Validate input
   if (sequence == NULL) {
       ereport(ERROR, (errmsg("K-mer sequence cannot be NULL")));
       return NULL;
   }

   if (!validate_kmer_sequence(sequence)) {
       ereport(ERROR, (errmsg("Invalid K-mer sequence: must contain only A, T, C, G and be at most 32 nucleotides long")));
       return NULL;
   }

   // Encode the K-mer sequence into the 64-bit bit_sequence
   kmer->bit_sequence = encode_kmer(sequence, length);

   return kmer;
}

/*
 * String representation of a K-mer
 */
static char *kmer_to_str(const Kmer *kmer)
{
    if (kmer->length == 0) {
        return pstrdup("");  // Return an empty string if the k-mer has no nucleotides
    }
    return decode_kmer(kmer->bit_sequence, kmer->length);  // Decode the k-mer and return the result
}

PG_FUNCTION_INFO_V1(kmer_in);
Datum
kmer_in(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);  // Get the input string
    Kmer *kmer = kmer_make(str);       // Use kmer_make to create the encoded Kmer object

    PG_RETURN_POINTER(kmer);          // Return the Kmer as a Datum
}

PG_FUNCTION_INFO_V1(kmer_out);
Datum
kmer_out(PG_FUNCTION_ARGS)
{
    Kmer *kmer = (Kmer *) PG_GETARG_POINTER(0);  // Get the Kmer object
    char *result = kmer_to_str(kmer);            // Convert the Kmer to a string
    PG_FREE_IF_COPY(kmer, 0);                    // Free memory if needed
    PG_RETURN_CSTRING(result);                   // Return the string
}

/*
 * This function is supposed to take in an existing K-mer sequence and return a new K-mer sequence with the same values!
 * An existing sequence means it's a binary encoded sequence in String format; just Postgres things!
 */
PG_FUNCTION_INFO_V1(kmer_recv);
Datum
kmer_recv(PG_FUNCTION_ARGS)
{
    StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);

    // Read the length of the K-mer
    int length = pq_getmsgint(buf, sizeof(int));

    // Allocate memory for the Kmer object
    Kmer *kmer = (Kmer *) palloc0(sizeof(Kmer));
    kmer->length = length;

    if (length <= 0 || length > 32) {
        ereport(ERROR, (errmsg("Invalid K-mer length: must be between 1 and 32")));
    }

    // Read the encoded bit_sequence from the message buffer
    kmer->bit_sequence = pq_getmsgint64(buf);

    PG_RETURN_POINTER(kmer);
}

/*
 * Does the same but in reverse, takes a DNA/K-mer sequence and encodes it into a binary format
 *
 * We can't use kmer_to_str here because
 * kmer_to_str converts the binary data back to a string format which is not what we want here!
 */
PG_FUNCTION_INFO_V1(kmer_send);
Datum
kmer_send(PG_FUNCTION_ARGS)
{
    Kmer *kmer = (Kmer *) PG_GETARG_POINTER(0);
    StringInfoData buf;

    pq_begintypsend(&buf);

    // Serialize the length of the K-mer
    pq_sendint(&buf, kmer->length, sizeof(int));

    // Serialize the 64-bit bit_sequence
    pq_sendint64(&buf, kmer->bit_sequence);

    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

/*
 * Convert a text object to a Kmer object
 */
PG_FUNCTION_INFO_V1(kmer_cast_from_text);
Datum
kmer_cast_from_text(PG_FUNCTION_ARGS)
{
    text *txt = PG_GETARG_TEXT_P(0);  // Get the input text
    char *str = DatumGetCString(DirectFunctionCall1(textout, PointerGetDatum(txt)));  // Convert to C string
    Kmer *kmer = kmer_make(str);  // Encode the string as a Kmer
    PG_RETURN_POINTER(kmer);  // Return the Kmer object
}

/*
 * Convert a Kmer object to a text object
 */
PG_FUNCTION_INFO_V1(kmer_cast_to_text);
Datum
kmer_cast_to_text(PG_FUNCTION_ARGS)
{
    Kmer *kmer = (Kmer *) PG_GETARG_POINTER(0);  // Get the Kmer object
    text *out = (text *) DirectFunctionCall1(textin,
                    PointerGetDatum(kmer_to_str(kmer)));  // Convert the Kmer to a string and then to text
    PG_FREE_IF_COPY(kmer, 0);  // Free memory if the Kmer was copied
    PG_RETURN_TEXT_P(out);  // Return the text object
}

/*
 * Get the input C string and create a Kmer object
 *
 * These are somehow required by Postgres to create a new Kmer object
 */
PG_FUNCTION_INFO_V1(kmer_constructor);
Datum
kmer_constructor(PG_FUNCTION_ARGS)
{
    char *sequence = PG_GETARG_CSTRING(0);
    PG_RETURN_POINTER(kmer_make(sequence));
}

PG_FUNCTION_INFO_V1(kmer_to_string);
Datum
kmer_to_string(PG_FUNCTION_ARGS)
{
    Kmer *kmer = (Kmer *) PG_GETARG_POINTER(0); // Get the Kmer object
    char *result = decode_kmer(kmer->bit_sequence, kmer->length);  // Decode the bit_sequence into a string
    result[kmer->length] = '\0';  // Null-terminate the string
    PG_FREE_IF_COPY(kmer, 0);  // Free memory if the Kmer was copied
    PG_RETURN_CSTRING(result);  // Return the decoded string
}

/**
 * Compare two Kmers for equality
 *
 * We compare the lengths and the bit sequences
 */
static bool kmer_eq_internal(Kmer *kmer1, Kmer *kmer2)
{
    // Check if lengths are different
    if (kmer1->length != kmer2->length) {
        return false;
    }

    // Now compare the bit sequences
    if (kmer1->bit_sequence != kmer2->bit_sequence) {
        return false;  // Sequences are different
    }

    return true;  // Both lengths and sequences are equal
}

PG_FUNCTION_INFO_V1(kmer_eq);
Datum
kmer_eq(PG_FUNCTION_ARGS)
{
    Kmer *kmer1 = (Kmer *) PG_GETARG_POINTER(0);  // Get the first Kmer object
    Kmer *kmer2 = (Kmer *) PG_GETARG_POINTER(1);  // Get the second Kmer object
    bool result = kmer_eq_internal(kmer1, kmer2); // Compare the two Kmers
    PG_FREE_IF_COPY(kmer1, 0);
    PG_FREE_IF_COPY(kmer2, 1);
    PG_RETURN_BOOL(result);  // Return the result
}

PG_FUNCTION_INFO_V1(kmer_length);
Datum
kmer_length(PG_FUNCTION_ARGS)
{
    Kmer *kmer = (Kmer *) PG_GETARG_POINTER(0);  // Get the Kmer object
    int length = kmer->length;  // Retrieve the length
    PG_FREE_IF_COPY(kmer, 0);
    PG_RETURN_INT32(length);  // Return the length as an integer
}

PG_FUNCTION_INFO_V1(kmer_ne);
Datum
kmer_ne(PG_FUNCTION_ARGS)
{
    Kmer *kmer1 = (Kmer *) PG_GETARG_POINTER(0);  // Get the first Kmer object
    Kmer *kmer2 = (Kmer *) PG_GETARG_POINTER(1);  // Get the second Kmer object

    bool result = !kmer_eq_internal(kmer1, kmer2);  // Negate the result of kmer_eq_internal

    PG_FREE_IF_COPY(kmer1, 0);
    PG_FREE_IF_COPY(kmer2, 1);
    PG_RETURN_BOOL(result);  // Return the result
}

PG_FUNCTION_INFO_V1(kmer_hash);
Datum
kmer_hash(PG_FUNCTION_ARGS)
{
    Kmer *input = (Kmer *) PG_GETARG_POINTER(0);  // Get the Kmer input
    const unsigned char *data;
    uint32 hash;

    // Use hash_any to hash the bit_sequence field (8 bytes)
    data = (unsigned char *) &input->bit_sequence;
    hash = hash_any(data, sizeof(input->bit_sequence));  // 64-bit input, 32-bit output

    PG_RETURN_UINT32(hash);  // Return the hash as uint32
}

/*
 * This is a set returning function that generates all possible k-mers from a given DNA sequence
 *
 * We don't return all kmers at once, we return them one by one, this is why we use SRF_RETURN_NEXT
 * Reference: https://www.postgresql.org/docs/current/xfunc-c.html#XFUNC-C-RETURN-SET
 */
PG_FUNCTION_INFO_V1(generate_kmers);
Datum
generate_kmers(PG_FUNCTION_ARGS)
{
    FuncCallContext *funcctx;
    MemoryContext oldcontext;

    // Declare state and vars
    struct {
        Dna *dna;
        int k;
    } *state; // We use this to store the DNA sequence and k value, seems counter-intuitive but there's no other way to store state in a SRF
    Dna *dna; // Just used to receive the DNA sequence (it's a pointer so we don't copy the whole thing)
    int k; // Input k value
    int current_index; // Current kmer index we're generating, part of the state
    char *kmer_sequence; // Even though we store kmers as binary, we have them as strings in the intermediate stage
    Kmer *kmer; // The Kmer object we will return

    // First call initialization
    if (SRF_IS_FIRSTCALL())
    {
        funcctx = SRF_FIRSTCALL_INIT();
        oldcontext = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

        // Extract arguments
        dna = (Dna *) PG_GETARG_VARLENA_P(0); // We know the first argument is a DNA sequence (and not just text)
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

    // Next, we get the state from the function context and run it as if we were in
    state = funcctx->user_fctx;
    // Current kmer index, we maintain "state" this way - the number of times we/they have called the function
    current_index = funcctx->call_cntr;

    if (current_index < funcctx->max_calls)
    {
        dna = state->dna;
        k = state->k;

        // Allocate memory for the kmer we will return
        kmer_sequence = palloc(k + 1);
        kmer_sequence[k] = '\0';  // Null-terminate the kmer, since it's a string

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
                case 0x0: kmer_sequence[i] = 'A'; break;
                case 0x1: kmer_sequence[i] = 'T'; break;
                case 0x2: kmer_sequence[i] = 'C'; break;
                case 0x3: kmer_sequence[i] = 'G'; break;
                default:
                    ereport(ERROR, (errmsg("Unexpected bit pattern in DNA sequence - this shouldn't ever happen!")));
            }
        }

        // Encode the decoded k-mer as a Kmer object
        kmer = kmer_make(kmer_sequence);
        pfree(kmer_sequence);  // Free the temporary k-mer string

        // Return just this kmer
        SRF_RETURN_NEXT(funcctx, PointerGetDatum(kmer));
    }
    else
    {
        // Finish and cleanup
        pfree(state);
        SRF_RETURN_DONE(funcctx);
    }
}
/*
 * Basically just checks if the prefix is the same as the first n nucleotides of the kmer
 */
PG_FUNCTION_INFO_V1(starts_with);
Datum
starts_with(PG_FUNCTION_ARGS)
{
    int prefix_bits;
    uint64_t mask;
    bool result;

    Kmer *kmer = (Kmer *) PG_GETARG_POINTER(0);  // Get the prefix Kmer
    Kmer *prefix = (Kmer *) PG_GETARG_POINTER(1);    // Get the target Kmer

    // Prefix length must not exceed Kmer length
    if (prefix->length > kmer->length) {
        ereport(ERROR, (errmsg("Prefix length cannot exceed kmer length")));
    }

    // Calculate the number of bits to compare
    prefix_bits = prefix->length * 2;

    // Check if the first prefix_bits match in both Kmers
    mask = ((uint64_t)1 << prefix_bits) - 1;  // Mask for prefix_bits
    result = (prefix->bit_sequence == (kmer->bit_sequence & mask));

    PG_RETURN_BOOL(result);
}

/********************************************************************************************
* Qkmer functions
********************************************************************************************/

/*
 * Enforces that the qkmer pattern is valid and contains only IUPAC nucleotide codes
 * Also, the pattern must not be empty and must not exceed 32 characters
 */
static bool validate_qkmer_pattern(const char *pattern) {
    if (pattern == NULL || *pattern == '\0') {
        ereport(ERROR, (errmsg("qkmer pattern cannot be empty")));
        return false;
    }

    // Check length
    if (strlen(pattern) > 32) {
        ereport(ERROR, (errmsg("Qkmer pattern length cannot exceed 32 characters")));
        return false;
    }

    for (const char *p = pattern; *p; p++) {
        switch (*p) {
            case 'A': case 'T': case 'C': case 'G': case 'U': case 'W': case 'S': case 'M': case 'K':
            case 'R': case 'Y': case 'B': case 'D': case 'H': case 'V': case 'N':
                break;
            default:
                ereport(ERROR, (errmsg("Invalid character in qkmer pattern: %c", *p)));
                return false;
        }
    }

    return true;
}

/**
 * Encoding function for Q-kmers
 *
 * Store the q-kmer in a char array, each nucleotide takes 1 byte
 * There's also potential to use shorter header with SET_VARSIZE_SHORT, but it's harder to work with somehow
 */
static Qkmer *qkmer_make(const char *sequence)
{
    int length = strlen(sequence);
    Size qkmer_size;
    Qkmer *qkmer;

    // Validate sequence characters
    if (!validate_qkmer_pattern(sequence)) {
        ereport(ERROR, (errmsg("Invalid Qkmer sequence: must contain valid IUPAC nucleotide codes")));
        return NULL;
    }

    qkmer_size = offsetof(Qkmer, sequence) + length + 1;  // +1 for null terminator
    qkmer = (Qkmer *) palloc0(qkmer_size);

    SET_VARSIZE(qkmer, qkmer_size);

    // Copy the sequence
    strncpy(qkmer->sequence, sequence, length);
    qkmer->sequence[length] = '\0'; // Null-terminate the char sequence

    return qkmer;
}

PG_FUNCTION_INFO_V1(qkmer_in);
Datum
qkmer_in(PG_FUNCTION_ARGS)
{
    char *str = PG_GETARG_CSTRING(0);
    Qkmer *qkmer = qkmer_make(str);
    PG_RETURN_QKMER_P(qkmer);
}

PG_FUNCTION_INFO_V1(qkmer_out);
Datum
qkmer_out(PG_FUNCTION_ARGS)
{
    Qkmer *qkmer = PG_GETARG_QKMER_P(0);
    char *result = pstrdup(qkmer->sequence); // Duplicate the string for output
    PG_RETURN_CSTRING(result);
}

PG_FUNCTION_INFO_V1(qkmer_length);
Datum
qkmer_length(PG_FUNCTION_ARGS)
{
    Qkmer *qkmer = PG_GETARG_QKMER_P(0);
    int length = strlen(qkmer->sequence); // Compute length
    PG_RETURN_INT32(length);
}

PG_FUNCTION_INFO_V1(qkmer_eq);
Datum
qkmer_eq(PG_FUNCTION_ARGS)
{
    Qkmer *qkmer1 = PG_GETARG_QKMER_P(0);
    Qkmer *qkmer2 = PG_GETARG_QKMER_P(1);
    bool result = strcmp(qkmer1->sequence, qkmer2->sequence) == 0;
    PG_RETURN_BOOL(result);
}

PG_FUNCTION_INFO_V1(qkmer_recv);
Datum
qkmer_recv(PG_FUNCTION_ARGS)
{
    Qkmer *qkmer;
    StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);

    // Read the length of the sequence
    int length = pq_getmsgint(buf, sizeof(int));

    // Allocate a temporary buffer for the sequence
    char *sequence = (char *) palloc(length + 1);
    pq_copymsgbytes(buf, sequence, length);
    sequence[length] = '\0'; // Null-terminate

    // Use qkmer_make to create the Qkmer
    qkmer = qkmer_make(sequence);
    pfree(sequence); // Free the temp buffer!

    PG_RETURN_QKMER_P(qkmer);
}

PG_FUNCTION_INFO_V1(qkmer_send);
Datum
qkmer_send(PG_FUNCTION_ARGS)
{
    Qkmer *qkmer = PG_GETARG_QKMER_P(0);
    StringInfoData buf;

    // Get the length of the Qkmer sequence
    int length = strlen(qkmer->sequence);

    // Start constructing the binary representation
    pq_begintypsend(&buf);

    // Write the length of the sequence as a 4-byte integer
    pq_sendint(&buf, length, sizeof(int));

    // Write the sequence as raw bytes (excluding the null terminator)
    pq_sendbytes(&buf, qkmer->sequence, length);

    // Return the serialized bytea
    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

/*
* Uses the output of generate_kmers to generate qkmers that match a given pattern
*
* List of IUPAC nucleotide codes we are using:
* https://en.wikipedia.org/wiki/Nucleic_acid_notation
* Symbol   Bases Represented
* A        A
* C        C
* G        G
* T        T
* U        U
* W        A, T
* S        C, G
* M        A, C
* K        G, T
* R        A, G
* Y        C, T
* B        C, G, T
* D        A, G, T
* H        A, C, T
* V        A, C, G
* N        A, C, G, T
*/
static bool nucleotide_matches(char nucleotide, char iupac) {
    switch (iupac) {
        case 'A': return nucleotide == 'A';
        case 'T': return nucleotide == 'T';
        case 'C': return nucleotide == 'C';
        case 'G': return nucleotide == 'G';
        case 'U': return nucleotide == 'U';  // Uracil (RNA equivalent of T)
        case 'W': return nucleotide == 'A' || nucleotide == 'T'; // Weak: A or T
        case 'S': return nucleotide == 'C' || nucleotide == 'G'; // Strong: C or G
        case 'M': return nucleotide == 'A' || nucleotide == 'C'; // Amino: A or C
        case 'K': return nucleotide == 'G' || nucleotide == 'T'; // Keto: G or T
        case 'R': return nucleotide == 'A' || nucleotide == 'G'; // Purine: A or G
        case 'Y': return nucleotide == 'C' || nucleotide == 'T'; // Pyrimidine: C or T
        case 'B': return nucleotide == 'C' || nucleotide == 'G' || nucleotide == 'T'; // Not A: C, G, or T
        case 'D': return nucleotide == 'A' || nucleotide == 'G' || nucleotide == 'T'; // Not C: A, G, or T
        case 'H': return nucleotide == 'A' || nucleotide == 'C' || nucleotide == 'T'; // Not G: A, C, or T
        case 'V': return nucleotide == 'A' || nucleotide == 'C' || nucleotide == 'G'; // Not T: A, C, or G
        case 'N': return true; // Any nucleotide: A, C, G, or T, assume that these are the only valid nucleotides (which we have enforced)
        default:
            ereport(ERROR, (errmsg("Invalid character in pattern: %c!", iupac)));
            return false;
    }
}

/*
 * Just checks if a kmer contains a given qkmer pattern using nucleotide_matches()
 */
PG_FUNCTION_INFO_V1(contains);
Datum
contains(PG_FUNCTION_ARGS)
{
    // Get args
    Qkmer *qkmer = (Qkmer *) PG_GETARG_POINTER(0); // Qkmer query pattern
    Kmer *kmer = (Kmer *) PG_GETARG_POINTER(1);    // kmer to match, iteratively generated from generate_kmers()
    int qkmer_length, kmer_length;

    // Validations for both types is already done in their respective types! We wouldn't even be here if they weren't valid!

    qkmer_length = strlen(qkmer->sequence);
    kmer_length = kmer->length;

    // Fail if lengths do not match
    if (qkmer_length != kmer_length) {
        ereport(ERROR, (errmsg("Qkmer pattern and kmer lengths do not match")));
    }

    // Compare the pattern with the encoded kmer bit sequence
    for (int i = 0; i < qkmer_length; i++) {
        // Decode the nucleotide at position i from the Kmer bit_sequence
        int bit_offset = i * 2; // Offset in the 64-bit bit_sequence
        uint64_t nucleotide_bits = (kmer->bit_sequence >> bit_offset) & 0x3; // Extract 2 bits for the nucleotide
        char nucleotide;

        // Map the decoded bits to a nucleotide character
        switch (nucleotide_bits) {
            case 0x0: nucleotide = 'A'; break;
            case 0x1: nucleotide = 'T'; break;
            case 0x2: nucleotide = 'C'; break;
            case 0x3: nucleotide = 'G'; break;
            default:
                ereport(ERROR, (errmsg("Unexpected bit pattern in Kmer")));
        }

        // Compare the nucleotide with the pattern using nucleotide_matches()
        if (!nucleotide_matches(nucleotide, qkmer->sequence[i])) {
            PG_RETURN_BOOL(false); // Return false if there's a mismatch
        }
    }

    // If we get here, the kmer matches the pattern!
    PG_RETURN_BOOL(true);
}

/********************************************************************************************
* SpGIST Indexing Functions
*
* These functions are used to create a SpGIST index on the Kmer type
* Modelled from: https://doxygen.postgresql.org/spgtextproc_8c_source.html
********************************************************************************************/

// Helper macro for forming Kmer datums
#define FORM_KMER_DATUM(kmer_ptr) KmerPGetDatum(kmer_ptr)

/*
Struct for sorting values in picksplit
 */
typedef struct spgNodePtr
{
    Datum       d;
    int         i;
    int16       c;
} spgNodePtr;

/*
* Find the length of the common prefix of a and b
*/
static int
commonPrefix(const char *a, const char *b, int lena, int lenb)
{
    int         i = 0;

    while (i < lena && i < lenb && *a == *b)
    {
        a++;
        b++;
        i++;
    }

    return i;
}

/*
* Binary search an array of int16 datums for a match to c
*
* On success, *i gets the match location; on failure, it gets where to insert
*/
static bool
searchChar(Datum *nodeLabels, int nNodes, int16 c, int *i)
{
    int StopLow = 0,
    StopHigh = nNodes;

    while (StopLow < StopHigh)
    {
        int StopMiddle = (StopLow + StopHigh) >> 1;
        int16 middle = DatumGetInt16(nodeLabels[StopMiddle]);

        if (c < middle)
            StopHigh = StopMiddle;
        else if (c > middle)
            StopLow = StopMiddle + 1;
        else
        {
            *i = StopMiddle;
            return true;
        }
    }

    *i = StopHigh;
    return false;
}

static int pg_cmp_s16(int16 a, int16 b)
{
    return (int32) a - (int32) b;
}

/*
qsort comparator to sort spgNodePtr structs by "c"
 */
static int
cmpNodePtr(const void *a, const void *b)
{
    const spgNodePtr *aa = (const spgNodePtr *) a;
    const spgNodePtr *bb = (const spgNodePtr *) b;

    return pg_cmp_s16(aa->c, bb->c);
}

/*
 A function to test and get the Oid of the Postgres type kmer
 Required for SP-GiST implementation (and is called from Sql)
*/
PG_FUNCTION_INFO_V1(get_oid);
Datum
get_oid(PG_FUNCTION_ARGS){
    Oid int4_oid = typenameTypeId(NULL, makeTypeName("kmer"));
    PG_RETURN_INT32(int4_oid);
}

PG_FUNCTION_INFO_V1(spgist_kmer_config);
Datum
spgist_kmer_config(PG_FUNCTION_ARGS)
{
    spgConfigIn  *cfgin  = (spgConfigIn *) PG_GETARG_POINTER(0);
    spgConfigOut *cfgout = (spgConfigOut *) PG_GETARG_POINTER(1);

    Oid KMEROID = typenameTypeId(NULL, makeTypeName("kmer"));
    elog(INFO, "spgist_kmer_config: Retrieved KMEROID = %u", KMEROID);

    cfgout->prefixType = KMEROID;
    cfgout->leafType = KMEROID;
    cfgout->labelType = INT2OID; // TODO CHECK THIS! Essentially, the type of the node labels since the sequence is binary, this should work!

    cfgout->canReturnData = true; // A node may contain already the value that we want without it being a leaf
    cfgout->longValuesOK = false;   // Suffixing will shorten long values

    PG_RETURN_VOID();
}

PG_FUNCTION_INFO_V1(spgist_kmer_choose);
Datum
spgist_kmer_choose(PG_FUNCTION_ARGS)
{
    spgChooseIn *in = (spgChooseIn *) PG_GETARG_POINTER(0);
    spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);

    Kmer *input_kmer = DatumGetKmerP(in->datum);
    const char *input_sequence = kmer_to_str(input_kmer);
    int input_length = input_kmer->length;

    Kmer *prefix_kmer = NULL;
    const char *prefix_sequence = NULL;
    int prefix_length = 0;
    int common_length = 0;
    int next_char = -1;

    elog(INFO, "spgist_kmer_choose: Input sequence = %s", input_sequence);
    elog(INFO, "spgist_kmer_choose: Input length = %d", input_length);
    elog(INFO, "spgist_kmer_choose: Level = %d", in->level);
    elog(INFO, "spgist_kmer_choose: Has prefix = %d", in->hasPrefix);

    if (in->hasPrefix)
    {
        prefix_kmer = DatumGetKmerP(in->prefixDatum);
        prefix_sequence = kmer_to_str(prefix_kmer);
        prefix_length = prefix_kmer->length;

        // Find common prefix length
        while (common_length < prefix_length &&
               common_length < (input_length - in->level) &&
               prefix_sequence[common_length] == input_sequence[in->level + common_length])
        {
            common_length++;
        }

        if (common_length == prefix_length)
        {
            next_char = (in->level + common_length < input_length)
                          ? input_sequence[in->level + common_length]
                          : -1;
            // Validate next_char if necessary
            if (next_char != 'A' && next_char != 'T' && next_char != 'C' && next_char != 'G' && next_char != -1)
            {
                ereport(ERROR, (errmsg("Invalid next character in K-mer sequence: %c", next_char)));
            }

        }
        else
        {
            // Must split prefix
            out->resultType = spgSplitTuple;
            if (common_length == 0)
            {
                out->result.splitTuple.prefixHasPrefix = false;
            }
            else
            {
                // Shorten the prefix sequence and create a new Kmer object
                char *truncated_prefix = palloc0(common_length + 1); // +1 for null terminator
                memcpy(truncated_prefix, prefix_sequence, common_length);
                truncated_prefix[common_length] = '\0';

                Kmer *new_prefix = kmer_make(truncated_prefix);
                elog(INFO, "spgist_kmer_choose: New prefix = %s", truncated_prefix);

                out->result.splitTuple.prefixHasPrefix = true;
                out->result.splitTuple.prefixPrefixDatum = FORM_KMER_DATUM(new_prefix);

                // Free the allocated truncated string
                pfree(truncated_prefix);
            }

            out->result.splitTuple.prefixNNodes = 1;
            out->result.splitTuple.prefixNodeLabels = (Datum *) palloc(sizeof(Datum));
            out->result.splitTuple.prefixNodeLabels[0] = Int16GetDatum(*(unsigned char *)(prefix_sequence + common_length));

            // Handle postfix splitting
            out->result.splitTuple.childNodeN = 0; // Child index for the current node

            if (prefix_length - common_length == 1)
            {
                out->result.splitTuple.postfixHasPrefix = false;
            }
            else
            {
                // Shorten the postfix sequence and create a new Kmer object
                char *postfix_prefix = palloc0(prefix_length - common_length);
                memcpy(postfix_prefix, prefix_sequence + common_length + 1, prefix_length - common_length - 1);
                postfix_prefix[prefix_length - common_length - 1] = '\0';

                Kmer *new_postfix = kmer_make(postfix_prefix);
                elog(INFO, "spgist_kmer_choose: New postfix = %s", postfix_prefix);

                out->result.splitTuple.postfixHasPrefix = true;
                out->result.splitTuple.postfixPrefixDatum = FORM_KMER_DATUM(new_postfix);

                // Free the allocated postfix string
                pfree(postfix_prefix);
            }

            PG_RETURN_VOID();
        }
    }
    else if (in->level < input_length)
    {
        next_char = input_sequence[in->level];
        elog(INFO, "spgist_kmer_choose: Next character = %c", next_char);
    }
    else
    {
        next_char = -1; // No more characters left, this is the initialisation but why not be explicit
    }

    // Look up next_char in the node label array
    int idx;
    if (searchChar(in->nodeLabels, in->nNodes, next_char, &idx))
    {
        // Match found: descend to the existing node
        int level_add = common_length;
        if (next_char >= 0)
            level_add++;

        out->resultType = spgMatchNode;
        out->result.matchNode.nodeN = idx;
        out->result.matchNode.levelAdd = level_add;
        elog(INFO, "spgist_kmer_choose: Matched node %d with level add %d", idx, level_add);

        // Construct restDatum for the unmatched suffix
        if (input_length - in->level - level_add > 0)
        {
            char *rest_sequence = palloc0(input_length - in->level - level_add + 1);
            memcpy(rest_sequence, input_sequence + in->level + level_add,
                   input_length - in->level - level_add);
            rest_sequence[input_length - in->level - level_add] = '\0';

            Kmer *rest_kmer = kmer_make(rest_sequence);
            out->result.matchNode.restDatum = FORM_KMER_DATUM(rest_kmer);
            elog(INFO, "spgist_kmer_choose: Rest sequence = %s", rest_sequence);

            pfree(rest_sequence);
        }
        else
        {
            out->result.matchNode.restDatum = FORM_KMER_DATUM(kmer_make("")); // TODO CHECK THIS! Empty string is usually not a good idea
            elog(INFO, "spgist_kmer_choose: Empty rest sequence");
        }
    }
    else if (in->allTheSame)
    {
        // Handle the all-the-same case by splitting the tuple
        Kmer *current_kmer = DatumGetKmerP(in->datum);
        const char *current_sequence = kmer_to_str(current_kmer);
        int current_length = current_kmer->length;
        elog(INFO, "spgist_kmer_choose: All the same case, current sequence = %s", current_sequence);

        Kmer *prefix_kmer = DatumGetKmerP(in->prefixDatum);
        const char *prefix_sequence = kmer_to_str(prefix_kmer);
        int prefix_length = prefix_kmer->length;
        int common_length = 0;
        while (common_length < prefix_length && common_length < current_length &&
               prefix_sequence[common_length] == current_sequence[common_length])
        {
            common_length++;
        }
        elog(INFO, "spgist_kmer_choose: Common length = %d", common_length);

        out->resultType = spgSplitTuple;

        // Shorten the prefix if necessary
        if (common_length > 0)
        {
            char *new_prefix = palloc0(common_length + 1); // +1 for null terminator
            memcpy(new_prefix, prefix_sequence, common_length);
            Kmer *shortened_prefix = kmer_make(new_prefix);
            elog(INFO, "spgist_kmer_choose: Shortened prefix = %s", new_prefix);

            out->result.splitTuple.prefixHasPrefix = true;
            out->result.splitTuple.prefixPrefixDatum = FORM_KMER_DATUM(shortened_prefix);

            pfree(new_prefix);
        }
        else
        {
            out->result.splitTuple.prefixHasPrefix = false;
        }

        // Assign the existing nodes to the split
        out->result.splitTuple.prefixNNodes = in->nNodes;
        out->result.splitTuple.prefixNodeLabels = palloc(sizeof(Datum) * in->nNodes);

        for (int i = 0; i < in->nNodes; i++)
        {
            out->result.splitTuple.prefixNodeLabels[i] = in->nodeLabels[i];
            elog(INFO, "spgist_kmer_choose: Adding node %d with label %d", i, DatumGetInt16(in->nodeLabels[i]));
        }

        // Handle the new value as a separate child node
        if (current_length > common_length)
        {
            char *new_postfix = palloc0(current_length - common_length + 1);
            memcpy(new_postfix, current_sequence + common_length, current_length - common_length);
            Kmer *postfix_kmer = kmer_make(new_postfix);

            out->result.splitTuple.postfixHasPrefix = true;
            out->result.splitTuple.postfixPrefixDatum = FORM_KMER_DATUM(postfix_kmer);

            pfree(new_postfix);
        }
        else
        {
            out->result.splitTuple.postfixHasPrefix = false;
        }
    }
    else
    {
        // Add a node for the not-previously-seen next_char value
        out->resultType = spgAddNode;
        out->result.addNode.nodeLabel = Int16GetDatum(next_char);
        out->result.addNode.nodeN = idx;
        //elog(INFO, "spgist_kmer_choose: Adding node %d with label %d", idx, next_char);
    }

    PG_RETURN_VOID();
}

PG_FUNCTION_INFO_V1(spgist_kmer_picksplit);
Datum
spgist_kmer_picksplit(PG_FUNCTION_ARGS)
{
    spgPickSplitIn *in = (spgPickSplitIn *) PG_GETARG_POINTER(0);
    spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);

    Kmer *first_kmer = DatumGetKmerP(in->datums[0]);
    char *first_sequence = kmer_to_str(first_kmer);
    int common_length = first_kmer->length; // Start with the full length of the first Kmer

    elog(INFO, "spgist_kmer_picksplit: First sequence = %s", first_sequence);
    elog(INFO, "spgist_kmer_picksplit: First length = %d", common_length);
    elog(INFO, "spgist_kmer_picksplit: Number of tuples = %d", in->nTuples);

    // Find the longest common prefix
    for (int i = 1; i < in->nTuples; i++)
    {
        Kmer *kmer = DatumGetKmerP(in->datums[i]);
        char *sequence = kmer_to_str(kmer);
        int len = 0;

        while (len < common_length && len < kmer->length && first_sequence[len] == sequence[len])
        {
            len++;
        }

        if (len < common_length)
        {
            common_length = len;
        }
    }

    // Enforce maximum prefix length FIXME: Change 10 to some other value, 32?
    common_length = Min(common_length, 10);
    elog(INFO, "spgist_kmer_picksplit: Common length = %d", common_length);

    // Construct the prefix if it exists
    if (common_length > 0)
    {
        char *truncated_prefix = palloc0(common_length + 1);
        memcpy(truncated_prefix, first_sequence, common_length);
        truncated_prefix[common_length] = '\0';

        Kmer *new_prefix = kmer_make(truncated_prefix);
        out->hasPrefix = true;
        out->prefixDatum = FORM_KMER_DATUM(new_prefix);
        elog(INFO, "spgist_kmer_picksplit: Out Prefix = %s", truncated_prefix);

        pfree(truncated_prefix);
    }
    else
    {
        out->hasPrefix = false;
    }

    spgNodePtr *nodes = palloc(sizeof(spgNodePtr) * in->nTuples);

    for (int i = 0; i < in->nTuples; i++)
    {
        Kmer *kmer = DatumGetKmerP(in->datums[i]);
        const char *sequence = kmer_to_str(kmer);
        elog(INFO, "spgist_kmer_picksplit: Considering sequence = %s", sequence);

        if (common_length < kmer->length)
        {
            nodes[i].c = *(unsigned char *) (sequence + common_length);
        }
        else
        {
            nodes[i].c = -1; // Entire string matches the prefix
        }

        nodes[i].i = i;
        nodes[i].d = in->datums[i];
    }

    // Sort nodes by their next character so we can group them, also allows binary search
    qsort(nodes, in->nTuples, sizeof(*nodes), cmpNodePtr);

    // Initialize output structures
    out->nNodes = 0;
    out->nodeLabels = palloc(sizeof(Datum) * in->nTuples);
    out->mapTuplesToNodes = palloc(sizeof(int) * in->nTuples);
    out->leafTupleDatums = palloc(sizeof(Datum) * in->nTuples);

    //elog(INFO, "spgist_kmer_picksplit: Now grouping tuples into nodes");

    // Group tuples into nodes
    for (int i = 0; i < in->nTuples; i++)
    {
        if (i == 0 || nodes[i].c != nodes[i - 1].c)
        {
            out->nodeLabels[out->nNodes] = Int16GetDatum(nodes[i].c);
            out->nNodes++;
        }

        Kmer *suffix = NULL;
        const char *sequence = kmer_to_str(DatumGetKmerP(nodes[i].d));
        elog(INFO, "spgist_kmer_picksplit: Considering sequence = %s", sequence);
        if (common_length < strlen(sequence))
        {
            char *suffix_sequence = palloc0(strlen(sequence) - common_length + 1);
            memcpy(suffix_sequence, sequence + common_length, strlen(sequence) - common_length);
            suffix_sequence[strlen(sequence) - common_length] = '\0';
            elog(INFO, "spgist_kmer_picksplit: Suffix = %s", suffix_sequence);

            suffix = kmer_make(suffix_sequence);
            pfree(suffix_sequence);
        }
        else
        {
            elog(INFO, "spgist_kmer_picksplit: Empty suffix");
            suffix = kmer_make("");
        }

        out->leafTupleDatums[nodes[i].i] = FORM_KMER_DATUM(suffix);
        out->mapTuplesToNodes[nodes[i].i] = out->nNodes - 1;
    }
    elog(INFO, "spgist_kmer_picksplit: out->nNodes = %d", out->nNodes);

    pfree(nodes);

    PG_RETURN_VOID();
}
//
//PG_FUNCTION_INFO_V1(spgist_kmer_inner_consistent);
//Datum
//spgist_kmer_inner_consistent(PG_FUNCTION_ARGS)
//{
//    spgInnerConsistentIn *in = (spgInnerConsistentIn *) PG_GETARG_POINTER(0);
//    spgInnerConsistentOut *out = (spgInnerConsistentOut *) PG_GETARG_POINTER(1);
//
//    Kmer *reconstructedKmer = NULL;
//    char *reconstructedSequence = NULL;
//    int reconstructedLength = in->level;
//    Kmer *prefixKmer = NULL;
//    char *prefixSequence = NULL;
//    int prefixLength = 0;
//    int maxReconstrLen;
//    char *reconstrSequence;
//    int i;
//    elog(INFO, "spgist_kmer_inner_consistent: Level = %d", in->level);
//    elog(INFO, "spgist_kmer_inner_consistent: Number of nodes = %d", in->nNodes);
//    elog(INFO, "spgist_kmer_inner_consistent: Number of keys = %d", in->nkeys);
//
//    /* Reconstruct the sequence from the parent value and prefix */
//    if (in->reconstructedValue != (Datum) 0)
//    {
//        reconstructedKmer = (Kmer *) DatumGetKmerP(in->reconstructedValue);
//        //elog(INFO, "spgist_kmer_inner_consistent: value of in->reconstructedValue = %d", reconstructedKmer->bit_sequence);
//        //elog(INFO, "spgist_kmer_inner_consistent: value of in->reconstructedValue = %d", reconstructedKmer->length);
//        reconstructedSequence = kmer_to_str(reconstructedKmer);
//        reconstructedLength = reconstructedKmer->length;
//        elog(INFO, "spgist_kmer_inner_consistent: Reconstructed sequence = %s", reconstructedSequence);
//    }
//
//    if (in->hasPrefix)
//    {
//        prefixKmer = DatumGetKmerP(in->prefixDatum);
//        prefixSequence = kmer_to_str(prefixKmer);
//        prefixLength = prefixKmer->length;
//        elog(INFO, "spgist_kmer_inner_consistent: Prefix sequence = %s", prefixSequence);
//    }
//
//    maxReconstrLen = reconstructedLength + prefixLength + 1; /* +1 for node label */
//    reconstrSequence = palloc0(maxReconstrLen);
//
//    if (reconstructedSequence)
//    {
//        memcpy(reconstrSequence, reconstructedSequence, reconstructedLength);
//    }
//
//    if (prefixSequence)
//    {
//        memcpy(reconstrSequence + reconstructedLength, prefixSequence, prefixLength);
//    }
//
//    elog(INFO, "spgist_kmer_inner_consistent: Reconstructed sequence 1 = %s", reconstrSequence);
//    elog(INFO, "spgist_kmer_inner_consistent: prefixSequence = %s", prefixSequence);
//
//    /* Prepare the output arrays */
//    out->nodeNumbers = (int *) palloc(sizeof(int) * in->nNodes);
//    out->levelAdds = (int *) palloc(sizeof(int) * in->nNodes);
//    out->reconstructedValues = (Datum *) palloc(sizeof(Datum) * in->nNodes);
//    out->nNodes = 0;
//
//    /* Process each child node */
//    for (i = 0; i < in->nNodes; i++)
//    {
//        int16 nodeLabel = DatumGetInt16(in->nodeLabels[i]);
//        int thisLen = reconstructedLength + prefixLength;
//        bool isConsistent = true;
//        int j;
//
//        elog(INFO, "spgist_kmer_inner_consistent: Node label = %d", nodeLabel);
//        //elog(INFO, "spgist_kmer_inner_consistent: thisLen = %d", thisLen);
//        elog(INFO, "spgist_kmer_inner_consistent: reconstructedSequence 0 = %s", reconstructedSequence);
//
//        /* Append the node label if it's not a dummy value */
//        if (nodeLabel >= 0)
//        {
//            reconstrSequence[thisLen] = (char) nodeLabel;
//            elog(INFO, "spgist_kmer_inner_consistent: Node label = %c", reconstrSequence[thisLen]);
//            thisLen++;
//            reconstrSequence[thisLen] = '\0'; // Add this line to null-terminate the string
//
//            //elog(INFO, "spgist_kmer_inner_consistent: thisLen = %d", thisLen);
//        }
//        elog(INFO, "spgist_kmer_inner_consistent: Reconstructed sequence 2 = %s", reconstrSequence);
//
//        /* Check consistency with the query conditions */
//        for (j = 0; j < in->nkeys; j++)
//        {
//            StrategyNumber strategy = in->scankeys[j].sk_strategy;
//            Kmer *queryKmer = DatumGetKmerP(in->scankeys[j].sk_argument);
//            const char *querySequence = kmer_to_str(queryKmer);
//            int queryLength = queryKmer->length;
//            int cmpResult = memcmp(reconstrSequence, querySequence, Min(thisLen, queryLength));
//
//            elog(INFO, "spgist_kmer_inner_consistent: Reconstructed sequence = %s", reconstrSequence);
//            elog(INFO, "spgist_kmer_inner_consistent: Query sequence = %s", querySequence);
//            elog(INFO, "spgist_kmer_inner_consistent: Comparison result = %d", cmpResult);
//
//            //elog(INFO, "spgist_kmer_inner_consistent: Strategy = %d", strategy);
//
//            switch (strategy)
//            {
//                case BTLessStrategyNumber:
//                case BTLessEqualStrategyNumber:
//                    if (cmpResult > 0)
//                        isConsistent = false;
//                    break;
//
//                case BTEqualStrategyNumber:
//                    if (cmpResult != 0 || queryLength != thisLen)
//                        isConsistent = false;
//                    break;
//
//                case BTGreaterEqualStrategyNumber:
//                case BTGreaterStrategyNumber:
//                    if (cmpResult < 0)
//                        isConsistent = false;
//                    break;
//
//                case RTPrefixStrategyNumber:
//                    if (cmpResult != 0)
//                        isConsistent = false;
//                    break;
//
//                default:
//                    elog(ERROR, "unrecognized strategy number: %d",
//                         strategy);
//                    break;
//            }
//
//            if (!isConsistent) // No need to check further
//                break;
//        }
//
//        // Add consistent nodes to the output
//        if (isConsistent)
//        {
//            out->nodeNumbers[out->nNodes] = i;
//            out->levelAdds[out->nNodes] = thisLen - in->level;
//
//            Kmer *newReconstructedKmer = kmer_make(reconstrSequence);
//            newReconstructedKmer->length = thisLen;
//            elog(INFO, "spgist_kmer_inner_consistent: New reconstructed sequence = %s", reconstrSequence);
//
//            out->reconstructedValues[out->nNodes] =
//                FORM_KMER_DATUM(newReconstructedKmer);
//
//            out->nNodes++;
//        }
//    }
//
//    // Clean up
//    pfree(reconstrSequence);
//
//    PG_RETURN_VOID();
//}
//
//PG_FUNCTION_INFO_V1(spgist_kmer_leaf_consistent);
//Datum
//spgist_kmer_leaf_consistent(PG_FUNCTION_ARGS)
//{
//    spgLeafConsistentIn *in = (spgLeafConsistentIn *) PG_GETARG_POINTER(0);
//    spgLeafConsistentOut *out = (spgLeafConsistentOut *) PG_GETARG_POINTER(1);
//    int level = in->level;
//    Kmer *leafKmer = DatumGetKmerP(in->leafDatum);
//    Kmer *reconstrKmer = NULL;
//    char *fullSequence = NULL;
//    int fullLength;
//    bool result = true;
//
//    /* All tests are exact */
//    out->recheck = false;
//
//    /* Retrieve reconstructed value from parent level, if any */
//    if (DatumGetPointer(in->reconstructedValue))
//    {
//        reconstrKmer = DatumGetKmerP(in->reconstructedValue);
//        elog(INFO, "spgist_kmer_leaf_consistent: Reconstructed sequence = %s", kmer_to_str(reconstrKmer));
//    }
//
//    /* Ensure that reconstructed Kmer's length matches the level */
//    Assert(reconstrKmer == NULL ? level == 0 : reconstrKmer->length == level);
//
//    /* Reconstruct the full sequence represented by this leaf tuple */
//    fullLength = level + leafKmer->length;
//    if (leafKmer->length == 0 && level > 0)
//    {
//        /* No leaf value, use only the reconstructed sequence */
//        fullSequence = kmer_to_str(reconstrKmer);
//        out->leafValue = in->reconstructedValue;
//    }
//    else
//    {
//        /* Combine reconstructed value and leaf value */
//        fullSequence = palloc(fullLength + 1); /* +1 for null-terminator */
//        if (level > 0)
//        {
//            memcpy(fullSequence, kmer_to_str(reconstrKmer), level);
//        }
//        if (leafKmer->length > 0)
//        {
//            memcpy(fullSequence + level, kmer_to_str(leafKmer), leafKmer->length);
//        }
//        fullSequence[fullLength] = '\0'; /* Null-terminate the combined sequence */
//
//        /* Create a new Kmer for the reconstructed full sequence */
//        Kmer *fullKmer = kmer_make(fullSequence);
//        out->leafValue = FORM_KMER_DATUM(fullKmer);
//    }
//
//    /* Perform the required comparisons */
//    for (int i = 0; i < in->nkeys; i++)
//    {
//        StrategyNumber strategy = in->scankeys[i].sk_strategy;
//        Kmer *queryKmer = DatumGetKmerP(in->scankeys[i].sk_argument);
//        const char *querySequence = kmer_to_str(queryKmer);
//        int queryLength = queryKmer->length;
//        int cmpResult;
//
//        if (strategy == RTPrefixStrategyNumber)
//        {
//            /*
//             * If level >= query length, then the reconstructed value must begin
//             * with the query string, so we don't need to check it again.
//             */
//            result = (level >= queryLength) ||
//                     (strncmp(fullSequence, querySequence, queryLength) == 0);
//
//            if (!result)
//                break; /* No need to consider further conditions */
//            continue;
//        }
//
//        /* Perform comparison based on strategy */
//        elog(INFO, "spgist_kmer_leaf_consistent: Full sequence = %s", fullSequence);
//        elog(INFO, "spgist_kmer_leaf_consistent: Query sequence = %s", querySequence);
//        cmpResult = strncmp(fullSequence, querySequence, Min(fullLength, queryLength));
//        elog(INFO, "spgist_kmer_leaf_consistent: Comparison result = %d", cmpResult);
//        if (cmpResult == 0)
//        {
//            if (queryLength > fullLength)
//                cmpResult = -1;
//            else if (queryLength < fullLength)
//                cmpResult = 1;
//        }
//
//        switch (strategy)
//        {
//            case BTLessStrategyNumber:
//                result = (cmpResult < 0);
//                break;
//            case BTLessEqualStrategyNumber:
//                result = (cmpResult <= 0);
//                break;
//            case BTEqualStrategyNumber:
//                elog(INFO, "spgist_kmer_leaf_consistent: REACHED BTEqualStrategyNumber for %s", querySequence);
//                result = (cmpResult == 0);
//                break;
//            case BTGreaterEqualStrategyNumber:
//                result = (cmpResult >= 0);
//                break;
//            case BTGreaterStrategyNumber:
//                result = (cmpResult > 0);
//                break;
//            default:
//                elog(ERROR, "unrecognized strategy number: %d",
//                     strategy);
//                result = false;
//                break;
//        }
//
//        if (!result)
//            break; /* No need to consider further conditions */
//    }
//
//    /* Return the final comparison result */
//    PG_RETURN_BOOL(result);
//}
PG_FUNCTION_INFO_V1(spgist_kmer_inner_consistent);
Datum
spgist_kmer_inner_consistent(PG_FUNCTION_ARGS)
{
    spgInnerConsistentIn *in = (spgInnerConsistentIn *) PG_GETARG_POINTER(0);
    spgInnerConsistentOut *out = (spgInnerConsistentOut *) PG_GETARG_POINTER(1);

    // Pointers for sequence reconstruction
    Kmer *reconstructedKmer = NULL;
    char *reconstructedSequence = NULL;
    Kmer *prefixKmer = NULL;
    char *prefixSequence = NULL;
    int reconstructedLength = 0;
    int prefixLength = 0;

    // Reconstruct the sequence from parent value and prefix
    if (in->reconstructedValue != (Datum) 0)
    {
        reconstructedKmer = DatumGetKmerP(in->reconstructedValue);
        reconstructedSequence = kmer_to_str(reconstructedKmer);
        reconstructedLength = reconstructedKmer->length;
        elog(INFO, "spgist_kmer_inner_consistent: Reconstructed sequence = %s", reconstructedSequence);
    }

    if (in->hasPrefix)
    {
        prefixKmer = DatumGetKmerP(in->prefixDatum);
        prefixSequence = kmer_to_str(prefixKmer);
        prefixLength = prefixKmer->length;
        elog(INFO, "spgist_kmer_inner_consistent: Prefix sequence = %s", prefixSequence);
        elog(INFO, "spgist_kmer_inner_consistent: Prefix length = %d", prefixLength);
        elog(INFO, "spgist_kmer_inner_consistent: Reconstructed length = %d", reconstructedLength);
    }

    // Calculate the total maximum possible reconstructed length
    int maxReconstructedLength = reconstructedLength + prefixLength + 1;
    char *fullReconstructedSequence = palloc0(maxReconstructedLength + 1);
    elog(INFO, "spgist_kmer_inner_consistent: Max reconstructed length = %d", maxReconstructedLength);
    elog(INFO, "spgist_kmer_inner_consistent: Full reconstructed sequence = %s", fullReconstructedSequence);

    // Safely copy reconstructed and prefix sequences
    if (reconstructedSequence)
        memcpy(fullReconstructedSequence, reconstructedSequence, reconstructedLength);
    if (prefixSequence)
        memcpy(fullReconstructedSequence + reconstructedLength, prefixSequence, prefixLength);

    // Prepare output arrays with dynamic allocation
    out->nodeNumbers = palloc(sizeof(int) * in->nNodes);
    out->levelAdds = palloc(sizeof(int) * in->nNodes);
    out->reconstructedValues = palloc(sizeof(Datum) * in->nNodes);
    out->nNodes = 0;

    // Process each child node for consistency
    for (int i = 0; i < in->nNodes; i++)
    {
        int16 nodeLabel = DatumGetInt16(in->nodeLabels[i]);
        int currentLength = reconstructedLength + prefixLength;
        bool isConsistent = true;

        // Append node label if valid
        if (nodeLabel >= 0)
        {
            fullReconstructedSequence[currentLength] = (char)nodeLabel;
            currentLength++;
            fullReconstructedSequence[currentLength] = '\0';
            elog(INFO, "spgist_kmer_inner_consistent: Node label = %c", fullReconstructedSequence[currentLength]);
            elog(INFO, "spgist_kmer_inner_consistent: Full reconstructed sequence = %s", fullReconstructedSequence);
        }

        // Check consistency against all scan keys
        for (int j = 0; j < in->nkeys; j++)
        {
            StrategyNumber strategy = in->scankeys[j].sk_strategy;
            Kmer *queryKmer = DatumGetKmerP(in->scankeys[j].sk_argument);
            const char *querySequence = kmer_to_str(queryKmer);
            int queryLength = queryKmer->length;
            elog(INFO, "spgist_kmer_inner_consistent: Query sequence = %s", querySequence);
            elog(INFO, "spgist_kmer_inner_consistent: Query length = %d", queryLength);
            elog(INFO, "spgist_kmer_inner_consistent: Current length = %d", currentLength);
            elog(INFO, "spgist_kmer_inner_consistent: Full reconstructed sequence = %s", fullReconstructedSequence);

            int comparisonResult = strncmp(fullReconstructedSequence, querySequence,
                                           Min(currentLength, queryLength));

            elog(INFO, "spgist_kmer_inner_consistent: Comparison result = %d", comparisonResult);
            // Adjust comparison for equal-length prefix
            if (comparisonResult == 0 && currentLength != queryLength)
                comparisonResult = (currentLength < queryLength) ? -1 : 1;
            elog(INFO, "spgist_kmer_inner_consistent: Comparison result = %d", comparisonResult);

            elog(INFO, "spgist_kmer_inner_consistent: Strategy = %d", strategy);
            switch (strategy)
            {
                case BTLessStrategyNumber:
                case BTLessEqualStrategyNumber:
                    isConsistent = (comparisonResult <= 0);
                    break;
                case BTEqualStrategyNumber:
                    elog(INFO, "spgist_kmer_inner_consistent: REACHED BTEqualStrategyNumber for %s", querySequence);
                    isConsistent = (comparisonResult == 0);
                    break;
                case BTGreaterEqualStrategyNumber:
                case BTGreaterStrategyNumber:
                    isConsistent = (comparisonResult >= 0);
                    break;
                case RTPrefixStrategyNumber:
                    isConsistent = (comparisonResult == 0 &&
                                    currentLength >= queryLength);
                    break;
                default:
                    elog(ERROR, "Unrecognized strategy number: %d", strategy);
                    isConsistent = false;
                    break;
            }

            if (!isConsistent)
                break;
        }

        // Add consistent nodes to output
        if (isConsistent)
        {
            Kmer *reconstructedNodeKmer = kmer_make(fullReconstructedSequence);
            reconstructedNodeKmer->length = currentLength;

            out->nodeNumbers[out->nNodes] = i;
            out->levelAdds[out->nNodes] = currentLength - in->level;
            out->reconstructedValues[out->nNodes] = FORM_KMER_DATUM(reconstructedNodeKmer);
            out->nNodes++;
        }
    }

    // Free temporary sequence
    pfree(fullReconstructedSequence);

    PG_RETURN_VOID();
}

PG_FUNCTION_INFO_V1(spgist_kmer_leaf_consistent);
Datum
spgist_kmer_leaf_consistent(PG_FUNCTION_ARGS)
{
    spgLeafConsistentIn *in = (spgLeafConsistentIn *) PG_GETARG_POINTER(0);
    spgLeafConsistentOut *out = (spgLeafConsistentOut *) PG_GETARG_POINTER(1);

    int level = in->level;
    Kmer *leafKmer = DatumGetKmerP(in->leafDatum);
    Kmer *reconstructedKmer = NULL;
    char *fullSequence = NULL;
    int fullLength;

    // No need for rescanning
    out->recheck = false;

    // Retrieve reconstructed value from parent level
    if (DatumGetPointer(in->reconstructedValue))
    {
        reconstructedKmer = DatumGetKmerP(in->reconstructedValue);
        Assert(reconstructedKmer->length == level);
    }

    // Reconstruct full sequence
    fullLength = level + leafKmer->length;
    fullSequence = palloc0(fullLength + 1);

    if (reconstructedKmer)
        memcpy(fullSequence, kmer_to_str(reconstructedKmer), level);

    if (leafKmer->length > 0)
        memcpy(fullSequence + level, kmer_to_str(leafKmer), leafKmer->length);

    // Create full kmer for leaf value
    Kmer *fullKmer = kmer_make(fullSequence);
    out->leafValue = FORM_KMER_DATUM(fullKmer);

    // Perform consistency checks
    bool result = true;
    for (int i = 0; i < in->nkeys; i++)
    {
        StrategyNumber strategy = in->scankeys[i].sk_strategy;
        Kmer *queryKmer = DatumGetKmerP(in->scankeys[i].sk_argument);
        const char *querySequence = kmer_to_str(queryKmer);
        int queryLength = queryKmer->length;

        int comparisonResult = strncmp(fullSequence, querySequence,
                                       Min(fullLength, queryLength));

        // Adjust for equal prefix but different total lengths
        if (comparisonResult == 0 && fullLength != queryLength)
            comparisonResult = (fullLength < queryLength) ? -1 : 1;

        switch (strategy)
        {
            case RTPrefixStrategyNumber:
                result = (strncmp(fullSequence, querySequence, queryLength) == 0);
                break;
            case BTLessStrategyNumber:
                result = (comparisonResult < 0);
                break;
            case BTLessEqualStrategyNumber:
                result = (comparisonResult <= 0);
                break;
            case BTEqualStrategyNumber:
                result = (comparisonResult == 0);
                break;
            case BTGreaterEqualStrategyNumber:
                result = (comparisonResult >= 0);
                break;
            case BTGreaterStrategyNumber:
                result = (comparisonResult > 0);
                break;
            default:
                elog(ERROR, "Unrecognized strategy number: %d", strategy);
                result = false;
                break;
        }

        if (!result)
            break;
    }

    // Free temporary sequence
    pfree(fullSequence);

    PG_RETURN_BOOL(result);
}