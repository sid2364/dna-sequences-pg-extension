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
#include "utils/pg_locale.h"
#include "utils/datum.h"
#include "utils/fmgrprotos.h"
#include "mb/pg_wchar.h"
#include "utils/sortsupport.h"
#include "nodes/nodes.h"
#include "nodes/makefuncs.h"
#include "nodes/parsenodes.h"

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

#define DatumGetKmerP(X)  ((Dna *) DatumGetPointer(X)) // We convert the datum pointer into a dna pointer
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

/**********************************************************************************************************************
* DNA functions
**********************************************************************************************************************/

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
  Dna *dna = PG_GETARG_VARLENA_P(0);
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
  Dna *dna  = PG_GETARG_VARLENA_P(0);
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
    Dna *dna = PG_GETARG_VARLENA_P(0);
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
  Dna *dna1 = PG_GETARG_VARLENA_P(0);
  Dna *dna2 = PG_GETARG_VARLENA_P(1);
  bool result = dna_eq_internal(dna1, dna2);
  PG_FREE_IF_COPY(dna1, 0);
  PG_FREE_IF_COPY(dna2, 1);
  PG_RETURN_BOOL(result);
}

PG_FUNCTION_INFO_V1(length);
Datum
length(PG_FUNCTION_ARGS)
{
    Dna *dna = PG_GETARG_VARLENA_P(0);
    uint64_t length = dna->length;  // Directly get the length field
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_INT32(length);
}

PG_FUNCTION_INFO_V1(dna_ne);
Datum
dna_ne(PG_FUNCTION_ARGS)
{
    Dna *dna1 = PG_GETARG_VARLENA_P(0);
    Dna *dna2 = PG_GETARG_VARLENA_P(1);
    bool result = !dna_eq_internal(dna1, dna2);  // ~ the result of dna_eq_internal, viola!
    PG_FREE_IF_COPY(dna1, 0);
    PG_FREE_IF_COPY(dna2, 1);
    PG_RETURN_BOOL(result);
}


/**********************************************************************************************************************
* Kmer functions
**********************************************************************************************************************/

/**
 * Encoding function for K-mers
 *
 * Store the k-mer in a single 64-bit variable, each nucleotide takes 2 bits
 * Encoding 'A' as 00, 'T' as 01, 'C' as 10, and 'G' as 11
*/
static uint64_t encode_kmer(const char *sequence, int length) {
    uint64_t bit_sequence = 0;  // Single 64-bit variable to store the k-mer!

    if (length <= 0 || length > 32) { // Literally cannot store a k-mer longer than 32 nucleotides
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
    elog(INFO, "Length: %d", length);
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

    if (sequence == NULL || *sequence == '\0') {
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
   //elog(INFO, "Debug: Creating Kmer for sequence '%s' with length %d", sequence, length);

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
        dna = PG_GETARG_VARLENA_P(0); // We know the first argument is a DNA sequence (and not just text)
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

/**********************************************************************************************************************
* Qkmer functions
**********************************************************************************************************************/

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

/**********************************************************************************************************************
* SP-GiST functions
**********************************************************************************************************************/

 /*
 Struct for sorting values in picksplit

 Taken from https://doxygen.postgresql.org/spgtextproc_8c_source.html
 */
typedef struct spgNodePtr
 {
     Datum       d;
     int         i;
     int16       c;
 } spgNodePtr;

typedef struct {
        char nextChar;
        int originalIndex;
        Datum originalDatum;
    } KmerNode;

/*
Find the length of the common prefix of a and b

Taken from https://doxygen.postgresql.org/spgtextproc_8c_source.html
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
 A function to test and get the Oid of the Postgres type kmer

 Required for SP-GiST implementation (and is called from Sql)
*/
PG_FUNCTION_INFO_V1(get_oid);
Datum
get_oid(PG_FUNCTION_ARGS){
    Oid int4_oid = typenameTypeId(NULL, makeTypeName("kmer"));
    PG_RETURN_INT32(int4_oid);
}

static int cmpKmerNode(const void *a, const void *b)
{
    KmerNode *nodeA = (KmerNode *) a;
    KmerNode *nodeB = (KmerNode *) b;
    return (int) (nodeA->nextChar - nodeB->nextChar);
}

/*
This function is used to find the common prefix between two strings
*/
//char* commonPrefix(char* str1, char* str2, int minLength) {
//    static char prefix[33]; // We already know that the size of a kmer is 32, we add one more place so that we can put the '\0' symbol at the end of it
//    int i = 0;
//    while (str1[i] == str2[i] && i < minLength) {
//        prefix[i] = str1[i];
//        i++;
//    }
//    prefix[i] = '\0';
//    return prefix;
//}

/*
Modelled from https://doxygen.postgresql.org/spgtextproc_8c_source.html
*/
static int cmpNodePtr(const void *a, const void *b)
{
    const spgNodePtr *aa = (const spgNodePtr *) a;
    const spgNodePtr *bb = (const spgNodePtr *) b;

    return (int) (aa->c - bb->c);
}

/*

Returns static information about the index implementation, including the data type OIDs of the prefix and node label data types
*/
PG_FUNCTION_INFO_V1(spgist_kmer_config);
Datum
spgist_kmer_config(PG_FUNCTION_ARGS)
{
    spgConfigIn  *cfgin  = (spgConfigIn *) PG_GETARG_POINTER(0);
    spgConfigOut *cfgout = (spgConfigOut *) PG_GETARG_POINTER(1);

    Oid KMEROID = typenameTypeId(NULL, makeTypeName("kmer"));
    elog(INFO, "Debug: Retrieved KMEROID = %u", KMEROID);

    cfgout->prefixType = KMEROID;
    cfgout->leafType = KMEROID;
    cfgout->labelType = VOIDOID;

    cfgout->canReturnData = true; // A node may contain already the value that we want without it being a leaf
    cfgout->longValuesOK = true;   // ¯\_(ツ)_/¯ (change back to false later)

    PG_RETURN_VOID();
}

/*
This function decides how to explore the tree when we have the spgChooseIn object containing the value to index

Determines which child node an inserted value should go to, or signals the creation of a new child node if none make sense
*/
//PG_FUNCTION_INFO_V1(spgist_dna_choose);
//Datum
//spgist_dna_choose(PG_FUNCTION_ARGS)
//{
//    spgChooseIn  *in  = (spgChooseIn *) PG_GETARG_POINTER(0);
//    spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);
//
//    Kmer *datum = (Kmer *) DatumGetPointer(in->datum); // Get the Kmer value to index
//    char *datum_str = kmer_to_str(datum); // Convert Kmer to string for easier comparison
//    elog(INFO, "Debug: Indexing kmer: %s", datum_str);
//
//    char *prefix = NULL; // The prefix or text contained in the node (for now)
//    if (in->hasPrefix) {
//        Kmer *prefix_kmer = (Kmer *) DatumGetPointer(in->prefixDatum);
//        prefix = kmer_to_str(prefix_kmer);
//        elog(INFO, "Debug: Node prefix: %s", prefix);
//    }
//
//    char *restDatum = datum_str; // Initialize restDatum with the full datum_str
//    // If there is a prefix, check if the datum_str starts with it
//    if (prefix != NULL) {
//        size_t prefix_len = strlen(prefix);
//        if (strncmp(datum_str, prefix, prefix_len) == 0) {
//            // restDatum is the part of the datum_str after the prefix
//            restDatum = datum_str + prefix_len;
//        } else {
//            // Prefix does not match, need to add a new node
//            elog(INFO, "Debug: Prefix mismatch. Adding new node.");
//            out->resultType = spgAddNode;
//            out->result.addNode.nodeLabel = CStringGetTextDatum(restDatum); // The entire restDatum
//            out->result.addNode.nodeN = in->nNodes; // New node index
//            PG_RETURN_VOID(); // Et voila!
//        }
//    }
//
//    // If restDatum is empty, we've reached a leaf node
//    if (*restDatum == '\0') {
//        elog(INFO, "Debug: Reached a leaf node.");
//        out->resultType = spgMatchNode; // We have a match!
//        PG_RETURN_VOID();
//    }
//
//    // Get the next character after the prefix
//    char nextChar[2]; // It's 2 because we need to store the restDatum[0] and the '\0' character
//    elog(INFO, "Debug: Next character after prefix: %s", nextChar);
//
//    nextChar[0] = *restDatum;
//    nextChar[1] = '\0';
//
//    int i;
//    bool matchFound = false;
//
//    // Search for a child node that matches the next character
//    for (i = 0; i < in->nNodes; i++) {
//        char *nodeLabel = TextDatumGetCString(in->nodeLabels[i]);
//        if (strcmp(nodeLabel, nextChar) == 0) {
//            // Found a matching child node
//            elog(INFO, "Debug: Found matching child node for %s", nextChar);
//
//            out->resultType = spgMatchNode; // We have a match!
//            out->result.matchNode.nodeN = i; // Node index
//            out->result.matchNode.levelAdd = 1;
//            out->result.matchNode.restDatum = CStringGetTextDatum(restDatum + 1); // Skip the already matched character
//            matchFound = true;
//            break;
//        }
//    }
//
//    if (!matchFound) { // Should be always False if we reach this point, but let's be explicit
//        // No matching child node, need to add a new node
//        elog(INFO, "Debug: No matching node found. Adding new node.");
//
//        out->resultType = spgAddNode;
//        out->result.addNode.nodeLabel = CStringGetTextDatum(nextChar); // We try to add the rest of the entire datum since we did not find any corresponding prefix
//        out->result.addNode.nodeN = in->nNodes; // Found that this needs to be set to the number of nodes + 1!
//    }
//    elog(INFO, "Debug: Returning from spgist_dna_choose");
//
//    PG_RETURN_VOID();
//}

 static bool
 searchChar(Datum *nodeLabels, int nNodes, int16 c, int *i)
 {
     int         StopLow = 0,
                 StopHigh = nNodes;

     while (StopLow < StopHigh)
     {
         int         StopMiddle = (StopLow + StopHigh) >> 1;
         int16       middle = DatumGetInt16(nodeLabels[StopMiddle]);

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

 Datum
 spg_kmer_choose(PG_FUNCTION_ARGS)
 {
     spgChooseIn *in = (spgChooseIn *) PG_GETARG_POINTER(0);
     spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);
     Kmer       *inText = DatumGetKmerP(in->datum);
     char       *inStr = VARDATA_ANY(inText);
     int         inSize = VARSIZE_ANY_EXHDR(inText);
     char       *prefixStr = NULL;
     int         prefixSize = 0;
     int         commonLen = 0;
     int16       nodeChar = 0;
     int         i = 0;

     /* Check for prefix match, set nodeChar to first byte after prefix */
     if (in->hasPrefix)
     {
         text       *prefixText = DatumGetTextPP(in->prefixDatum);

         prefixStr = VARDATA_ANY(prefixText);
         prefixSize = VARSIZE_ANY_EXHDR(prefixText);

         commonLen = commonPrefix(inStr + in->level,
                                  prefixStr,
                                  inSize - in->level,
                                  prefixSize);

         if (commonLen == prefixSize)
         {
             if (inSize - in->level > commonLen)
                 nodeChar = *(unsigned char *) (inStr + in->level + commonLen);
             else
                 nodeChar = -1;
         }
         else
         {
             /* Must split tuple because incoming value doesn't match prefix */
             out->resultType = spgSplitTuple;

             if (commonLen == 0)
             {
                 out->result.splitTuple.prefixHasPrefix = false;
             }
             else
             {
                 out->result.splitTuple.prefixHasPrefix = true;
                 out->result.splitTuple.prefixPrefixDatum =
                     formTextDatum(prefixStr, commonLen);
             }
             out->result.splitTuple.prefixNNodes = 1;
             out->result.splitTuple.prefixNodeLabels =
                 (Datum *) palloc(sizeof(Datum));
             out->result.splitTuple.prefixNodeLabels[0] =
                 Int16GetDatum(*(unsigned char *) (prefixStr + commonLen));

             out->result.splitTuple.childNodeN = 0;

             if (prefixSize - commonLen == 1)
             {
                 out->result.splitTuple.postfixHasPrefix = false;
             }
             else
             {
                 out->result.splitTuple.postfixHasPrefix = true;
                 out->result.splitTuple.postfixPrefixDatum =
                     formTextDatum(prefixStr + commonLen + 1,
                                   prefixSize - commonLen - 1);
             }

             PG_RETURN_VOID();
         }
     }
     else if (inSize > in->level)
     {
         nodeChar = *(unsigned char *) (inStr + in->level);
     }
     else
     {
         nodeChar = -1;
     }

     /* Look up nodeChar in the node label array */
     if (searchChar(in->nodeLabels, in->nNodes, nodeChar, &i))
     {
         /*
          * Descend to existing node.  (If in->allTheSame, the core code will
          * ignore our nodeN specification here, but that's OK.  We still have
          * to provide the correct levelAdd and restDatum values, and those are
          * the same regardless of which node gets chosen by core.)
          */
         int         levelAdd;

         out->resultType = spgMatchNode;
         out->result.matchNode.nodeN = i;
         levelAdd = commonLen;
         if (nodeChar >= 0)
             levelAdd++;
         out->result.matchNode.levelAdd = levelAdd;
         if (inSize - in->level - levelAdd > 0)
             out->result.matchNode.restDatum =
                 formTextDatum(inStr + in->level + levelAdd,
                               inSize - in->level - levelAdd);
         else
             out->result.matchNode.restDatum =
                 formTextDatum(NULL, 0);
     }
     else if (in->allTheSame)
     {
         /*
          * Can't use AddNode action, so split the tuple.  The upper tuple has
          * the same prefix as before and uses a dummy node label -2 for the
          * lower tuple.  The lower tuple has no prefix and the same node
          * labels as the original tuple.
          *
          * Note: it might seem tempting to shorten the upper tuple's prefix,
          * if it has one, then use its last byte as label for the lower tuple.
          * But that doesn't win since we know the incoming value matches the
          * whole prefix: we'd just end up splitting the lower tuple again.
          */
         out->resultType = spgSplitTuple;
         out->result.splitTuple.prefixHasPrefix = in->hasPrefix;
         out->result.splitTuple.prefixPrefixDatum = in->prefixDatum;
         out->result.splitTuple.prefixNNodes = 1;
         out->result.splitTuple.prefixNodeLabels = (Datum *) palloc(sizeof(Datum));
         out->result.splitTuple.prefixNodeLabels[0] = Int16GetDatum(-2);
         out->result.splitTuple.childNodeN = 0;
         out->result.splitTuple.postfixHasPrefix = false;
     }
     else
     {
         /* Add a node for the not-previously-seen nodeChar value */
         out->resultType = spgAddNode;
         out->result.addNode.nodeLabel = Int16GetDatum(nodeChar);
         out->result.addNode.nodeN = i;
     }

     PG_RETURN_VOID();
 }

//
//PG_FUNCTION_INFO_V1(spgist_dna_choose);
//Datum
//spgist_dna_choose(PG_FUNCTION_ARGS)
//{
//    spgChooseIn  *in  = (spgChooseIn *) PG_GETARG_POINTER(0);
//    spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);
//
//    // Decode the incoming k-mer using kmer_make
//    Kmer *datum = (Kmer *) DatumGetPointer(in->datum);
//    char *datum_str = kmer_to_str(datum);
//    elog(INFO, "Debug: Indexing kmer: %s", datum_str);
//
//    Kmer *prefix_kmer = NULL;
//    char *prefix_str = NULL;
//    if (in->hasPrefix) {
//        prefix_kmer = (Kmer *) DatumGetPointer(in->prefixDatum);
//        prefix_str = kmer_to_str(prefix_kmer);
//        elog(INFO, "Debug: Node prefix: %s", prefix_str);
//    }
//
//    // Determine the "rest" of the k-mer after the prefix
//    const char *restDatum = datum_str;
//    if (prefix_str != NULL) {
//        size_t prefix_len = strlen(prefix_str);
//        elog(INFO, "Debug: Prefix length: %zu", prefix_len);
//
//        if (strncmp(datum_str, prefix_str, prefix_len) == 0) {
//            // If the prefix matches, restDatum is the remainder
//            restDatum = datum_str + prefix_len;
//            elog(INFO, "Debug: Rest of datum after prefix: %s", restDatum);
//        } else {
//            // Prefix mismatch: split the node
//            elog(INFO, "Debug: Prefix mismatch. Splitting node.");
//            out->resultType = spgSplitTuple;
//
//            // Set the new prefix for the split node
//            if (prefix_len > 0) {
//                out->result.splitTuple.prefixHasPrefix = true;
//                out->result.splitTuple.prefixPrefixDatum = PointerGetDatum(kmer_make(prefix_str));
//            } else {
//                out->result.splitTuple.prefixHasPrefix = false;
//            }
//
//            // Set child node label to the first non-matching character
//            out->result.splitTuple.prefixNNodes = 1;
//            char split_label = datum_str[prefix_len];
//            out->result.splitTuple.prefixNodeLabels = (Datum *) palloc(sizeof(Datum));
//            out->result.splitTuple.prefixNodeLabels[0] = PointerGetDatum(kmer_make(&split_label));
//
//            // Handle the "postfix" part of the split node
//            if (strlen(prefix_str) > prefix_len + 1) {
//                out->result.splitTuple.postfixHasPrefix = true;
//                out->result.splitTuple.postfixPrefixDatum = PointerGetDatum(kmer_make(prefix_str + prefix_len + 1));
//            } else {
//                out->result.splitTuple.postfixHasPrefix = false;
//            }
//
//            out->result.splitTuple.childNodeN = 0;
//            PG_RETURN_VOID();
//        }
//    }
//
//    // If the rest of the k-mer is empty, we've reached a leaf node
//    if (*restDatum == '\0') {
//        elog(INFO, "Debug: Reached a leaf node.");
//        out->resultType = spgMatchNode;
//        PG_RETURN_VOID();
//    }
//
//    // Determine the next character in the "rest" of the k-mer
//    char nextChar[2] = {*restDatum, '\0'};
//    elog(INFO, "Debug: Next character after prefix: %s", nextChar);
//
//    // Search for a child node matching the next character
//    int i;
//    bool matchFound = false;
//    for (i = 0; i < in->nNodes; i++) {
//        elog(INFO, "Debug: Comparing node label %d", i);
//
//        // Log the raw value of in->nodeLabels
//        elog(INFO, "Debug: Raw in->nodeLabels[%d] = %p", i, in->nodeLabels[i]);
//
//        // Validate that the nodeLabel is not null or invalid before dereferencing
//        if (in->nodeLabels[i] == 0) {
//            elog(ERROR, "Debug: Null node label encountered at index %d", i);
//            continue;
//        }
//
//        Kmer *nodeKmer = (Kmer *) DatumGetPointer(in->nodeLabels[i]);
//        elog(INFO, "Debug: Node label (as Kmer): %p", nodeKmer);
//        //elog(INFO, "Debug: Node label (as string): %s", kmer_to_str(nodeKmer));
//        if (nodeKmer == NULL) {
//            elog(ERROR, "Debug: Invalid Kmer pointer for node label at index %d", i);
//            continue;
//        }
//        elog(INFO, "Debug: Node label length: %d", nodeKmer->length);
//        if (nodeKmer->length < 0 || nodeKmer->length > 32) {
//            elog(INFO, "Debug: Invalid Kmer length %d at index %d. Skipping node.",
//                 nodeKmer->length, i);
//            continue;
//        }
//
//        // Use kmer_to_str safely to convert the Kmer to a string
//        char *nodeLabelStr = kmer_to_str(nodeKmer);
//        if (nodeLabelStr == NULL) {
//            elog(ERROR, "Debug: Failed to convert Kmer to string at index %d", i);
//            continue;
//        }
//
//        // Log the converted node label
//        elog(INFO, "Debug: Node label (as string): %s", nodeLabelStr);
//
//        if (strcmp(nodeLabelStr, nextChar) == 0) {
//            matchFound = true;
//            break;
//        }
//    }
//
//    if (matchFound) {
//        // Found a matching child node
//        elog(INFO, "Debug: Found matching child node for %s", nextChar);
//        out->resultType = spgMatchNode;
//        out->result.matchNode.nodeN = i;
//        out->result.matchNode.levelAdd = 1;
//        out->result.matchNode.restDatum = PointerGetDatum(kmer_make(restDatum + 1));
//    } else {
//        // No matching child node: add a new node
//        elog(INFO, "Debug: No matching node found. Adding new node for %s", nextChar);
//        out->resultType = spgAddNode;
//        out->result.addNode.nodeLabel = PointerGetDatum(kmer_make(nextChar));
//        out->result.addNode.nodeN = in->nNodes; // New node index
//    }
//
//    PG_RETURN_VOID();
//}

/*
Decides how to create a new inner tuple over a set of leaf tuples.

Partitions a set of input values into subsets for new child nodes and assigns prefixes and labels for the inner node
*/
//PG_FUNCTION_INFO_V1(kmer_picksplit);
//Datum kmer_picksplit(PG_FUNCTION_ARGS) {
//    spgPickSplitIn *in = (spgPickSplitIn *) PG_GETARG_POINTER(0);
//    spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);
//
//    int i;
//    Kmer **data = (Kmer **) palloc(in->nTuples * sizeof(Kmer *)); // Will contain the content of the different leaves
//    char **kmer_strings = (char **) palloc(in->nTuples * sizeof(char *)); // Will contain the string representation (!) of the different leaves
//    int minLength = INT_MAX; // Needs to be large enough to be replaced by the first value :) Could also be 33, but let's be dramatic
//
//    elog(INFO, "Debug: Number of tuples: %d", in->nTuples);
//
//    // We store the content of the different tuples
//    for (i = 0; i < in->nTuples; i++) {
//        data[i] = (Kmer *) DatumGetPointer(in->datums[i]);
//        kmer_strings[i] = kmer_to_str(data[i]);
//        elog(INFO, "Debug: Leaf tuple %d: %s", i, kmer_strings[i]);
//
//        // We find the size of the smallest size of the different tuples
//        // Because the prefix of the inner node will never have a size greater than the size of the smallest tuple
//        int len = strlen(kmer_strings[i]);
//        if (len < minLength) {
//            minLength = len;
//        }
//    }
//
//    // Find the common prefix with kmer_strings[]; reinventing the wheel a bit since we already have commonPrefix()
//    char *prefix = pstrdup(kmer_strings[0]);
//    prefix[minLength] = '\0'; // Ensure prefix is not longer than the shortest string
//
//    for (i = 1; i < in->nTuples; i++) {
//        int j = 0;
//        while (j < minLength && prefix[j] == kmer_strings[i][j]) {
//            j++;
//        }
//        prefix[j] = '\0';
//        minLength = j;
//        if (minLength == 0) {
//            elog(INFO, "Debug: No common prefix found for this tuple with value %s at index %d", kmer_strings[i], i);
//            elog(INFO, "Debug: Prefix is %s", prefix); // This will be a NULL string
//            break;
//        }
//    }
//
//    if (minLength == 0) {
//        // No common prefix
//        out->hasPrefix = false;
//        out->prefixDatum = (Datum) 0;
//        // Debug: Print all tuples to understand why no common prefix was found
//        elog(INFO, "Debug: No common prefix found for %d tuples. Tuples are:", in->nTuples);
//    } else {
//        out->hasPrefix = true;
//        out->prefixDatum = CStringGetTextDatum(prefix);
//        elog(INFO, "Debug: Common prefix: %s", prefix);
//    }
//
//    // Collect unique next characters after the prefix
//    char *nextChars = (char *) palloc(in->nTuples * sizeof(char));
//    int nUniqueChars = 0;
//    char *uniqueChars = (char *) palloc(256 * sizeof(char)); // Max 4 unique characters (TODO change this back to 4 later)
//    int charToNodeIndex[256]; // Mapping from character to node index
//
//    memset(charToNodeIndex, -1, sizeof(charToNodeIndex)); // Initialize to -1 so we can check if a character has been seen before
//
//    for (i = 0; i < in->nTuples; i++) {
//        char nextChar = (kmer_strings[i][minLength] == '\0') ? '\0' : kmer_strings[i][minLength];
//
//        if (charToNodeIndex[(unsigned char)nextChar] == -1) {
//            // New unique character found
//            if (nUniqueChars >= 4) {
//                elog(ERROR, "Too many unique characters for nodes.");
//            }
//            uniqueChars[nUniqueChars] = nextChar;
//            charToNodeIndex[(unsigned char)nextChar] = nUniqueChars;
//            nUniqueChars++;
//        }
//        nextChars[i] = nextChar;
//    }
//
//    out->nNodes = nUniqueChars;
//    elog(INFO, "Debug: Number of child nodes: %d", nUniqueChars);
//
//    out->nodeLabels = (Datum *) palloc(nUniqueChars * sizeof(Datum));
//    for (i = 0; i < nUniqueChars; i++) {
//        char label[2] = { uniqueChars[i], '\0' }; // Explicit af
//        out->nodeLabels[i] = CStringGetTextDatum(label);
//        elog(INFO, "Debug: Node label: %s", label);
//    }
//
//    // Map tuples to nodes
//    out->mapTuplesToNodes = (int *) palloc(in->nTuples * sizeof(int));
//    for (i = 0; i < in->nTuples; i++) {
//        char nextChar = nextChars[i];
//        elog(INFO, "Debug: Tuple %d next char: %c", i, nextChar);
//        int nodeIndex = charToNodeIndex[(unsigned char) nextChar];
//        elog(INFO, "Debug: Node index: %d", nodeIndex);
//        if (nodeIndex == -1) {
//            elog(ERROR, "Unexpected error mapping tuple %d to node.", i);
//        }
//        out->mapTuplesToNodes[i] = nodeIndex;
//        elog(INFO, "Debug: Tuple %d mapped to node %d", i, nodeIndex);
//    }
//
//    // Clean up allocated memory
//    pfree(uniqueChars);
//    pfree(nextChars);
//    pfree(data);
//    pfree(kmer_strings);
//
//    PG_RETURN_VOID();
//}
PG_FUNCTION_INFO_V1(kmer_picksplit);
Datum
kmer_picksplit(PG_FUNCTION_ARGS)
{
    spgPickSplitIn *in = (spgPickSplitIn *) PG_GETARG_POINTER(0);
    spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);

    elog(INFO, "Debug: Number of tuples: %d", in->nTuples);

    // Extract Kmer pointers
    Kmer **data = (Kmer **) palloc(in->nTuples * sizeof(Kmer *));
    char **kmer_strings = (char **) palloc(in->nTuples * sizeof(char *));

    elog(INFO, "Debug: Number of tuples: %d", in->nTuples);

    // Decode k-mers into strings
    int i, minLength = INT_MAX;
    for (i = 0; i < in->nTuples; i++) {
        data[i] = (Kmer *) DatumGetPointer(in->datums[i]);
        kmer_strings[i] = kmer_to_str(data[i]);
        int len = strlen(kmer_strings[i]);
        if (len < minLength) {
            minLength = len;
        }
    }

    // Find the longest common prefix
    int commonLen = minLength;
    for (i = 1; i < in->nTuples && commonLen > 0; i++) {
        int tmp = commonPrefix(kmer_strings[0], kmer_strings[i], strlen(kmer_strings[0]), strlen(kmer_strings[i]));
        if (tmp < commonLen) {
            commonLen = tmp;
        }
    }

    // Set the prefix
    if (commonLen > 0) {
        out->hasPrefix = true;
        char *prefix = palloc(commonLen + 1);
        strncpy(prefix, kmer_strings[0], commonLen);
        prefix[commonLen] = '\0';
        out->prefixDatum = PointerGetDatum(kmer_make(prefix));
    } else {
        out->hasPrefix = false;
    }

    // Extract node labels (next character after the prefix)
    spgNodePtr *nodes = (spgNodePtr *) palloc(sizeof(spgNodePtr) * in->nTuples);
    for (i = 0; i < in->nTuples; i++) {
        if (commonLen < strlen(kmer_strings[i])) {
            nodes[i].c = kmer_strings[i][commonLen];
        } else {
            nodes[i].c = '\0';
        }
        nodes[i].i = i;
        nodes[i].d = in->datums[i];
    }

    // Sort nodes by label values
    qsort(nodes, in->nTuples, sizeof(spgNodePtr), cmpNodePtr);

    // Populate output
    out->nNodes = 0;
    out->nodeLabels = (Datum *) palloc(sizeof(Datum) * in->nTuples);
    out->mapTuplesToNodes = (int *) palloc(sizeof(int) * in->nTuples);
    out->leafTupleDatums = (Datum *) palloc(sizeof(Datum) * in->nTuples);

    for (i = 0; i < in->nTuples; i++) {
        if (i == 0 || nodes[i].c != nodes[i - 1].c) {
            char label[2] = {nodes[i].c, '\0'};
            out->nodeLabels[out->nNodes] = PointerGetDatum(kmer_make(label));
            out->nNodes++;
        }

        // Determine remaining k-mer after the prefix and the node label
        char *remaining = kmer_strings[nodes[i].i] + commonLen + 1;
        out->leafTupleDatums[nodes[i].i] = PointerGetDatum(kmer_make(remaining));
        out->mapTuplesToNodes[nodes[i].i] = out->nNodes - 1;
    }

    PG_RETURN_VOID();
}



/*
Returns set of nodes (branches) to follow during tree search.

The function is called when we are at an inner node and we need to decide which child nodes to explore.
*/
PG_FUNCTION_INFO_V1(inner_consistent);
Datum inner_consistent(PG_FUNCTION_ARGS) {
    spgInnerConsistentIn *in = (spgInnerConsistentIn *) PG_GETARG_POINTER(0);
    spgInnerConsistentOut *out = (spgInnerConsistentOut *) PG_GETARG_POINTER(1);

    out->nNodes = 0;
    out->nodeNumbers = palloc(in->nNodes * sizeof(int)); // I give huge space for the worst case

    for (int i = 0; i < in->nNodes; i++) {
      char *label = TextDatumGetCString(in->nodeLabels[i]);
      for (int j = 0; j < in->nkeys; j++ ){
        ScanKey key =  &in->scankeys[j]; // If I am correct, the array ScanKey contains a ScanKey which can be used here
        Kmer *query_kmer = (Kmer *) DatumGetPointer(key->sk_argument);
        char *query_str = kmer_to_str(query_kmer);

        if ((key->sk_strategy == BTEqualStrategyNumber) && // Strategy is the equal to operator
                strncmp(query_str, label, strlen(label)) == 0) {
          out->nodeNumbers[out->nNodes] = i; // If our condition is fulfilled then we add the node to explore
          out->nNodes++;
        }
      }
    }
    PG_RETURN_VOID();
}

/*
Checks if a value stored in a leaf node satisfies the query conditions and determines if additional rechecks are needed.
*/
PG_FUNCTION_INFO_V1(leaf_consistent);
Datum
leaf_consistent(PG_FUNCTION_ARGS)
{
    spgLeafConsistentIn *in = (spgLeafConsistentIn *) PG_GETARG_POINTER(0);
    spgLeafConsistentOut *out = (spgLeafConsistentOut *) PG_GETARG_POINTER(1);
    bool result = true;

    out->recheck = false;
    out->recheckDistances = false;
    out->distances = NULL;

    Kmer *leaf_kmer = (Kmer *) DatumGetPointer(in->leafDatum);

    for (int i = 0; i < in->nkeys; i++) {   // We are basically doing the same thing as in the inner_consistent function
        ScanKey key = &in->scankeys[i];
        Kmer *query_kmer = (Kmer *) DatumGetPointer(key->sk_argument);

        if (key->sk_strategy == BTEqualStrategyNumber) {
            if (!kmer_eq_internal(leaf_kmer, query_kmer)) {
                result = false;
                break;
            }
        } // TODO: If we implement comparison wo Qkmer, we might need to add cases here for that straegy
    }

    PG_RETURN_BOOL(result);
}