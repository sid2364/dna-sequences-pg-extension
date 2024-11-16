#include "postgres.h"
#include "utils/varlena.h"
#include "varatt.h" // For SET_VARSIZE!
#include "fmgr.h"
#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"
#include <stddef.h>  // Include for offsetof

#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdint.h>

PG_MODULE_MAGIC;

typedef struct __attribute__((packed)) Qkmer
{
    char vl_len_[4];
    int32 length;                      
    uint64_t bit_sequence[FLEXIBLE_ARRAY_MEMBER]; 
} Qkmer;


#define DatumGetQKmerP(X)  ((Qkmer *) DatumGetPointer(X))
#define QKmerPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_QKMER_P(n) DatumGetQKmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_QKMER_P(x) return QKmerPGetDatum(x)


static void encode_qkmer(const char *sequence, uint64_t *bit_sequence, int length) {
    for (int i = 0; i < length; i++) {
        int offset = (i * 4) % 64;     // 4 bits for each symbol
        int index = i / 16;            // The size is larger than for the adn, so we must divide the index by 16

        switch (sequence[i]) {
            case 'A': /* 0000 */ break;
            case 'C': bit_sequence[index] |= ((uint64_t)0x1 << offset); break;
            case 'G': bit_sequence[index] |= ((uint64_t)0x2 << offset); break;
            case 'T': bit_sequence[index] |= ((uint64_t)0x3 << offset); break;
            case 'W': bit_sequence[index] |= ((uint64_t)0x4 << offset); break;
            case 'S': bit_sequence[index] |= ((uint64_t)0x5 << offset); break;
            case 'M': bit_sequence[index] |= ((uint64_t)0x6 << offset); break;
            case 'K': bit_sequence[index] |= ((uint64_t)0x7 << offset); break;
            case 'R': bit_sequence[index] |= ((uint64_t)0x8 << offset); break;
            case 'Y': bit_sequence[index] |= ((uint64_t)0x9 << offset); break;
            case 'B': bit_sequence[index] |= ((uint64_t)0xA << offset); break;
            case 'D': bit_sequence[index] |= ((uint64_t)0xB << offset); break;
            case 'H': bit_sequence[index] |= ((uint64_t)0xC << offset); break;
            case 'V': bit_sequence[index] |= ((uint64_t)0xD << offset); break;
            case 'N': bit_sequence[index] |= ((uint64_t)0xE << offset); break;
            case 'U': bit_sequence[index] |= ((uint64_t)0xF << offset); break;
            default:
                ereport(ERROR, (errmsg("Invalid character in qkmer sequence: %c", sequence[i])));
        }
    }
}

static char* decode_qkmer(const uint64_t *bit_sequence, int length) {
    char *sequence = palloc0(length + 1);  

    for (int i = 0; i < length; i++) {
        int offset = (i * 4) % 64; // 4 bits for each symbol
        int index = i / 16;    // The size is larger than for the adn, so we must divide the index by 16
        uint64_t bits = (bit_sequence[index] >> offset) & 0xF;

        switch (bits) {
            case 0x0: sequence[i] = 'A'; break;
            case 0x1: sequence[i] = 'C'; break;
            case 0x2: sequence[i] = 'G'; break;
            case 0x3: sequence[i] = 'T'; break;
            case 0x4: sequence[i] = 'W'; break;
            case 0x5: sequence[i] = 'S'; break;
            case 0x6: sequence[i] = 'M'; break;
            case 0x7: sequence[i] = 'K'; break;
            case 0x8: sequence[i] = 'R'; break;
            case 0x9: sequence[i] = 'Y'; break;
            case 0xA: sequence[i] = 'B'; break;
            case 0xB: sequence[i] = 'D'; break;
            case 0xC: sequence[i] = 'H'; break;
            case 0xD: sequence[i] = 'V'; break;
            case 0xE: sequence[i] = 'N'; break;
            case 0xF: sequence[i] = 'U'; break;
        }
    }

    return sequence;
}

static bool validate_qkmer_sequence(const char *sequence) {
    if (sequence == NULL || *sequence == '\0') {
        ereport(ERROR, (errmsg("QKmer sequence cannot be empty")));
        return false;
    }
    for (const char *p = sequence; *p; p++) {
        if (strchr("ACGTWSMKRYBDHVNU", *p) == NULL) { // A simpler comparaison
            ereport(ERROR, (errmsg("Invalid character in QKmer sequence: %c", *p)));
            return false;
        }
    }
    return true;
}

static Qkmer *qkmer_make(const char *sequence)
{
    int length = strlen(sequence);

    int num_bits = length * 4;  // 4 bits per symbol
    int bit_length = (num_bits + 63) / 64;  // Number of 64-bit chunks we need, rest will be padded with zeros
    Size qkmer_size = offsetof(Qkmer, bit_sequence) + bit_length * sizeof(uint64_t);

    
    Qkmer *qkmer = (Qkmer *) palloc0(qkmer_size);

    if (sequence != NULL) {
        if (!validate_qkmer_sequence(sequence)) {
            ereport(ERROR, (errmsg("Invalid qkmer sequence: must contain only valid symbols (A, T, C, G, W, S, M, K, R, Y, B, D, H, V, N)")));
            return NULL;
        }
    }

    SET_VARSIZE(qkmer, qkmer_size); 
    qkmer->length = length;

    // Encode the qkmer sequence
    encode_qkmer(sequence, qkmer->bit_sequence, length);

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



PG_FUNCTION_INFO_V1(qkmer_recv);
Datum
qkmer_recv(PG_FUNCTION_ARGS)
{
    StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);

    int length = pq_getmsgint(buf, sizeof(int));
    int num_bits = length * 4;  // Each qkmer symbol is 4 bits
    int bit_length = (num_bits + 63) / 64;  // Number of 64-bit chunks needed, rounding up

    Size qkmer_size = offsetof(Qkmer, bit_sequence) + bit_length * sizeof(uint64_t);

    Qkmer *qkmer = (Qkmer *) palloc0(qkmer_size);
    SET_VARSIZE(qkmer, qkmer_size);
    qkmer->length = length;

    for (int i = 0; i < bit_length; i++) {
        qkmer->bit_sequence[i] = pq_getmsgint64(buf);
    }

    PG_RETURN_POINTER(qkmer);
}

PG_FUNCTION_INFO_V1(qkmer_send);
Datum
qkmer_send(PG_FUNCTION_ARGS)
{
    Qkmer *qkmer = (Qkmer *) PG_GETARG_POINTER(0);
    StringInfoData buf;

    int bit_length = (qkmer->length * 4 + 63) / 64;

    pq_begintypsend(&buf);
    pq_sendint(&buf, qkmer->length, sizeof(int));

    for (int i = 0; i < bit_length; i++) {
        pq_sendint64(&buf, qkmer->bit_sequence[i]);
    }
    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}


PG_FUNCTION_INFO_V1(qkmer_cast_from_text);
Datum
qkmer_cast_from_text(PG_FUNCTION_ARGS)
{
    text *txt = PG_GETARG_TEXT_P(0);
    char *str = DatumGetCString(DirectFunctionCall1(textout, PointerGetDatum(txt)));
    Qkmer *qkmer = qkmer_make(str);  
    PG_RETURN_POINTER(qkmer);
}


PG_FUNCTION_INFO_V1(qkmer_cast_to_text);
Datum
qkmer_cast_to_text(PG_FUNCTION_ARGS)
{
    Qkmer *qkmer = (Qkmer *) PG_GETARG_QKMER_P(0);
    char *result = decode_qkmer(qkmer->bit_sequence, qkmer->length); 
    text *out = (text *) DirectFunctionCall1(textin, PointerGetDatum(result));
    PG_FREE_IF_COPY(qkmer, 0);
    PG_RETURN_TEXT_P(out);
}
/*****************************************************************************/
PG_FUNCTION_INFO_V1(qkmer_constructor);
Datum
qkmer_constructor(PG_FUNCTION_ARGS)
{
    char *sequence = PG_GETARG_CSTRING(0);
    PG_RETURN_POINTER(qkmer_make(sequence));
}/*****************************************************************************/

static char * qkmer_to_str(const Qkmer *qkmer)
{
    if (qkmer->length == 0) {
        return pstrdup("");
    }
    return decode_qkmer(qkmer->bit_sequence, qkmer->length);
}

PG_FUNCTION_INFO_V1(qkmer_out);
Datum
qkmer_out(PG_FUNCTION_ARGS)
{
    Qkmer *qkmer = PG_GETARG_QKMER_P(0);
    char *result = qkmer_to_str(qkmer);
    PG_FREE_IF_COPY(qkmer, 0);
    PG_RETURN_CSTRING(result);
}