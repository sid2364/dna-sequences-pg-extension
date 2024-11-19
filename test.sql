SELECT equals(dna('ATCG'), dna('ATGCG'));

SELECT equals(dna('ATCG'), dna('ATCG'));
SELECT equals(dna('ATCG'), dna('GTCA'));

SELECT dna_ne(dna('ATCG'), dna('ATCG'));

SELECT length(dna('ATCG')); -- Should be 4
SELECT pg_column_size(dna('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'));
-- 64 characters/bytes, should become 24!
-- Calculation =>
-- 16 for the bit sequence of 64 bytes (2 bits per byte/nucleotide base, so 128 bits total, or 16 bytes)
-- + 4 for the int length
-- + 4 for the varlena header VARHDRSZ

-- Test creating a kmer
SELECT kmer('ACGTA');

-- Test equality function
SELECT kmer_equals(kmer('ACGTA'), kmer('ACGTA'));  -- Should return true

-- Test inequality function
SELECT kmer_ne(kmer('ACGTA'), kmer('TGCAA'));  -- Should return true

-- Test length function
SELECT kmer_length(kmer('ACGTA'));  -- Should return 5

-- Test casting kmer to text
SELECT kmer_cast_to_text(kmer('ACGTA'));  -- Should return 'ACGTA'

-- Test casting text to kmer
SELECT kmer_cast_from_text('ACGTA'::text);  -- Should return kmer('ACGTA')
