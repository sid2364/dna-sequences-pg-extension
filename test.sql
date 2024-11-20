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

SELECT generate_kmers('ATCGTAGCGT', 3); -- Should return 8 kmers / non-uniques!
-- Above is the same as SELECT generate_kmers(dna('ATCGTAGCGT'), 3);

SELECT k.kmer FROM generate_kmers('ACGTACGT', 6) AS k(kmer) WHERE k.kmer = 'ACGTAC';
--  kmer
----------
-- ACGTAC
--(1 row)

SELECT k.kmer FROM generate_kmers('ACTGACGTACC', 3) AS k(kmer) WHERE k.kmer ^@ 'AC';
-- kmer
--------
-- ACT
-- ACG
-- ACC
--(3 rows)

SELECT k.kmer FROM generate_kmers('ACGTACGCACGT', 6) AS k(kmer) WHERE 'DNMSRN' @> k.kmer ;
--  kmer
----------
-- GTACGC
-- GCACGT
--(2 rows)

-- K-mer counting!
SELECT k.kmer, count(*) FROM generate_kmers('ATCGATCGATCGATCGACG', 5) AS k(kmer) GROUP BY k.kmer ORDER BY count(*) DESC;
-- kmer  | count
---------+-------
-- ATCGA |     4
-- CGATC |     3
-- GATCG |     3
-- TCGAT |     3
-- TCGAC |     1
-- CGACG |     1
--(6 rows)

-- Just creating a table with some test data
CREATE TABLE k AS SELECT kmer FROM generate_kmers('ACGTACGTACGT', 6) AS k(kmer);

-- Counting total, distinct and unique k-mers in a table
WITH kmers AS (
SELECT k.kmer, count(*)
FROM generate_kmers('ACGTACGTACGTAG', 5) AS k(kmer)
GROUP BY k.kmer
)
SELECT sum(count) AS total_count,
count(*) AS distinct_count,
count(*) FILTER (WHERE count = 1) AS unique_count
FROM kmers;
-- total_count | distinct_count | unique_count
---------------+----------------+--------------
--          10 |              5 |            1
--(1 row)

