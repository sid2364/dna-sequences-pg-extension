SELECT equals(dna('ATCG'), dna('ATGCG'));
-- equals
----------
-- f
--(1 row)

SELECT equals(dna('ATCG'), dna('ATCG'));
-- equals
----------
-- t
--(1 row)

SELECT equals(dna('ATCG'), dna('GTCA'));
--equals
----------
-- f
--(1 row)


SELECT dna_ne(dna('ATCG'), dna('ATCG'));
-- dna_ne
----------
-- f
--(1 row)

SELECT length(dna('ATCG')); -- Should be 4
-- length
----------
--      4
--(1 row)

SELECT pg_column_size(dna('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'));
-- 64 characters/bytes, should become 24!
-- Calculation =>
-- 16 for the bit sequence of 64 bytes (2 bits per byte/nucleotide base, so 128 bits total, or 16 bytes)
-- + 4 for the int length
-- + 4 for the varlena header VARHDRSZ

-- pg_column_size
------------------
--             24

SELECT length(dna('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'));


SELECT generate_kmers('ATCGTAGCGT', 3); -- Should return 8 kmers / non-uniques!
-- Above is the same as SELECT generate_kmers(dna('ATCGTAGCGT'), 3);
-- generate_kmers
------------------
-- ATC
-- TCG
-- CGT
-- GTA
-- TAG
-- AGC
-- GCG
-- CGT
--(8 rows)


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

SELECT qkmer('ATCGUWSMKRYBDHVN');
--------------------
-- ATCGUWSMKRYBDHVN
--(1 row)

SELECT equals(qkmer('KRYBDHVN'), qkmer('KRYBDHVN'));
-- equals
----------
-- t
--(1 row)

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

-- Just creating a table with some test data
-- CREATE TABLE k AS SELECT kmer FROM generate_kmers('ACGTACGTACGT', 6) AS k(kmer);
DROP TABLE IF EXISTS dna_sequences;
CREATE TABLE dna_sequences (
    id SERIAL PRIMARY KEY,
    sequence dna
);

COPY dna_sequences (sequence)
FROM '/tmp/random_nucleotides.txt'
WITH (FORMAT text);

SELECT id, pg_column_size(sequence) AS size
FROM dna_sequences;

SELECT id, length(sequence) AS length, pg_column_size(sequence) AS size FROM dna_sequences;

--- Now count kmers on the massive table
WITH kmers AS (
    SELECT k.kmer, COUNT(*) AS count
    FROM dna_sequences d,
         LATERAL generate_kmers(d.sequence, 5) AS k(kmer)
    GROUP BY k.kmer
)
SELECT
    SUM(count) AS total_count,
    COUNT(*) AS distinct_count,
    COUNT(*) FILTER (WHERE count = 1) AS unique_count
FROM kmers;

-- test this: whenever you access dna, use VARDATA(dna) and then get use that instead of the pointer directly
-- also check PG_DETOAST_DATUM_COPY(dna)