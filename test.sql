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

SELECT k.kmer FROM generate_kmers('ACGTACGCACGT', 6) AS k(kmer) WHERE 'DNMSRN' @> k.kmer; -- Contains
SELECT k.kmer FROM generate_kmers('ACGTACGCACGT', 6) AS k(kmer) WHERE contains('DNMSRN', k.kmer); -- Same as above
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
FROM generate_kmers('ACGTACGTACGTAG', 8) AS k(kmer)
GROUP BY k.kmer
)
SELECT sum(count) AS total_count,
count(*) AS distinct_count,
count(*) FILTER (WHERE count = 1) AS unique_count
FROM kmers;
-- total_count | distinct_count | unique_count
---------------+----------------+--------------
--           7 |              5 |            3
--(1 row)

-- Create a table we can store DNA sequences in
DROP TABLE IF EXISTS dna_sequences;
CREATE TABLE dna_sequences (
    id SERIAL PRIMARY KEY,
    sequence dna
);

COPY dna_sequences (sequence)
FROM '/tmp/random_nucleotides.txt'
WITH (FORMAT text);

-- Compare the length of actual DNA sequence and the size of the stored sequence once casted to DNA (75% reduction in size!)
SELECT id, length(sequence) AS length, pg_column_size(sequence) AS size FROM dna_sequences;
-- id | length | size
------+--------+-------
--  1 | 100000 | 25012
--(1 row)

--- Now count kmers on the massive table
WITH kmers AS (
    SELECT k.kmer, COUNT(*) AS count
    FROM dna_sequences d,
         LATERAL generate_kmers(d.sequence, 10) AS k(kmer)
    GROUP BY k.kmer
)
SELECT
    SUM(count) AS total_count,
    COUNT(*) AS distinct_count,
    COUNT(*) FILTER (WHERE count = 1) AS unique_count
FROM kmers;
-- total_count | distinct_count | unique_count
---------------+----------------+--------------
--       99991 |          95320 |        90764
--(1 row)

-- Testing the SpGIST index
-- Create a table with k-mers which we can index once it's populated
DROP TABLE IF EXISTS kmer_data_t;
DO $$
BEGIN
    CREATE TABLE IF NOT EXISTS kmer_data_t ( -- Shouldn't exist if we DROPped it, but why not
        id SERIAL PRIMARY KEY,
        kmer_sequence kmer
    );
END $$;

-- Populate kmer_data_t with k-mers generated from MASSIVE sequence in dna_sequences
DO $$
DECLARE
    kmer_record RECORD;
BEGIN
    FOR kmer_record IN
        SELECT k.kmer
        FROM dna_sequences d,
             LATERAL generate_kmers(d.sequence, 5) AS k(kmer) -- Generate k-mers of length 5
    LOOP
        INSERT INTO kmer_data_t (kmer_sequence) VALUES (kmer_record.kmer);
    END LOOP;
END $$;

--SET client_min_messages = INFO;
--SET log_min_messages = INFO;

VACUUM ANALYZE kmer_data_t;

------ First check without the index
EXPLAIN ANALYZE
SELECT * FROM kmer_data_t WHERE kmer_sequence = 'ATCGT';
--                                                  QUERY PLAN
-----------------------------------------------------------------------------------------------------------------
-- Seq Scan on kmer_data_t  (cost=0.00..1886.95 rows=49998 width=20) (actual time=0.009..4.168 rows=116 loops=1)
--   Filter: (kmer_sequence = 'ATCGT'::kmer)
--   Rows Removed by Filter: 99880
-- Planning Time: 0.084 ms
-- Execution Time: 4.181 ms
--(5 rows)
--


-- Create the SP-GiST index
CREATE INDEX spgist_kmer_idx
ON kmer_data_t USING spgist (kmer_sequence spgist_kmer_ops);

-- Disable sequential scan and test the index
SET enable_seqscan = OFF;
EXPLAIN ANALYZE
SELECT * FROM kmer_data_t WHERE kmer_sequence = 'ATCGT';
--                                                           QUERY PLAN
-----------------------------------------------------------------------------------------------------------------------------------
-- Bitmap Heap Scan on kmer_data_t  (cost=1551.76..2813.74 rows=49998 width=20) (actual time=0.046..0.138 rows=102 loops=1)
--   Recheck Cond: (kmer_sequence = 'ATCGT'::kmer)
--   Heap Blocks: exact=96
--   ->  Bitmap Index Scan on spgist_kmer_idx  (cost=0.00..1539.27 rows=49998 width=0) (actual time=0.031..0.031 rows=102 loops=1)
--         Index Cond: (kmer_sequence = 'ATCGT'::kmer)
-- Planning Time: 0.111 ms
-- Execution Time: 0.185 ms
--(7 rows)


-- Enable sequential scan and check starts_with operator
SET enable_seqscan = ON;
EXPLAIN ANALYZE
SELECT * FROM kmer_data_t WHERE kmer_sequence ^@ 'ACTG';
--                                                  QUERY PLAN
-----------------------------------------------------------------------------------------------------------------
-- Seq Scan on kmer_data_t  (cost=0.00..1886.95 rows=49998 width=20) (actual time=0.013..4.620 rows=399 loops=1)
--   Filter: (kmer_sequence ^@ 'ACTG'::kmer)
--   Rows Removed by Filter: 99597
-- Planning Time: 0.032 ms
-- Execution Time: 4.648 ms
--(5 rows)


-- Disable sequential scan and check starts_with operator with index
SET enable_seqscan = OFF;
EXPLAIN ANALYZE
SELECT * FROM kmer_data_t WHERE kmer_sequence ^@ 'ACTG';
--                                                           QUERY PLAN
-----------------------------------------------------------------------------------------------------------------------------------
-- Bitmap Heap Scan on kmer_data_t  (cost=1551.76..2813.74 rows=49998 width=20) (actual time=0.094..0.298 rows=388 loops=1)
--   Recheck Cond: (kmer_sequence ^@ 'ACTG'::kmer)
--   Heap Blocks: exact=282
--   ->  Bitmap Index Scan on spgist_kmer_idx  (cost=0.00..1539.27 rows=49998 width=0) (actual time=0.069..0.069 rows=388 loops=1)
--         Index Cond: (kmer_sequence ^@ 'ACTG'::kmer)
-- Planning Time: 0.037 ms
-- Execution Time: 0.316 ms
--(7 rows)

-- Check qkmer query
SET enable_seqscan = OFF;
EXPLAIN ANALYZE
SELECT * FROM kmer_data_t WHERE 'MRKYN' @> kmer_sequence;
--                                                            QUERY PLAN
-------------------------------------------------------------------------------------------------------------------------------------
-- Seq Scan on kmer_data_t  (cost=10000000000.00..10000001886.95 rows=49998 width=20) (actual time=50.231..54.008 rows=6315 loops=1)
--   Filter: ('MRKYN'::qkmer @> kmer_sequence)
--   Rows Removed by Filter: 93681
-- Planning Time: 0.012 ms
-- JIT:
--   Functions: 2
--   Options: Inlining true, Optimization true, Expressions true, Deforming true
--   Timing: Generation 0.091 ms, Inlining 36.008 ms, Optimization 9.219 ms, Emission 4.989 ms, Total 50.306 ms
-- Execution Time: 64.507 ms
--(9 rows)


-- Check the index usage
SELECT * FROM pg_stat_user_indexes WHERE indexrelname = 'spgist_kmer_idx';
-- relid  | indexrelid | schemaname |   relname   |  indexrelname   | idx_scan |        last_idx_scan         | idx_tup_read | idx_tup_fetch
----------+------------+------------+-------------+-----------------+----------+------------------------------+--------------+---------------
-- 610030 |     610036 | public     | kmer_data_t | spgist_kmer_idx |        2 | 2024-12-03 12:12:56.42183+01 |          490 |             0
--(1 row)

