# DNA Sequence PostgreSQL Extension

## Overview

The **DNA Sequence PostgreSQL Extension** provides support for DNA sequence and k-mer analysis directly in PostgreSQL. This extension introduces custom types and functions to efficiently store, retrieve, and manipulate DNA sequences and perform k-mer-based operations.

## Features

- **Custom Types**:
    - `dna`: Stores DNA sequences in a compressed binary format.
    - `kmer`: Represents fixed-length sub-sequences (k-mers) of DNA, optimized for performance.
    - `qkmer` (informal): Facilitates querying k-mers with patterns (IUPAC codes).

- **Efficient Storage**:
    - DNA sequences are encoded using 2 bits per nucleotide (`A`, `T`, `C`, `G`).
    - We use a variable-length encoding for DNA sequences to minimize storage requirements with PostgreSQL's varlena format.
    - K-mers use a fixed 64-bit representation for sequences up to 32 nucleotides. The assumption is that most (useful) k-mers are 32 nucleotides or less.

- **Advanced Query Support**:
    - Generate k-mers from DNA sequences.
    - Query k-mers with IUPAC patterns.
    - Perform k-mer counting, filtering, and aggregation.

## Operators
- **Custom Operators**:
    - DNA and k-mer equality (`=`), inequality (`<>`), and prefix matching (`^@`).
    - Pattern matching with `@>` for `qkmer`.

## Installation

1. Clone the repository:
   ```bash
   git clone git@github.com:sid2364/dna-sequences-pg-extension.git
   cd dna-extension

2. Build and install the extension:
   ```bash
   make
   sudo make install
   ```
3. Load the extension in PostgreSQL:
   ```sql
   CREATE EXTENSION dna;
   ```
   
## Usage
### Bit-Packed DNA Encoding
Below is an example of the size of a DNA sequence stored in the `dna` type:
```sql

SELECT pg_column_size(dna('ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG'));

-- pg_column_size
------------------
--             24
```
The calculation is as follows:
- 16 bytes for the bit sequence of 64 bytes (2 bits per byte/nucleotide base, so 128 bits total, or 16 bytes)
- 4 bytes for the integer length
- 4 bytes for the varlena header `VARHDRSZ`
- Total: 24 bytes

K-mers are stored in a similar way, with a fixed 64-bit representation (a single `uint64`).
### K-mer Generation
```sql
SELECT generate_kmers('ATCGTAGCGT', 3); -- Should return 8 kmers / non-uniques!

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
```
The `generate_kmers()` function generates all possible k-mers of a given length from a DNA sequence. The function returns a set of k-mers, which can be used in queries or further processing. These k-mers are not unique, so the same k-mer can appear multiple times in the result set (which is also useful for some applications).

### K-mer Querying
```sql
SELECT k.kmer FROM generate_kmers('ACTGACGTACC', 3) AS k(kmer) WHERE k.kmer ^@ 'AC';
-- kmer
--------
-- ACT
-- ACG
-- ACC
--(3 rows)
```
This shows the `starts_with()` or `^@` operator in action, which is used to query k-mers that start with a given pattern.

### Contains Query
```sql
SELECT k.kmer FROM generate_kmers('ACGTACGCACGT', 6) AS k(kmer) WHERE 'DNMSRN' @> k.kmer ;
--  kmer
----------
-- GTACGC
-- GCACGT
--(2 rows)
```
This shows the `contains()` or `@>` operator in action, which is used to query k-mers that are contained within a given pattern.

### K-mer Counting
```sql
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
```
This shows the counting of k-mers in a DNA sequence. The `generate_kmers()` function is used to generate all possible k-mers of a given length from a DNA sequence. The result set is then grouped by k-mer and counted. The `GROUP BY` clause requires a hash function to be defined for the `kmer` type (apart from the equality operator).

### More k-mer - Total, Distinct, Unique
```sql
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
```
Here, we calculate the total, distinct, and unique k-mer counts in a DNA sequence. The `generate_kmers()` function is used to generate all possible k-mers of a given length from a DNA sequence. The result set is then grouped by k-mer and counted. The `WITH` clause is used to create a temporary table `kmers` that contains the k-mer and its count. The final query calculates the total count, distinct count, and unique count of k-mers in the DNA sequence. This is specially useful for k-mer analysis!

### DNA Sequence Data
- https://www.ncbi.nlm.nih.gov/nuccore/HQ287898.1
- https://www.ncbi.nlm.nih.gov/genbank/samplerecord/#SequenceLengthA
