\echo Use "CREATE EXTENSION dna" to load this file. \quit

--Input/Output

CREATE OR REPLACE FUNCTION dna_in(cstring)
  RETURNS dna
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION dna_out(dna)
  RETURNS cstring
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION dna_recv(internal)
  RETURNS dna
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION dna_send(dna)
  RETURNS bytea
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE TYPE dna (
  internallength = variable,
  input          = dna_in,
  output         = dna_out,
  receive        = dna_recv,
  send           = dna_send,
  alignment      = int,
  storage        = extended -- Required for passing in large strings, otherwise postgres doesn't allow more than 8kb
);

CREATE OR REPLACE FUNCTION dna(text)
  RETURNS dna
  AS 'MODULE_PATHNAME', 'dna_cast_from_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(dna)
  RETURNS text
  AS 'MODULE_PATHNAME', 'dna_cast_to_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE CAST (text as dna) WITH FUNCTION dna(text) AS IMPLICIT;
CREATE CAST (dna as text) WITH FUNCTION text(dna);


--Constructor

CREATE FUNCTION dna_construct(text)
  RETURNS dna
  AS 'MODULE_PATHNAME', 'dna_constructor'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

--Operators

CREATE FUNCTION equals(dna, dna)
  RETURNS boolean
  AS 'MODULE_PATHNAME', 'equals'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION dna_ne(dna, dna)
  RETURNS boolean
  AS 'MODULE_PATHNAME', 'dna_ne'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;


CREATE OPERATOR ~= (
  LEFTARG = dna, RIGHTARG = dna,
  PROCEDURE = equals,
  COMMUTATOR = ~=, NEGATOR = <>
);

CREATE OPERATOR = (
    LEFTARG = dna, RIGHTARG = dna,
    PROCEDURE = equals,
    COMMUTATOR = ~=, NEGATOR = <>
);

CREATE OPERATOR <> (
  LEFTARG = dna, RIGHTARG = dna,
  PROCEDURE = dna_ne,
  COMMUTATOR = <>, NEGATOR = ~=
);

--Functions

 CREATE FUNCTION length(dna)
  RETURNS bigint
  AS 'MODULE_PATHNAME', 'length'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

--K-mers

--Input/Output functions
CREATE OR REPLACE FUNCTION kmer_in(cstring)
  RETURNS kmer
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_out(kmer)
  RETURNS cstring
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_recv(internal)
  RETURNS kmer
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_send(kmer)
  RETURNS bytea
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

--Type definition
 -- Fixed size: 4 bytes for length + 8 bytes for bit_sequence (total 12 bytes, aligned to 16 bytes)
CREATE TYPE kmer (
  internallength = 16,
  input          = kmer_in,
  output         = kmer_out,
  receive        = kmer_recv,
  send           = kmer_send,
  alignment      = int
);

--Casting Functions
CREATE OR REPLACE FUNCTION kmer(text)
  RETURNS kmer
  AS 'MODULE_PATHNAME', 'kmer_cast_from_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(kmer)
  RETURNS text
  AS 'MODULE_PATHNAME', 'kmer_cast_to_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE CAST (text AS kmer) WITH FUNCTION kmer(text) AS IMPLICIT;
CREATE CAST (kmer AS text) WITH FUNCTION text(kmer);

--Constructor
CREATE FUNCTION kmer_construct(text)
  RETURNS kmer
  AS 'MODULE_PATHNAME', 'kmer_constructor'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

--Operators
CREATE FUNCTION kmer_eq(kmer, kmer)
  RETURNS boolean
  AS 'MODULE_PATHNAME', 'kmer_eq'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION kmer_ne(kmer, kmer)
  RETURNS boolean
  AS 'MODULE_PATHNAME', 'kmer_ne'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OPERATOR ~= (
  LEFTARG = kmer,
  RIGHTARG = kmer,
  PROCEDURE = kmer_eq,
  COMMUTATOR = ~=, NEGATOR = <>
);

CREATE OPERATOR = (
    LEFTARG = kmer,
    RIGHTARG = kmer,
    PROCEDURE = kmer_eq,
    COMMUTATOR = ~=, NEGATOR = <>
);

CREATE OPERATOR <> (
  LEFTARG = kmer, RIGHTARG = kmer,
  PROCEDURE = kmer_ne,
  COMMUTATOR = <>, NEGATOR = ~=
);

--Length
CREATE FUNCTION length(kmer)
  RETURNS int
  AS 'MODULE_PATHNAME', 'kmer_length'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

--Generate K-mers from DNA
--First arg is cast from string to DNA with "dna_cast_from_text" directly
--Returns a set of k-mers (of type kmer!)
CREATE FUNCTION generate_kmers(dna dna, k int)
RETURNS SETOF kmer
AS 'MODULE_PATHNAME', 'generate_kmers'
LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION starts_with(kmer, kmer) RETURNS boolean
AS 'MODULE_PATHNAME', 'starts_with'
LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR ^@ (
    LEFTARG = kmer,
    RIGHTARG = kmer,
    PROCEDURE = starts_with
);

--For defining the hash function for the kmer type
CREATE FUNCTION kmer_hash(kmer)
    RETURNS INTEGER
    AS 'MODULE_PATHNAME', 'kmer_hash'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OPERATOR CLASS kmer_hash_ops
DEFAULT FOR TYPE kmer USING HASH AS
    OPERATOR 1 = (kmer, kmer),
    FUNCTION 1 kmer_hash(kmer);

-- Qkmer type
CREATE FUNCTION qkmer_in(cstring) RETURNS qkmer
    AS 'MODULE_PATHNAME', 'qkmer_in'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION qkmer_out(qkmer) RETURNS cstring
    AS 'MODULE_PATHNAME', 'qkmer_out'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION qkmer_recv(internal) RETURNS qkmer
    AS 'MODULE_PATHNAME', 'qkmer_recv'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION qkmer_send(qkmer) RETURNS bytea
    AS 'MODULE_PATHNAME', 'qkmer_send'
    LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE qkmer (
    INTERNALLENGTH = variable,
    INPUT = qkmer_in,
    OUTPUT = qkmer_out,
    RECEIVE = qkmer_recv,
    SEND = qkmer_send,
    ALIGNMENT = int
);

CREATE FUNCTION length(qkmer) RETURNS int
    AS 'MODULE_PATHNAME', 'qkmer_length'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION equals(qkmer, qkmer) RETURNS boolean
    AS 'MODULE_PATHNAME', 'qkmer_eq'
     LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION qkmer(text)
  RETURNS qkmer
  AS 'MODULE_PATHNAME', 'qkmer_cast_from_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(qkmer)
    RETURNS text
    AS 'MODULE_PATHNAME', 'qkmer_cast_to_text'
    LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE CAST (text AS qkmer) WITH FUNCTION qkmer(text) AS IMPLICIT;
CREATE CAST (qkmer AS text) WITH FUNCTION text(qkmer);

CREATE OPERATOR = (
    LEFTARG = qkmer,
    RIGHTARG = qkmer,
    PROCEDURE = equals
);

--For Qkmer pattern search
CREATE FUNCTION contains(qkmer, kmer) RETURNS boolean
AS 'MODULE_PATHNAME', 'contains'
LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR @> (
    LEFTARG = qkmer,
    RIGHTARG = kmer,
    PROCEDURE = contains
);

-- For the SpGiST index
CREATE FUNCTION get_oid(kmer)
  RETURNS int
  AS 'MODULE_PATHNAME', 'get_oid'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION spgist_kmer_config(internal, internal) RETURNS void
    AS 'MODULE_PATHNAME', 'spgist_kmer_config'
    LANGUAGE C IMMUTABLE;

CREATE FUNCTION spgist_kmer_choose(internal, internal) RETURNS void
    AS 'MODULE_PATHNAME', 'spgist_kmer_choose'
    LANGUAGE C IMMUTABLE;

CREATE FUNCTION spgist_kmer_picksplit(internal, internal) RETURNS void
    AS 'MODULE_PATHNAME', 'spgist_kmer_picksplit'
    LANGUAGE C IMMUTABLE;

CREATE FUNCTION spgist_kmer_inner_consistent(internal, internal) RETURNS void
    AS 'MODULE_PATHNAME', 'spgist_kmer_inner_consistent'
    LANGUAGE C IMMUTABLE;

CREATE FUNCTION spgist_kmer_leaf_consistent(internal, internal) RETURNS bool
    AS 'MODULE_PATHNAME', 'spgist_kmer_leaf_consistent'
    LANGUAGE C IMMUTABLE;

CREATE OPERATOR CLASS spgist_kmer_ops
DEFAULT FOR TYPE kmer USING spgist AS
    OPERATOR 1 = (kmer, kmer),
    OPERATOR 2 ^@ (kmer, kmer),
    OPERATOR 3 @> (qkmer, kmer), -- DOES NOT WORK
    FUNCTION 1 spgist_kmer_config(internal, internal),
    FUNCTION 2 spgist_kmer_choose(internal, internal),
    FUNCTION 3 spgist_kmer_picksplit(internal, internal),
    FUNCTION 4 spgist_kmer_inner_consistent(internal, internal),
    FUNCTION 5 spgist_kmer_leaf_consistent(internal, internal),
    STORAGE kmer;

