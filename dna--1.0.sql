\echo Use "CREATE EXTENSION dna" to load this file. \quit

/******************************************************************************
 * Input/Output
 ******************************************************************************/

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

CREATE TYPE kmer (
  internallength = 100,
  input          = kmer_in,
  output         = kmer_out,
  receive        = kmer_recv,
  send           = kmer_send,
  alignment      = char
);

CREATE OR REPLACE FUNCTION kmer(text)
  RETURNS kmer
  AS 'MODULE_PATHNAME', 'kmer_cast_from_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(kmer)
  RETURNS text
  AS 'MODULE_PATHNAME', 'kmer_cast_to_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE CAST (text as kmer) WITH FUNCTION kmer(text) AS IMPLICIT;
CREATE CAST (kmer as text) WITH FUNCTION text(kmer);

CREATE FUNCTION kmer_construct(text)
  RETURNS kmer
  AS 'MODULE_PATHNAME', 'kmer_constructor'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION kmer_eq(kmer, kmer)
  RETURNS boolean
  AS 'MODULE_PATHNAME', 'kmer_eq'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION kmer_ne(kmer, kmer)
  RETURNS boolean
  AS 'MODULE_PATHNAME', 'kmer_ne'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;


CREATE OPERATOR ~= (
  LEFTARG = kmer, RIGHTARG = kmer,
  PROCEDURE = kmer_eq,
  COMMUTATOR = ~=, NEGATOR = <>
);
CREATE OPERATOR <> (
  LEFTARG = kmer, RIGHTARG = kmer,
  PROCEDURE = kmer_ne,
  COMMUTATOR = <>, NEGATOR = ~=
);


CREATE FUNCTION kmer_dist(kmer, kmer)
  RETURNS double precision
  AS 'MODULE_PATHNAME', 'kmer_dist'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OPERATOR <-> (
  LEFTARG = kmer, RIGHTARG = kmer,
  PROCEDURE = kmer_dist,
  COMMUTATOR = <->
);

CREATE FUNCTION length(kmer)
  RETURNS int 
  AS 'MODULE_PATHNAME', 'kmer_length'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

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
  internallength = 100,
  input          = dna_in,
  output         = dna_out,
  receive        = dna_recv,
  send           = dna_send,
  alignment      = char
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


/******************************************************************************
 * Constructor
 ******************************************************************************/

CREATE FUNCTION dna_construct(text)
  RETURNS dna
  AS 'MODULE_PATHNAME', 'dna_constructor'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

/******************************************************************************
 * Operators
 ******************************************************************************/

CREATE FUNCTION dna_eq(dna, dna)
  RETURNS boolean
  AS 'MODULE_PATHNAME', 'dna_eq'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE FUNCTION dna_ne(dna, dna)
  RETURNS boolean
  AS 'MODULE_PATHNAME', 'dna_ne'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;


CREATE OPERATOR ~= (
  LEFTARG = dna, RIGHTARG = dna,
  PROCEDURE = dna_eq,
  COMMUTATOR = ~=, NEGATOR = <>
);
CREATE OPERATOR <> (
  LEFTARG = dna, RIGHTARG = dna,
  PROCEDURE = dna_ne,
  COMMUTATOR = <>, NEGATOR = ~=
);


CREATE FUNCTION dna_dist(dna, dna)
  RETURNS double precision
  AS 'MODULE_PATHNAME', 'dna_dist'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OPERATOR <-> (
  LEFTARG = dna, RIGHTARG = dna,
  PROCEDURE = dna_dist,
  COMMUTATOR = <->
);

/******************************************************************************
 * Functions
 ******************************************************************************/

CREATE FUNCTION length(dna)
  RETURNS int 
  AS 'MODULE_PATHNAME', 'dna_length'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;
