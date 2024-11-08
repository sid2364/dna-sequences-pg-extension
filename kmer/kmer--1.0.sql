\echo Use "CREATE EXTENSION kmer" to load this file. \quit

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


/******************************************************************************
 * Constructor
 ******************************************************************************/

CREATE FUNCTION kmer_construct(text)
  RETURNS kmer
  AS 'MODULE_PATHNAME', 'kmer_constructor'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

/******************************************************************************
 * Operators
 ******************************************************************************/

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

/******************************************************************************
 * Functions
 ******************************************************************************/

CREATE FUNCTION length(kmer)
  RETURNS int 
  AS 'MODULE_PATHNAME', 'kmer_length'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;
