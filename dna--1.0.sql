\echo Use "CREATE EXTENSION dna" to load this file. \quit

/******************************************************************************
 * Input/Output
 ******************************************************************************/

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
  alignment      = int
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

/*
CREATE FUNCTION dna_dist(dna, dna)
  RETURNS double precision
  AS 'MODULE_PATHNAME', 'dna_dist'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OPERATOR <-> (
  LEFTARG = dna, RIGHTARG = dna,
  PROCEDURE = dna_dist,
  COMMUTATOR = <->
);
*/

/******************************************************************************
 * Functions
 ******************************************************************************/

 CREATE FUNCTION length(dna)
  RETURNS int 
  AS 'MODULE_PATHNAME', 'length'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

/******************************************************************************
* K-mers
******************************************************************************/

/* First arg is cast from string to DNA with "dna_cast_from_text" directly */
CREATE FUNCTION generate_kmers(dna dna, k int)
RETURNS SETOF text
AS 'MODULE_PATHNAME', 'generate_kmers'
LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmer_equals(text, text) RETURNS boolean
AS 'MODULE_PATHNAME', 'kmer_equals'
LANGUAGE C IMMUTABLE STRICT;

CREATE OPERATOR = (
    LEFTARG = text,
    RIGHTARG = text,
    PROCEDURE = kmer_equals
); /* This is different from dna equals! */

/******************************************************************************
* For Qkmer
******************************************************************************/

/*
 CREATE OR REPLACE FUNCTION qkmer_in(cstring)
  RETURNS qkmer
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION qkmer_out(qkmer)
  RETURNS cstring
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION qkmer_recv(internal)
  RETURNS qkmer
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION qkmer_send(qkmer)
  RETURNS bytea
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

  CREATE TYPE qkmer (
  internallength = variable,
  input          = qkmer_in,
  output         = qkmer_out,
  receive        = qkmer_recv,
  send           = qkmer_send,
  alignment      = int
);

CREATE FUNCTION qkmer_construct(text)
  RETURNS qkmer
  AS 'MODULE_PATHNAME', 'qkmer_constructor'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;


CREATE OR REPLACE FUNCTION qkmer(text)
  RETURNS qkmer
  AS 'MODULE_PATHNAME', 'qkmer_cast_from_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(qkmer)
  RETURNS text
  AS 'MODULE_PATHNAME', 'qkmer_cast_to_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE CAST (text as qkmer) WITH FUNCTION qkmer(text) AS IMPLICIT;
CREATE CAST (qkmer as text) WITH FUNCTION text(qkmer);

*/