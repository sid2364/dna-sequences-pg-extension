\echo Use "CREATE EXTENSION dna" to load this file.
\quit

/******************************************************************************
 * Input/Output for DNATypes
 ******************************************************************************/

-- DNA Type Input and Output Functions
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

-- K-mer Type Input and Output Functions
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


/******************************************************************************
 * DNA Type
 ******************************************************************************/

CREATE TYPE dna (
  internallength = variable,
  input          = dna_in,
  output         = dna_out,
  receive        = dna_recv,
  send           = dna_send,
  alignment      = int
);

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
 * K-mer Type
 ******************************************************************************/
-- Création du type kmer
CREATE TYPE kmer AS (
    length integer,
    bit_sequence bytea
);

-- Fonction pour créer un kmer à partir d'une séquence
CREATE FUNCTION kmer_in(cstring) RETURNS kmer AS 'MODULE_PATHNAME', 'kmer_in' LANGUAGE C IMMUTABLE STRICT;

-- Fonction pour convertir un kmer en une chaîne de caractères
CREATE FUNCTION kmer_out(kmer) RETURNS cstring AS 'MODULE_PATHNAME', 'kmer_out' LANGUAGE C IMMUTABLE STRICT;

-- Fonction pour convertir un kmer en une chaîne
CREATE FUNCTION kmer_to_string(kmer) RETURNS cstring AS 'MODULE_PATHNAME', 'kmer_to_string' LANGUAGE C IMMUTABLE STRICT;

-- Opérations de comparaison pour kmer
CREATE FUNCTION kmer_eq(kmer, kmer) RETURNS boolean AS 'MODULE_PATHNAME', 'equals' LANGUAGE C IMMUTABLE STRICT;

-- Fonction pour obtenir la longueur d'un kmer
CREATE FUNCTION kmer_length(kmer) RETURNS integer AS 'MODULE_PATHNAME', 'length' LANGUAGE C IMMUTABLE STRICT;


/******************************************************************************
 * For Qkmer Type (Commented Out)
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

-- End of the extension setup
