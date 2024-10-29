\echo Use "CREATE EXTENSION adn" to load this file. \quit

/******************************************************************************
 * Input/Output
 ******************************************************************************/

CREATE OR REPLACE FUNCTION adn_in(cstring)
  RETURNS adn
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION adn_out(adn)
  RETURNS cstring
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION adn_recv(internal)
  RETURNS adn
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION adn_send(adn)
  RETURNS bytea
  AS 'MODULE_PATHNAME'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE TYPE adn (
  internallength = 100,
  input          = adn_in,
  output         = adn_out,
  receive        = adn_recv,
  send           = adn_send,
  alignment      = char
);

CREATE OR REPLACE FUNCTION adn(text)
  RETURNS adn
  AS 'MODULE_PATHNAME', 'adn_cast_from_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(adn)
  RETURNS text
  AS 'MODULE_PATHNAME', 'adn_cast_to_text'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE CAST (text as adn) WITH FUNCTION adn(text) AS IMPLICIT;
CREATE CAST (adn as text) WITH FUNCTION text(adn);


/******************************************************************************
 * Constructor
 ******************************************************************************/

CREATE FUNCTION adn_construct(text)
  RETURNS adn
  AS 'MODULE_PATHNAME', 'adn_constructor'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

/******************************************************************************
 * Operators
 ******************************************************************************/

CREATE FUNCTION adn_eq(adn, adn)
  RETURNS boolean
  AS 'MODULE_PATHNAME', 'adn_eq'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;
CREATE FUNCTION adn_ne(adn, adn)
  RETURNS boolean
  AS 'MODULE_PATHNAME', 'adn_ne'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;


CREATE OPERATOR ~= (
  LEFTARG = adn, RIGHTARG = adn,
  PROCEDURE = adn_eq,
  COMMUTATOR = ~=, NEGATOR = <>
);
CREATE OPERATOR <> (
  LEFTARG = adn, RIGHTARG = adn,
  PROCEDURE = adn_ne,
  COMMUTATOR = <>, NEGATOR = ~=
);


CREATE FUNCTION adn_dist(adn, adn)
  RETURNS double precision
  AS 'MODULE_PATHNAME', 'adn_dist'
  LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OPERATOR <-> (
  LEFTARG = adn, RIGHTARG = adn,
  PROCEDURE = adn_dist,
  COMMUTATOR = <->
);
