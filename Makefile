EXTENSION   = dna
MODULES 	= dna
DATA        = dna--1.0.sql dna.control

PG_CONFIG   ?= pg_config
PGXS        := $(shell $(PG_CONFIG) --pgxs)
INCLUDE_DIR := $(shell $(PG_CONFIG) --includedir-server)

CFLAGS += -I$(INCLUDE_DIR)

include $(PGXS)
