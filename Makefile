EXTENSION   = dna
MODULES 	= dna qkmer
DATA        = dna--1.0.sql dna.control
OBJS        = qkmer.o dna.o

PG_CONFIG   ?= pg_config
PGXS        := $(shell $(PG_CONFIG) --pgxs)
INCLUDE_DIR := $(shell $(PG_CONFIG) --includedir-server)

CFLAGS += -I$(INCLUDE_DIR)

include $(PGXS)
