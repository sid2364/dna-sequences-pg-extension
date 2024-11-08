#!/bin/bash

sudo systemctl status postgresql

make
sudo make install

sudo -u postgres psql -d dna -c "DROP EXTENSION IF EXISTS dna CASCADE;"
sudo -u postgres psql -d dna -c "CREATE EXTENSION dna;"

sudo -u postgres psql kmer
sudo -u postgres psql dna
