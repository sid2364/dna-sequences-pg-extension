#!/bin/bash

sudo systemctl status postgresql

sudo make install

sudo -u postgres psql -c "CREATE DATABASE dna" # Too complicated to check if it already exists, just "ignore" the error

sudo -u postgres psql -d dna -c "DROP EXTENSION IF EXISTS dna CASCADE;"
sudo -u postgres psql -d dna -c "CREATE EXTENSION dna;"

sudo -u postgres psql -d dna -f test.sql # Run the tests

sudo -u postgres psql dna