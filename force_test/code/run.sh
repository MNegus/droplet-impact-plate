#!/bin/bash

# Deletes the previous directory and test files
rm -r force
rm *.s
rm *.s.d
rm *.tests 
rm *.deps
rm *.tst

# Runs the code using the Makefile
make force.tst 
