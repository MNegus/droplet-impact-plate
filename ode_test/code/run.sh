#!/bin/bash

# Deletes the previous directory and test files
rm -r ode
rm *.s
rm *.s.d
rm *.tests 
rm *.deps
rm *.tst

# Runs the code using the Makefile
make ode.tst 
