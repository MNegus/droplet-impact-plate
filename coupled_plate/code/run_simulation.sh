#!/bin/bash

# run_simulation.sh
# Bash script to run the simulation. It removes any previous outputs,
# runs the code and then moves the output into the parent directory
# The first input is the number of threads to run the simulation with

# Sets the number of OpenMP threads
export OMP_NUM_THREADS=$1

# Deletes the previous directory and test files
rm -r cantilever
rm *.s
rm *.s.d
rm *.tests 
rm *.deps
rm *.tst

# Runs the code using the Makefile
make cantilever.tst 

# Remove any existing raw data in the parent directory
rm -r ../raw_data 

# Moves the new data to the parent directory
mv cantilever ../raw_data

