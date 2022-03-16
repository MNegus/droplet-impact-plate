#!/bin/bash

# run_simulation.sh
# Bash script to run a simulation. It takes two inputs:
# Input 1: Name of the C file (without the .C extension)
# Input 2: Number of threads to run the simulation on (default 1)
# It removes any previous outputs, runs the code and then moves the output into 
# a directory one level up called "raw_data"

# Saves script name, which will also be the name of the directory that the 
#output gets saved for (crucially this does not contain the .c extension)
script_name=$1

# Sets the number of OpenMP threads
export OMP_NUM_THREADS=$2

# Deletes the previous directory and test files
rm -r ${script_name}
rm *.s
rm *.s.d
rm *.tests 
rm *.deps
rm *.tst

# Runs the code using the Makefile
make ${script_name}.tst 

# Remove any existing raw data in the parent directory
rm -r ../raw_data 

# Moves the new data to the parent directory
mv ${script_name} ../raw_data

