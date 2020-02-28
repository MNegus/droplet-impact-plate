#!/bin/bash
# Bash file to run the simulation

# Set the number of OMP threads
export OMP_NUM_THREADS=4

# Runs the simulation
make droplet_impact.tst
