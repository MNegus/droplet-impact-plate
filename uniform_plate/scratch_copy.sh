#!/bin/bash

# scratch_copy.sh
# Script to copy the code and data cleaning directories over to the scratch 
# directory.
# The first input is the location of the scratch directory, and the
# second is the desired name for the sub-directory inside the scratch

SCRATCH_DIR=$1 # Scratch directory address
SUB_DIR_NAME=$2 # Name of the sub-directory

# Makes the sub-directory
mkdir ${SCRATCH_DIR}/${SUB_DIR_NAME}

# Copies the code directory over to the scratch
cp -r code ${SCRATCH_DIR}/${SUB_DIR_NAME}

# Copies the data cleaning directory over to the scratch
cp -r data_cleaning ${SCRATCH_DIR}/${SUB_DIR_NAME}