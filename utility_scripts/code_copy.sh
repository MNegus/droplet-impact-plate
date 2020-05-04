#!/bin/bash

# code_copy.sh
# Script to copy the code into the directory it is to be run in (i.e. the 
# /scratch directory)
# The first input is the local directory of the code (e.g. solid_wall, 
# coupled_plate etc), the second directory is the parent destination directory
# and the third input is the sub-directory name (i.e. the name of the directory
# the code will be stored in within the scratch)

LOCAL_DIR=$1 # Local directory of the code
DEST_DIR=$2 # Destination of the parent directory
SUB_DIR_NAME=$3 # Name of the sub-directory

# Makes the sub-directory
mkdir ${DEST_DIR}/${SUB_DIR_NAME}

# Copies the code directory over to the destination
cp -r ${LOCAL_DIR}/code ${DEST_DIR}/${SUB_DIR_NAME}

# Copies the run script over to the destination
cp run_simulation.sh ${DEST_DIR}/${SUB_DIR_NAME}/code

# Copies the Makefile over to the destination
cp Makefile ${DEST_DIR}/${SUB_DIR_NAME}/code
