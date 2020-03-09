#!/bin/bash

# init_remote.sh
# Script to initialise the directory structure and run scripts in a remote 
# machine by copying over the code directory

CODE_DIR=$1 # Directory where the original code scripts are
SCRATCH_NAME=$2 # Address of the scratch directory
DEST_DIR_NAME=$3 # Desired name of the directory 

# Makes the directories in the scratch folder
mkdir ${SCRATCH_NAME} 
mkdir ${SCRATCH_NAME}/${DEST_DIR_NAME}
mkdir ${SCRATCH_NAME}/${DEST_DIR_NAME}/raw_data

# Copies over the code
cp -r ${CODE_DIR} ${SCRATCH_NAME}/${DEST_DIR_NAME}
