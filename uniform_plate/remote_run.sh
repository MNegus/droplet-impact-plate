#!/bin/bash

# remote_run.sh
# An overall run script for use on a remote machine. It copies the source code
# from the home directory into scratch, then goes into the remote scratch, runs 
# the code and copies the output back into the local machine's scratch.

###############################################################################
# DEFINITIONS
################################################################################

# Home directory code location
CODE_DIR=~/repos/plate-impact/uniform_plate/code 

# Name of the scratch in the remote computer
REMOTE_SCRATCH_NAME=/scratch/negus

# Name of the directory where the source code is to be copied to
DEST_DIR_NAME=uniform_plate 

# Name of the local machine to copy the results back to
LOCAL_MACHINE=swamp-beast

# Name of the local scratch directory to copy the data into
LOCAL_SCRATCH_NAME=/scratch/uniform_plate


################################################################################
# EXECUTION
################################################################################

# Runs the script for copying over the source code
source ~/repos/plate-impact/remote_scripts/init_remote.sh \
    ${CODE_DIR} ${REMOTE_SCRATCH_NAME} ${DEST_DIR_NAME}

# Moves into the remote source code location
cd ${SCRATCH_NAME}/${DEST_DIR_NAME}/code

# Runs the Basilisk script
./run.sh 

# Copies the data into the local scratch
REMOTE_RAW_DATA=${REMOTE_SCRATCH_NAME}/${DEST_DIR_NAME}/raw_data
source ~/repos/plate-impact/post_processing/data_copy.sh \
    ${LOCAL_MACHINE} ${REMOTE_RAW_DATA} ${LOCAL_SCRATCH_NAME}

# Cleans the data
ssh ${LOCAL_MACHINE} "source ~/repos/plate-impact/uniform_plate/data_cleaning/data_clean.sh"

