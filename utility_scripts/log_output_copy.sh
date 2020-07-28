#!/bin/bash

# log_output_copy.sh

REMOTE_MACHINE=$1 # Name of compute machine e.g. nocturne, vanisher, ...
REMOTE_PARENT_DIR=$2 # Parent directory in compute machine e.g. /scratch
SUB_DIR_NAME=$3 # Sub-directory where data is stored (i.e. run_name)
LOCAL_PARENT_DIR=$4 # Parent directory in the local machine

# Full directory name on the remote machine
remote_dir=${REMOTE_MACHINE}:${REMOTE_PARENT_DIR}/${SUB_DIR_NAME}

# Full directory name on the local machine
local_dir=${LOCAL_PARENT_DIR}/${SUB_DIR_NAME}

# Makes the directory on the office machine
mkdir ${local_dir}

# Makes the raw data directory on the office machine
mkdir ${local_dir}/raw_data

# Secure copies the log file and code file
scp -r ${remote_dir}/code ${local_dir}
scp -r ${remote_dir}/raw_data/log ${local_dir}/raw_data

# Cleans the log file
./output_clean.sh ${local_dir}


