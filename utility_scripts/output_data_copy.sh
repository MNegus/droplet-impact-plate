#!/bin/bash

# output_data_copy.sh

# Copies output data (and the corresponding parameters file) from a remote
# machine to a different location (typically a long-term storage location like a
# hard-drive). Script designed to run from the local machine where the data is 
# to be copied to, and uses scp to copy from a remote machine

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

# Secure copies the parameters file and raw data
scp -r ${remote_dir}/{code/parameters.h,raw_data} ${local_dir}
