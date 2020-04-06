#!/bin/bash

#data_copy.sh

# Copies results from the scratch of a compute machine into the office machine's
# scratch folder

COMPUTE_MACHINE=$1 # Name of compute machine e.g. nocturne, vanisher, ...
COMPUTE_SCRATCH=$2 # Directory of scratch in compute machine e.g. /scratch/negus
SUB_DIR_NAME=$3 # Sub-directory where data is stored e.g. plate_vel_0.1 
OFFICE_SCRATCH=$4 # Directory of scratch in office machine e.g. /scratch/uniform_plate

# Full directory name on the compute machine
compute_dir=${COMPUTE_MACHINE}:${COMPUTE_SCRATCH}/${SUB_DIR_NAME}

# Full directory name on the office machine
office_dir=${OFFICE_SCRATCH}/${SUB_DIR_NAME}

# Makes the directory on the office machine
mkdir ${office_dir}

# Secure copies the parameters file and raw data
scp -r ${compute_dir}/{code/parameters.h,raw_data} ${office_dir}
