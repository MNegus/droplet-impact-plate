#!/bin/bash

# data_copy.sh
# Copies data from a remote machine into the default scratch

LOCAL_MACHINE_NAME=$1
REMOTE_SCRATCH_DIR=$2
LOCAL_SCRATCH_DIR=$3

scp -r ${REMOTE_SCRATCH_DIR}/ ${LOCAL_MACHINE_NAME}:${LOCAL_SCRATCH_DIR}/