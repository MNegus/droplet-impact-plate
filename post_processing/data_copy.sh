#!/bin/bash

# data_copy.sh
# Copies data from a remote machine into the default scratch

HOME_MACHINE_NAME=$1
REMOTE_SCRATCH_DIR=$2
HOME_SCRATCH_DIR=$3

scp -r ${REMOTE_SCRATCH_DIR}/ ${HOME_MACHINE_NAME}:${HOME_SCRATCH_DIR}/