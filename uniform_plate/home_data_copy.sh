#!/bin/bash

# home_data_copy.sh 

# Copies results from the scratch of the office machine into the home machine's
# external drive

OFFICE_DIR=$1
HOME_DIR=$2

scp office-machine:${OFFICE_DIR} ${HOME_DIR}