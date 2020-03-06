#!/bin/bash

# Script to copy code over from the home directory to the scratch of one of the 
# remote machines. This is to be run inside the remote machine

HOME_MACHINE=swamp-beast # Default machine name
HOME_DIR=$1 # Home directory location
DEST_DIR=$2 # Destination directory in the remote machine

# Uses rsync to copy over the folder
rsync -aHSv ${HOME_MACHINE}:${HOME_DIR}/. ${DEST_DIR}/.
