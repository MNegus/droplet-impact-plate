#!/bin/bash

code_dir=$1
dest_dir=$2


for MAXLEVEL in 10 11 12 13 14
do
    # Changes the parameters file
    sed -i "/MAXLEVEL/c\const int MAXLEVEL = $MAXLEVEL; // Maximum refinement level" parameters.h

    # Copies over the code 
    ./code_copy.sh $code_dir $dest_dir max_level_$MAXLEVEL
done
