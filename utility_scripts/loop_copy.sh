#!/bin/bash

code_dir=$1
dest_dir=$2


for AXISYMMETRIC in 0 1 
do
    # Change the axisymmetric value in parameters file
    sed -i "/#define AXISYMMETRIC/c\#define AXISYMMETRIC $AXISYMMETRIC" parameters.h

    # Make the respective directory
    mkdir $dest_dir/axi_$AXISYMMETRIC
    
    for MAXLEVEL in 10 11 12 13 14
    do
        # Changes the parameters file
        sed -i "/MAXLEVEL/c\const int MAXLEVEL = $MAXLEVEL; // Maximum refinement level" parameters.h

        # Copies over the code 
        ./code_copy.sh $code_dir $dest_dir/axi_$AXISYMMETRIC max_level_$MAXLEVEL
    done
done