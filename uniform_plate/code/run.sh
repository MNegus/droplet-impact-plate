#!/bin/bash

RAW_DATA_DIR=/scratch/negus/uniform_plate/raw_data

export OMP_NUM_THREADS=8

# Used to set a special directory name, if wanted
export special_dir_name="oscillation"

# for PLATE_VEL in 0 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5
for PLATE_VEL in 0.1
do
    # Deletes the previous directory and test files
    rm -r uniform_plate
    rm *.s
    rm *.s.d
    rm *.tests 
    rm *.deps
    rm *.tst

    make uniform_plate.tst PLATE_VEL=-$PLATE_VEL

    if [ -z "$special_dir_name" ]
    then
        export DATA_DIRECTORY=${RAW_DATA_DIR}/plate_vel_${PLATE_VEL}
    else
        export DATA_DIRECTORY=${RAW_DATA_DIR}/plate_vel_${PLATE_VEL}_${special_dir_name}
    fi

    rm -r $DATA_DIRECTORY
    mv uniform_plate $DATA_DIRECTORY
    cp parameters.h $DATA_DIRECTORY/
done

