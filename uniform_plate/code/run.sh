#!/bin/bash

RAW_DATA_DIR=/scratch/negus/uniform_plate/raw_data

export OMP_NUM_THREADS=8

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

    export DATA_DIRECTORY=${RAW_DATA_DIR}/plate_vel_${PLATE_VEL}

    rm -r $DATA_DIRECTORY
    mv uniform_plate $DATA_DIRECTORY
done

