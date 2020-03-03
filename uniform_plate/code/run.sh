#!/bin/bash

export OMP_NUM_THREADS=8

# for PLATE_VEL in 0.2 0.3 0.4 0.5
for PLATE_VEL in 0.2
do
    # Deletes the previous directory and test files
    rm -r uniform_plate
    rm *.s
    rm *.s.d
    rm *.tests 
    rm *.deps
    rm *.tst

    make uniform_plate.tst PLATE_VEL=-$PLATE_VEL

    export DATA_DIRECTORY=plate_vel_${PLATE_VEL}

    rm -r $DATA_DIRECTORY
    mv uniform_plate $DATA_DIRECTORY
done

