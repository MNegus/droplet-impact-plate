#!/bin/bash

export OMP_NUM_THREADS=8

for PLATE_ALIGNMENT in 0 4 2 -1
do
    # Deletes the previous directory if it exists
    rm -r plate_impact
    rm *.s
    rm *.s.d
    rm *.tests 
    rm *.deps
    rm *.tst


    make plate_impact.tst PLATE_ALIGNMENT=$PLATE_ALIGNMENT
    
    if [ $PLATE_ALIGNMENT == "-1" ]; then
        export DATA_DIRECTORY=wall_data
    else
        export DATA_DIRECTORY=plate_alignment_${PLATE_ALIGNMENT}
    fi

    rm -r $DATA_DIRECTORY
    mv plate_impact $DATA_DIRECTORY

done


