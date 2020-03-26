#!/bin/bash


RAWDATADIRECTORY=/scratch/uniform_plate/raw_data


# for PLATEVEL in 0 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5
for PLATEVEL in 0.1
do
    coderun=plate_vel_${PLATEVEL}_${special_dir_name}
    CLEANEDDIRECTORY=/scratch/uniform_plate/cleaned_data/$coderun
    source ~/repos/plate-impact/post_processing/data_cleaning/simulation_output_clean.sh $RAWDATADIRECTORY $coderun $CLEANEDDIRECTORY
done

