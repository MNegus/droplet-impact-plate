#!/bin/bash


RAWDATADIRECTORY=/scratch/Uniform_plate/raw_data

for coderun in plate_vel_0.1 plate_vel_0.2 plate_vel_0.3 plate_vel_0.4 plate_vel_0.5
do
    CLEANEDDIRECTORY=/scratch/Uniform_plate/cleaned_data/$coderun
    source ~/repos/plate-impact/post_processing/data_cleaning/simulation_output_clean.sh $RAWDATADIRECTORY $coderun $CLEANEDDIRECTORY
done

