#!/bin/bash

################################################################################
# Cleans all of the data for the Plate_validation simulations by iteratively 
# running the simulation_output_clean.sh script
################################################################################

RAWDATADIRECTORY=/scratch/Plate_validation_data/raw_data

for coderun in plate_alignment_0 plate_alignment_2 plate_alignment_4 wall_data
do
    ./simulation_output_clean.sh $RAWDATADIRECTORY $coderun
done