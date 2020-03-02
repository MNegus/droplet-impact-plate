#!/bin/bash





RAWDATADIRECTORY=~/repos/plate-impact/uniform_plate/code/
coderun=uniform_plate
CLEANEDDIRECTORY=/scratch/Uniform_plate/cleaned_data/$coderun

source ~/repos/plate-impact/post_processing/data_cleaning/simulation_output_clean.sh $RAWDATADIRECTORY $coderun $CLEANEDDIRECTORY

