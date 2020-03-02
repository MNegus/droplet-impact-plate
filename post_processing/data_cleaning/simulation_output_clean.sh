#!/bin/bash

################################################################################
# This script cleans data outputted from droplet impact simulations. 
# Specifically, simulations where the volume of the droplet is outputted in time
# in the log file. Also where there are "plate_output_n.txt" files, which have 
# their first line listing the value of time when they were made, and the rest
# of the lines in column format giving the values of y, h, p, u_x and u_y along 
# the plate (e.g. y = ..., h = ..., p = ..., u_x = ...., u_y = ...). It converts
# these plate_output files into CSV style text files in order to make 
# quantitative analysis easier.
################################################################################

# Name of the directory where the raw data is stored
RAWDATADIRECTORY=$1
coderun=$2
CLEANEDDIRECTORY=$3

DATADIRECTORY=$RAWDATADIRECTORY/$coderun

echo Cleaning $coderun

# Creates the directories to store cleaned data
# CLEANEDDIRECTORY=/scratch/Plate_validation_data/cleaned_data/$coderun
mkdir $CLEANEDDIRECTORY
mkdir $CLEANEDDIRECTORY/plate_outputs

################################################################################
# Cleans the log file for the volume outputs
################################################################################
# We are interested in the lines which start with "t = ", and can remove all 
# all other lines
sed -n '/^t = /p' $DATADIRECTORY/log > $CLEANEDDIRECTORY/volumes.txt

# Removes all characters apart from numbers, commas and full stops. Therefore 
# first column is t, the second is volume
sed -i 's/[^0-9,\.]//g' $CLEANEDDIRECTORY/volumes.txt

echo Cleaned volumes file


################################################################################
# Cleans the plate output files
################################################################################
# DATADIRECTORY will contain a certain number of files called 
# plate_output_n.txt, where n goes from 1 up to some number. The first line 
# of plate_output_n.txt will be of the form "t = ...", which is the time the 
# output was made. The rest of the lines are "y = ..., h = ..., etc".

# Counts the number of these files in DATADIRECTORY
NOFILES=$(ls $DATADIRECTORY/ | grep plate_output_ | wc -l)

# Creates a file to store the times recorded by each file
TIMESFILE=$CLEANEDDIRECTORY/plate_outputs/times.txt
> $TIMESFILE

# Loops over all of the plate output files
for (( filenum=1; filenum<=$NOFILES; filenum++))
do
    # Name of raw data file
    DATAFILE=$DATADIRECTORY/plate_output_$filenum.txt
    
    # Extract the time from the first line
    TIME=$(head -n 1 $DATAFILE | sed 's/[^0-9\.-]//g')

    # Appends the time to the time file
    echo $filenum, $TIME >> $TIMESFILE
    
    # Creates the file to output cleaned data to
    OUTPUTFILE=$CLEANEDDIRECTORY/plate_outputs/output_$filenum.txt

    # Clean the data file by deleting first line and non-essential characters. 
    # Note we need to include the letter e, for the scientific notation numbers
    # Resulting columns: y, h, p, u_x, u_y
    sed 1d $DATAFILE > $OUTPUTFILE
    sed -i 's/[^0-9,\.e-]//g' $OUTPUTFILE

    echo Cleaned plate_output_$filenum.txt 
done

