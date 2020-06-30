#!/bin/bash

# output_clean.sh

# This script cleans data outputted from droplet impact simulations. 
# Specifically, simulations where quantities are outputted at each timestep in 
# the log file (such as pressure, force and position of the plate).
# Also where there are "plate_output_n.txt" files, which have 
# their first line listing the value of time when they were made, and the rest
# of the lines in column format giving the values of quantities such as radial
# coordinate (y), vertical coordinate (x), pressure (p) and other quantities in
# the format y = ..., x = ..., p = ..., etc. It converts
# these plate_output files into CSV style text files in order to make 
# quantitative analysis easier.
# NOTE: For the plate_output files, it assumes the output is in scientific 
# notation (e.g. 1.083e-7), so the "e" character is not removed. Hence care is
# needed if there are "e" characters that are not a part of the numerical output

# Parent directory where the raw data is
PARENT_DIR=$1

# Directory where raw data is stored
RAW_DATA_DIR=${PARENT_DIR}/raw_data

# Directory where the cleaned data is stored
CLEANED_DATA_DIR=${PARENT_DIR}/cleaned_data

# Creates the directories to store cleaned data
mkdir ${CLEANED_DATA_DIR}
mkdir ${CLEANED_DATA_DIR}/plate_outputs

################################################################################
# Cleans the log file 
################################################################################
# We are interested in the lines which start with "t = ", and can remove all 
# all other lines
sed -n '/^t = /p' ${RAW_DATA_DIR}/log > ${CLEANED_DATA_DIR}/output.txt

# Removes readable qualifies
sed -e "s/ds_dt = //g" -i ${CLEANED_DATA_DIR}/output.txt
sed -e "s/d2s_dt2 = //g" -i ${CLEANED_DATA_DIR}/output.txt
sed -e "s/t = //g" -i ${CLEANED_DATA_DIR}/output.txt
sed -e "s/v = //g" -i ${CLEANED_DATA_DIR}/output.txt
sed -e "s/F = //g" -i ${CLEANED_DATA_DIR}/output.txt
sed -e "s/F_avg = //g" -i ${CLEANED_DATA_DIR}/output.txt
sed -e "s/f_avg = //g" -i ${CLEANED_DATA_DIR}/output.txt
sed -e "s/s = //g" -i ${CLEANED_DATA_DIR}/output.txt
sed -e "s/bubble_no = //g" -i ${CLEANED_DATA_DIR}/output.txt
sed -e "s/drop_no = //g" -i ${CLEANED_DATA_DIR}/output.txt
sed -e "s/force_term = //g" -i ${CLEANED_DATA_DIR}/output.txt
sed -e "s/avg = //g" -i ${CLEANED_DATA_DIR}/output.txt
sed -e "s/std = //g" -i ${CLEANED_DATA_DIR}/output.txt

echo Cleaned log file

# Moves the log file into the main directory, for more easy access
cp ${RAW_DATA_DIR}/log ${PARENT_DIR}

################################################################################
# Moves the videos into a separate directory
################################################################################

# Makes a directory for the movies
MOVIE_DIR=${PARENT_DIR}/movies
mkdir ${MOVIE_DIR}

# Moves all mp4 files from the raw_data directory into the movies directory
mv ${RAW_DATA_DIR}/*.mp4 ${MOVIE_DIR}

################################################################################
# Moves the gfs files into a separate directory
################################################################################

# Makes a directory for the gfs files
GFS_DIR=${PARENT_DIR}/gfs_files
mkdir ${GFS_DIR}

# Moves all the gfs files from the raw_data direcotry into the movies directory
mv ${RAW_DATA_DIR}/*.gfs ${GFS_DIR}

################################################################################
# Cleans the plate output files
################################################################################
# The raw data directory will contain a certain number of files called 
# plate_output_n.txt, where n goes from 1 up to some number. The first line 
# of plate_output_n.txt will be of the form "t = ...", which is the time the 
# output was made. The rest of the lines are "y = ..., h = ..., etc".

# Counts the number of these files in the raw data directory
NOFILES=$(ls ${RAW_DATA_DIR}/ | grep plate_output_ | wc -l)

# Creates a file to store the times recorded by each file
TIMESFILE=${CLEANED_DATA_DIR}/plate_outputs/times.txt
> $TIMESFILE

# Loops over all of the plate output files
for (( filenum=1; filenum<=$NOFILES; filenum++))
do
    # Name of raw data file
    DATAFILE=${RAW_DATA_DIR}/plate_output_$filenum.txt
    
    # Extract the time from the first line
    TIME=$(head -n 1 $DATAFILE | sed 's/[^0-9\.-]//g')

    # Appends the time to the time file
    echo $filenum, $TIME >> $TIMESFILE
    
    # Creates the file to output cleaned data to
    OUTPUTFILE=${CLEANED_DATA_DIR}/plate_outputs/output_$filenum.txt

    # Clean the data file by deleting first line and non-essential characters. 
    # Note we need to include the letter e, for the scientific notation numbers
    # Resulting columns: y, h, p, u_x, u_y
    sed 1d $DATAFILE > $OUTPUTFILE
    sed -i 's/[^0-9,\.e-]//g' $OUTPUTFILE

    echo Cleaned plate_output_$filenum.txt 
done

