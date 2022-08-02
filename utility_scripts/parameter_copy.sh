#!/bin/bash

code_dir=$1
dest_dir=$2


# Change the axisymmetric value in parameters file
sed -i "/#define AXISYMMETRIC/c\#define AXISYMMETRIC 1" parameters.h


####################################
# ALPHA varying case
####################################

# Make the ALPHA varying directory
mkdir $dest_dir/ALPHA_varying

# Set BETA and GAMMA
sed -i "/const double BETA/c\const double BETA = 0.; // Damping term" parameters.h
sed -i "/const double GAMMA/c\const double GAMMA = 0.; // Elastic term" parameters.h

for ALPHA in 100.0 20.0 10.0 5.0 2.0 1.0
do
    # Edit ALPHA
    sed -i "/const double ALPHA/c\const double ALPHA = $ALPHA; // Mass term" parameters.h

    # Copy over the code
    ./code_copy.sh $code_dir $dest_dir/ALPHA_varying ALPHA_$ALPHA
done


####################################
# BETA varying case
####################################

# Make the BETA varying directory
mkdir $dest_dir/BETA_varying

# Set ALPHA and GAMMA
sed -i "/const double ALPHA/c\const double ALPHA = 2.0; // Mass term" parameters.h
sed -i "/const double GAMMA/c\const double GAMMA = 100.0; // Elastic term" parameters.h

for BETA in 0.0 7.07 28.28 141.42
do
    # Edit BETA
    sed -i "/const double BETA/c\const double BETA = $BETA; // Damping term" parameters.h

    # Copy over the code
    ./code_copy.sh $code_dir $dest_dir/BETA_varying BETA_$BETA
done


####################################
# GAMMA varying case
####################################

# Make the GAMMA varying directory
mkdir $dest_dir/GAMMA_varying

# Set ALPHA and BETA
sed -i "/const double ALPHA/c\const double ALPHA = 2.0; // Mass term" parameters.h
sed -i "/const double BETA/c\const double BETA = 0.0; // Damping term" parameters.h

for GAMMA in 0.0 10.0 100.0 500.0 1000.0
do
    # Edit GAMMA
    sed -i "/const double GAMMA/c\const double GAMMA = $GAMMA; // Elastic term" parameters.h

    # Copy over the code
    ./code_copy.sh $code_dir $dest_dir/GAMMA_varying GAMMA_$GAMMA
done