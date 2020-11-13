#!/bin/bash

MAINDIR=/home/michael/Documents/supplementary_material/dns_videos
STATDIR=$MAINDIR/stationary_videos
MOVDIR=$MAINDIR/moving_videos
RESDIR=$MAINDIR/combined_dns_videos

# Pressure videos
for COEFF in 1.0 1.1 1.2 1.3 1.4
do
    ./dns_snapshot.sh $STATDIR/pressure_$COEFF.mp4 $MOVDIR/pressure_$COEFF.mp4 \
        $RESDIR/combined_pressure_$COEFF.mp4
done

# Velocity videos
for COEFF in 1 2 3
do
    ./dns_snapshot.sh $STATDIR/horizontal_vel_$COEFF.mp4 \
        $MOVDIR/horizontal_vel_$COEFF.mp4 \
        $RESDIR/combined_horizontal_vel_$COEFF.mp4
    
    ./dns_snapshot.sh $STATDIR/vertical_vel_$COEFF.mp4 \
        $MOVDIR/vertical_vel_$COEFF.mp4 \
        $RESDIR/combined_vertical_vel_$COEFF.mp4
    
done