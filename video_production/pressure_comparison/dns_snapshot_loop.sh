#!/bin/bash

MAINDIR=/media/michael/newarre/presentation_data
STATDIR=$MAINDIR/stationary_videos
MOVDIR=$MAINDIR/moving_videos
RESDIR=$MAINDIR/combined_dns_videos

for COEFF in 1.0 1.1 1.2 1.3 1.4
do
    ./dns_snapshot.sh $STATDIR/pressure_$COEFF.mp4 $MOVDIR/pressure_$COEFF.mp4 \
        $RESDIR/combined_pressure_$COEFF.mp4
done