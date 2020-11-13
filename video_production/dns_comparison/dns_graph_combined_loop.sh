#!/bin/bash

MAINDIR=/home/michael/Documents/supplementary_material/dns_videos
OVERLAY=$MAINDIR/pressure_overlays/overlay_height_257_max_3.avi
HEIGHT=250


for COEFF in 1.0 1.1 1.2 1.3 1.4
do
    DNS=$MAINDIR/combined_dns_videos/combined_pressure_$COEFF.mp4
    RES=$MAINDIR/dns_with_graphs/pressure_${COEFF}_overlay_height_257_max_3.mp4
    ./dns_graph_combined.sh $OVERLAY $DNS $RES $HEIGHT
done