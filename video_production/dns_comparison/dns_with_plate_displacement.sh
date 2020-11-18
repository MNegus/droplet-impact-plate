#!/bin/bash

# Script to combine the video of the DNS pressure alongside a moving graph
# outputted from MATLAB. The MATLAB graph will be in the video
# pressure_overlay.avi, which needs to be padded and resolution changed in order
# to fit underneath the DNS video combined_pressure.mp4. The final video will be
# called pressure_with_graph.mp4

PLATEORIG=/home/michael/Documents/supplementary_material/plate_displacement_animation/original_animation.avi
DNS=/home/michael/Documents/supplementary_material/pressure_with_graph.mp4
RESULT=/home/michael/Documents/supplementary_material/dns_videos/combined_with_displacement/pressure.mp4
HEIGHT=250

# Crops the displacement video
ffmpeg -y -i $PLATEORIG -filter:v "crop=1728:$HEIGHT:146:$HEIGHT" cropped_plate.mp4

# Rescale the cropped pressure video
ffmpeg -y -i cropped_plate.mp4 -vf scale=2048:-1 plate_anim.mp4

# # Combine videos
ffmpeg -y -i $DNS -i plate_anim.mp4 -filter_complex vstack \
   $RESULT

# Remove videos
rm cropped_plate.mp4