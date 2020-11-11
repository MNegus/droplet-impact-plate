#!/bin/bash

# Script to combine two videos outputted from DNS of the pressure
# Assumes in the current directory there are the videos:
#
# stationary.mp4 (Video from stationary plate case)
# moving.mp4 (Video from moving plate case)
#
# Subsequently produces dns_pressure.mp4, a video with the
# stationary case on the left and moving on the right

# Flips the stationary pressure
ffmpeg -i stationary.mp4 -vf hflip -c:a copy stationary_mirrored.mp4

# Combines mirrored stationary and moving 
ffmpeg -i stationary_mirrored.mp4 -i moving.mp4 -filter_complex hstack \
    dns_pressure.mp4

# Removes the mirrored video
rm stationary_mirrored.mp4
