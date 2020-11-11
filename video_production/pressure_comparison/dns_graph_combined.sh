#!/bin/bash

# Script to combine the video of the DNS pressure alongside a moving graph
# outputted from MATLAB. The MATLAB graph will be in the video
# pressure_overlay.avi, which needs to be padded and resolution changed in order
# to fit underneath the DNS video combined_pressure.mp4. The final video will be
# called pressure_with_graph.mp4

# # Adjusts pressure resolution
# ffmpeg -y -i pressure_overlay.avi -vf scale=2048:-1 scaled_pressure.mp4

# # Crops the scaled pressure video
# ffmpeg -y -i scaled_pressure.mp4 -filter:v "crop=1675:512:190:512" cropped_pressure.mp4

# # Rescale the cropped pressure video
# ffmpeg -y -i cropped_pressure.mp4 -vf scale=2048:-1 graph_pressure.mp4

# Pads the DNS video
ffmpeg -y -i dns_pressure.mp4 -vf "pad=width=2162:height=1024:x=100:y=0:color=white" padded_dns.mp4

# Rescales DNS video
ffmpeg -y -i padded_dns.mp4 -vf scale=2048:-1 scaled_dns.mp4

# Combine videos
ffmpeg -y -i scaled_dns.mp4 -i graph_pressure.mp4 -filter_complex vstack \
    dns_with_graph_pressure.mp4


