# data_analysis

MATLAB scripts for analysing the data output from the simulations

## Interfacial analysis 
The *interface_analysis* directory contains scripts to analyse the interfacial position data outputted from the simulations using output_facets. These are all designed to be used in conjunction with eachother to create plots of the interfaces. 
* **read_interface_points.m**: Reads the interface files and converts them into a pair of MATLAB arrays for the start and end points of each line segment
* **line_segment_plot.m**: Takes the arrays for the line segments outputted from *read_interface_points.m* and plots the interface as a series of line segments
* **coarsen_interface.m**: Takes in the arrays for the line segments and coarsens them by halving the number of line segments. Used for highly refined simulations which result in the interface taking a long time to plot using *line_segment_plot.m*
