# wall_code
Code for droplet impact onto a solid wall which is at the bottom of the computational domain. 

## Features
* Axisymmetric implementation of a droplet falling downwards towards the "left" computational boundary
* No-slip and impermeability conditions at the left boundary result in it acting as a solid wall
* Can specify physical properties of the droplet (e.g. density ratios and initial velocity) and computational properties (e.g. maximum refinement levels)

## Outputs:
* Pressure along the cell directly above the wall
* Videos of the droplet volume fraction and vorticity in time
* GFS files to be viewed using gfsview2D
* Location of droplet interface using output_facets 
* Area of the droplet phase (for checking mass convergence)

## Running
1. Specify desired physical and computational properties in the parameters.h 
2. Specify the number of OpenMP threads to use in run.sh
3. Execute run.sh to being the simulation