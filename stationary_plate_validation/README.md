# stationary_plate_validation

Scripts for validating the accuracy of the method for creating a solid plate 
using a volume fraction comparing the stationary plate case with the solid wall
code. 

## Features
* Implements impact onto a plate, where the plate is defined to be a lower 
section of the domain (i.e. x < PLATE_HEIGHT)
* *run.sh* conducts 4 different simulations, with the plate being aligned to 
cells in different ways:
    * *PLATE_ALIGNMENT = 2*: The top of the plate cuts a cell in half (1/2)
    * *PLATE_ALIGNMENT = 4*: The top of the plate cuts a cell in quarter (1/4)
    * *PLATE_ALIGNMENT = 0*: The top of the plate is aligned with the edge of a 
    cell
    * *PLATE_ALIGNMENT = -1*: There is no plate and instead the droplet impacts 
    the solid wall of the domain
* With these four simulations, the initial height the droplet is above the 
plate/wall is kept the same, so that in theory they impact at the same time. All
other physical/computational quantities are also kept the same
* The aim is to compare the outputs from the four simulations to compare how 
they differ. We ideally want the results to be the same.

## Outputs
* **Plate output**: At regular time intervals, the pressure and velocities along 
the plate are outputted to files called *plate_output_n.txt*, where *n* is an 
integer. 
    * The *plate_output_n.txt* files are of the following format
    ``` 
    t = 0.01
    y = 1.2, h = 0, p = 1.05e-4, u_x = 1.078e-7, u_y = 1.034e0 
    y = 1.2, h = 1, p = 8.67e-5, u_x = 1.223e-9 , u_y = 6.053e-1 
    y = 1.2, h = 2, ...
    ...
    ```
    * The first line of the *plate_output_n.txt* files give the computational time that the file was created
    * The rest of the lines are in column format, giving the value of *y* (which is the radial coordinate), *h* (the number of cells above the plate), the pressure *p*, vertical velocity *u_x* and the radial velocity *u_y*
    * These files are in the format to be cleaned using *simulation_output_clean.sh* in the *post_processing* directory
* **Volume**: The volume of the droplet phase in time is outputted into the log file in the format of 
    ```
    t = 0, volume = 4.18879
    t = 0.001, volume = 4.18879
    ...
    ```
* **GFS files**: GFS files are outputted regularly in time in the form *gfs_output_n.txt* where *n* is an integer
* **Video**: Videos of the volume fractions in time are created in the *images* event.

## Execution
The *run.sh* scripts runs the simulation 4 times, with the different options for plate/wall alignment. It will create 4 separate directories which stores this raw data in

## Data cleaning
The *data_cleaning* directory contains a script *data_clean.sh*, which cleans the data output from the 4 data directories. This script requires access to the *data_cleaning* subdirectory of the *post_processing* directory. 

## Data analysis
The *data_analysis* directory contains *MATLAB* scripts for analysing the cleaned data outputted from the simulations.
* **Interfacial overlay**: The *interface_overlay.m* overlays the interfacial locations of the four simulations (with the plate height normalised to be at *z = 0* in all cases). The function has two options for the input variable *plot_type*:
    * *plot_type = macro*: This creates a video of the "macroscopic" view of the impact, where the entire droplet is visible. It overlays the four lines outputted from the simulations in time and plots them on the same graph. 
    * *plot_type = jet*: This creates a video which follows along the formation of the jet. The "camera view" moves to the right with the speed of the turnover point, and zooms out as the jet grows so all of it is visible. This video is used to compare the differences in the four simulations at this small scale. Additionally, the theoretical predictions of the interfacial location using Wagner theory is plotted on the same graph. 
* **Volume conservation**: In the *data_analysis.m* script, the ratio of the droplet volume to the initial volume in time is plotted for each simulation run.
* **Pressure along plate**: In the *data_analysis.m* script, the pressure measured along the plate (for *h = 0, 1, 2*) for each simulation is plotted as a function of time, and compared against the theoretical prediction of Wagner theory. 

## Figures and Videos
The *Figures* directory stores any figures outputted from the data analysis, and similarly with the *Videos* directory, which also contains the video outputs directly from the simulation.