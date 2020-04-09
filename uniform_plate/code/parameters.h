/* parameters.h
Header file for the parameters to feed into the simulations for droplet impact*/

/* Physical constants */
const double RHO_R = 0.003; // Density ratio
const double MU_R = 0.002; // Viscosity ratio
const double REYNOLDS = 1000.0; // Reynolds number
const double WEBER = 10000.0; // Weber number
const double FR = 10.1; // Froude number
const double PLATE_VEL = -0.1; // Velocity of the plate
const double DROP_VEL = -1.0; // Initial velocity of the droplet 
const double DROP_RADIUS = 1.0; // Radius of droplet
const double DROP_CENTRE = 2.125; // Initial centre position of droplet
const double INITIAL_PLATE_TOP = 1.0; // Initial top of plate 
const double PLATE_WIDTH = 3.0; // Width of plate (horizontal direction)
const double PLATE_THICKNESS = 0.15; // Thickness (vertical direction)

/* Computational constants */
const int MINLEVEL = 4; // Minimum refinement level 
const int MAXLEVEL = 12; // Maximum refinement level
const double BOX_WIDTH = 4.0; // Width of the computational box
const double START_OUTPUT_TIME = 0.0; // Time to start outputs
const double END_OUTPUT_TIME = 2.0; // Time to end outputs
const double GFS_OUTPUT_TIMESTEP = 1e-2; // Time between gfs outputs
const double PLATE_OUTPUT_TIMESTEP = 1e-3; // Time between plate outputs
const double INTERFACE_OUTPUT_TIMESTEP = 1e-3; // Time between interface outputs
const double HARD_MAX_TIME = 2.0; // Hard maximum time 
const int FILTERED = 0; // Set to zero to prevent filtering
const double POISSON_TOLERANCE = 1e-3; // Tolerance of Poisson solver (default 1e-3)
