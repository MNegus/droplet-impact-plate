/* parameters.h
Header file for the parameters to feed into the simulations for droplet impact
*/

/* Fluid properties */
const double RHO_R = 0.00120; // Density ratio
const double MU_R = 0.0183; // Viscosity ratio
const double REYNOLDS = 4990; // Reynolds number
const double WEBER = 342; // Weber number
const double FR = 50.5 // Froude number

/* Droplet definition */
const double DROP_VEL = -1.0; // Initial velocity of the droplet 
const double DROP_RADIUS = 1.0; // Radius of droplet
const double INITIAL_DROP_HEIGHT = 0.125; // Initial gap between drop and plate

/* Plate definition */
const double INITIAL_PLATE_TOP = 1.0; // Initial top of plate 
const double PLATE_WIDTH = 3.0; // Width of plate (horizontal direction)
const double PLATE_THICKNESS = 0.15; // Thickness (vertical direction)

/* Prescribed plate parameters (only used when plate velocity is prescribed) */
const double PLATE_VEL = 0; // Velocity of the plate

/* Plate ODE terms (only used when plate is coupled). Corresponds to 
ALPHA s''(t) + BETA s'(t) + GAMMA s(t) = F(t) */
const double ALPHA = 0.01; // Mass term
const double BETA = 0; // Damping term
const double GAMMA = 1; // Elastic term

/* Computational constants */
const int MOVIES = 1; // Boolean for producing movies
const int MINLEVEL = 4; // Minimum refinement level 
const int MAXLEVEL = 10; // Maximum refinement level
const double BOX_WIDTH = 4.0; // Width of the computational box
const double START_OUTPUT_TIME = 0.0; // Time to start outputs
const double END_OUTPUT_TIME = 2.0; // Time to end outputs
const double GFS_OUTPUT_TIMESTEP = 1e-2; // Time between gfs outputs
const double PLATE_OUTPUT_TIMESTEP = 1e-3; // Time between plate outputs
const double INTERFACE_OUTPUT_TIMESTEP = 1e-3; // Time between interface outputs
const double HARD_MAX_TIME = 2.0; // Hard maximum time 
