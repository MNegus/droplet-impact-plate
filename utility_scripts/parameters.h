/* parameters.h
Header file for the parameters to feed into the simulations for droplet impact. 

Most should be left as they are, but some can be edited to change either the 
computational or physical options
*/

/* Dimensional fluid properties */
const double R = 1.0e-3; // Radius of droplet (metres)
const double V = 5.0; // Velocity of droplet (metres per second)
const double RHO_L = 998.0; // Density of liquid (kilograms per metre cubed)
const double RHO_G = 1.23; // Density of gas (kilograms per metre cubed)
const double MU_L = 1.0e-3; // Viscosity of liquid (Pascal seconds)
const double MU_G = 1.81e-5; // Viscosity of gas (Pascal seconds)
const double SIGMA = 7.29e-2; // Surface tension (Newtons per metre)

/* Dimensionless droplet definition */
const double DROP_VEL = -1.0; // Initial velocity of the droplet 
const double DROP_RADIUS = 1.0; // Radius of droplet
const double INITIAL_DROP_HEIGHT = 0.125; // Initial gap between drop and plate

/* Plate definition */
const double PLATE_WIDTH = 2.0; // Width of plate (horizontal direction)

/* Plate ODE terms (only used when plate is coupled). Corresponds to 
ALPHA s''(t) + BETA s'(t) + GAMMA s(t) = F(t) */
const double ALPHA = 2.; // Mass term
const double BETA = 0.; // Damping term
const double GAMMA = 500.; // Elastic term

/* Override: Constant acceleration. Use this to set the acceleration of the 
plate to be constant (decoupling the plate motion from the droplet) */
const int CONST_ACC = 0; // Set to 1 to specify a constant acceleration
const double PLATE_ACC = 0.; // Constnat acceleration of the plate

/* Computational options. I.e. ones relating to the computational setup 
and not related to the actual physics */
// General
const double HARD_MAX_TIME = 0.8; // Hard maximum time (end time may be shorter)
const double BOX_WIDTH = 3.0; // Width of the computational box
const double FORCE_DELAY_TIME = 0.01; // Delay time before force is applied on plate
// Refinement options
const int MINLEVEL = 4; // Minimum refinement level 
const int MAXLEVEL = 13; // Maximum refinement level
const int PLATE_REFINE_NO = 4; // Number of max refinement cells above plate
const double DROP_REFINED_WIDTH = 0.04; // width of refined region around droplet
// Output options
const int MOVIES = 1; // Set 1 to produce movies
const double START_OUTPUT_TIME = 0.0; // Time to start outputs
const double END_OUTPUT_TIME = 2.0; // Time to end outputs
const double GFS_OUTPUT_TIMESTEP = 1e-2; // Time between gfs outputs
const double PLATE_OUTPUT_TIMESTEP = 1e-3; // Time between plate outputs
const double LOG_OUTPUT_TIMESTEP = 1e-4; // Time between log outputs
const double INTERFACE_OUTPUT_TIMESTEP = 1e-3; // Time between interface outputs
// Removal options
const double REMOVAL_DELAY = 0.02; // Time after pinch-off to start removal
const int REMOVE_ENTRAPMENT = 0; // If 1, completely remove entrapped air
// Peak detect options
const int PEAK_DETECT = 1; // If 1, remove peaks in the force
const int PEAK_LAG = 4; // Lag used in peak detection
const double PEAK_THRESHOLD = 4.0; // Number of std devs away from mean
const double PEAK_INFLUENCE = 0.1; // Influence weighting from peak data
const double PEAK_DELAY = 0.135; // Delay before peak detection starts


