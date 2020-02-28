/* parameters.h
Header file for the parameters to feed into the simulations for droplet impact*/

/* Physical constants */
const double RHO_R = 0.003; // Density ratio
const double MU_R = 0.002; // Viscosity ratio
const double REYNOLDS = 1000.0; // Reynolds number
const double WEBER = 10000.0; // Weber number
const double FR = 10.1; // Froude number
const double DROP_VEL = -1.0; // Initial velocity of the droplet 
const double DROP_RADIUS = 1.0; // Radius of droplet
const double INTITIAL_DROP_HEIGHT = 0.56; // Distance of drop from surface

/* Computational constants */
const int MINLEVEL = 4; // Minimum refinement level 
const int MAXLEVEL = 10; // Maximum refinement level
const double BOX_WIDTH = 7.5; // Width of the computational box

