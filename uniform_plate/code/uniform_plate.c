/* uniform_plate.c
    An axisymmetric droplet falling towards an impermeable plate which can move
    in the vertical direction at a prescribed velocity, PLATE_VEL, which is 
    assumed to be negative

    We run the simulation until the turnover point reaches the droplet radius. 
    According to Wagner theory, for DROP_VEL = 1 and DROP_RADIUS = 1, the
    turnover point is at
        r = d(t) = sqrt(3 * (1 + PLATE_VEL) * (t - IMPACT_TIME))
    where IMPACT_TIME is the time at which the droplet hits the plate.
    So for d = DROP_RADIUS = 1 when t = IMPACT_TIME + 1 / (3 * (1 + PLATE_VEL))
*/

#include "parameters.h" // Include all defined parameters
#include "axi.h" // Axisymmetric coordinates
#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h" // Implements two-phase flow
#include "view.h" // Creating movies using bview
#include "tension.h" // Adds forces due to surface tension
#include "tag.h" // For removing small droplets
#include <omp.h> // For openMP parallel

/* Computational constants derived from parameters */
double MIN_CELL_SIZE; // Size of the smallest cell
double DROP_REFINED_WIDTH; // Width of the refined area around the droplet
double PLATE_REFINED_WIDTH; // Width of the refined area around the plate
double IMPACT_TIME; // Theoretical time of impact
double MAX_TIME; // Maximum time to run the simulation for
double INTERPOLATE_DISTANCE; // Distance above plate to read the pressure

/* Boundary conditions */
// Outflow through boundaries far from impact
u.n[right] = neumann(0.); 
u.n[top] = neumann(0.);
u.n[left] = neumann(0.);

// Zero pressure far from impact
p[right] = dirichlet(0.); 
p[top] = dirichlet(0.); 

/* Fields */
scalar plate[]; // VOF field for the plate

/* Global variables */
double plate_position; // x-position of the top of the plate at time t
int gfs_output_no = 1; // Records how many GFS files have been outputted
int plate_output_no = 1; // Records how many plate data files there have been
int interface_output_no = 1; // Records how many interface files there have been
double start_wall_time;
double end_wall_time;

/* Function declarations */
double plate_region(double xp, double yp); // Defines VOF field for plate
double distance_from_plate (double xp, double yp); // Gives distance from plate 
double force_on_plate (); // Calculates the force on the plate via integration

int main() {
/* Main function for running the simulation */

    /* If filtering is set in the parameters file, then define the FILTERED 
    macro */
    if (FILTERED) {
        #define FILTERED 1
    }

    /* Changes the Poisson tolerence to what was defined in parameters.h */
    TOLERANCE = POISSON_TOLERANCE;
    NITERMAX = POISSON_NITERMAX;

    /* Create the computational domain */
    init_grid(1 << MINLEVEL); // Create grid according to the minimum level
    size(BOX_WIDTH); // Size of the domain

    /* Set physical constants */
    rho1 = 1.; // Density of water phase
    rho2 = RHO_R; // Density of air phase
    mu1 = 1. / REYNOLDS; // Viscosity of water phase
    mu2 = mu1 * MU_R; // Viscosity of air phase
    f.sigma = 1. / WEBER; // Surface tension at interface

    /* Derived constants */
    MIN_CELL_SIZE = BOX_WIDTH / pow(2, MAXLEVEL); // Size of the smallest cell
    PLATE_REFINED_WIDTH = 0.3 * PLATE_THICKNESS; // Refined region around plate
    DROP_REFINED_WIDTH = 0.05; // Refined region around droplet
    IMPACT_TIME = (DROP_CENTRE - DROP_RADIUS - INITIAL_PLATE_TOP) \
        / (PLATE_VEL - DROP_VEL);
    
    INTERPOLATE_DISTANCE = MIN_CELL_SIZE; // Distance above plate to read pressure

    // Maximum time is shortly after the Wagner theory prediction of the 
    // turnover point reaching the radius of the droplet
    double wagner_max_time = 1.5 * (IMPACT_TIME + 1. / (3. * (1. + PLATE_VEL)));
    MAX_TIME = min(HARD_MAX_TIME, wagner_max_time); 

    /* Run the simulation */
    run();
}

event init (t = 0) {
/* Initialises the computational domain */

    // Records the wall time
    start_wall_time = omp_get_wtime();

    plate_position = INITIAL_PLATE_TOP; // Initialises top of plate position

    /* Refines around the droplet */
    refine(sq(x - DROP_CENTRE) + sq(y) < sq(DROP_RADIUS + DROP_REFINED_WIDTH) \
        && sq(x - DROP_CENTRE) + sq(y) > sq(DROP_RADIUS - DROP_REFINED_WIDTH) \
        && level < MAXLEVEL);
    
    /* Initialises the droplet volume fraction */
    fraction(f, -sq(x - DROP_CENTRE) - sq(y) + sq(DROP_RADIUS));

     /* Refines around the plate */
    refine(distance_from_plate(x, y) < PLATE_REFINED_WIDTH \
        && level < MAXLEVEL - 2);
    
    /* Initialises the plate volume fraction */
    fraction(plate, 1 - plate_region(x, y)); 

    /* Initialises the velocities */
    foreach () {
        u.x[] = DROP_VEL * f[]; // Droplet velocity
        u.x[] = plate[] * PLATE_VEL + (1 - plate[]) * u.x[]; // Plate velocity
    }
    boundary ((scalar *){u});
}


event moving_plate (i++) {
    /* Updates plate position */
    plate_position = INITIAL_PLATE_TOP + PLATE_VEL * t;

     /* Refines around the plate */
    refine(distance_from_plate(x, y) < PLATE_REFINED_WIDTH \
        && level < MAXLEVEL - 2);

    /* Definition of the plate volume fraction at time t */
    fraction(plate, 1 - plate_region(x, y));
    
    /* Alters the velocity of the fluid depending on if is touching the plate */
    foreach() {
        // Sets vertical velocity to match the plate at the interface
        u.x[] = plate[] * PLATE_VEL + (1 - plate[]) * u.x[];

        // Sets horizontal velocity inside plate to be zero
        u.y[] = (1. - plate[]) * u.y[]; 
    }
    boundary ((scalar *){u}); // Redefine boundary conditions for u

    /* EXPERIMENTAL: Sets the pressure */
    /* We try a fix of by setting the pressure in the mixed cell (where the 
    plate cuts through) to be equal to the pressure in the cell above. This is 
    kind of like a retroactive Neumann condition on the pressure */
    // foreach() {
    //     /* Finds the cells which the plate cuts through. If the plate is along a
    //     cell boundary, then nothing is done. More thoroughly:
    //         1) For a cell to be an interfacial cell, then it MUST be true that 
    //         the cell above has plate[] == 0 and the cell below has plate[] == 1
    //         2) If the plate cuts through a cell, then point 1) then for each y, 
    //         point 1) is only true for one cell
    //         3) However if the plate is lying on the boundary of a cell, then 
    //         there will be two cells where this is true for each y. In this case,
    //         plate[] == 0 or 1. We want to do nothing in these cases, so we
    //         ensure that plate[] > 0 and plate[] < 1.
    //         4) We also only test the top part of the plate, away from the curved
    //         edge, so ensure y < PLATE_WIDTH
    //     */
    //     if ((plate[1, 0] == 0) && (plate[-1, 0] == 1) && (plate[] > 0) \
    //         && (plate[] < 1) && (y < PLATE_WIDTH)) {
    //         /* Sets pressure to be equal to pressure in above cell */
    //         p[] = p[1, 0]; 
    //      }
    // }
}

event refinement (i++) {

    // Adapts with respect to velocities and volume fractions 
    adapt_wavelet ({u.x, u.y, f, plate}, (double[]){1e-2, 1e-2, 1e-2, 1e-2}, 
        minlevel = MINLEVEL, maxlevel = MAXLEVEL);

    // Refines region above plate
    refine(y < PLATE_WIDTH && x >= plate_position \
        && x < plate_position + 0.5 * PLATE_REFINED_WIDTH && level < MAXLEVEL);
}

event gravity (i++) {
/* Adds acceleration due to gravity in the vertical direction */
    face vector av = a; // Acceleration at each face
    foreach_face(x) av.x[] += - 1/sq(FR); // Adds acceleration due to gravity
}

event small_droplet_removal (i++) {
/* Removes any small droplets or bubbles that have formed, that are smaller than
    a specific size */
    // Removes droplets of diameter 5 cells or less
    remove_droplets(f, 5);

    // Removes bubbles of diameter 5 cells or less
    remove_droplets(f, 5, 1e-4, true); 
}

event output_volume (t += 0.001) {
/* Outputs the volume of the liquid phase to the log file */
    fprintf(stderr, "t = %g, v = %g, F = %g\n", t, 2 * pi * statsf(f).sum, force_on_plate());
}

event output_interface (t += INTERFACE_OUTPUT_TIMESTEP) {
/* Outputs the interface locations of the droplet */
    if ((t >= START_OUTPUT_TIME) && (t <= END_OUTPUT_TIME)) {
    // Creates text file to save output to
    char interface_filename[80];
    sprintf(interface_filename, "interface_%d.txt", interface_output_no);
    FILE *interface_file = fopen(interface_filename, "w");

    // Outputs the interface locations and closes the file
    output_facets(f, interface_file);
    fclose(interface_file);

    interface_output_no++;
    }
}

event output_values_along_plate (t += PLATE_OUTPUT_TIMESTEP) {
/* Outputs the values of pressure and velocities along the plate */

    if ((t >= START_OUTPUT_TIME) && (t <= END_OUTPUT_TIME)) {
    // Creates the file for outputting data
    char plate_output_filename[80];
    sprintf(plate_output_filename, "plate_output_%d.txt", plate_output_no);
    FILE *plate_output_file = fopen(plate_output_filename, "w");

    // Adds the time to the first line of the file
    fprintf(plate_output_file, "t = %g\n", t);

    // /* Iterates over the cells on the boundary of the plate */
    // foreach() {
    //     // Determines if current cell is along the boundary of the plate
    //     if ((plate[1, 0] == 0) && (plate[-1, 0] == 1) && (y < PLATE_WIDTH)) {
    //         /* Outputs the values along the plate and the two cells above it */
    //         for (int h = 0; h <= 2; h++) {
    //             fprintf(plate_output_file, \
    //                 "y = %g, h = %d, p = %g, u_x = %g, u_y = %g\n",
    //                 y, h, p[h, 0], u.x[h, 0], u.y[h, 0]);
    //         }
    //     }
    // }

    /* Uses interpolation to read pressure just above the plate */
    foreach() {
        // Determines if current cell is along the boundary of the plate
        if ((plate[1, 0] == 0) && (plate[-1, 0] == 1) && (y < PLATE_WIDTH)) {
            // x position to interpolate from
            double interpolate_read_pos = plate_position + INTERPOLATE_DISTANCE;

            // Interpolated pressure
            double pressure = interpolate(p, interpolate_read_pos, y);

            // Interpolated velocities
            double ux = interpolate(u.x, interpolate_read_pos, y);
            double uy = interpolate(u.y, interpolate_read_pos, y);

            fprintf(plate_output_file, \
                    "y = %g, x = %g, p = %g, u_x = %g, u_y = %g\n",
                    y, interpolate_read_pos, pressure, ux, uy);
        }
    }

    fclose(plate_output_file);

    plate_output_no++; // Increments output number
    }
}

event gfs_output (t += GFS_OUTPUT_TIMESTEP) {
/* Saves a gfs file */
    if ((t >= START_OUTPUT_TIME) && (t <= END_OUTPUT_TIME)) {
    char gfs_filename[80];
    sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
    output_gfs(file = gfs_filename);

    gfs_output_no++;
    }
}

// event images (t += 0.005) {
// /* Produces movies and images using bview */

//     // Creates a string with the time to put on the plots
//     char time_str[80];
//     sprintf(time_str, "t = %g\n", t);

//     // Set up bview box
//     view (width = 512, height = 512, fov = 30, ty = -0.5, \
//         quat = {0, 0, -0.707, 0.707});

//     // Movie of the volume fraction of the droplet
//     clear();
//     draw_vof("f", lw = 2);
//     draw_vof("plate", lw = 2);
//     squares("f", linear = true);
//     squares("plate", linear = true);
//     mirror ({0,1}) {
//         draw_vof("f", lw = 2);
//         squares("f", linear = true);
//         draw_vof("plate", lw = 2);
//         squares("plate", linear = true);
//     }
//     draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
//     save ("tracer.mp4");
// }

event end (t = MAX_TIME) {
    end_wall_time = omp_get_wtime();
    fprintf(stderr, "Finished after %g seconds\n", end_wall_time - start_wall_time);
}


double distance_from_plate (double xp, double yp) {
/* Returns the distance the point (xp, yp) is from the interface of the plate */

    // Midpoint of the plate in the x axis
    double plate_midpoint = plate_position - 0.5 * PLATE_THICKNESS;

    /* If y is less than PLATE_WIDTH, then the closest interface point will be 
        the straight edges before the curved end point. Else the closest point 
        will be along the curved end */
    if (yp < PLATE_WIDTH) {
        if (xp > plate_midpoint) {
            return fabs(xp - (plate_midpoint + 0.5 * PLATE_THICKNESS));
        } else {
            return fabs(xp - (plate_midpoint - 0.5 * PLATE_THICKNESS));
        }
    } else
    {
        return fabs(sqrt(sq(xp - plate_midpoint) + sq(yp - PLATE_WIDTH)) \
            - 0.5 * PLATE_THICKNESS);
    }
}

double plate_region (double xp, double yp) {
/* Return is < 1 if (xp, yp) is within the plate region and a 
    > 1 if it is outside. The plate is a rectangular region with a 
    semi-circular end */

    // Midpoint of the plate in the x axis
    double plate_midpoint = plate_position - 0.5 * PLATE_THICKNESS;

    /* Deviations from the plate region, which are greater than 1 if xp or yp is
        outside the plate region */
    // Vertical deviation from the plate
    double xp_deviation = fabs(xp - plate_midpoint) / (0.5 * PLATE_THICKNESS);
    // Deviation from the circle at the end
    double circle_deviation = \
        sqrt((sq(xp - plate_midpoint) + sq(yp - PLATE_WIDTH)) \
            / sq(0.5 * PLATE_THICKNESS));
    // Horizontal deviation from either the circle or the flat parts
    double yp_deviation = min(yp / PLATE_WIDTH, circle_deviation);

    // Returns the greatest value of the deviation
    return max(xp_deviation, yp_deviation);
}

double force_on_plate() {
/* Calculates the force that the fluid is exerting on the top of the plate */
    double force = 0; // Initialse force to be zero
    int found_cell = 0;
    foreach() {
        if ((plate[1, 0] == 0) && (plate[-1, 0] == 1) \
            && (x >= plate_position - 0.5 * Delta) && y < (PLATE_WIDTH)) {
                found_cell = 1;
                double u_x_deriv = (u.x[2, 0] - u.x[1, 0]) / Delta;

                double ff = f[1, 0];
                double avg_mu = ff * (mu1 - mu2) + mu2;
                // Work out how to take into account variable visosity
                force += y * y * Delta * (p[1, 0] - 2 * avg_mu * u_x_deriv);
                // force += y * y * Delta * (p[1, 0]);
        }
    }
    if (found_cell == 0) {
        fprintf(stderr, "No cell found\n");
    }
    force = 2 * pi * force;
    return force;
}