/* impact_coupled.c
    An axisymmetric droplet falling towards an impermeable plate which can move
    in the vertical direction depending on the ODE:
    alpha s''(t) + beta s'(t) + gamma s(t) = F(t)
    where the plate is at z = INITIAL_PLATE_TOP + s(t) at time t, 
    s(0) = s'(0) = 0, F(t) is the force on the plate and alpha, beta, gamma are
    parameters.

    We run the simulation until the turnover point reaches the droplet radius. 
    According to Wagner theory, for DROP_VEL = 1 and DROP_RADIUS = 1, the
    turnover point is at
        r = d(t) = sqrt(3 * (t - IMPACT_TIME - s(t)))
    where IMPACT_TIME is the time at which the droplet hits the plate.
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
double DROP_CENTRE; // Initial centre of the droplet
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
double start_wall_time; // Time the simulation was started
double end_wall_time; // Time the simulation finished
int gfs_output_no = 1; // Records how many GFS files have been outputted
int plate_output_no = 1; // Records how many plate data files there have been
int interface_output_no = 1; // Records how many interface files there have been
// Stores time the interface was outputted
char interface_time_filename[80] = "interface_times.txt"; 

/* Plate position variables, which is at x = INITIAL_PLATE_TOP - s(t) */
double s_previous, s_current, s_next; // Values of s 
double ds_dt; // First time derivative of s
double d2s_dt2; // Second time derivative of s
double plate_vel; // Velocity of plate, equal to -s'(t)
double plate_position; // x position of plate = INITIAL_PLATE_TOP - s(t)

/* Function declarations */
double plate_region(double xp, double yp); // Defines VOF field for plate
double distance_from_plate (double xp, double yp); // Gives distance from plate 
double force_on_plate (); // Calculates the force on the plate via integration


int main() {
/* Main function for running the simulation */

    /* Create the computational domain */
    init_grid(1 << MINLEVEL); // Create grid according to the minimum level
    size(BOX_WIDTH); // Size of the domain

    /* Set physical constants */
    rho1 = 1.; // Density of water phase
    rho2 = RHO_R; // Density of air phase
    mu1 = 1. / REYNOLDS; // Viscosity of water phase
    mu2 = mu1 * MU_R; // Viscosity of air phase
    f.sigma = 1. / WEBER; // Surface tension at interface
    
    /* Initialise plate parameters */
    s_current = 0.; // Initialise displacement as 0
    s_previous = 0.; // Ensures initial velocity is zero
    plate_vel = 0.; // Initial zero plate velocity
    plate_position = INITIAL_PLATE_TOP; // Initialises top of plate position

    /* Derived constants */
    MIN_CELL_SIZE = BOX_WIDTH / pow(2, MAXLEVEL); // Size of the smallest cell
    PLATE_REFINED_WIDTH = 0.3 * PLATE_THICKNESS; // Refined region around plate
    DROP_REFINED_WIDTH = 0.05; // Refined region around droplet
    DROP_CENTRE = INITIAL_PLATE_TOP + INITIAL_DROP_HEIGHT + DROP_RADIUS;
    IMPACT_TIME = INITIAL_DROP_HEIGHT / (-DROP_VEL);
    INTERPOLATE_DISTANCE = MIN_CELL_SIZE; // Distance above plate to read pressure

    // Maximum time is shortly after the Wagner theory prediction of the 
    // turnover point reaching the radius of the droplet
    double wagner_max_time = 1.5 * (IMPACT_TIME + 1. / 3.);
    MAX_TIME = min(HARD_MAX_TIME, wagner_max_time); 

    /* Initialises interface time file */
    FILE* interface_time_file = fopen(interface_time_filename, "w");
    fclose(interface_time_file);

    /* Run the simulation */
    run();
}


event init (t = 0) {
/* Initialises the computational domain */

    // Records the wall time
    start_wall_time = omp_get_wtime();

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
        u.x[] = plate[] * plate_vel + (1 - plate[]) * u.x[]; // Plate velocity
    }
    boundary ((scalar *){u});
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

        // Appends the interface time file with the time and plate position (0)
        FILE *interface_time_file = fopen(interface_time_filename, "a");
        fprintf(interface_time_file, "%d, %g, %g\n", \
            interface_output_no, t, plate_position);
        fclose(interface_time_file);

        interface_output_no++;
    }
}


event output_data (t += PLATE_OUTPUT_TIMESTEP) {
/* Outputs data about the flow */
    if ((t >= START_OUTPUT_TIME) && (t <= END_OUTPUT_TIME)) {
        // Creates the file for outputting data along the plate
        char plate_output_filename[80];
        sprintf(plate_output_filename, "plate_output_%d.txt", plate_output_no);
        FILE *plate_output_file = fopen(plate_output_filename, "w");

        // Adds the time to the first line of the plate output file
        fprintf(plate_output_file, "t = %.4f\n", t);

        // Initialises the force variable
        double force = 0.; 

        // Iterates over the cells above the plate
        foreach(reduction(+:force)) {
            // Identifies the cells along the plate
            if ((plate[1, 0] == 0) && (plate[-1, 0] == 1) && (y < PLATE_WIDTH)) {
                    
                    /* Force calculation */
                    // Viscosity average in the cell above the plate
                    double avg_mu = f[1, 0] * (mu1 - mu2) + mu2;

                    // Viscous stress in the cell above the plate
                    double viscous_stress = \
                        - 2 * avg_mu * (u.x[2, 0] - u.x[1, 0]) / Delta;

                    // Adds the contribution to the force using trapeze rule
                    force += y * Delta * (p[1, 0] + viscous_stress);

                    /* Plate output */
                    fprintf(plate_output_file, \
                        "y = %.8f, x = %.8f, p = %.8f, strss = %.8f\n",\
                        y, x, p[1, 0], viscous_stress);
            }
        }

        // Close plate output file
        fclose(plate_output_file);
        plate_output_no++; // Increments output number

        // Integrates force about the angular part
        force = 2 * pi * force; 

        /* Outputs data to log file */
        fprintf(stderr, \
            "t = %.8f, v = %.8f, F = %g, s = %g, ds_dt = %g, d2s_dt2 = %g\n", \
            t, 2 * pi * statsf(f).sum, force, s_current, ds_dt, d2s_dt2);
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

event images (t += 0.005) {
 /* Produces movies and images using bview */ 

    // Creates a string with the time to put on the plots
    char time_str[80];
    sprintf(time_str, "t = %g\n", t);

    // Set up bview box
    view (width = 512, height = 512, fov = 30, ty = -0.5, \
        quat = {0, 0, -0.707, 0.707});

    /* Movie of the volume fraction of the droplet */
    clear();
    draw_vof("f", lw = 2);
    draw_vof("plate", lw = 2);
    squares("f", linear = true);
    squares("plate", linear = true);
    mirror ({0,1}) {
        draw_vof("f", lw = 2);
        squares("f", linear = true);
        draw_vof("plate", lw = 2);
        squares("plate", linear = true);
    }
    draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
    save ("tracer.mp4");

    /* Movie of the vertical velocity */
    clear();
    draw_vof("f", lw = 2);
    draw_vof("plate", lw = 2);
    squares("u.x", linear = false);
    mirror ({0,1}) {
        draw_vof("f", lw = 2);
        squares("u.x", linear = false);
        draw_vof("plate", lw = 2);
    }
    draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
    save ("vertical_vel.mp4");


    /* Movie of the horizontal velocity */
    clear();
    draw_vof("f", lw = 2);
    draw_vof("plate", lw = 2);
    squares("u.y", linear = false);
    mirror ({0,1}) {
        draw_vof("f", lw = 2);
        squares("u.y", linear = false);
        draw_vof("plate", lw = 2);
    }
    draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
    save ("horizontal_vel.mp4");

    /* Movie of the pressure */
    clear();
    draw_vof("f", lw = 2);
    draw_vof("plate", lw = 2);
    squares("p", linear = false);
    mirror ({0,1}) {
        draw_vof("f", lw = 2);
        squares("p", linear = false);
        draw_vof("plate", lw = 2);
    }
    draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
    save ("pressure.mp4");
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


event moving_plate (i++) {
/* Moves the plate as a function of the force on it */

    /* Calculate the force on the plate by integrating using trapeze rule*/
    double force = 0.; // Initialise to be zero

    // Iterates over the cells above the plate
    foreach(reduction(+:force)) {
        // Identifies the cells along the plate
        if ((plate[1, 0] == 0) && (plate[-1, 0] == 1) && (y < PLATE_WIDTH)) {
                // Viscosity average in the cell above the plate
                double avg_mu = f[1, 0] * (mu1 - mu2) + mu2;

                // Viscous stress in the cell above the plate
                double viscous_stress = \
                    - 2 * avg_mu * (u.x[2, 0] - u.x[1, 0]) / Delta;

                // Adds the contribution to the force using trapeze rule
                force += y * Delta * (p[1, 0] + viscous_stress);
        }
    }

    force = 2 * pi * force; // Integrates about the angular part

    /* Solves the ODE for the updated plate position and acceleration */
    s_next = (dt * dt * force \
        + (2. * ALPHA - dt * dt * GAMMA) * s_current \
        - (ALPHA - dt * BETA / 2.) * s_previous) \
        / (ALPHA + dt * BETA / 2.);

    /* Updates values of s */
    ds_dt = (s_next - s_previous) / (2. * dt);
    d2s_dt2 = (s_next - 2 * s_current + s_previous) / (dt * dt);
    s_previous = s_current;
    s_current = s_next; 

    /* Updates plate position and velocity */
    plate_position = INITIAL_PLATE_TOP - s_current;
    plate_vel = -ds_dt;

    /* Redefine the plate volume fraction */
    fraction(plate, 1 - plate_region(x, y));
    
    /* Alters the velocity of the fluid depending on if is touching the plate */
    foreach() {
        // Sets vertical velocity to match the plate at the interface
        u.x[] = plate[] * plate_vel + (1 - plate[]) * u.x[];

        // Sets horizontal velocity inside plate to be zero
        u.y[] = (1. - plate[]) * u.y[]; 
    }
    boundary ((scalar *){u}); // Redefine boundary conditions for u
}


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