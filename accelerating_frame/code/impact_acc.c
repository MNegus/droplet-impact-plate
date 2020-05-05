/* impact_acc.c
    A Basilisk script to model the impact of a droplet of water impacting onto a
    moving plate. The domain is set to be in an accelerating frame with the 
    plate, so an additional body force is added.
*/

#include "parameters.h" // Includes all defined parameters
#include "axi.h" // Axisymmetric coordinates
#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h" // Implements two-phase flow
#include "view.h" // Creating movies using bview
#include "tension.h" // Surface tension of droplet
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

/* Global variables */
double start_wall_time; // Time the simulation was started
double end_wall_time; // Time the simulation finished
int gfs_output_no = 1; // Records how many GFS files have been outputted
int plate_output_no = 1; // Records how many plate data files there have been
int interface_output_no = 1; // Records how many interface files there have been
// Stores time the interface was outputted
char interface_time_filename[80] = "interface_times.txt"; 

/* Plate position variables */
double s_previous, s_current, s_next;
double d2s_dt2; // Second time derivative of s

/* Boundary conditions */
// Conditions for entry from above
u.n[right] = neumann(0.); // Allows outflow through boundary
p[right] = dirichlet(0.); // 0 pressure far from surface

// Conditions far from the droplet in the radial direction
u.n[top] = neumann(0.); // Allows outflow through boundary
p[top] = dirichlet(0.); // 0 pressure far from surface

// Conditions on surface
u.n[left] = dirichlet(0.); // No flow through surface
u.t[left] = dirichlet(0.); // No slip at surface

int main() {

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
    DROP_REFINED_WIDTH = 0.05; // Refined region around droplet
    PLATE_REFINED_WIDTH = 0.3 * PLATE_THICKNESS; // Refined region around plate
    DROP_CENTRE = INITIAL_DROP_HEIGHT + DROP_RADIUS;
    IMPACT_TIME = INITIAL_DROP_HEIGHT / (-DROP_VEL);
    INTERPOLATE_DISTANCE = MIN_CELL_SIZE; // Distance above plate to read pressure

    /* Maximum time is shortly after Wagner theory would predict the turnover 
    point reaches the radius of the droplet */
    double wagner_max_time = 2.0 * (IMPACT_TIME + 1. / 3.);
    MAX_TIME = min(HARD_MAX_TIME, wagner_max_time);

    /* Initialises interface time file */
    FILE* interface_time_file = fopen(interface_time_filename, "w");
    fclose(interface_time_file);

    run(); // Runs the simulation
}


event init(t = 0) {
/* Initialises the flow as a spherical droplet falling downwards */

    // Records the wall time
    start_wall_time = omp_get_wtime();

    /* Refines around the droplet */
    refine(sq(x - DROP_CENTRE) + sq(y) < sq(DROP_RADIUS + DROP_REFINED_WIDTH) \
        && sq(x - DROP_CENTRE) + sq(y) > sq(DROP_RADIUS - DROP_REFINED_WIDTH) \
        && level < MAXLEVEL);
    
    /* Initialises the droplet volume fraction */
    fraction(f, -sq(x - DROP_CENTRE) - sq(y) + sq(DROP_RADIUS));

    /* Initialise the droplet velocity downwards */
    foreach() {
        u.x[] = DROP_VEL * f[];
    }

    /* Initialise the plate position, velocity and acceleration to be zero */
    s_previous = 0.;
    s_current = 0.;
    d2s_dt2 = 0.;
}

event output_stats (t += PLATE_OUTPUT_TIMESTEP) {
/* Outputs the stats about the flow*/
    if ((t >= START_OUTPUT_TIME) && (t <= END_OUTPUT_TIME)) {
        // Creates the file for outputting data along the plate
        char plate_output_filename[80];
        sprintf(plate_output_filename, "plate_output_%d.txt", plate_output_no);
        FILE *plate_output_file = fopen(plate_output_filename, "w");

        // Adds the time to the first line of the file
        fprintf(plate_output_file, "t = %g\n", t);

        // Initialises the force variable
        double force = 0.; 

        // Loops over the left hand boundary (the interface)
        foreach_boundary(left, reduction(+:force)) {
            /* Force calculation */
            // Viscosity average in the cell above the plate
            double avg_mu = f[1, 0] * (mu1 - mu2) + mu2;

            // Viscous stress in the cell above the plate
            double viscous_stress = - 2 * avg_mu * (u.x[2, 0] - u.x[1, 0]) / Delta;

            // Adds the contribution to the force using trapeze rule
            force += y * Delta * (p[1, 0] + viscous_stress);

            /* Plate output */
            fprintf(plate_output_file, "y = %g, p = %g, strss = %g\n",\
                 y, p[1, 0], viscous_stress);
        }

        // Close plate output file
        fclose(plate_output_file);
        plate_output_no++; // Increments output number

        // Integrates force about the angular part
        force = 2 * pi * force; 

        /* Outputs to force and volume to log file */
        fprintf(stderr, "t = %g, v = %g, F = %g, s = %g, d2s_dt2 = %g\n", \
            t, 2 * pi * statsf(f).sum, force, s_current, d2s_dt2);
    }
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
            interface_output_no, t, 0.);
        fclose(interface_time_file);

        interface_output_no++;
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


event images (t += 0.01) {
// Produces movies and images using bview 

    // Voriticty calculation
    scalar omega[];
    vorticity (u, omega);

    // Set up bview box
    view (width = 512, height = 512, fov = 20, ty = -0.5, quat = {0, 0, -0.707, 0.707});

    // Movie of the volume fraction of the droplet
    clear();
    draw_vof("f", lw = 2);
    squares("f", linear = true);
    mirror ({0,1}) {
        draw_vof("f", lw = 2);
        squares("f", linear = true);
    }
    save ("tracer.mp4");

    // Movie of the vorticity
    clear();
    squares("omega", linear = true);
    draw_vof("f", lw = 2);
    mirror ({0,1}) {
        draw_vof("f", lw = 2);
        squares("omega", linear = true);
    }
    save ("vort.mp4");
}


event refinement (i++) {
/* Refines the grid where appropriate */

    /* Adapts with respect to velocities and volume fraction */
    adapt_wavelet ({u.x, u.y, f}, (double[]){1e-2, 1e-2, 1e-2},
        minlevel = MINLEVEL, maxlevel = MAXLEVEL);
    
    /* Refines above the plate */
    refine((y < PLATE_WIDTH) && (x < 0.5 * PLATE_REFINED_WIDTH) \
        && level < MAXLEVEL);
}


event acceleration (i++) {
/* Adds acceleration due to gravity and the moving plate at each time step */
    face vector av = a; // Acceleration at each face

    /* Adds acceleration due to gravity and the plate */
    foreach_face(x){
        av.x[] += d2s_dt2 - 1./sq(FR);
    }
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

    // Iterates over the solid boundary
    foreach_boundary(left, reduction(+:force)) {
        // Viscosity average in the cell above the plate
        double avg_mu = f[1, 0] * (mu1 - mu2) + mu2;

        // Viscous stress in the cell above the plate
        double viscous_stress = - 2 * avg_mu * (u.x[2, 0] - u.x[1, 0]) / Delta;

        // Adds the contribution to the force using trapeze rule
        force += y * Delta * (p[1, 0] + viscous_stress);
    }

    force = 2 * pi * force; // Integrates about the angular part

    /* Solves the ODE for the updated plate position and acceleration */
    s_next = (dt * dt * force \
        + (2. * ALPHA - dt * dt * GAMMA) * s_current \
        - (ALPHA - dt * BETA / 2.) * s_previous) \
        / (ALPHA + dt * BETA / 2.);

    // Updates second derivative
    d2s_dt2 = (s_next - 2 * s_current + s_previous) / (dt * dt);

    // Updates current and previous values of s
    s_previous = s_current;
    s_current = s_next; 

    // Redefine boundary conditions for u
    boundary ((scalar *){u}); 
}


event end (t = MAX_TIME) {
    end_wall_time = omp_get_wtime();
    fprintf(stderr, "Finished after %g seconds\n", end_wall_time - start_wall_time);
}
