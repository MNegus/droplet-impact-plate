/* impact_wall.c
    An axisymmetric droplet falling towards an impermeable, rigid surface lying 
    along the boundary at z = 0. 

    Run until the turnover point approximately reaches the initial droplet 
    radius. 

    Author: Michael Negus
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
double DROP_CENTRE; // Initial centre position of the droplet
double IMPACT_TIME; // Theoretical time of impact
double MAX_TIME; // Maximum time to run the simulation for

/* Global variables */
double start_wall_time; // Time the simulation was started
double end_wall_time; // Time the simulation finished
int gfs_output_no = 1; // Records how many GFS files have been outputted
int plate_output_no = 1; // Records how many plate data files there have been
int interface_output_no = 1; // Records how many interface files there have been
// Stores time the interface was outputted
char interface_time_filename[80] = "interface_times.txt"; 

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

    /* Derived constants */
    MIN_CELL_SIZE = BOX_WIDTH / pow(2, MAXLEVEL); // Size of the smallest cell
    PLATE_REFINED_WIDTH = 0.3 * PLATE_THICKNESS; // Refined region around plate
    DROP_REFINED_WIDTH = 0.05; // Refined region around droplet
    DROP_CENTRE = INITIAL_DROP_HEIGHT + DROP_RADIUS; // Initial centre of drop
    IMPACT_TIME = INITIAL_DROP_HEIGHT / (-DROP_VEL); // Theoretical impact time

    /* Maximum time is shortly after Wagner theory would predict the turnover 
    point reaches the radius of the droplet */
    double wagner_max_time = 1.5 * (IMPACT_TIME + 1. / 3.);
    MAX_TIME = min(HARD_MAX_TIME, wagner_max_time);

    /* Initialises interface time file */
    FILE* interface_time_file = fopen(interface_time_filename, "w");
    fclose(interface_time_file);

    /* Runs the simulation */
    run(); 
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
    boundary ((scalar *){u});
}


event refinement (i++) {
/* Adaptive grid refinement */

    // Adapts with respect to velocities and volume fraction 
    adapt_wavelet ({u.x, u.y, f}, (double[]){1e-2, 1e-2, 1e-2},
        minlevel = MINLEVEL, maxlevel = MAXLEVEL);
    
    /* Refines above the plate */
    refine((y < PLATE_WIDTH) && (x < 0.5 * PLATE_REFINED_WIDTH) \
        && level < MAXLEVEL);
}


event gravity (i++) {
/* Adds acceleration due to gravity in the vertical direction */
    face vector av = a; // Acceleration at each face
    foreach_face(x) av.x[] -= 1./sq(FR); // Adds acceleration due to gravity
}


event small_droplet_removal (i++) {
/* Removes any small droplets or bubbles that have formed, that are smaller than
    a specific size */
    // Removes droplets of diameter 5 cells or less
    remove_droplets(f, 5);
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
        foreach_boundary(left,reduction(+:force)) {
            if (y < PLATE_WIDTH) {
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
            "t = %.8f, v = %.8f, F = %.8f\n", t, 2 * pi * statsf(f).sum, force);
    }
}


event output_interface (t += INTERFACE_OUTPUT_TIMESTEP) {
/* Outputs the interface locations of the droplet */
    if ((t >= START_OUTPUT_TIME) && (t <= END_OUTPUT_TIME)) {
        // Creates text file to output interface positions to
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

        // Increments output number
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


event movies (t += 0.005) {
/* Produces movies using bview */ 
    if (MOVIES) {
        // Creates a string with the time to put on the plots
        char time_str[80];
        sprintf(time_str, "t = %g\n", t);

        // Set up bview box
        view (width = 512, height = 512, fov = 30, ty = -0.5, \
            quat = {0, 0, -0.707, 0.707});

        /* Movie of the volume fraction of the droplet */
        clear();
        draw_vof("f", lw = 2);
        squares("f", linear = true);
        mirror ({0,1}) {
            draw_vof("f", lw = 2);
            squares("f", linear = true);
        }
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("tracer.mp4");

        /* Movie of the vertical velocity */
        clear();
        draw_vof("f", lw = 2);
        squares("u.x", linear = false);
        mirror ({0,1}) {
            draw_vof("f", lw = 2);
            squares("u.x", linear = false);
        }
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("vertical_vel.mp4");


        /* Movie of the horizontal velocity */
        clear();
        draw_vof("f", lw = 2);
        squares("u.y", linear = false);
        mirror ({0,1}) {
            draw_vof("f", lw = 2);
            squares("u.y", linear = false);
        }
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("horizontal_vel.mp4");

        /* Movie of the pressure */
        clear();
        draw_vof("f", lw = 2);
        squares("p", linear = false);
        mirror ({0,1}) {
            draw_vof("f", lw = 2);
            squares("p", linear = false);
        }
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("pressure.mp4");
    }
}


event end (t = MAX_TIME) {
/* Ends the simulation */ 

    end_wall_time = omp_get_wtime(); // Records the time of finish

    fprintf(stderr, "Finished after %g seconds\n", \
        end_wall_time - start_wall_time);
}

