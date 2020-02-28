/* plate_impact.c

*/

#include "parameters.h" // Include all defined parameters
#include "axi.h" // Axisymmetric coordinates
#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h" // Implements two-phase flow
#include "view.h" // Creating movies using bview
#include "tension.h" // Adds forces due to surface tension
#include "tag.h" // For removing small droplets
#include <omp.h> // For openMP parallel

/* Physical constants derived from parameters */
double DROP_CENTRE; // Initial position of the centre of the droplet
double PLATE_HEIGHT; // Height of the plate

/* Computational constants derived from parameters */
double MIN_CELL_SIZE; // Size of the smallest cell
double DROP_REFINED_WIDTH; // Width of the refined area around the droplet
double PLATE_REFINED_WIDTH; // Width of the refined area around the plate

/* Boundary conditions */
// Outflow through boundaries far from impact
u.n[right] = neumann(0.); 
u.n[top] = neumann(0.);

// Zero pressure far from impact
p[right] = dirichlet(0.); 
p[top] = dirichlet(0.); 

// No-slip or permeability at the wall
u.t[left] = dirichlet(0.);
u.n[left] = dirichlet(0.);

/* Fields */
scalar plate[]; // VOF field for the plate

/* Data output constants */
int gfs_output_no = 1; // Records how many GFS files have been outputted
int plate_output_no = 1; // Records how many plate data files there have been
int interface_output_no = 1; // Records how many interface files there have been

int main() {
/* Main function for running the simulation */

    /* Create the computational domain */
    init_grid(1 << MINLEVEL); // Create grid according to the minimum level
    size(BOX_WIDTH); // Size of the domain

    /* Set physical constants */
    rho1 = 1.; // Density of water phase
    rho2 = RHO_R; // Density of+= 0.05 air phase
    mu1 = 1. / REYNOLDS; // Viscosity of water phase
    mu2 = mu1 * MU_R; // Viscosity of air phase
    f.sigma = 1. / WEBER; // Surface tension at interface

    /* Derived constants */
    MIN_CELL_SIZE = BOX_WIDTH / pow(2, MAXLEVEL); // Size of the smallest cell
    DROP_REFINED_WIDTH = 0.05; 
    PLATE_REFINED_WIDTH = 20 * MIN_CELL_SIZE;

    fprintf(stderr, "Plate alignment = %d\n", PLATE_ALIGNMENT);
    if (PLATE_ALIGNMENT == -1) {
        PLATE_HEIGHT = 0;
    } else if (PLATE_ALIGNMENT == 0) {
        PLATE_HEIGHT = 1000. * MIN_CELL_SIZE;
    } else {
        PLATE_HEIGHT = (1000. + 1. / PLATE_ALIGNMENT) * MIN_CELL_SIZE;
    }
    
    DROP_CENTRE = PLATE_HEIGHT + INTITIAL_DROP_HEIGHT + DROP_RADIUS;
    run();
}

event init (t = 0) {
/* Initialises the simulation */

    /* Initial refinement */
    // Refine around the plate
    // refine(x < PLATE_HEIGHT + PLATE_REFINED_WIDTH \
    //     && x > PLATE_HEIGHT - 0.2 * PLATE_REFINED_WIDTH && level < MAXLEVEL);

    // Refine around the droplet
    if (PLATE_ALIGNMENT == -1) {
        refine(sq(x - DROP_CENTRE) + sq(y) < sq(DROP_RADIUS +  DROP_REFINED_WIDTH) \
            && sq(x - DROP_CENTRE) + sq(y) > sq(DROP_RADIUS - DROP_REFINED_WIDTH) \
            && level < MAXLEVEL);
    } else {
        refine((((sq(x - DROP_CENTRE) + sq(y) < sq(DROP_RADIUS + DROP_REFINED_WIDTH)) \
            && (sq(x - DROP_CENTRE) + sq(y) > sq(DROP_RADIUS - DROP_REFINED_WIDTH))) \
        || ((x < PLATE_HEIGHT + PLATE_REFINED_WIDTH) \
            && (x > PLATE_HEIGHT - 0.2 * PLATE_REFINED_WIDTH))) \
        && level < MAXLEVEL);
    }
    refine(x < PLATE_HEIGHT + PLATE_REFINED_WIDTH \
        && x > PLATE_HEIGHT - 0.2 * PLATE_REFINED_WIDTH && level < MAXLEVEL);

    

    /* Initialise volume fractions */
    fraction(f, - sq(x - DROP_CENTRE) - sq(y) + sq(DROP_RADIUS)); // Droplet
    fraction(plate, PLATE_HEIGHT - x); // Plate

    /* Initialise droplet velocity */
    foreach() {
        u.x[] = DROP_VEL * f[];
    }
}


event refinement (i++) {
/* Refines the grid where appropriate */
    
    // Adapts with respect to velocities and volume fractions 
    adapt_wavelet ({u.x, u.y, f}, (double[]){1e-2, 1e-2, 1e-4}, 
        minlevel = MINLEVEL, maxlevel = MAXLEVEL);
    
    /* Refines around the boundaries of the plate */
    refine(x < PLATE_HEIGHT + PLATE_REFINED_WIDTH \
        && x > PLATE_HEIGHT - 0.2 * PLATE_REFINED_WIDTH \
        && level < MAXLEVEL);

}

event plate_definition (i++) {
/* Redefines the plate phase at each timestep in order to keep it stationary */

    /* Definition of the plate volume fraction */
    fraction(plate, PLATE_HEIGHT - x);

    /* Alters the velocity of the fluid depending on if is touching the plate */
    foreach() {
        // Stephane's trick to make plate velocity zero
        u.x[] = (1 - plate[]) * u.x[];

        // Sets horizontal velocity inside plate to be zero
        u.y[] = (1. - plate[]) * u.y[]; 
    }
    boundary ((scalar *){u}); // Redefine boundary conditions for u
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

event output_values_along_plate (t += 0.001) {
/* Outputs the values of pressure and velocities along the plate */

    // Creates the file for outputting data
    char plate_output_filename[80];
    sprintf(plate_output_filename, "plate_output_%d.txt", plate_output_no);
    FILE *plate_output_file = fopen(plate_output_filename, "w");

    // Adds the time to the first line of the file
    fprintf(plate_output_file, "t = %g\n", t);

    if (PLATE_ALIGNMENT == -1) {
        /* Iterates over the computational boundary */
        foreach_boundary(left) {
            /* Outputs the values along the plate and the two cells above it */
                for (int h = 0; h <= 2; h++) {
                    fprintf(plate_output_file, \
                        "y = %g, h = %d, p = %g, u_x = %g, u_y = %g\n",
                        y, h, p[h, 0], u.x[h, 0], u.y[h, 0]);
                }
        }
    } else if (PLATE_ALIGNMENT == 0) {
        /* Iterates over the cells with the plate on their lower boundary */
        foreach() {
            // Determines if current cell is along the boundary of the plate
            if ((plate[1, 0] == 0) && (plate[-1, 0] == 1) && x > PLATE_HEIGHT) {
                /* Outputs the values along the plate and the two cells above it */
                for (int h = 0; h <= 2; h++) {
                    fprintf(plate_output_file, \
                        "y = %g, h = %d, p = %g, u_x = %g, u_y = %g\n",
                        y, h, p[h, 0], u.x[h, 0], u.y[h, 0]);
                }
            }
        }
    } {
        /* Iterates over the cells on the boundary of the plate */
        foreach() {
            // Determines if current cell is along the boundary of the plate
            if ((plate[1, 0] == 0) && (plate[-1, 0] == 1)) {

                /* Outputs the values along the plate and the two cells above it */
                for (int h = 0; h <= 2; h++) {
                    fprintf(plate_output_file, \
                        "y = %g, h = %d, p = %g, u_x = %g, u_y = %g\n",
                        y, h, p[h, 0], u.x[h, 0], u.y[h, 0]);
                }
            }
        }
    }

    fclose(plate_output_file);

    plate_output_no++; // Increments output number
}


event output_volume (t += 0.001) {
/* Outputs the volume of the liquid phase to the log file */
    fprintf(stderr, "t = %g, volume = %g\n", t, 2 * pi * statsf(f).sum);
}


event output_interface (t += 0.001) {
/* Outputs the interface locations of the droplet */

    // Creates text file to save output to
    char interface_filename[80];
    sprintf(interface_filename, "interface_%d.txt", interface_output_no);
    FILE *interface_file = fopen(interface_filename, "w");

    // Outputs the interface locations and closes the file
    output_facets(f, interface_file);
    fclose(interface_file);

    interface_output_no++;
}

event gfs_output (t += 0.01) {
/* Saves a gfs file */

    char gfs_filename[80];
    sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
    output_gfs(file = gfs_filename);

    gfs_output_no++;
}

event images (t += 0.005) {
/* Produces movies and images using bview */

    // // Voriticty calculation
    // scalar omega[];
    // vorticity (u, omega);

    // Set up bview box
    view (width = 512, height = 512, fov = 20, ty = -0.5, \
        quat = {0, 0, -0.707, 0.707});

    // Movie of the volume fraction of the droplet
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
    save ("tracer.mp4");
}

event end (t = 0.75) {
    fprintf(stderr, "Finished simulation");
}
