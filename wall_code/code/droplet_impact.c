/* droplet_impact.c
    A Basilisk script to model the impact of a droplet of water impacting onto a
    rigid surface. 

    Set up to compare to the case where the rigid surface is defined using a 
    volume fraction
*/

#include "parameters.h" // Includes all defined parameters
#include "axi.h" // Axisymmetric coordinates
#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h" // Implements two-phase flow
#include "view.h" // Creating movies using bview
#include "tension.h" // Surface tension of droplet
#include "tag.h" // For removing small droplets
#include <omp.h> // For openMP parallel

/* Physical constants */
double DROP_CENTRE; // Initial position of the centre of the droplet

/* Computational constants */
double DROP_REFINED_WIDTH; // Width of refined areas around droplet
double REFINED_BOX_HEIGHT; // Height of the refined areas around substrate
double REFINED_BOX_WIDTH; // Width of the refined area around the substrate
double MIN_CELL_SIZE; // Size of the cells at the most refined

/* Outputting constants */
int gfs_output_no = 1; // Used to count how many gfsview outputs there have been
int interface_output_no = 1; // How many interface outputs there have been
int pressure_output_no = 1; // How many pressure outputs there have been
// Stores times the interface was outputted
char interface_time_filename[] = "interface_times.txt"; 
char droplet_area_filename[] = "area.txt";


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

    /* Derived constants */
    DROP_CENTRE = INTITIAL_DROP_HEIGHT + DROP_RADIUS; // Centre of droplet
    DROP_REFINED_WIDTH = 0.1 * DROP_RADIUS; // Width of refined around droplet
    REFINED_BOX_HEIGHT = 0.02 * DROP_RADIUS; // Height of refined box around wall
    REFINED_BOX_WIDTH = BOX_WIDTH / 4.; // Width of refined box around wall
    MIN_CELL_SIZE = BOX_WIDTH / pow(2, MAXLEVEL); // Size of the smallest cell

    init_grid(1 << MINLEVEL); // Create grid according to the minimum level

    size(BOX_WIDTH); // Size of the domain

    /* Set physical constants */
    rho1 = 1.; // Density of water phase
    rho2 = RHO_R; // Density of air phase
    mu1 = 1. / REYNOLDS; // Viscosity of water phase
    mu2 = mu1 * MU_R; // Viscosity of air phase
    f.sigma = 1. / WEBER; // Surface tension at interface

    run(); // Runs the simulation
}


event init(t = 0) {
/* Initialises the flow as a spherical droplet falling downwards */

    /* Initial refinement */
    // Refines region around the droplet
    refine(sq(x - DROP_CENTRE) + sq(y) < sq(DROP_RADIUS + DROP_REFINED_WIDTH) \
        && sq(x - DROP_CENTRE) + sq(y) > sq (DROP_RADIUS - DROP_REFINED_WIDTH) \
        && level < MAXLEVEL);

    // Refine around the substrate
    refine(x < REFINED_BOX_HEIGHT && y < REFINED_BOX_WIDTH && level < MAXLEVEL);

    /* Initialise droplet volume fraction */ 
    fraction(f, - sq(x - DROP_CENTRE) - sq(y) + sq(DROP_RADIUS)); 
    

    /* Initialise the droplet velocity downwards */
    foreach() {
        u.x[] = -f[];
    }
    
     /* Initialise the file for outputting the times at which interfaces are 
        saved */
    FILE *interface_time_file = fopen(interface_time_filename, "w");
    fclose(interface_time_file);
}


event acceleration (i++) {
/* Adds acceleration due to gravity at each time step */
    face vector av = a; // Acceleration at each face
    foreach_face(x) av.x[] -= 1./sq(FR); // Adds acceleration due to gravity
}


event refinement (i++) {
/* Refines the grid where appropriate */

    // Adapt the mesh for rapid changes in the free surface, velocity and 
    // vorticity.
    adapt_wavelet ({u.x, u.y, f}, (double[]){1e-2, 1e-2, 1e-2},
        maxlevel = MAXLEVEL, minlevel = MINLEVEL);
    
    // Refine around the substrate
    refine(x < REFINED_BOX_HEIGHT && y < REFINED_BOX_WIDTH && level < MAXLEVEL);
}


// event small_droplet_removal (i++) {
// /* Removes any small droplets that have formed, that are smaller than a specific    
//     size */
//     remove_droplets(f, 20); // Removes droplets of diameter 6 cells or less
// }


event pressure_output (t += 0.01) {
/* Outputs the pressure along the surface of the plate */

    // Creates and opens the file to output pressure into
    char pressure_filename[80];
    sprintf(pressure_filename, "pressure_%d.txt", pressure_output_no);
    FILE *pressure_file = fopen(pressure_filename, "w");

    foreach() {
        // If x < MIN_CELL_SIZE, then we are currently at the cell directly 
        // above the wall
        if (x < MIN_CELL_SIZE ) {
            // Outputs the value of pressure to the file
            fprintf(pressure_file, "%g %g %g\n", y, x, p[]);
        }
    }

    fclose(pressure_file); // Closes the pressure file

    pressure_output_no++; // Increments the number of outputs there have been
}

event images (t += 0.01; t <= 2.) {
/* Produces movies and images using bview */

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


event gfs_output (t += 0.01) {
    // Saves a gfs file

    char gfs_filename[80];
    sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
    output_gfs(file = gfs_filename);

    gfs_output_no++;
}


event interface_output (t += 0.01) {
/* Outputs the location of the interfaces of the droplet and the plate */

    /* Filenames for the outputs */
    char droplet_interface_filename[80];
    sprintf(droplet_interface_filename, \
        "droplet_interface_%d.txt", interface_output_no);
    FILE *droplet_interface_file = fopen(droplet_interface_filename, "w");

    
    /* Appends the timestep number and the time for this specific output_no, so 
    that in post-processing we know what these are for each interface file */
    FILE *interface_time_file = fopen(interface_time_filename, "a");
    fprintf(interface_time_file, "%d %d %g\n", interface_output_no, i, t); 
    fclose(interface_time_file);

    /* Outputs the interface locations */
    output_facets(f, droplet_interface_file); // Droplet output
    fclose(droplet_interface_file);

    /* Outputs the area of the droplet */
    FILE *area_file = fopen(droplet_area_filename, "a");
    fprintf(area_file, "%g %g\n", t, 2 * pi * statsf(f).sum);
    fclose(area_file);

    interface_output_no++; // Increments number of outputs
}



