/* impact_acc.c
    A Basilisk script to model the impact of a droplet of water impacting onto a
    moving plate. The domain is set to be in an accelerating frame with the 
    plate, so an additional body force is added.
*/

// RC This is commonly used for large viscosity contrasts
#define FILTERED
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

#include "parameters.h" // Includes all defined parameters
#include "axi.h" // Axisymmetric coordinates
#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h" // Implements two-phase flow
#include "view.h" // Creating movies using bview
#include "tension.h" // Surface tension of droplet
#include "tag.h" // For removing small droplets
#include <omp.h> // For openMP parallel

#include "contact.h" // RC for imposing contact angle on surface


/* Computational constants derived from parameters */
double MIN_CELL_SIZE; // Size of the smallest cell
double DROP_REFINED_WIDTH; // Width of the refined area around the droplet
double PLATE_REFINED_WIDTH; // Width of the refined area around the plate
double DROP_CENTRE; // Initial centre of the droplet
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
double pinch_off_time = 0.; // Time pinch-off of the entrapped bubble occurs
double drop_thresh = 1e-4; // Remove droplets threshold

/* Force interpolation */
double * forces_array; // Forces of the previous timesteps
double * times_array; // Times of the previous timesteps
double current_force; // Current force on plate
double force_term; // Force term used in the ODE 

/* Plate position variables */
double s_previous = 0.; // Value of s at previous timestep
double s_current = 0.; // Value of s at current timestep
double s_next; // Values of s at next timestep
double ds_dt; // First time derivative of s
double d2s_dt2; // Second time derivative of s

FILE * fp_stats; // RC
char interp_stats_filename[80] = "interp_stats.txt";

vector h[];  // RC
double theta0 = 90;  // RC

/* Boundary conditions */
// Conditions for entry from above
u.n[right] = neumann(0.); // Free flow condition
p[right] = dirichlet(0.); // 0 pressure far from surface

// Conditions far from the droplet in the radial direction
u.n[top] = neumann(0.); // Allows outflow through boundary
u.t[top] = dirichlet(0.); // Stationary vertical flow
p[top] = dirichlet(0.); // 0 pressure far from surface

// Conditions on surface
u.n[left] = dirichlet(0.); // No flow through surface
u.t[left] = dirichlet(0.); // No slip at surface
h.t[left] = contact_angle (theta0*pi/180.); // RC contact angle

void remove_droplets_region(struct RemoveDroplets p,\
        double ignore_region_x_limit, double ignore_region_y_limit);
        
int main() {
/* Main function to set up the simulation */

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
    DROP_REFINED_WIDTH = 0.04; // Refined region around droplet
    PLATE_REFINED_WIDTH \
        = PLATE_REFINE_NO * MIN_CELL_SIZE; // Refined region width around plate
    DROP_CENTRE = INITIAL_DROP_HEIGHT + DROP_RADIUS;
    IMPACT_TIME = INITIAL_DROP_HEIGHT / (-DROP_VEL);

    /* Maximum time is shortly after Wagner theory would predict the turnover 
    point reaches the radius of the droplet */
    double wagner_max_time = 2.0 * (IMPACT_TIME + 1. / 3.);
    MAX_TIME = min(HARD_MAX_TIME, wagner_max_time);

    /* Allocates memory for the force and times arrays */
    if (INTERPOLATE) {
        forces_array = malloc(INTERP_NO * sizeof(double));
        times_array = malloc(INTERP_NO * sizeof(double));
    }

    /* Initialises interface time file */
    FILE* interface_time_file = fopen(interface_time_filename, "w");
    fclose(interface_time_file);

    /* Initialises interp stats file */
    FILE * interp_stats_file = fopen(interp_stats_filename, "w");
    fclose(interp_stats_file);

    // RC
    {
      char name[200];
      sprintf(name, "logstats.dat");
      fp_stats = fopen(name, "w");
    }

    // RC these give the solver a helping hand (at some extra cost)
    DT = 1.0e-4;
    NITERMIN = 1; // default 1
    NITERMAX = 400; // default 100
    TOLERANCE = 1e-5; // default 1e-3

    run();

    fclose(fp_stats); // RC
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
/* Refines the grid where appropriate */

    /* Adapts with respect to velocities and volume fraction */
    adapt_wavelet ({u.x, u.y, f}, (double[]){1e-2, 1e-2, 1e-6}, // RC changed the interfacial tolerance to something more stringent
        minlevel = MINLEVEL, maxlevel = MAXLEVEL);
    
    /* Refines above the plate */
    refine((y < PLATE_WIDTH) && (x <= PLATE_REFINED_WIDTH) \
        && level < MAXLEVEL);

    /* Refines box around the origin for entrapped bubble */
    // RC this one shouldn't be needed if all else goes well
    // refine((y < 0.05 * DROP_RADIUS) && (x < 0.05 * DROP_RADIUS) \
    //    && level < MAXLEVEL);
}


event moving_plate (i++) {
/* Moves the plate as a function of the force on it */

    /* Calculate the force on the plate by integrating using trapeze rule*/
    current_force = 0.; // Initialise to be zero

    // Iterates over the solid boundary
    foreach_boundary(left, reduction(+:current_force)) {
        if (y < PLATE_WIDTH) {
            // Viscosity average in the cell above the plate
            double avg_mu = f[] * (mu1 - mu2) + mu2;

            // Viscous stress in the cell above the plate
            double viscous_stress = \
                - 2 * avg_mu * (u.x[1, 0] - u.x[]) / Delta;

            // Adds the contribution to the force using trapeze rule
            current_force += y * Delta * (p[] + viscous_stress);
        }
    }

    current_force = 2 * pi * current_force; // Integrates about the angular part

    /* Force interpolation. If the current value of force has deviated from the
    previous value by too much, then we instead interpolate the last two values
    to calculate an approximation for this timestep */
    
    if (INTERPOLATE) {
        // If the interpolate parameter is specified
        if (INTERP_NO != 2) {
            // Current version only works for INTERP_NO = 2
            fprintf(stderr, "INTERP_NO must equal 2\n");
            exit(1);
        }

        if (i < INTERP_NO) {
            /* For the first INTERP_NO timesteps, we populate the force and 
            times arrays with the current force and time values. The forcing 
            term is set to be equal to the current force */ 

            // Populate arrays
            forces_array[i] = current_force;
            times_array[i] = t;

            // Specify forcing term for ODE
            force_term = current_force;

        } else {
            // Else, all the values of the force and times array are populated
            if (t < INTERP_DELAY)  {
                /* Beforen t == INTERP_DELAY, no interpolation is done and the 
                force term in the ODE is taken to be the current force  */
                force_term = current_force;
            } else {
                /* We make a prediction for the force at the current timestep
                using the last INTERP_NO timesteps. If the difference between
                the current force and this prediction is greater than the
                threshold, then we take the forcing term in the ODE to be equal
                to the prediction instead */


                // ONLY WORKS FOR INTERP_NO = 2
                // Calculate predicted value of the force
                double gradient = (forces_array[1] - forces_array[0]) \
                        / (times_array[1] - times_array[0]);
                double f_predict \
                    = forces_array[1] + gradient * (t -  times_array[1]);
                
                // Check difference between current force and predicted value
                if ((current_force < f_predict + INTERP_THRESHOLD) \
                    && (current_force > f_predict - INTERP_THRESHOLD)) {
                    // If the force is within the set threshold of the
                    // prediction
                    force_term = current_force;
                } else {
                    // Else the force is outside the range
                    force_term = f_predict;

                    // Output stats on interpolation
                    FILE *interp_stats_file = fopen(interp_stats_filename, "a");
                    fprintf(interp_stats_file, \
                        "t = %g, F = %g, F0 = %g, t0 = %g, F1 = %g, t1 = %g, F_int = %g\n", \
                        t, current_force, forces_array[0], times_array[0], \
                        forces_array[1], times_array[1], force_term);
                    fprintf(interp_stats_file, "lower_thresh = %g, upper_thresh = %g\n",\
                    f_predict - INTERP_THRESHOLD, f_predict + INTERP_THRESHOLD);
                    fprintf(interp_stats_file, "\n");
                    fclose(interp_stats_file);
                }
            }

            // Populate arrays
            #pragma omp critical
            for (int j = 0; j < INTERP_NO - 1; j++) {
                forces_array[j] = forces_array[j + 1];
                times_array[j] = times_array[j + 1];
                
            }
            forces_array[INTERP_NO - 1] = force_term;
            times_array[INTERP_NO - 1] = t;
        }
    } else {
        force_term = current_force;
    }


    if (t < FORCE_DELAY_TIME) force_term = 0;

    /* Solves the ODE for the updated plate position and acceleration */
    s_next = (dt * dt * force_term \
        + (2. * ALPHA - dt * dt * GAMMA) * s_current \
        - (ALPHA - dt * BETA / 2.) * s_previous) \
        / (ALPHA + dt * BETA / 2.);

    /* Updates values of s and its derivatives */
    ds_dt = (s_next - s_previous) / (2. * dt);
    d2s_dt2 = (s_next - 2 * s_current + s_previous) / (dt * dt);
    s_previous = s_current;
    s_current = s_next; 

    /* OVERRIDE: CONSTANT ACCELERATION */
    if (CONST_ACC) {
        if (t < IMPACT_TIME) {
            d2s_dt2 = 0.;
            ds_dt = 0.;
            s_current = 0.;
        } else {
            d2s_dt2 = PLATE_ACC;
            ds_dt = d2s_dt2 * (t - IMPACT_TIME);
            s_current = 0.5 * d2s_dt2 * (t - IMPACT_TIME) * (t - IMPACT_TIME);
        }
    }

    /* Updates velocity BC's */
    u.t[top] = dirichlet(ds_dt);
    u.n[left] = y < PLATE_WIDTH ? dirichlet(0.) : dirichlet(ds_dt);

    boundary ((scalar *){u}); // Redefine boundary conditions for u
    // boundary((scalar *){p}); // Redefine pressure BC
}


event acceleration (i++) {
/* Adds acceleration due to gravity and the moving plate at each time step */
    face vector av = a; // Acceleration at each face

    /* Adds acceleration due to gravity and the plate */
    foreach_face(x){
        av.x[] += d2s_dt2 - 1./sq(FR);
    }
}


event small_droplet_removal (i += 5) { // RC doing this every n iterations or every 1e-3/4 in time is rather important (to avoid choking the run)
// Removes any small droplets or bubbles that have formed, that are smaller than
//    a specific size

    // Minimum diameter (in cells) a droplet/bubble has to be, else it will be 
    // removed
    int drop_min_cell_width = 16;

    // Region to ignore
    double ignore_region_x_limit = 0.02; //0.01;
    double ignore_region_y_limit = 0.02; //2.2 * 0.05;
    
    // Counts the number of bubbles there are
    scalar bubbles[];
    foreach() {
        bubbles[] = 1. - f[] > drop_thresh;
    }
    int bubble_no = tag(bubbles);

    // Determines if we are before or after the pinch-off time
    if (pinch_off_time == 0.) {
        // The first time the bubble number is above 1, we define it to be the 
        // pinch off time
        if (bubble_no > 1) {
            pinch_off_time = t;
        }
    } else if (t >= pinch_off_time + REMOVAL_DELAY) {
        // If we are a certain time after the pinch-off time, remove drops and 
        // bubbles below the specified minimum size

        struct RemoveDroplets remove_struct;
        remove_struct.f = f;
        remove_struct.minsize = drop_min_cell_width;
        remove_struct.threshold = drop_thresh;
        remove_struct.bubbles = false;

        // Remove droplets
        remove_droplets_region(remove_struct, ignore_region_x_limit, ignore_region_y_limit);

        // Remove bubbles
        remove_struct.bubbles = true;
        remove_droplets_region(remove_struct, ignore_region_x_limit, ignore_region_y_limit);

        // Completely removes bubble if specified
        if (REMOVE_ENTRAPMENT) {
            foreach(){ 
                if (x < 0.01 && y < 2 * 0.05) {
                    f[] = 1.;
                }
            }
        }
    }
}

event output_data (t += PLATE_OUTPUT_TIMESTEP) {
/* Outputs data about the flow*/
    if ((t >= START_OUTPUT_TIME) && (t <= END_OUTPUT_TIME)) {
        // Creates the file for outputting data along the plate
        char plate_output_filename[80];
        sprintf(plate_output_filename, "plate_output_%d.txt", plate_output_no);
        FILE *plate_output_file = fopen(plate_output_filename, "w");

        // Adds the time to the first line of the file
        fprintf(plate_output_file, "t = %g\n", t);

        // Loops over the left hand boundary (the interface)
        foreach_boundary(left) {
            /* Force calculation */
            // Viscosity average in the cell above the plate
            double avg_mu = f[] * (mu1 - mu2) + mu2;

            // Viscous stress in the cell above the plate
            double viscous_stress = - 2 * avg_mu * (u.x[1, 0] - u.x[]) / Delta;

            /* Plate output */
            fprintf(plate_output_file, "y = %g, x = %g, p = %g, strss = %g\n",\
                 y, x, p[], viscous_stress);
        }

        // Close plate output file
        fclose(plate_output_file);
        plate_output_no++; // Increments output number

        // /* Counts the number of droplets and bubbles there are */
        // scalar drop_field[];
        // scalar bubble_field[];
        // foreach() {
        //     drop_field[] = f[] > drop_thresh;
        //     bubble_field[] = 1. - f[] > drop_thresh;
        // }
        // int curr_drop_no = tag(drop_field);
        // int curr_bubble_no = tag(bubble_field);

        /* Outputs data to log file */
        fprintf(stderr, \
            "t = %.8f, v = %.8f, F = %.8f, force_term = %.8f, s = %g, ds_dt = %g, d2s_dt2 = %g\n", \
            t, 2 * pi * statsf(f).sum, current_force, force_term, s_current, ds_dt, d2s_dt2);
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


event movies (t += 0.001) {
/* Produces movies using bview */ 
    if (MOVIES) {
        // Creates a string with the time to put on the plots
        char time_str[80];
        sprintf(time_str, "t = %g\n", t);

	// RC Changed so that the view is zoomed in on the target region
	// Can leave general view as well, but for debugging this is more informative

        // Set up bview box
        view (width = 1024, height = 1024, fov = 5.0, ty = -0.1, \
            quat = {0, 0, -0.707, 0.707});

        /* Movie of the volume fraction of the droplet */
        clear();
        draw_vof("f", lw = 2);
        squares("f", linear = true, spread = -1, linear = true, map = cool_warm); // RC - minor changes here and beyond
        mirror ({0,1}) {
            draw_vof("f", lw = 2);
            squares("f", linear = true, spread = -1, linear = true, map = cool_warm);
        }
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("tracer.mp4");

        /* Movie of the vertical velocity */
        clear();
        draw_vof("f", lw = 2);
        squares("u.x", linear = false, spread = -1, linear = true, map = cool_warm);
        mirror ({0,1}) {
            draw_vof("f", lw = 2);
            squares("u.x", linear = false, spread = -1, linear = true, map = cool_warm);
        }
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("vertical_vel.mp4");


        /* Movie of the horizontal velocity */
        clear();
        draw_vof("f", lw = 2);
        squares("u.y", linear = false, spread = -1, linear = true, map = cool_warm);
        mirror ({0,1}) {
            draw_vof("f", lw = 2);
            squares("u.y", linear = false, spread = -1, linear = true, map = cool_warm);
        }
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("horizontal_vel.mp4");

        /* Movie of the pressure */
        clear();
        draw_vof("f", lw = 2);
        squares("p", linear = false, spread = -1, linear = true, map = cool_warm);
        mirror ({0,1}) {
            draw_vof("f", lw = 2);
            squares("p", linear = false, spread = -1, linear = true, map = cool_warm);
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

/* Alternative remove_droplets definitions */
void remove_droplets_region(struct RemoveDroplets p,\
        double ignore_region_x_limit, double ignore_region_y_limit) 
{
    scalar d[], f = p.f;
    double threshold = p.threshold ? p.threshold : 1e-4;
    foreach()
    d[] = (p.bubbles ? 1. - f[] : f[]) > threshold;
    int n = tag (d), size[n], keep_tags[n];

    for (int i = 0; i < n; i++) {
        size[i] = 0;
        keep_tags[i] = 1;
    }
    foreach_leaf() {
        if (d[] > 0) {
            int j = ((int) d[]) - 1;
            size[j]++;
            if ((x < ignore_region_x_limit) && (y < ignore_region_y_limit)) {
                keep_tags[j] = 0;
            }
        }
    }
    #if _MPI
    MPI_Allreduce (MPI_IN_PLACE, size, n, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    #endif
    int minsize = pow (p.minsize ? p.minsize : 3, dimension);
    foreach() {
        int j = ((int) d[]) - 1;
        if (d[] > 0 && size[j] < minsize && keep_tags[j] == 1)
            f[] = p.bubbles;
    }
    boundary ({f});
}

// RC

event logstats (t += 0.01) {

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}
