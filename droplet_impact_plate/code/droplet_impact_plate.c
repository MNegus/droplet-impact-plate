/* droplet_impact_plate.c
    A Basilisk script to model the impact of a droplet of water impacting onto a
    moving plate. The domain is set to be in an accelerating frame with the 
    plate, so an additional body force is added. 
*/

// Filtering for large viscosity ratios
#define FILTERED

// Filtered viscosity field
#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))

#include "parameters.h" // Includes all defined parameters
#if AXISYMMETRIC
#include "axi.h" // Axisymmetric coordinates
#endif
#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h" // Implements two-phase flow
#include "view.h" // Creating movies using bview
#include "tension.h" // Surface tension of droplet
#include "tag.h" // For removing small droplets
#include "contact.h" // For imposing contact angle on the surface
#include <omp.h> // For openMP parallel

/* Physical constants */
double REYNOLDS; // Reynolds number of liquid
double WEBER; // Weber number of liquid
double FROUDE; // Froude number of liquid
double RHO_R; // Density ratio
double MU_R; // Viscosity ratio

/* Computational constants derived from parameters */
double MIN_CELL_SIZE; // Size of the smallest cell
double PLATE_REFINED_WIDTH; // Width of the refined area around the plate
double DROP_CENTRE; // Initial centre of the droplet
double IMPACT_TIME; // Theoretical time of impact
double MAX_TIME; // Maximum time to run the simulation for

/* Global variables */
double start_wall_time; // Time the simulation was started
double end_wall_time; // Time the simulation finished
int gfs_output_no = 0; // Records how many GFS files have been outputted
int plate_output_no = 0; // Records how many plate data files there have been
int interface_output_no = 0; // Records how many interface files there have been
double pinch_off_time = 0.; // Time pinch-off of the entrapped bubble occurs
double drop_thresh = 1e-4; // Remove droplets threshold
double bubble_area = 0.; // Area of entrapped bubble
char interface_time_filename[80] \
    = "interface_times.txt"; // Stores the time the interface was outputted

/* Force peak detection */
double * filtered_forces; // Filtered forces of previous timesteps
double current_force; // Current force on plate
double force_term; // Force term used in the ODE 
double avgFilter; // Average of force over the last PEAK_LAG timesteps
double stdFilter; // Standard deviation of force over last PEAK_LAG timesteps
int peak_no = 0; // Number of times we have done peak detection
double previous_avg = 0; // Value of avgFilter in previous timestep
double previous_std = 0; // Value of stdFilter in previous timestep

/* Plate position variables */
double s_previous = 0.; // Value of s at previous timestep
double s_current = 0.; // Value of s at current timestep
double s_next; // Values of s at next timestep
double ds_dt; // First time derivative of s
double d2s_dt2; // Second time derivative of s

/* Stats output */
FILE * fp_stats; 
char interp_stats_filename[80] = "interp_stats.txt";

/* Contact angle variables */ 
vector h[]; // Height function
double theta0 = 90; // Contact angle in degrees

/* Bubble counting */
scalar bubbles[]; // Tag field for bubbles

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

// Conditions on line of symmetry
#if !(AXISYMMETRIC)
u.n[bottom] = dirichlet(0.);
#endif

// Function for doing peak detection
void peak_detect(double current_force);

// Function for removing droplets away from a specific region
void remove_droplets_region(struct RemoveDroplets p,\
        double ignore_region_x_limit, double ignore_region_y_limit);


int main() {
/* Main function to set up the simulation */

    /* Create the computational domain */
    init_grid(1 << MINLEVEL); // Create grid according to the minimum level
    size(BOX_WIDTH); // Size of the domain

    /* Determine physical constants */
    REYNOLDS = RHO_L * V * R / MU_L; // Reynolds number of liquid
    WEBER = RHO_L * V * V * R / SIGMA; // Weber number of liquid
    FROUDE = V / sqrt(9.81 * R); // Froude number of liquid
    RHO_R = RHO_G / RHO_L; // Density ratio
    MU_R = MU_G / MU_L; // Viscosity ratio

    /* Set VOF constants */
    rho1 = 1.; // Density of water phase
    rho2 = RHO_R; // Density of air phase
    mu1 = 1. / REYNOLDS; // Viscosity of water phase
    mu2 = mu1 * MU_R; // Viscosity of air phase
    f.sigma = 1. / WEBER; // Surface tension at interface

    /* Derived constants */
    MIN_CELL_SIZE = BOX_WIDTH / pow(2, MAXLEVEL); // Size of the smallest cell
    // PLATE_REFINED_WIDTH \
    //     = PLATE_REFINE_NO * MIN_CELL_SIZE; // Refined region width around plate
    PLATE_REFINED_WIDTH = 0.01;
    DROP_CENTRE = INITIAL_DROP_HEIGHT + DROP_RADIUS;
    IMPACT_TIME = INITIAL_DROP_HEIGHT / (-DROP_VEL);

    /* Maximum time is shortly after Wagner theory would predict the turnover 
    point reaches the radius of the droplet */
    #if AXISYMMETRIC
    double wagner_max_time = 2.0 * (IMPACT_TIME + 1. / 3.);
    #else
    double wagner_max_time = 2.0 * (IMPACT_TIME + 1. / 4.);
    #endif
    // MAX_TIME = min(HARD_MAX_TIME, wagner_max_time);
    MAX_TIME = HARD_MAX_TIME;

    /* Allocates memory for the force and times arrays */
    if (PEAK_DETECT) {
        filtered_forces = malloc(PEAK_LAG * sizeof(double));
    }

    /* Initialises interface time file */
    FILE* interface_time_file = fopen(interface_time_filename, "w");
    fclose(interface_time_file);

    /* Initialises interp stats file */
    FILE * interp_stats_file = fopen(interp_stats_filename, "w");
    fclose(interp_stats_file);

    /* Open stats file */
    char name[200];
    sprintf(name, "logstats.dat");
    fp_stats = fopen(name, "w");

    /* Poisson solver constants */
    DT = 1.0e-4; // Minimum timestep
    NITERMIN = 1; // Min number of iterations (default 1)
    NITERMAX = 300; // Max number of iterations (default 100)
    TOLERANCE = 1e-5; // Possion solver tolerance (default 1e-3)

    // Run the simulation
    run();

    // Close stats file
    fclose(fp_stats);
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
    adapt_wavelet ({u.x, u.y, f}, (double[]){1e-3, 1e-3, 1e-6}, 
        minlevel = MINLEVEL, maxlevel = MAXLEVEL);
    
    /* Refines above the plate */
    refine((y < PLATE_WIDTH) && (x <= PLATE_REFINED_WIDTH) \
        && level < MAXLEVEL);
}


event moving_plate (t += 1e-4) {
/* Moves the plate as a function of the force on it */

    /* Calculate the force on the plate by integrating using trapezoidal rule */
    current_force = 0.; // Initialise to be zero

    // Iterates over the solid boundary
    foreach_boundary(left, reduction(+:current_force)) {
        if (y < PLATE_WIDTH) {
            // Viscosity average in the cell above the plate
            double avg_mu = f[] * (mu1 - mu2) + mu2;

            // Viscous stress in the cell above the plate
            double viscous_stress = \
                - 2 * avg_mu * (u.x[1, 0] - u.x[]) / Delta;

            // Adds the contribution to the force using trapeze rule, depending 
            // on if we are in the axisymmetric setting or not
            #if AXISYMMETRIC
            current_force += y * Delta * (p[] + viscous_stress);
            #else
            current_force += Delta * (p[] + viscous_stress);
            #endif
        }
    }

    // Integrates about angular part in axisymmetric, or doubles in 2D to take 
    // into account the other side of the plate

    #if AXISYMMETRIC
    current_force = 2 * pi * current_force; 
    #else
    current_force = 2 * current_force; 
    #endif
    
    // If the peak detect parameter is satisfied
    if (PEAK_DETECT) {
        peak_detect(current_force);
    } else {
        // Without peak detect, we set the force term to the current force
        force_term = current_force;
    }

    // If before force delay time, we set the force term to be zero
    if (t < FORCE_DELAY_TIME) force_term = 0;

    /* Solves the ODE for the updated plate position and acceleration using 
    a second-order explicit finite difference scheme */
    s_next = (DT * DT * force_term \
        + (2. * ALPHA - DT * DT * GAMMA) * s_current \
        - (ALPHA - DT * BETA / 2.) * s_previous) \
        / (ALPHA + DT * BETA / 2.);

    /* Updates values of s and its derivatives */
    ds_dt = (s_next - s_previous) / (2. * DT);
    d2s_dt2 = (s_next - 2 * s_current + s_previous) / (DT * DT);
    s_previous = s_current;
    s_current = s_next; 

    
    if (CONST_ACC) {
        /* If CONST_ACC is set, then we override this and set the variables to 
        their pre-defined values */
        if (t < IMPACT_TIME) {
            d2s_dt2 = 0.;
            ds_dt = 0.;
            s_current = 0.;
        } else {
            d2s_dt2 = PLATE_ACC; // PLATE_ACC is the pre-defined constant acc.
            ds_dt = d2s_dt2 * (t - IMPACT_TIME);
            s_current = 0.5 * d2s_dt2 * (t - IMPACT_TIME) * (t - IMPACT_TIME);
        }
    } else if (IMPOSED) {
        /* Else if IMPOSED is set, then we define the plate motion to be an 
        imposed sinusoidal curve */
        if (t < IMPACT_TIME) {
            d2s_dt2 = 0.;
            ds_dt = 0.;
            s_current = 0.;
        } else {
            double tShift = t - IMPACT_TIME;
            double k = 12.0;
            d2s_dt2 = IMPOSED_COEFF * (2 + sq(k) * cos(k * tShift));
            ds_dt = IMPOSED_COEFF * (2 * tShift + k * sin(k * tShift));
            s_current = IMPOSED_COEFF * (1 - cos(k * tShift) + sq(tShift));
        }
    }

    /* Updates velocity boundary conditions */
    u.t[top] = dirichlet(ds_dt);
    u.n[left] = y < PLATE_WIDTH ? dirichlet(0.) : dirichlet(ds_dt);

    boundary ((scalar *){u}); // Redefine boundary conditions for u
}


event acceleration (i++) {
/* Adds acceleration due to gravity and the moving plate at each time step */
    face vector av = a; // Acceleration at each face

    /* Adds acceleration due to gravity and the plate */
    foreach_face(x){
        av.x[] += d2s_dt2 - 1./sq(FROUDE);
    }
}


event small_droplet_removal (t += 1e-4) { 
/* Removes any small droplets or bubbles that have formed, that are smaller than
 a specific size. Uses the remove_droplets_region code to leave the area near 
 the point of impact alone in order to properly resolve the entrapped bubble */

    // Minimum diameter (in cells) a droplet/bubble has to be, else it will be 
    // removed
    int drop_min_cell_width = 36;
    int bubble_min_cell_width = 8;

    // Region to ignore
    double ignore_region_x_limit = 0.1; 
    double ignore_region_y_limit = 0.1; 
    
    // Counts the number of bubbles there are using the tag function
    foreach() {
        bubbles[] = 1. - f[] > drop_thresh;
    }
    int bubble_no = tag(bubbles);

    // Determines if we are before or after the pinch-off time
    if (pinch_off_time == 0.) {
        /* The first time the bubble number is above 1, we define it to be the 
        pinch off time */
        if (bubble_no > 1) {
            pinch_off_time = t;
        }
    } else {
        /* Determine area of entrapped bubble */

        // We assume that entrapped bubbles are any bubble which is not the #
        // surrounding air, so we determine the tag of the surrounding air, and 
        // then add up the volume of all of the cells which do not have that 
        // tag

        // We initialise the tag to be 0, which is the tag the
        // bulk droplet will have, so that if the bubble is not found, the 
        // resulting area will be huge and easy to identify that we're wrong.
        int air_tag = 0; 

        // Serial foreach loop along the right boundary to find the tag of the 
        // air, which in theory should end after one iteration
        foreach_boundary(right, serial) {
            if (f[] == 0.) {
                air_tag = bubbles[];
                break;
            }
        }

        // Determine the area of all of the entrapped air cells, which will have
        // a tag not equal to 0 (which will be liquid) or air_tag, which is the
        // tag of the surrounding air
        bubble_area = 0.;
        foreach(reduction(+:bubble_area)) {
            if ((bubbles[] > 0) && (bubbles[] != air_tag)) {
                bubble_area += (1. - f[]) * dv();
            }
        }

        /* After the removal delay, remove drops and bubbles as necessary */    
        if (t >= pinch_off_time + REMOVAL_DELAY) {
            // Set up RemoveDroplets struct
            struct RemoveDroplets remove_struct;
            remove_struct.f = f;
            remove_struct.minsize = drop_min_cell_width;
            remove_struct.threshold = drop_thresh;
            remove_struct.bubbles = false;

            if (t < 0.3) {
                // Remove droplets outside of the specified region
                remove_droplets_region(remove_struct, ignore_region_x_limit, \
                    ignore_region_y_limit);

                // Remove bubbles outside of the specified region
                remove_struct.bubbles = true;
                remove_struct.minsize = bubble_min_cell_width;
                remove_droplets_region(remove_struct, ignore_region_x_limit, \
                    ignore_region_y_limit);
            } else {
                remove_droplets_region(remove_struct, 0, 0);
                
                remove_struct.bubbles = true;
                remove_struct.minsize = bubble_min_cell_width;
                remove_droplets_region(remove_struct, 0, 0);
            }

            // Remove the entrapped bubble if specified
            if (REMOVE_ENTRAPMENT) {
                foreach() { 
                    if (x < 0.01 && y < 2 * 0.05) {
                        f[] = 1.;
                    }
                }
            }
        }
    }
        
}


event output_plate (t += PLATE_OUTPUT_TIMESTEP) {
/* Outputs data along the plate */

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
    }
}


event output_log (t += LOG_OUTPUT_TIMESTEP) {
/* Outputs data about the general flow */
    if ((t >= START_OUTPUT_TIME) && (t <= END_OUTPUT_TIME)) {
        /* Outputs data to log file */
        fprintf(stderr, \
            "t = %.4f, F = %.6f, force_term = %.6f, avg = %.6f, std = %.6f, s = %g, ds_dt = %g, d2s_dt2 = %g, bubble_area = %.7f\n", \
            t, current_force, force_term, previous_avg, previous_std, \
            s_current, ds_dt, d2s_dt2, bubble_area);
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
        // Output gfs file
        char gfs_filename[80];
        sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
        output_gfs(file = gfs_filename);

        // Output fields
        char field_filename[80];
        sprintf(field_filename, "field_output_%d.txt", gfs_output_no);
        FILE *field_file = fopen(field_filename, "w");

        int N_output = (int) floor(pow(2, MAXLEVEL) * 2. / 6.);
        output_field ({p,f,u}, field_file, N_output, box = {{0,0},{2.5,2.5}});

        fclose(field_file);

        gfs_output_no++;
    }
}


event movies (t += 1e-3) {
/* Produces movies using bview */ 
    if (MOVIES) {
        // Creates a string with the time to put on the plots
        char time_str[80];
        sprintf(time_str, "t = %g\n", t);

        /* Zoomed out view */
        // Set up bview box
        view (width = 1024, height = 1024, fov = 9, ty = -0.235, tx = -0.235, \
            quat = {0, 0, -0.707, 0.707});

        /* Movie of the volume fraction of the droplet */
        clear();
        mirror({0, 1}) {
            draw_vof("f", lw = 2);
            squares("f", linear = true, spread = -1, linear = true, map = cool_warm); 
        }
        draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);
        save ("tracer.mp4");

        /* Pressure video, scaled by the stationary Wagner maximum */
        for (int pcoeff = 0; pcoeff <= 4; pcoeff++) {
            clear();
            mirror({0, 1}) {
                draw_vof("f", lw = 2);
                if (t <= IMPACT_TIME) {
                    squares("p", min = 0, linear = false, spread = -1, linear = true, map = cool_warm);
                } else {
                    double wagner_pmax = 3 / (8 * (t - IMPACT_TIME));
                    squares("p", min = 0, max = wagner_pmax * (1 + 0.25 * pcoeff), linear = false, spread = -1, linear = true, map = cool_warm);
                }
            }
            char pressure_vid_filename[80];
            sprintf(pressure_vid_filename, "pressure_%d.mp4", pcoeff);
            save (pressure_vid_filename);
        }

        /* Velocity videos. Aim is for each velocity component and norm, produce multiple videos with
        different (fixed) colour maps */
        for (int velmax = 1; velmax <= 3; velmax++) {
            /* Movie of vertical velocity */
            int velmax = 2;
            clear();
            mirror({0, 1}) {
                draw_vof("f", lw = 2);
                squares("u.x", min = -velmax, max = velmax, linear = false, spread = -1, linear = true, map = cool_warm);
            }
            // draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);

            char vertical_vid_filename[80];
            sprintf(vertical_vid_filename, "vertical_vel_%d.mp4", velmax);
            save (vertical_vid_filename);

            /* Movie of horizontal velocity */
            velmax = 3;
            clear();
            mirror({0, 1}) {
                draw_vof("f", lw = 2);
                squares("u.y", min = -velmax, max = velmax, linear = false, spread = -1, linear = true, map = cool_warm);
            }
            // draw_string(time_str, pos=1, lc= { 0, 0, 0 }, lw=2);

            char horizontal_vid_filename[80];
            sprintf(horizontal_vid_filename, "horizontal_vel_%d.mp4", velmax);
            save (horizontal_vid_filename);
        }
    }
}


event logstats (t += 0.01) {
/* Event to regularly output relevant statistics */

    timing s = timer_timing (perf.gt, i, perf.tnc, NULL);
 
    // i, timestep, no of cells, real time elapsed, cpu time
    fprintf(fp_stats, "i: %i t: %g dt: %g #Cells: %ld Wall clock time (s): %g CPU time (s): %g \n", \
        i, t, dt, grid->n, perf.t, s.cpu);
    fflush(fp_stats);
}


event end (t = MAX_TIME) {
/* Ends the simulation */ 

    end_wall_time = omp_get_wtime(); // Records the time of finish

    fprintf(stderr, "Finished after %g seconds\n", \
        end_wall_time - start_wall_time);

    if (PEAK_DETECT) {
        free(filtered_forces);
    }
}

/* Peak detect algorithm */
void peak_detect(double current_force) {
    /* Peak detection. We attempt to use a peak detection algorithm to check if 
    the current force value is as expected, or if it has peaked to a 
    non-desirable value. 
    Details of the algorithm are found at:
    https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data 
    */
   
    /* For the first PEAK_LAG timesteps, we populate the filtered force
        array with the current values of the force. The forcing term is set
        to the current force */ 
    if (peak_no < PEAK_LAG) {

        // Populate array
        filtered_forces[peak_no] = current_force;

        // Specify forcing term for ODE
        force_term = current_force;

        peak_no++;

    } else {
        double new_filtered; // New filtered value of force

        /* Before the specified delay, no peak detection happens. This is to
        ensure the problem has regularised */
        if (t < PEAK_DELAY)  {
            force_term = current_force;
            new_filtered = force_term;
        /* Else, we initiate the peak detection algorithm */
        } else {
            double previous_force = filtered_forces[PEAK_LAG - 1];
            double diff_from_previous \
                = (current_force - previous_force) / previous_force;

            if ((current_force <= 0) || fabs(diff_from_previous) > 0.25) {
                /* If current force is less than or equal to zero, or 
                differs from the previous force by more than 25%, then 
                completely ignore */
                force_term = filtered_forces[PEAK_LAG - 1];
                new_filtered = force_term;

                // Output the force data
                FILE * interp_stats_file = fopen(interp_stats_filename, "a");
                fprintf(interp_stats_file, "t = %g, F = %g, avgFilter = %g, stdFilter = %g, force_term = %g\n", \
                    t, current_force, avgFilter, stdFilter, force_term);
                fclose(interp_stats_file);
            } else if (fabs(current_force - avgFilter) > PEAK_THRESHOLD * stdFilter) {
                /* If current force deviates from the mean more than 
                PEAK_THRESHOLD number of standard deviations, then take 
                force term to be an influenced value */

                force_term = avgFilter;
                    
                new_filtered = PEAK_INFLUENCE * current_force \
                    + (1 - PEAK_INFLUENCE) * filtered_forces[PEAK_LAG - 1];

                // Output the force data
                FILE * interp_stats_file = fopen(interp_stats_filename, "a");
                fprintf(interp_stats_file, "t = %g, F = %g, avgFilter = %g, stdFilter = %g, force_term = %g\n", \
                    t, current_force, avgFilter, stdFilter, force_term);
                fclose(interp_stats_file);
            } else {
                /* Else current_force is kept */
                force_term = current_force;
                new_filtered = force_term;
            }
        }

        // Re-populate filtered forces array
        #pragma omp critical
        for (int j = 0; j < PEAK_LAG - 1; j++) {
            filtered_forces[j] = filtered_forces[j + 1];
        }
        filtered_forces[PEAK_LAG - 1] = new_filtered;
    }

    // Update average and standard deviation of filtered forces
    if (peak_no >= PEAK_LAG - 1) {
        previous_avg = avgFilter;
        previous_std = stdFilter;

        // Average
        avgFilter = 0;
        #pragma omp for reduction(+:avgFilter)
        for (int j = 0; j < PEAK_LAG; j++) {
            avgFilter = avgFilter + filtered_forces[j];
        }
        avgFilter = avgFilter / ((double) PEAK_LAG);

        // Standard deviation
        stdFilter = 0;
        #pragma omp for reduction(+:stdFilter)
        for (int j = 0; j < PEAK_LAG; j++) {
            stdFilter += (filtered_forces[j] - avgFilter) \
                * (filtered_forces[j] - avgFilter);
        }
        stdFilter = sqrt(stdFilter / ((double) PEAK_LAG));
    }
}


/* Alternative remove_droplets definition */
void remove_droplets_region(struct RemoveDroplets p,\
        double ignore_region_x_limit, double ignore_region_y_limit) {
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
