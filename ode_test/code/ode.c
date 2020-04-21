/* ode.c
    Script to test the ODE solving, within a Basilisk framework.

    Solving the ODE:
        ALPHA s''(t) + BETA s'(t) + GAMMA s(t) = F(t)
    where ALPHA, BETA, GAMMA are defined in parameters.h and F(t) is defined in
    the function force, with initial conditions s(0) = s'(0) = 0.
*/

#include "parameters.h"
#include "navier-stokes/centered.h"

/* Global variables */
double ALPHA, BETA, GAMMA; // Coefficients in front of terms in ODE 
double omega; // Oscillation period used in analytic solution
double delta_t = 0.1 / pow(2, RUN_NO); // Timestep value

// Global variables for s values at each timestep
double previous_s;
double current_s;
double next_s;

// Output file name
FILE * output_file;

/* Function declarations */
double force(double tt); // Forcing term on RHS
double s_analytic(double tt); // Analytic solution

int main() {

    /* Set values of ALPHA, BETA and GAMMA based on the parameters given */
    ALPHA = 1.;
    BETA = 2 * nu;
    GAMMA = omega0 * omega0;

    /* Defines omega in terms of omega0 and nu */
    omega = sqrt(omega0 * omega0 - nu * nu);

    fprintf(stderr, "delta_t = %g\n", delta_t);
    run();
}

event init(t = 0) {
/* Initialises the problem */

    /* Discrete form of the initial condition s(0) = s'(0) = 0 */
    previous_s = 0.;
    current_s = 0;

    // Opens the output file
    char output_filename[80];
    sprintf(output_filename, "output_%d.txt", RUN_NO);
    output_file = fopen(output_filename, "w");
}


event output (t += delta_t) {
    // Outputs current time, current s and analytic s into the output file
    fprintf(output_file, "%.8f, %.8f, %.8f\n", t, current_s, s_analytic(t));
}

event ode_solving (t += delta_t) {
    // Solves for next s, central difference for beta term
    next_s = (delta_t * delta_t * force(t) \
        + (2. * ALPHA - delta_t * delta_t * GAMMA) * current_s \
        - (ALPHA - delta_t * BETA / 2.) * previous_s) \
        / (ALPHA + delta_t * BETA / 2.);
    
    // Redefines current_s and previous_s
    previous_s = current_s;
    current_s = next_s;
}

event end(t = MAX_TIME) {
/* Event when max time is reached */
    // Closes the output file
    fclose(output_file);
}

double force(double tt) {
/* Force function */
    return omega0 * omega0 * q;
}

double s_analytic(double tt) {
/* Analytic solution to the equation
s''(t) + 2 nu s'(t) + omega0^2 s(t) = omega0^2 * q
*/
    return \
    q * (1 - exp(-nu * tt) * (cos(omega * t) + (nu / omega) * sin(omega * t)));
}