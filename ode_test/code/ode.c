/* ode.c
    Script to test the ODE solving, within a Basilisk framework.
    Solving the ODE
    ALPHA s''(t) + BETA s'(t) + GAMMA s(t) = F(t)
    where ALPHA, BETA, GAMMA are defined in parameters.h and F(t) is defined in
    the function force
*/

#include "parameters.h"
#include "navier-stokes/centered.h"

// Function declaration for force
double force(double tt);

// Global variables for s values at each timestep
double previous_s;
double current_s;
double next_s;

// double ALPHA, BETA, GAMMA;

int main() {

    // // Sinusoidal forcing parameters
    // ALPHA = 1.;
    // BETA = nu;
    // GAMMA = omega0 * omega0;
    run();
}

event init(t = 0) {
    // Initial conditions for s
    previous_s = 0.;
    current_s = 0;
}

event output (t += Delta) {
    // Outputs current time and current s into the log file
    fprintf(stderr, "%g, %g\n", t, current_s);
}

event ode_solving (t += Delta) {
    // Solves for next s, central difference for beta term
    next_s = (Delta * Delta * force(t) \
        + (2. * ALPHA - Delta * Delta * GAMMA) * current_s \
        - (ALPHA - Delta * BETA / 2.) * previous_s) \
        / (ALPHA + Delta * BETA / 2.);


    // next_s = (Delta * Delta * force(t) - ALPHA * previous_s \
    //     - (Delta * Delta * GAMMA - Delta * BETA - 2 * ALPHA) * current_s) \
    //         / (ALPHA + Delta * BETA);
    
    // Redefines current_s and previous_s
    previous_s = current_s;
    current_s = next_s;
}

event end(t = MAX_TIME) {
    
}


double force(double tt) {
/* Force function */
    double omega = 1.2;
    // return sin(omega * t);
    // return omega0 * omega0 * F0 * sin(omega * tt);
    // return tt * tt;
    // return exp(sin(omega * tt));
    return exp(sin(omega * pow(tt, 1.5)));
}