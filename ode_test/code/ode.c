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

int main() {
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
    // Solves for next s
    next_s = (Delta * Delta * force(t) - ALPHA * previous_s \
        - (Delta * Delta * GAMMA - Delta * BETA - 2 * ALPHA) * current_s) \
            / (ALPHA + Delta * BETA);
    
    // Redefines current_s and previous_s
    previous_s = current_s;
    current_s = next_s;
}

event end(t = MAX_TIME) {
    
}


double force(double tt) {
/* Force function */
    double a = 1;
    return a;
}