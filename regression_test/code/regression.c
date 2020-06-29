#include "parameters.h" // Includes all defined parameters
#include "axi.h" // Axisymmetric coordinates
#include "navier-stokes/centered.h" // To solve the Navier-Stokes
#include "two-phase.h" // Implements two-phase flow
// #include <omp.h> // For openMP parallel
#include <gsl/gsl_fit.h>


double * forces_array;
double * times_array;

double gradient = 0.5;
double intercept = 4.9;

int main () {

    forces_array = malloc(INTERP_NO * sizeof(double));
    times_array = malloc(INTERP_NO * sizeof(double));

    DT = 1e-2;

    run();

    free(forces_array);
    free(times_array);
}

event force(i++) {
    if (i < INTERP_NO) {
        forces_array[i] = gradient * t + intercept;
        times_array[i] = t;
    } else {
        double c0, c1, cov00, cov01, cov11, sumsq;
        gsl_fit_linear ( times_array, 1, forces_array, 1, INTERP_NO, &c0, &c1, \
            &cov00, &cov01, &cov11, &sumsq);

        fprintf(stderr, "t = %g, c0 = %g, c1 = %g\n", t, c0, c1);
        
        double current_force = gradient * t + intercept;
        double interp_force = c0 + c1 * t;
        fprintf(stderr, "f = %g, f_interp = %g", current_force, interp_force);
        fprintf(stderr, "\n");

        // #pragma omp critical
        for (int j = 0; j < INTERP_NO - 1; j++) {
            forces_array[j] = forces_array[j + 1];
            times_array[j] = times_array[j + 1];
        }
        forces_array[INTERP_NO - 1] = current_force;
        times_array[INTERP_NO - 1] = t;
    }
}

event end(t = 1) {
    fprintf(stderr, "Ended\n");
}