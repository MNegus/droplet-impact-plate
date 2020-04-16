/* force.c
    Basilisk script to test the force calculation. We set a uniform grid with a 
    prescribed pressure across the domain which only depends on the radial 
    coordinate. We calculate the force on the bottom of the domain by 
    approximating a surface integral in cylindrical coordinates using 
    trapazoidal rule.
*/

#include "navier-stokes/centered.h"
#include "axi.h"
#include "parameters.h"

double force(); // C function for calculating the force
double pressure(double rr, double tt); // Prescribed pressure 

int gfs_output_no = 1;
int plate_output_no = 1;

int main() {

    init_grid(1 << LEVEL);
    size(BOX_WIDTH);

    run();
}

event update_pressure(t += TIMESTEP) {
/* Updates the pressure across the domain */
    foreach() {
        p[] = pressure(y, t);
    }
}

// event gfs_output (t += TIMESTEP) {
//     char gfs_filename[80];
//     sprintf(gfs_filename, "gfs_output_%d.gfs", gfs_output_no);
//     output_gfs(file = gfs_filename);
//     gfs_output_no++;
// }

event output_force (t += TIMESTEP) {
    double analytic_force = \
        2 * pi * (exp(t) - (BOX_WIDTH + 1) * exp(t - BOX_WIDTH));
    fprintf(stderr, "t = %g, F = %g, F0 = %g\n", t, force(), analytic_force);
}

event end(t = MAX_TIME) {

}

double pressure(double rr, double tt) {
/* Prescribed pressure across the domain as a function of the radial coordinate,
rr and the time, tt. */
    return exp(-(rr - tt));
}

double force() {
/* Calculates the force by performing a surface integral about bottom of the 
domain, which is the left hand boundary. Recall the vertical coordinate (z) is
x and the radial coordinate (r) is y. */
    double force_value = 0.; // Initialises to zero

    /* Plate output file */
    char plate_output_filename[80];
    sprintf(plate_output_filename, "plate_output_%d.txt", plate_output_no);
    FILE *plate_output_file = fopen(plate_output_filename, "w");

    foreach_boundary(left) {
        fprintf(plate_output_file, "Delta = %g, p[1, 0] = %g, x = %g, y = %g\n", \
            Delta, p[1, 0], x, y);
        force_value += Delta * y * p[1, 0];
    }
    fclose(plate_output_file);
    plate_output_no++;

    force_value = 2 * pi * force_value; // Integrates angularly

    return force_value;
}