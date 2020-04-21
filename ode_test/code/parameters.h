/* parameters.h
Header file for the ODE test to feed in relevant parameters */

/* Computational parameters */
const double dt0 = 1e-1; // First timestep
const double MAX_TIME = 5.; // Maximum time to run the ODE for

/* Oscillation parameters:
Parameters such that the ODE is:
s''(t) + 2 nu s'(t) + omega0^2 s(t) = omega0^2 q,
hence in terms of the solver, ALPHA = 1, BETA = 2 * nu, GAMMA = omega0^2 and 
F(t) = omega0^2 * q */
const double nu = 1.; // Damping term
const double omega0 = 10.; // Elastic term
const double q = 1.; // Forcing term
