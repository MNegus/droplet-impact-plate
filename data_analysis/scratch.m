addpath("forces");

% Parameters
alpha = 2;
beta = 0;
gamma = 0;
tmax = 0.2;
t0 = 1e-9

%% epsilon = 1
eps = 1;
tvals = linspace(t0, tmax + t0, 1e3);
% Solves for s
[tvals, s, sdot, sddot] = s_solution_alt(tvals, alpha, beta, gamma, eps)

% Determine forces
outer_F = outer_force(tvals, s, sdot, sddot, eps);
inner_F = inner_force(tvals, s, sdot, sddot, eps);
overlap_F = overlap_force(tvals, s, sdot, sddot, eps);
composite_F = composite_force(tvals, s, sdot, sddot, eps);

%% epsilon = 0.1
% eps = 0.1;
% scaled_tvals = linspace(t0, tmax + t0, 1e3) / eps^2;
% scaled_alpha = alpha / eps^2;
% scaled_beta = beta;
% scaled_gamma = eps^2 * gamma;
% % Solves for s
% [scaled_tvals, s_small, sdot_small, sddot_small] = s_solution_alt(scaled_tvals, scaled_alpha, scaled_beta, scaled_gamma, eps);
% t_small = eps^2 * t;
% 
% % Determine forces
% outer_F_small = outer_force(scaled_tvals, s_small, sdot_small, sddot_small, eps);
% inner_F_small = inner_force(scaled_tvals, s_small, sdot_small, sddot_small, eps);
% overlap_F_small = overlap_force(scaled_tvals, s_small, sdot_small, sddot_small, eps);
% composite_F_small = composite_force(scaled_tvals, s_small, sdot_small, sddot_small, eps);

%% Plots comparison
close all;
figure(1);
hold on;
plot(tvals, s);
% plot(eps^2 * scaled_tvals, eps^2 * s_small);
% plot(tvals, outer_F);
% plot(scaled_tvals, outer_F_small);


figure(2);

