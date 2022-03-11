%% analytical_pressure_comparison


% Adds analytical pressures to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));




%% Parameters
close all;
eps = 1;
t_max = 0.1;
t0 = 1e-9;
tvals = linspace(t0, t_max + t0, 1e3);

% Analysis directory
analysis_directory ...
    = "/home/michael/Documents/supplementary_material/composite_pressure_comparison";

%% Creates figure and plots stationary
figure(3);
hold on;
p0_stationary = sqrt(3 ./ tvals) / (pi * eps);
plot(tvals, p0_stationary, 'Linewidth', 1);

%% alpha = 0.1, beta = 0, gamma = 0
% Plate parameters
alpha = 0.1;
gamma = 0;
beta = 0;

% Solve for s
[tvals, ss, sdots, sddots] = s_solution_alt(tvals, alpha, beta, gamma, eps);
[ds, ddots, dddots, Js] = s_dependents(tvals, ss, sdots, sddots);

% s depdendents
p0_moving = (4 * ds / (3 * pi)) .* (2 * ddots.^2 + ds .* dddots) / eps;

plot(tvals, p0_moving, 'Linewidth', 1);

%% alpha = 0.1, beta = 0, gamma = 400
% Plate parameters
alpha = 0.1;
gamma = 400;
beta = 0;

% Solve for s
[tvals, ss, sdots, sddots] = s_solution_alt(tvals, alpha, beta, gamma, eps);
[ds, ddots, dddots, Js] = s_dependents(tvals, ss, sdots, sddots);

% s depdendents
p0_moving = (4 * ds / (3 * pi)) .* (2 * ddots.^2 + ds .* dddots) / eps;

plot(tvals, p0_moving, 'Linewidth', 1);

%% alpha = 0.1, beta = 12.6491, gamma = 400
% Plate parameters
alpha = 0.1;
gamma = 400;
beta = 2 * sqrt(alpha * gamma);

% Solve for s
[tvals, ss, sdots, sddots] = s_solution_alt(tvals, alpha, beta, gamma, eps);
[ds, ddots, dddots, Js] = s_dependents(tvals, ss, sdots, sddots);

% s depdendents
p0_moving = (4 * ds / (3 * pi)) .* (2 * ddots.^2 + ds .* dddots) / eps;

plot(tvals, p0_moving, 'Linewidth', 1);

%% Saves figure
xlim([1e-3, 1.01 * max(tvals)]);
ylabel("p(0, 0, t)");
xlabel("t");
grid on;
title("Pressure at origin")
plot_name = sprintf("%s/origin_pressure_comparison.png", analysis_directory);
legend(["Stationary plate", "alpha = 0.1, beta = 0, gamma = 0", "alpha = 0.1, beta = 0, gamma = 400", "alpha = 0.1, beta = 12.6491, gamma = 400"]);
pause(0.1);
ax = gca;
exportgraphics(ax, plot_name, 'resolution', 300);
