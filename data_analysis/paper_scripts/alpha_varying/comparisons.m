%% force_comparison.m
% Script to compare the forces on the plate for where the initial height
% has been varied. The times are shifted so the theoretical time of impact
% is always the same
%

% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions


% Parent directory where all the data is held
parent_directory = '/mnt/newarre/cantilever_paper_data/alpha_varying/';

% Directory to save the figure(s)
analysis_directory = "Analysis";
analysis_directory = strcat(parent_directory, analysis_directory);

% Range of alphas
alphas = [1, 2, 5, 10, 20, 50, 100, 1000];

% Defines arrays for all the values of alpha
data_directories = string(length(alphas));
legend_entries = string(length(alphas));

for k = 1 : length(alphas)
    alpha = alphas(k);
    
    data_directories(k) = [parent_directory, '/alpha_', num2str(alpha)];
    
    legend_entries(k) = ['$\alpha =$ ', num2str(alpha)] ;
end

% Stationary plate data
stationary_plate_directory = '/mnt/newarre/cantilever_paper_data/stationary_plate';


%% Parameters

% Physical parameters
rho_w = 998;
R0 = 1e-3;
U0 = 5;
T0 = R0 / U0;
Patm = 10^5;

% Dimensional functions
F_newton = @(F) rho_w * U0^2 * R0^2 * F;
t_milli = @(t) T0 * 1000 * t;
r_milli = @(r) R0 * 1000 * r;

% Value of epsilon
eps = 1;

% Plate parameters
alpha = 1000;
beta = 0;
gamma = 0;

% Initial drop height 
initial_drop_heights = 0.125;

% Impact time 
impact_time = initial_drop_heights;

% Maximum time
t_max = 0.5;

%% Wagner solution
% Stationary plate Wagner solution
wagner_t = linspace(1e-7, t_max - impact_time, 1e3);
s = zeros(size(wagner_t));
sdot = zeros(size(wagner_t));
sddot = zeros(size(wagner_t));

% Composite force solution
wagner_force = composite_force(wagner_t, s, sdot, sddot, eps);

%% Plotting
close all;

% Force comparison
figure(1);
hold on;

% s displacement comparison
figure(2);
hold on;

% sdot comparison
figure(3);
hold on;

% sddot comparison
figure(4);
hold on;

for k = 1 : length(data_directories)
   output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(k)));
   
   ts = output_mat(:, 1);
   Fs = output_mat(:, 3);
   ss = output_mat(:, 6);
   sdots = output_mat(:, 7);
   sddots = output_mat(:, 8);
   
   figure(1);
   plot(t_milli(ts), F_newton(Fs));
   
   figure(2);
   plot(t_milli(ts), r_milli(ss));
   
   figure(3);
   plot(t_milli(ts), sdots * U0);
   
   figure(4);
   plot(t_milli(ts), sddots * U0 / T0);
end

% Force comparison
figure(1);
output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", ...
    stationary_plate_directory));
ts = output_mat(:, 1);
Fs = output_mat(:, 3);
plot(t_milli(ts), F_newton(Fs));
plot(t_milli(wagner_t + impact_time), F_newton(wagner_force), ...
    'color', 0.5 * [1 1 1], 'Linewidth', 2, 'Linestyle', '--');
legend([legend_entries, "Stationary plate", "Analytical"], ...
    "Interpreter", "latex", "location", "northwest", "Fontsize", 12);


xlim([t_milli(-0.1) t_milli(t_max)]);
grid on;
xlabel("$t^*$ / ms", "Interpreter", "latex");
ylabel("$F^*(t^*)$ / N ", 'Interpreter', 'latex');

ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title(['Force on plate: $\beta =$ ', ...
    num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
    'Fontsize', 14);
print(gcf, sprintf("%s/force_comparison.png", analysis_directory), ...
    '-dpng', '-r300');

% s comparisons
figure(2);
% plot(t_milli(wagner_t + impact_time), r_milli(s), ...
%     'color', 0.5 * [1 1 1], 'Linewidth', 2, 'Linestyle', '--');
legend([legend_entries], ...
    "Interpreter", "latex", "location", "northwest", "Fontsize", 12);


xlim([0 t_milli(t_max)]);
grid on;
xlabel("$t^*$ / ms", "Interpreter", "latex");
ylabel("$s^*(t^*)$ / mm", 'Interpreter', 'latex');

ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title(['Plate displacement: $\beta =$ ', ...
    num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
    'Fontsize', 14);
print(gcf, sprintf("%s/plate_displacement_comparison.png", analysis_directory), ...
    '-dpng', '-r300');


% sdot comparisons
figure(3);
% plot(t_milli(wagner_t + impact_time), sdot * U0, ...
%     'color', 0.5 * [1 1 1], 'Linewidth', 2, 'Linestyle', '--');
legend([legend_entries], ...
    "Interpreter", "latex", "location", "northwest", "Fontsize", 12);


xlim([0 t_milli(t_max)]);
grid on;
xlabel("$t^*$ / ms", "Interpreter", "latex");
ylabel("$\dot{s}^*(t^*)$ / m/s", 'Interpreter', 'latex');

ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title(['Plate velocity: $\beta =$ ', ...
    num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
    'Fontsize', 14);
print(gcf, sprintf("%s/plate_velocity_comparison.png", analysis_directory), ...
    '-dpng', '-r300');

% sddot comparisons
figure(4);
% plot(t_milli(wagner_t + impact_time), sddot * U0 / T0, ...
%     'color', 0.5 * [1 1 1], 'Linewidth', 2, 'Linestyle', '--');
legend([legend_entries], ...
    "Interpreter", "latex", "location", "northwest", "Fontsize", 12);


xlim([t_milli(-0.05) t_milli(t_max)]);
grid on;
xlabel("$t^*$ / ms", "Interpreter", "latex");
ylabel("$\ddot{s}^*(t^*)$ / m/s$^2$", 'Interpreter', 'latex');

ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title(['Plate acceleration: $\beta =$ ', ...
    num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
    'Fontsize', 14);
print(gcf, sprintf("%s/plate_acceleration_comparison.png", analysis_directory), ...
    '-dpng', '-r300');