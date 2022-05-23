%% turnover_and_jet_comparison.m
% Compares the turnover point and jet thickness of a stationary plate run
% and a moving plate run

close all;

%% Data definitions
% Here we specify the locations where the plate output files are stored. We
% expect a different directory for each simulation result.

% Parent directories where all of the data is stored under (e.g. external
% hard drive location)
data_directory = "/home/michael/Documents/supplementary_material/turnover_and_jet";
% stationary_directory = "/media/michael/newarre/cantilever_paper_data/stationary_plate";
% moving_directory = "/media/michael/newarre/cantilever_paper_data/gamma_varying/gamma_500";

% Directory where the resulting videos are to be stored
results_directory = "/home/michael/Documents/supplementary_material/turnover_and_jet";

%% Parameters

% Plate parameters
alpha = 2;
beta = 0;
gamma = 500;
eps = 1;

impact_time = 0.125;
output_range = 125:575; % Range of outputs
timestep = 1e-3; % (Constant) timestep
tmax = timestep * max(output_range) - impact_time;


full_tvals = timestep * output_range - impact_time;

%% Analytical solutions
% Analytical time variables
analytical_tvals = timestep : timestep : tmax;

% Solves for the plate displacement
[tvals, s, sdot, sddot] = s_solution_alt(analytical_tvals, alpha, beta, gamma, eps);

% Moving plate solutions
[d_mov, ddot_mov, dddot_mov, J_mov] = s_dependents(tvals, s, sdot, sddot);

% Stationary plate solutions
[d_stat, ddot_stat, dddot_stat, J_stat] = s_dependents(tvals, 0, 0, 0);


%% Numerical solutions
stationary_file = dlmread(sprintf("%s/stationary_turnover_points.txt", ...
    data_directory));
moving_file = dlmread(sprintf("%s/moving_turnover_points.txt", ...
    data_directory));


% Stationary plate solutions
num_d_stat = stationary_file(output_range, 2);
num_ddot_stat = diff(num_d_stat) / timestep;
num_dddot_stat = diff(num_ddot_stat) / timestep;
num_J_stat = stationary_file(output_range, 3);

% Moving plate solutions
num_d_mov = moving_file(output_range, 2);
num_ddot_mov = diff(num_d_mov) / timestep;
num_dddot_mov = diff(num_ddot_mov) / timestep;
num_J_mov = moving_file(output_range, 3);

% %% Alternative
% gamma = 10;
% % Solves for the plate displacement
% [tvals, s, sdot, sddot] = s_solution_alt(analytical_tvals, alpha, beta, gamma, eps);
% [d_gamma_10, ddot_gamma_10, dddot_gamma_10, J_gamma_10] = s_dependents(tvals, s, sdot, sddot);
% 
% moving_file = dlmread(sprintf("%s/gamma_10_turnover_points.txt", ...
%     data_directory));
% num_d_mov = moving_file(output_range, 2);
% num_ddot_mov = diff(num_d_mov) / timestep;
% num_dddot_mov = diff(num_ddot_mov) / timestep;
% num_J_mov = moving_file(:, 3);


%% Plots
close all;

%% Turnover point position
close(figure(1));
figure(1);
hold on;
% scatter(full_tvals, num_d_stat, [], [0 0 0], 'o');
% scatter(full_tvals, num_d_mov, [], 0.5 * [1 1 1], 'o');
plot(full_tvals, num_d_stat, 'color', [0 0 0]);
plot(full_tvals, num_d_mov, 'color', 0.5 * [1 1 1]);
plot(tvals, d_stat, 'Linestyle', '--', 'color', [0 0 0], 'Linewidth', 2);
plot(tvals, d_mov, 'Linestyle', '--', 'color', 0.5 * [1 1 1], 'Linewidth', 2);
grid on;
xlabel("$t$", "Interpreter", "latex", 'Fontsize',30);
ylabel("$d(t)$", "Interpreter", "latex", 'Fontsize', 30);
ax = gca;
ax.FontSize = 16;
set(gca,'TickLabelInterpreter','latex');

L = legend(["$d_s(t)$, numerical", "$d_m(t)$, numerical", ...
    "$d_s(t)$, analytical", "$d_m(t)$, analytical"]);
set(L, 'Interpreter', 'latex');
set(L, 'FontSize', 15);
set(L, 'location', 'southeast');
plot_name = sprintf("%s/turnover_point_full.png", results_directory);
exportgraphics(ax, plot_name, 'resolution', 300);

%% Turnover point comparison
close(figure(2));
figure(2);
hold on;
analytical_values = (d_stat - d_mov) ./ d_stat;
scatter(full_tvals, (num_d_stat - num_d_mov) ./ num_d_stat, [], 0.5 * [1 1 1], 'o');
plot(tvals, analytical_values, 'Linestyle', '--', 'color', [0 0 0], 'Linewidth', 2);
% ylim(1.2 * [0, max(analytical_values)]);
ylim([0, 0.025]);
grid on;
xlabel("$t$", "Interpreter", "latex", 'Fontsize',30);
ylabel("$\overline{d_s(t) - d_m(t)}$", "Interpreter", "latex", 'Fontsize', 30);
ax = gca;
ax.FontSize = 16;
set(gca,'TickLabelInterpreter','latex');

L = legend(["Numerical", "Analytical"]);
set(L, 'Interpreter', 'latex');
set(L, 'FontSize', 15);
plot_name = sprintf("%s/turnover_point_comparison.png", results_directory);
exportgraphics(ax, plot_name, 'resolution', 300);


%% Turnover point velocity comparison
figure(3);
hold on;
analytical_values = (ddot_stat - ddot_mov) ./ ddot_stat;
plot(tvals, analytical_values, 'Linestyle', '--', 'color', [0 0 0], 'Linewidth', 2);
% scatter(full_tvals(1 : end - 1), (num_ddot_stat - num_ddot_mov) ./ num_ddot_stat, [], 0.5 * [1 1 1], 'o');
grid on;
ylim([-0.05, 0.05]);
yticks([-0.05, -0.025, 0, 0.025, 0.05]);
xlabel("$t$", "Interpreter", "latex", 'Fontsize',30);
ylabel("$\overline{v_s(t) - v_m(t)}$", "Interpreter", "latex", 'Fontsize', 30);
ax = gca;
ax.FontSize = 16;
set(gca,'TickLabelInterpreter','latex');
L = legend(["Analytical"]);
set(L, 'Interpreter', 'latex');
set(L, 'FontSize', 15);
plot_name = sprintf("%s/turnover_velocity_comparison.png", results_directory);
exportgraphics(ax, plot_name, 'resolution', 300);

%% Jet thickness comparison
% figure(4);
% hold on;
% % analytical_values = (J_stat - J_mov) ./ J_stat;
% analytical_values = (J_stat - J_mov);
% plot(tvals, analytical_values);
% 
% % scatter(full_tvals, (num_J_stat - num_J_mov) ./ num_J_stat, [], 0.5 * [1 1 1], 'o');
% 
% scatter(full_tvals, (num_J_stat - num_J_mov), [], 0.5 * [1 1 1], 'o');
% ylim(1.2 * [0, max(analytical_values)]);


%% Jet thickness
% figure(5);
% hold on;
% plot(tvals, J_stat);
% scatter(full_tvals, num_J_stat);
% plot(tvals, J_mov);
% scatter(full_tvals, num_J_mov);
% legend("Stationary, analytical", "Stationary, numerical", "Moving, analytical", ...
%     "Moving, numerical");

