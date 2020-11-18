%% force_comparison.m
% Script to compare the displacements s. 
% The times are shifted so the theoretical time of impact
% is always the same
%

% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions


% Parent directory where all the data is held
parent_directory = '/mnt/newarre/cantilever_paper_data/';

% Directory to save the figure(s)
analysis_directory = "special_case/Analysis";
analysis_directory = strcat(parent_directory, analysis_directory);

% Defines array with both directories in
data_directories = ["gamma_varying/gamma_500"];

% Concatenates arrays to include parent directory
for k = 1 : length(data_directories)
    data_directories(k) = strcat(parent_directory, data_directories(k)); 
end

%% Parameters

% Value of epsilon
eps = 1;

% Plate parameters
alpha = 2;
beta = 0;
gamma = 500;

% Initial drop height 
initial_drop_heights = 0.125;

% Impact time 
impact_time = initial_drop_heights;

% Maximum computational time
t_max = 0.8;

%% Wagner solution

% Plate displacement solution
[wagner_t, s, sdot, sddot] = s_solution(t_max - impact_time, alpha, beta, gamma, eps);

% Composite force solution
wagner_force = composite_force(wagner_t, s, sdot, sddot, eps);

%% Plotting
close all;

figure(1);
hold on;

% Read computational solution
output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(1)));
ts = output_mat(:, 1);
ss = output_mat(:, 6);

% Shift time so t = 0 happens at impact time
ts = ts - impact_time;

% Plot computational solution
h(1) = plot(ts, ss, 'color', 0.6 * [1 1 1], 'Linewidth', 2);

% Plot analytical solution
h(2) = plot(wagner_t, s, 'color', 0.6 * [1 1 1], 'Linewidth', 1.5, 'Linestyle', '--');

% Plots vertical lines at labeled points
h(5) = xline(0.015);
h(6) = xline(0.295);
h(7) = xline(0.535);
h(8) = xline(0.675);


legend(h([5, 6, 7, 8]), ...
    ["Stationary plate, analytical", "Stationary plate, numerical", ...
        "Moving plate, analytical", "Moving plate, numerical"], ...
    "Interpreter", "latex", "location", "northoutside", "Fontsize", 12, ...
    "Numcolumns", 2);

xlim([-impact_time t_max - impact_time]);
ylim([0 1.1 * max(s)]);
set(gca, 'XTick', -impact_time : impact_time : t_max - impact_time);
grid on;
xlabel("$t$", "Interpreter", "latex");
ylabel("$s(t)$", 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
set(gcf, 'Position',  [0, 0, 700, 250]);
plot_name = sprintf("%s/plate_displacement_comparison.png", analysis_directory);
pause(0.1);
exportgraphics(ax, plot_name, 'resolution', 300);
