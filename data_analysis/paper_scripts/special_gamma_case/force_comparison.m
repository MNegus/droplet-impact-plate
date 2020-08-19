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
data_directory = ["gamma_varying/gamma_500"];

% Concatenates arrays to include parent directory
data_directory = strcat(parent_directory, data_directory); 

% Stationary plate data
stationary_plate_directory = '/mnt/newarre/cantilever_paper_data/stationary_plate';


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

%% Stationary plate Wagner solution
stationary_force = composite_force(wagner_t, zeros(size(s)), ...
    zeros(size(s)), zeros(size(s)), eps);

%% Plotting
close all;

figure(1);
hold on;

% Stationary plate solution directory
stationary_output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", ...
    stationary_plate_directory));

% Reads times and forces from stationary plate solution directory
ts = stationary_output_mat(:, 1);
Fs = stationary_output_mat(:, 3);
ts = ts - impact_time; % Shift time so t = 0 happens at impact time


% Plots stationary plate solution
plot(ts, Fs, 'Linewidth', 2, 'color', 0.75 * [1 1 1]);


% Moving plate output matrix
output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directory(1)));

% Reads times and forces from data directory
ts = output_mat(:, 1);
Fs = output_mat(:, 3);
ts = ts - impact_time; % Shift time so t = 0 happens at impact time

% Plots moving plate solution
plot(ts, Fs, 'Linewidth', 2, 'color', 0.25 * [1 1 1]);


% Plots Wagner solutions
plot(wagner_t, stationary_force, 'Linewidth', 2, 'linestyle', '--', ...
    'color', 'black');
plot(wagner_t, wagner_force, 'Linewidth', 2, 'linestyle', ':', ...
    'color', 'black');

% legend(["Stationary computational", "Stationary analytical", ...
%     "Moving computational", "Moving analytical"], ...
%     "Interpreter", "latex", "location", "northwest", "Fontsize", 12);
xlim([-impact_time t_max - impact_time]);
grid on;
xlabel("$t$", "Interpreter", "latex");
ylabel("$F(t)$", 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 12;
ylim([0 5.1]);
set(gca, 'XTick', -impact_time : impact_time : t_max - impact_time);
set(gca, 'YTick', 0 : 1 : 5);
set(gca,'TickLabelInterpreter','latex');
set(gcf, 'Position',  [0, 0, 700, 200]);
ax = gca;
plot_name = sprintf("%s/force_comparison.png", analysis_directory);
pause(0.1);
exportgraphics(ax, plot_name, 'resolution', 300);
