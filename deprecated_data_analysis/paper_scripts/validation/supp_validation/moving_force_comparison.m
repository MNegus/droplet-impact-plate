%% force_comparison.m
% Script to compare the forces on the plate for where the initial height
% has been varied. The times are shifted so the theoretical time of impact
% is always the same
%

% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions


% Parent directory where all the data is held
parent_directory = '/media/negus/newarre/supplementary_material/validation/moving_validation/';

% Directory to save the figure(s)
analysis_directory = "Analysis";
analysis_directory = strcat(parent_directory, analysis_directory);

% Defines array with both directories in
data_directories = ["level_10", "level_11", "level_12", "level_13", "level_14"];

% Concatenates arrays to include parent directory
for k = 1 : length(data_directories)
    data_directories(k) = strcat(parent_directory, data_directories(k)); 
end

legend_entries = ["Level 10", "Level 11", "Level 12", "Level 13", "Level 14"];


%% Parameters

% Value of epsilon
eps = 1;

% Plate parameters
alpha = 2;
beta = 7.07;
gamma = 100;

% Impact time 
impact_time = 0.125;

% Maximum time
t_max = 0.8;

% Timestep
dt = 1e-4;

%% Wagner solution

% Plate displacement solution
[wagner_t, s, sdot, sddot] = s_solution(t_max - impact_time, alpha, beta, gamma, eps);
% wagner_t = dt : dt : t_max - impact_time;

% Composite force solution
wagner_force = composite_force(wagner_t, s, sdot, sddot, eps);
% wagner_force = composite_force(wagner_t, 0, 0, 0, eps);

%% Plotting
close all;

% Force comparison
figure(2);
hold on;


for k = 1 : length(data_directories)
   output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(k)));
   
   ts = output_mat(:, 1);
   Fs = output_mat(:, 3);
   sdots = output_mat(:, 7);
   
   figure(2);
   plot(ts, Fs, 'Linewidth', 2);
   
end

figure(2);
plot(wagner_t + impact_time, wagner_force, ...
    'color', 0.5 * [1 1 1], 'Linewidth', 2, 'Linestyle', '--');
legend([legend_entries, "Analytical"], ...
    "Interpreter", "latex", "location", "northwest", "Fontsize", 12);


xlim([0 t_max]);
grid on;
xlabel("$t$", "Interpreter", "latex");
ylabel("$F(t)$", 'Interpreter', 'latex');

ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title(['Force on plate: $\alpha =$ ', num2str(alpha), ', $\beta =$ ', ...
    num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
    'Fontsize', 14);
print(gcf, sprintf("%s/force_comparison.png", analysis_directory), ...
    '-dpng', '-r300');
savefig(sprintf("%s/force_comparison.fig", analysis_directory));