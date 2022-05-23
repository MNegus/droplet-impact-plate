%% force_comparison.m
% Script to compare the forces on the plate for different plate widths,
% plus comparing to the analytical solution
%

% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions

% Parent directory where all the data is held
parent_directory = '/mnt/newarre/cantilever_paper_data/validation_alpha_100/plate_width_validation/';


% Directory to save the figure(s)
analysis_directory = "Analysis";
analysis_directory = strcat(parent_directory, analysis_directory);

% Defines array with both directories in
data_directories = ["plate_width_2", "plate_width_3", "plate_width_4", "plate_width_5"];

% Concatenates arrays to include parent directory
for k = 1 : length(data_directories)
    data_directories(k) = strcat(parent_directory, data_directories(k)); 
end

legend_entries = ["Plate width = 2", "Plate width = 3", "Plate width = 4", "Plate width = 5"];


%% Parameters

% Value of epsilon
eps = 1;

% Plate parameters
alpha = 100;
beta = 0;
gamma = 0;

% Initial drop height
initial_drop_height = 0.125;

% Impact time
impact_time = initial_drop_height;

% Maximum time
t_max = 0.5;

%% Wagner solution

% Plate displacement solution
[wagner_t, s, sdot, sddot] = s_solution(t_max - impact_time, alpha, beta, gamma, eps);

% Composite force solution
wagner_force = composite_force(wagner_t, s, sdot, sddot, eps);

%% Plotting
close all;

% Force comparison
figure(1);
hold on;

% % s displacement comparison
% figure(2);
% hold on;

for k = 1 : length(data_directories)
   output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(k)));
   
   ts = output_mat(:, 1);
   Fs = output_mat(:, 3);
%    ss = output_mat(:, 6);
   
   figure(1);
   plot(ts, Fs);
    
%    figure(2);
%    plot(ts, ss);
end

figure(1);
plot(wagner_t + impact_time, wagner_force, ...
    'color', 0.5 * [1 1 1], 'Linewidth', 2, 'Linestyle', '--');
legend([legend_entries, "Analytical"], ...
    "Interpreter", "latex", "location", "southeast", "Fontsize", 12);


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


% figure(2);
% plot(wagner_t + impact_time, s);