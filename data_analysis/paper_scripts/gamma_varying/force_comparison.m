%% force_comparison.m
% Script to compare the forces on the plate for where the initial height
% has been varied. The times are shifted so the theoretical time of impact
% is always the same
%

% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions


% Parent directory where all the data is held
parent_directory = '/mnt/newarre/cantilever_paper_data/gamma_varying/';

% Directory to save the figure(s)
analysis_directory = "Analysis";
analysis_directory = strcat(parent_directory, analysis_directory);

% Range of alphas
gammas = [0, 10, 100, 500, 1000];

% Defines arrays for all the values of alpha
data_directories = string(length(gammas));
legend_entries = string(length(gammas));

for k = 1 : length(gammas)
    gamma = gammas(k);
    
    data_directories(k) = [parent_directory, '/gamma_', num2str(gamma)];
    
    legend_entries(k) = ['$\gamma =$ ', num2str(gamma)] ;
end

% Stationary plate data
stationary_plate_directory = '/mnt/newarre/cantilever_paper_data/stationary_plate';


%% Parameters

% Value of epsilon
eps = 1;

% Plate parameters
alpha = 2;
beta = 0;

% Initial drop height 
initial_drop_heights = 0.125;

% Impact time 
impact_time = initial_drop_heights;

% Maximum time
t_max = 0.8;

%% Wagner solution
% % Stationary plate Wagner solution
% wagner_t = linspace(1e-7, t_max - impact_time, 1e3);
% s = zeros(size(wagner_t));
% sdot = zeros(size(wagner_t));
% sddot = zeros(size(wagner_t));
% 
% % Composite force solution
% wagner_force = composite_force(wagner_t, s, sdot, sddot, eps);

%% Plotting force comparisons
close all;

figure(1);
hold on;

for k = 1 : length(data_directories)
   output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(k)));
   
   ts = output_mat(:, 1);
   Fs = output_mat(:, 3);
   
   figure(1);
   plot(ts, Fs, 'Displayname', legend_entries(k));
end

% % Analytical solution
% plot(wagner_t + impact_time, wagner_force, ...
%     'color', 0.5 * [1 1 1], 'Linewidth', 2, 'Linestyle', '--', ...
%     'Displayname', ['Analytical', newline, 'solution']);

% Stationary plate solution
output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", ...
    stationary_plate_directory));
ts = output_mat(:, 1);
Fs = output_mat(:, 3);
plot(ts, Fs, 'Linewidth', 2, 'color', 'black', ...
    'Displayname', ['Stationary', newline, 'plate']);


% x limits
xlim([0 t_max]);


% % Arrow for increasing alpha
% X = [0.4 0.5];
% Y = [0.2 0.6];
% annotation('arrow', X, Y);
% 
% % Arrow label
% txt = '$\alpha$';
% text(0.39, 2.86, txt, "Interpreter", "Latex", "Fontsize", 14);

legend("Interpreter", "latex", "location", "northwest", "Fontsize", 12);

grid on;
xlabel("$t$", "Interpreter", "latex");
ylabel("$F(t)$", 'Interpreter', 'latex');

ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title(['Force on plate: $\alpha =$ ', ...
    num2str(alpha), ', $\beta =$ ' num2str(beta)], "Interpreter", "latex", ...
    'Fontsize', 14);
print(gcf, sprintf("%s/force_comparison.png", analysis_directory), ...
    '-dpng', '-r300');