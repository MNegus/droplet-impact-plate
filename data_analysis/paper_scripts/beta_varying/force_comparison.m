%% force_comparison.m
% Script to compare the forces on the plate for where the initial height
% has been varied. The times are shifted so the theoretical time of impact
% is always the same
%

% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions


% Parent directory where all the data is held
parent_directory = '/mnt/newarre/cantilever_paper_data/damping_alpha_2_gamma_100/';

% Directory to save the figure(s)
analysis_directory = "Analysis";
analysis_directory = strcat(parent_directory, analysis_directory);

% Range of alphas
betas = [0, 7.07,  28.28, 141.42];

% Defines arrays for all the values of alpha
data_directories = string(length(betas));
legend_entries = string(length(betas));

for k = 1 : length(betas)
    beta = betas(k);
    
    data_directories(k) = [parent_directory, '/beta_', num2str(beta)];
    
    legend_entries(k) = ['$\beta =$ ', num2str(beta)] ;
end

% Stationary plate data
stationary_plate_directory = '/mnt/newarre/cantilever_paper_data/stationary_plate';


%% Parameters

% Value of epsilon
eps = 1;

% Plate parameters
alpha = 2;
beta = 0;
gamma = 100;

% Initial drop height 
initial_drop_heights = 0.125;

% Impact time 
impact_time = initial_drop_heights;

% Maximum time
t_max = 0.8;

%% Wagner solution
% Stationary plate Wagner solution
% wagner_t = linspace(1e-7, t_max - impact_time, 1e3);
% s = zeros(size(wagner_t));
% sdot = zeros(size(wagner_t));
% sddot = zeros(size(wagner_t));
% 
% % Composite force solution
% wagner_force = composite_force(wagner_t, s, sdot, sddot, eps);

%% Plotting
close all;

% Force comparison
figure(1);
hold on;


for k = 1 : length(data_directories)
   output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(k)));
   
   ts = output_mat(:, 1);
   Fs = output_mat(:, 3);
   ss = output_mat(:, 6);
   sdots = output_mat(:, 7);
   sddots = output_mat(:, 8);
   
   figure(1);
   plot(ts, Fs, 'Displayname', legend_entries(k));
   
end

% Force comparison
figure(1);
output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", ...
    stationary_plate_directory));
ts = output_mat(:, 1);
Fs = output_mat(:, 3);
plot(ts, Fs, 'Linewidth', 2, 'color', 'black', ...
    'Displayname', ['Stationary', newline, 'plate']);

% plot(wagner_t + impact_time, wagner_force, ...
%     'color', 0.5 * [1 1 1], 'Linewidth', 2, 'Linestyle', '--', ...
%     'Displayname', ['Analytical', newline, 'solution']);
legend("Interpreter", "latex", "location", "southeast", "Fontsize", 12);

% x limits
xlim([0 t_max]);


% Arrow for increasing beta
X = [0.5 0.4];
Y = [0.5 0.65];
annotation('arrow', X, Y);
grid on;
xlabel("$t$", "Interpreter", "latex");
ylabel("$F(t)$", 'Interpreter', 'latex');

% Arrow label
txt = '$\beta$';
text(0.31, 2.65, txt, "Interpreter", "Latex", "Fontsize", 14);



ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title(['Force on plate: $\alpha =$ ', ...
    num2str(alpha), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
    'Fontsize', 14);
print(gcf, sprintf("%s/force_comparison.png", analysis_directory), ...
    '-dpng', '-r300');

