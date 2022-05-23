% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions


% Parent directory where all the data is held
parent_directory = '/mnt/newarre/cantilever_paper_data/alpha_varying/';

% Directory to save the figure(s)
analysis_directory = "Analysis";
analysis_directory = strcat(parent_directory, analysis_directory);

% Range of alphas
alphas = 1;

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

% Value of epsilon
eps = 1;

% Plate parameters
beta = 0;
gamma = 0;

% Initial drop height 
initial_drop_heights = 0.125;

% Impact time 
impact_time = initial_drop_heights;

% Maximum time
t_max = 0.8;

%% Wagner solution

% Plate displacement solution
[wagner_t, s, sdot, sddot] = s_solution(t_max - impact_time, alpha, beta, gamma, eps);

% Composite force solution
wagner_force = composite_force(wagner_t, s, sdot, sddot, eps);


%% Plotting force comparisons
close all;

fig = figure(1);
hold on;


for k = 1 : length(data_directories)
   output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(k)));
   
   ts = output_mat(:, 1);
   ss = output_mat(:, 6);
   sdots = output_mat(:, 7);
   sddots = output_mat(:, 8);
   
   % Rescale ts with the impact time
   ts = ts - impact_time;
end

% Plots the numerical values
plot(ts, sddots);
plot(ts, (1 - sdots).^2 ./ (2 * (ts - ss)));


% Plots analytical values
plot(wagner_t, sddot);
plot(wagner_t, (1 - sdot).^2 ./ (2 * (wagner_t - s)));
ylim([0 2]);



