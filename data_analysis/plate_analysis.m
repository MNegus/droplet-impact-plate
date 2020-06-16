%% plate_analysis.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyses the plate position outputted by the simulations and compares
% them to the Wagner theory prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters

% Physical parameters
rho_w = 998;
R0 = 1e-3;
U0 = 10;
T0 = R0 / U0;
Patm = 10^5;

% Plate ODE coefficients
alpha = 1;
beta = 0;
gamma = 0;

% Initial drop height
initial_drop_height = 0.125;

impact_time = initial_drop_height;


%% Data definitions
% Here we specify the locations where the plate output files are stored. We
% expect a different directory for each simulation result.

% Parent directory where all of the data is stored under (e.g. external
% hard drive location)
parent_directory = "/mnt/newarre/alpha_1_gamma_1/";

% Directory where the resulting videos are to be stored
results_directory = sprintf("%s/Analysis", parent_directory);

% Individual directory names under the master directory. We assume that the
% plate output files are stored under
% master_directory/data_directory(k)/cleaned_data
data_directories = ["avg_no_2"];
no_dirs = length(data_directories); % Number of entries

% Adds the parent directory to the start of the data directories
for k = 1 : length(data_directories)
    data_directories(k) = strcat(parent_directory, data_directories(k)); 
end
% Readable names to label the plots for each of the data directories
legend_entries = ["Avg no = 2"];

%% s comparison (one data directory for now)
output_mat = dlmread(sprintf('%s/cleaned_data/output.txt',...
                        data_directories(1)));
comp_ts = output_mat(:, 1);
comp_ss = output_mat(:, 5);

t_max = max(comp_ts) - impact_time
[wagner_t, wagner_s, wagner_sdot, wagner_sddot] = s_solution(t_max, alpha, ...
    beta, gamma, eps);

figure(1);
hold on;
plot(comp_ts, comp_ss);
plot(wagner_t + impact_time, wagner_s);

