eps = 1;
alpha = 1;
beta = 0;
gamma = 40;
t_max = 1.0;
% impact_time = 0.125;

% %% Data definitions
% % Here we specify the locations where the plate output files are stored. We
% % expect a different directory for each simulation result.
% 
% % Parent directory where all of the data is stored under (e.g. external
% % hard drive location)
% parent_directory = "/mnt/newarre/coupled_entrapped_bubble/";
% 
% % Directory where the resulting videos are to be stored
% results_directory = sprintf("%s/Analysis", parent_directory);
% 
% % Individual directory names under the master directory. We assume that the
% % plate output files are stored under
% % master_directory/data_directory(k)/cleaned_data
% data_directories = ["alpha_2_gamma_4000_peak_detect"];
% no_dirs = length(data_directories); % Number of entries
% 
% % Adds the parent directory to the start of the data directories
% for k = 1 : length(data_directories)
%     data_directories(k) = strcat(parent_directory, data_directories(k)); 
% end
% % Readable names to label the plots for each of the data directories
% legend_entries = ["Peak detect"];

%% Solve for s

[t, s, sdot, sddot] = s_solution(t_max, alpha, beta, gamma, eps);

%% Plotting

figure(1);
hold on;
plot(t, s');

figure(2);
plot(t, sdot);

figure(3);
plot(t, sddot);

figure(4);
hold on
plot(t, outer_force(t, s, sdot, sddot, eps));
plot(t, composite_force(t, s, sdot, sddot, eps));
% 
