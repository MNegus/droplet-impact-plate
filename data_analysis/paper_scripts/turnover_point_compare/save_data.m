addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions

% Parent directory where all the data is held
parent_directory = '/mnt/newarre/cantilever_paper_data/alpha_varying/';

% Directory to save the figure(s)
turnover_directory = "turnover_points";
turnover_directory = strcat(parent_directory, turnover_directory);

% Range of alphas
alphas = [1, 2, 5, 10, 20, 100];

% Defines arrays for all the values of alpha
data_directories = string(length(alphas));
legend_entries = string(length(alphas));

for k = 1 : length(alphas)
    alpha = alphas(k);
    
    data_directories(k) = [parent_directory, '/alpha_', num2str(alpha)];
end

% Stationary plate data
stationary_plate_directory = '/mnt/newarre/cantilever_paper_data/stationary_plate';

%% Parameters
output_range = 1 : 800;
plate_position = 0;
plate_tol = 1e-3;
bubble_box_width = 5e-2;

%% Saves turnovers
% for k = [2, 3, 4, 5, 6]
    
%     data_dir = data_directories(k);
    data_dir = stationary_plate_directory;
    
    ds = turnover_points(output_range, data_dir, plate_position, ...
    plate_tol, bubble_box_width);

    output_matrix = zeros(length(output_range), 3);
    output_matrix(:, 1) = output_range;
    output_matrix(:, 2) = ds(:, 1);
    output_matrix(:, 3) = ds(:, 2);
    
%     dlmwrite(sprintf('%s/alpha_%d.txt', turnover_directory, alphas(k)), ...
%         output_matrix);
    dlmwrite(sprintf('%s/stationary.txt', turnover_directory), ...
        output_matrix);
% end