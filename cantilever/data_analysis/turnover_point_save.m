%% turnover_point_save.m
% Script goes through the interface files, extracts the turnover point 
% locations and saves them to a file

% Adds the interface analysis to path
addpath("~/repos/plate-impact/data_analysis/interface_analysis");

%% Data definitions
% Here we specify the location of the data and other parameters required to
% conduct the analysis

% Master directory where all the data is stored
master_directory = '/mnt/newarre/first_cantilever';

% Analysis directory to save the resulting file to
save_directory = sprintf('%s/data_analysis', master_directory);


%% Parameters
% Later we'll implement this being extracted from the parameters.h file
interface_output_timestep = 1e-3;

%% Data definitions

% Defines the range interface_n.txt files to analyse
interface_range = 1:626;

% Creates the matrix to store the data in
output_matrix = zeros(length(interface_range), 3);

% Loops over all the interface files
for analyse_no = 1 : length(interface_range)
    
    % Calculates which interface_n.txt file to look at
    interface_no = interface_range(analyse_no);
    
    % Calculates the time
    t = interface_no * interface_output_timestep;
    
    % Saves the time in the output matrix
    output_matrix(analyse_no, 1) = t;
    
    % Name of the interface file
    filename = sprintf("%s/raw_data/interface_%d.txt", ...
        master_directory, interface_no);
    
    % Reads the starting and ending points of the line segments
    transpose_coordinates = false; % Makes the r coordinate along the y axis
    [start_points, end_points] = ...
        read_interface_points(filename, transpose_coordinates);
    
%     line_segment_plot(start_points, end_points);
    
    % Finds unique values of x in all the points
    all_points = [start_points; end_points];
    [~, uniq_idxs, ~] = uniquetol(all_points(:, 1), 1e-4);
    uniq_points = all_points(uniq_idxs, :);
    
    % Sorts the resulting vector
    [~, sorted_idxs] = sort(uniq_points(:, 1));
    sorted_points = uniq_points(sorted_idxs, :);
    
    % Define continuous function of x
    all_points_fcn = @(x) interp1(sorted_points(:, 1), sorted_points(:, 2), x);
    
    % Decides range to search for local min in
    xs = linspace(min(all_points(:, 1)), 0.5 * max(all_points(:, 1)), 1e5);
    
    % Finds local min
    local_mins = islocalmin(all_points_fcn(xs));
    if nnz(local_mins) == 0
        turnover_z = 0;
        turnover_r = 0;
    else
        % Conclude the turnover point is the local min with the highest x value
        turnover_z = max(xs(local_mins));
        turnover_r = all_points_fcn(max(xs(local_mins)));
    end
    
    output_matrix(analyse_no, 2) = turnover_r;
    output_matrix(analyse_no, 3) = turnover_z;
    
    % Plot (optional for testing)
%     figure(1);
%     ys = all_points_fcn(xs);
%     title(sprintf("t = %g\n", t));
%     plot(xs, ys, turnover_z, turnover_r, 'r*');
end
writematrix(output_matrix, ...
    sprintf('%s/turnover_points.txt', save_directory));



