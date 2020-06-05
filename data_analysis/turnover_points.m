function ds = turnover_points(output_range, parent_dir, plate_position, ...
    plate_tol, bubble_box_width)
%TURNOVER_POINTS Reads simulation output to find the turnover points
%   Function to read the interface files, which are of the form
%   "interface_n.txt", where n is an integer in output_range and the files
%   are in data_dir. plate_position indicates the vertical position we
%   expect the plate to be in, which can be re-defined as the zero-vertical
%   component. 

% Scripts for analysing interface files
addpath("interfaces");

% Matrix to save turnover points in
ds = zeros(length(output_range), 2);

for k = 1 : length(output_range)
    % Name of interface file
    filename ...
        = sprintf('%s/raw_data/interface_%d.txt', ...
            parent_dir, output_range(k));
        
    % Reads in the start and end points of the line segments, with z along
    % the horizontal axis and r along the vertical
    [start_points, end_points] = read_interface_points(filename, false); 
    
    % Finds unique values of z in all the points
    all_points = [start_points; end_points];
    [~, uniq_idxs, ~] = uniquetol(all_points(:, 1), 1e-4);
    uniq_points = all_points(uniq_idxs, :);
    
    % Shifts z coordinate to be zero along the plate
    uniq_points(:, 1) = uniq_points(:, 1) - plate_position;
    
    % Removes all points that are below the plate, i.e. only keep points
    % that have a z coordinate above the plate tolerance
    uniq_points = uniq_points(uniq_points(:, 1) > plate_tol, :);
    
    % Removes all points inside the bubble box
    uniq_points = uniq_points((uniq_points(:, 1) > bubble_box_width) ...
        | (uniq_points(:, 2) > bubble_box_width), :);
    
    % Sorts the resulting vector in increasing z order
    [~, sorted_idxs] = sort(uniq_points(:, 1));
    sorted_points = uniq_points(sorted_idxs, :);
    
    % Define r as a continuous function of z using interpolation
    r_interp = @(z) interp1(sorted_points(:, 1), sorted_points(:, 2), z);
    
    % Vertical range to search for turnover point in. Corresponds to the
    % bottom of the droplet up to the middle, assuming that the impact is
    % not far along enough that the turnover point is above the vertical
    % mid-point of the droplet.
    zs = linspace(min(all_points(:, 1)), 0.5 * max(all_points(:, 1)), 1e5);
    
    % Finds local minima of the r_interp function
    local_mins = islocalmin(r_interp(zs));
    
    if nnz(local_mins) == 0
        % If no local minima are found, set the turnover point to just be
        % at r = 0, z = 0
        turnover_z = 0;
        turnover_r = 0;
    else  
        % Else, the turnover point is the local minimum with the highest z 
        % value
        turnover_z = max(zs(local_mins));
        turnover_r = r_interp(turnover_z);
    end
        
    % Saves r and z coordinate of the turnover point
    ds(k, 1) = turnover_r;
    ds(k, 2) = turnover_z + plate_position;
    
    % Plots the turnover point on a graph
    figure(1);
    plot(zs, r_interp(zs));
    hold on;
    scatter(ds(k, 2), ds(k, 1));
    hold off;
    title(num2str(k));
    
end

end