function interface_overlay(plot_type)
%% interface_overlay.m
% Overlays the droplet interface of various simulations for comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to create videos where the interface of the droplet in different
% simulations are overlayed for visual comparison. They could be the same
% type of simulation (e.g. solid wall, embedded plate, or even different
% types. The outputs are translated so the x axis of the video is along the
% surface of the plate, which allows for more meaningful comparisons of
% simulations where the plate is moving at different speeds. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Data definitions
% Here we specify the locations where the interface files are stored. We
% expect a different directory for each simulation result.

% Directory where the resulting videos are to be stored
results_directory = "/mnt/newarre/analysis_tests";

% Master directory where all of the data is stored under (e.g. external
% hard drive location)
master_directory = "/mnt/newarre/analysis_tests/";

% Individual directory names under the master directory. We assume that the
% interface files are stored under
% master_directory/data_directory(k)/raw_data
data_directories = ["wall_test", "stationary_prescribed", "coupled_test"];
no_dirs = length(data_directories); % Number of entries

% Entries in the legend
legend_entries = ["Wall", "Embedded", "$\alpha = 0.01, \gamma = 1$"];


% Concatenates the data_directories array with the master directory
for k = 1 : length(data_directories)
   data_directories(k) = strcat(master_directory, data_directories(k)); 
end

%% Plotting parameters
% Some parameters for customising the plots

% MAKE BASED ON TIME
% Range of plot numbers to record
plot_range = 115:600;

% Colors of the lines (assuming there are a maximum of 4 lines) 
line_colors = [[0    0.4470    0.7410]; 
    [0.8500    0.3250    0.0980];
    [0.9290    0.6940    0.1250];
    [0.4940    0.1840    0.5560]];

% Set true if the start and end points of the segments are to be plotted
include_points = false; 

% Specify x and y limits (empty = auto)
y_lower = -0.05
x_limits = [0, 2]; y_limits = [y_lower, 2 + y_lower]; 

% Width and height of figure (in pixels)
width=800; height=800;

% Video frame rate
frame_rate = 5;

%% Analytic solutions
% Defines anonymous functions for the analytic solutions to the interface
% location (i.e. Wagner theory results)

%% Reads interface_times.txt files 
% There should be one file per simulation which has three columns of
% numbers, the first is the "interface_output_no", which just counts how
% many times the interface has been outputted. The second column is the
% time it was outputted and the third the plate position at that time. The
% simulations should be synced such that if two interface output numbers
% are the same, then the time is the same. 

% Initialise array (sets to a max length)
interface_times = zeros(1000, 2 + no_dirs);

max_file_length = 0; % Current "maximum" length of a file

% Loops over the directories
for k = 1 : no_dirs
    % Read the impact_times file for this simulation
    read_mat = dlmread(...
        sprintf('%s/raw_data/interface_times.txt', data_directories(k))); 

    % Finds the length of the file
    file_length = length(read_mat(:, 1));
    if file_length > max_file_length
        % If the file is the longest so far, then make the first column
        % (plot_nos) and the second column (times) equal to this files
        % values
        interface_times(1:file_length, 1:2) = read_mat(:, 1:2);
        max_file_length = file_length;
    end
    
    % Reads the plate position from the third column
    interface_times(1:file_length, 2 + k) = read_mat(:, 3);
end

% Remove any zeros at the end
interface_times = interface_times(interface_times(:, 1) > 1, :);

%% Plotting
% Creates a video overlaying the interfaces of the different simulations

% Creates the figure
figure_no = 1; % Figure number
close(figure(figure_no)); % Closes existing plot
figure(figure_no); % Creates figure

% Resolution and aspect ratio of figure
set(gcf,'position',[10,10,width,height]);
set(gca,'DataAspectRatio',[1 1 1]);

% Create and open video writer
writerObj ...
    = VideoWriter(sprintf('%s/interface_overlay.avi', results_directory));
writerObj.FrameRate = frame_rate;
open(writerObj);

for plot_no = plot_range

    % Loop over the different simulation outputs
    for k = 1 : no_dirs
        % Notes the plate position at this timestep
        plate_position = interface_times(plot_no, 2 + k);
        
        % Name of the specific interface file
        interface_filename = sprintf("%s/raw_data/interface_%d.txt",...
            data_directories(k), plot_no);
        
        % Reads the interface file into a matrix. Column 1 gives the values
        % of z and column 2 gives the values of r
        interface_mat = dlmread(interface_filename);
        
        % Removes non-unique points up to a tolerance
        [~, unique_idxs, ~] = uniquetol(interface_mat(:, 1), 1e-5);
        unique_mat = interface_mat(unique_idxs, :);
        
        % Sorts the matrix in ascending values of z (column 1)
        [~, sorted_idxs] = sort(unique_mat(:, 1));
        sorted_mat = unique_mat(sorted_idxs, :);
        
        % Creates a line-of-best-fit of r as a function of z from the data
        z_query = linspace(min(sorted_mat(:, 1)), max(sorted_mat(:, 1)), 1e3);
        r_interp = interp1(sorted_mat(:, 1), sorted_mat(:, 2), z_query);
        
        % Shift the z values by the plate position
        z_query = z_query - plate_position;
        
        % Plot the line of best fit 
        figure(figure_no);
        plot(r_interp, z_query, 'color', line_colors(k, :), ...
            'Linewidth', 2);
        
        % Turns on the hold after the first plot
        if k == 1
            hold on;
        end
    end
    
    % Draws a horizontal line at z = 0 (where the plate is)
    yline(0, '--', 'color', 0.2 * [1 1 1]);
    
    hold off;
    
    %%%%%%%%%%%%%%%%%%%
    % Figure properties
    %%%%%%%%%%%%%%%%%%%
    % Title
    title(sprintf("Interfaces at $t = %.03f$", ...
        interface_times(plot_no, 2)), ...
        "Interpreter", "latex", 'Fontsize', 15);
    
    % Axes limits, labels and font sizes
    xlim(x_limits); ylim(y_limits); % Axes limits
    xlabel("$r$", "Interpreter", "latex", 'Fontsize', 30);
    ylabel("$z$", "Interpreter", "latex", 'Fontsize', 30);
    ax = gca;
    ax.FontSize = 16;
    set(gca,'TickLabelInterpreter','latex');
    
    % Set up the legend
    legend(legend_entries, 'Interpreter', 'latex', 'FontSize', 15, ...
        'Location', 'northeast');
    
    grid on; % Puts on a grid
    
    % Draws frame to video
    drawnow;
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
end
close(writerObj);


end