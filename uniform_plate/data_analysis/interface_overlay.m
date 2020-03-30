function interface_overlay(plot_type)
%% interface_overlay.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creates a video plot of the interfaces of the droplet in the
% Plate_validation simulations overlayed, in order to see where the
% simulations differ. If plot_type == "macro", then a video of the
% macroscopic view of the whole droplet is created. If plot_type == "jet",
% then then video will follow the jet tip, also plotting the Wagner
% solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Computational definitions
% Here we specify the parameters that were used for the simulation of
% relevance to this analysis

% Adds the interface analysis to path
addpath("~/repos/plate-impact/data_analysis/interface_analysis");

MAXLEVEL = 12; % Maximum refinement level of the simulation
BOX_WIDTH = 4.0; % Width of the computational box
min_cell_size = BOX_WIDTH / (2 ^ MAXLEVEL); % Size of smallest cell
impact_time = 0.2015; % Theoretical time of impact

% Base height for the plate, which is perturbed on the three plate runs
base_plate_height = 1000 * min_cell_size;

% Heights for full alignment, half alignment, quarter alignment and solid
% wall (where the plate height is zero)
plate_heights = [base_plate_height, ...
    base_plate_height + 0.5 * min_cell_size, ...
    base_plate_height + 0.25 * min_cell_size, 0];

%% Plotting parameters
%
% Colors of the lines
line_colors = [[0    0.4470    0.7410]; 
    [0.8500    0.3250    0.0980];
    [0.9290    0.6940    0.1250];
    [0.4940    0.1840    0.5560]];

%% Data definitions
% Here we specify the location of the data and other parameters required to
% conduct the analysis

% Master directory where all the data is stored
master_directory = '/mnt/newarre/oscillation_study/cleaned_data';

% Names of the individual directories where the data is stored
data_directories = ["uniform_plate"];

% Edits so the full directory address is given in data_directories
for k = 1:length(data_directories)
    data_directories(k) = strcat(master_directory, data_directories(k));
end

% Readable names to label the plots for each of the data directories
legend_entries = ["Plate"];

dir_arr = [1]; % Specifies which directories to plot

%% Function definitions
% Defines any anonymous functions needed

% Free surface location from the leading order outer solution of Wagner
% theory
outer_h = @(r, d, t) 0.5 * r.^2 - t ...
    + (2/pi) * (asin(d ./ r) .* (t - 0.5 * r.^2) ...
        + 0.5 * d * sqrt(r.^2 - d^2));

% Free surface location in the inner region, parametrised by the array of
% real numbers numbers sigma
function [r, h] = inner_h(sigma, d, t)
    J = 2 * t^(3/2) / (sqrt(3) * pi); % Jet thickness
    r = d + (J / pi) * (sigma - log(sigma) - 1);
    h = J * (1 + 4 * sqrt(sigma) / pi);
end

% Composite expansion
function [r, h] = composite_h(sigma, d, t)
    
    [r, inner_h_vals] = inner_h(sigma, d, t);
    outer_h_vals = outer_h(r, d, t);
    
    overlap_term = (2/pi) * sqrt(2 * (r - d)) * (d^(3/2) - t / sqrt(d)) ;
    
    h = inner_h_vals + outer_h_vals - overlap_term;
    
end

%% Plotting
% Creates videos of the macro-view of the simulation, i.e. zoomed out of
% the whole droplet. As the interface is highly refined, this plot utilises
% coarsening to reduce the amount of time it takes to plot the interface

% Sets up options for plot
figure_no = 1; % Figure number
close(figure(figure_no)); % Closes existing plot
figure(figure_no); % Creates figure

% Resolution and aspect ratio of figure
width=800;
height=800;
set(gcf,'position',[10,10,width,height])
set(gca,'DataAspectRatio',[1 1 1])

include_points = false; % Excludes data points from being plotted

if plot_type == "macro"
    x_limits = [0, 2]; % Specify x limits (empty = auto)
    y_limits = [0, 2]; % Specify y limits (empty = auto)
    
    % Coarsening options
    coarsen = true; % Coarsens the data set being plotted
    coarsen_number = 4; % Number of times to conduct coarsening
    coarsen_tol = 1e-2;
    remove_rest = true; % Removes points such as entrapped bubbles
    
    % Create the video writer object
    writerObj = VideoWriter('Videos/interface_overlay.avi');
elseif plot_type == "jet"
    % Coarsens if the number of points on a line is more than the threshold 
    coarsen = true; 
    coarsen_threshold = 300; 
    coarsen_tol = 1e-2;
    remove_rest = true;
    
    orig_plot_width = 0.1; % Width of plot
    
    % Create the video writer object
    writerObj = VideoWriter('Videos/jet_tip.avi');
end

writerObj.FrameRate = 5;
open(writerObj);

% Loops over plots for specified time interval
for plot_no = 200:500
    
    
    t = plot_no * 0.001; % Current time
    
    % Calculates turnover point location
    if t > impact_time
        d = sqrt(3 * (t - impact_time)); % From Wagner theory
    else
        d = 0; % Turnover point 0 before impact
    end
    
    if plot_type == "jet"
        plot_width = max([orig_plot_width, 0.4 * d]);
        if d < 0.25 * plot_width
            x_limits = [0, plot_width];
        else
            % Adjusts x and y limits to be centred on the turnover point
            x_limits = [d - 0.25 * plot_width, d + 0.75 * plot_width];
        end
        y_limits = [-0.005, plot_width - 0.005];
    end

    % Holds current figure
    hold on; 
    
    % Loops over the different simulation outputs
    for k = dir_arr
        % Name of file where interface data is stored
        filename = sprintf("%s/interface_%d.txt", data_directories(k), ...
            plot_no);
        
        % Reads the starting and ending points of the line segments
        [start_points, end_points] = ...
            read_interface_points(filename, true);
        
        % Adjusts the z components to reflect the varying plate heights
        start_points(:, 2) = start_points(:, 2) - plate_heights(k);
        end_points(:, 2) = end_points(:, 2) - plate_heights(k);
        
        % Finds the indexes of points which are in the limits
        keep_points = (start_points(:, 1) > x_limits(1)) ...
            & (start_points(:, 1) < x_limits(2)) ...
            & (start_points(:, 2) > y_limits(1)) ...
            & (start_points(:, 2) < y_limits(2));
        
        % Remove all points outside the limits
        start_points = start_points(keep_points, :);
        end_points = end_points(keep_points, :);

        % Coarsens the data
        if plot_type == "jet"
            if coarsen == true && length(start_points) > coarsen_threshold
                loop_no = 1;
                while length(start_points) > 2 * coarsen_threshold ...
                        && loop_no < 5
                    loop_no = loop_no + 1;
                    disp("hello")
                    % Coarsens from top
                    [~, first_seg] = max(start_points(:, 2));

                    % Conducts coarsening without removing bubbles
                    [start_points, end_points] ...
                        = coarsen_interface(start_points, end_points, ...
                            first_seg, 1, coarsen_tol, false);
                end
                % Coarsens from top
                [~, first_seg] = max(start_points(:, 2));

                % Conducts coarsening with removal, if specified
                [start_points, end_points] ...
                    = coarsen_interface(start_points, end_points, ... 
                        first_seg, 1, 1e-2, remove_rest);
            end
        end
        length(start_points)
        % Plots the line segments with associated colour and label
        line_segment_plot(start_points, end_points, include_points, ...
            x_limits, y_limits, figure_no, line_colors(k, :), legend_entries(k));
    
    end
    
    % Plots Wagner solution if appropriate
    if plot_type == "jet"
        % Outer solution
       rs = linspace(d, x_limits(2), 1e3);
       plot(rs, outer_h(rs, d, t - impact_time), 'color', 'black', ...
           'Linestyle', '--', 'Displayname', 'Outer Wagner');
       
       % Inner solution
       indices = linspace(-6, 3, 1e4);
       sigmas = 10.^indices;
       [inner_rs, inner_hs] = inner_h(sigmas, d, t - impact_time);
       plot(inner_rs, inner_hs, 'color', 'black', ...
           'Linestyle', ':', 'Displayname', 'Inner Wagner');
       
       % Composite solution
       [composite_rs, composite_hs] = composite_h(sigmas, d, t - impact_time);
       plot(composite_rs, composite_hs, 'color', 'black', ...
           'Displayname', 'Composite');
    end
    
    % Sets up figure properties
    set(gca,'DataAspectRatio',[1 1 1]); % Enforces aspect ratio
    
    % Title listing current time
    title(sprintf("Interfaces at $t = %.03f$. Impact time = %.03f\n", t, impact_time), "Interpreter", "latex", 'Fontsize', 15);
    
    % Axes labels and font sizes
    xlabel("$r$", "Interpreter", "latex", 'Fontsize', 30);
    ylabel("$z$", "Interpreter", "latex", 'Fontsize', 30);
    ax = gca;
    ax.FontSize = 16;
    set(gca,'TickLabelInterpreter','latex');
    grid on;
    
    % Sets up the legend
    L = legend;
    set(L, 'Interpreter', 'latex');
    set(L, 'FontSize', 15);
    set(L, 'Location', 'northeast');
    
    % Draws the figure and writes to the video writer object
    grid on;
    drawnow;
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    
    % Clears the frame
    clf;
end
close(writerObj); % Closes the object, saving the video

end
