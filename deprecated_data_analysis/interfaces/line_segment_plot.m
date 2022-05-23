function line_segment_plot(start_points, end_points, include_points, ...
    x_limits, y_limits, figure_no, color, line_name)
%line_segment_plot.m Creates a plot of a series of line segments given 
%their start and end points
%   Function creates a figure of a collection of line segments. If there
%   are N line segments, then start_points is an Nx2 matrix, with each row
%   giving the (x, y) coordinates of the first point in a line segment.
%   Similarly, end_points is an Nx2 matrix of the end points of the
%   segment. This function then plots all the line segments on the same
%   figure. 
%
%   start_points = Matrix of the starting points of the line segments
%   end_points = Matrix of the end points of the line segments
%   include_points = true if the start and end points need to be included
%   in the figure
%   x_limits = The x limits of the plot
%   y_limits = The y limits of the plot
%   figure_no = The number to assign the figure
    
% If figure_no is passed in, then the figure is given that specific
% number. Otherwise it is created without a number
    if exist('figure_no', 'var')
        figure(figure_no);
    else
        figure(1);
    end
    
    if ~exist('line_name', 'var')
        line_name = "Interface";
    end
    
    % Plots the individual line segments
    if exist('color', 'var')
        first = plot(...
            [start_points(1, 1)'; end_points(1, 1)'], ...
            [start_points(1, 2)'; end_points(1, 2)'], ...
            'color', color, 'Linewidth', 1, 'Displayname', line_name);
        interface_lines = plot(...
            [start_points(:, 1)'; end_points(:, 1)'], ...
            [start_points(:, 2)'; end_points(:, 2)'], ...
            'color', color, 'Linewidth', 1, 'HandleVisibility','off');
    else
        first = plot(...
            [start_points(1, 1)'; end_points(1, 1)'], ...
            [start_points(1, 2)'; end_points(1, 2)'], ...
            'color', 'b', 'Linewidth', 1, 'Displayname', line_name);
        interface_lines = plot(...
            [start_points(:, 1)'; end_points(:, 1)'], ...
            [start_points(:, 2)'; end_points(:, 2)'], ...
            'color', 'b', 'Linewidth', 1, 'HandleVisibility','off');
    end
    
    % If include_points is specified, then the start and end points are
    % plotted 
    if exist('include_points', 'var')
        if include_points 
            hold on;
            scatter(start_points(:, 1), start_points(:, 2), [], 'g');
            scatter(end_points(:, 1), end_points(:, 2), [], 'r'); 
        end
    end
    
    % Set the x and y limits if specified
    if exist('x_limits', 'var')
        if not(isempty(x_limits))
            xlim(x_limits);
        end
    end
    if exist('y_limits', 'var')
       if not(isempty(y_limits))
            ylim(y_limits);
        end 
    end
end

