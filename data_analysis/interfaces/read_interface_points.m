function [start_points, end_points] = ...
    read_interface_points(interface_filename, transpose_coordinates)
%read_interface_points.m Reads the points which define line segments from
%output_facets in a Basilisk script into two arrays
%   This function is designed to work with the output_facets function in
%   Basilisk, which outputs the location of the interface of a fluid into a
%   text file as a series of coupled points. If there are N line segments,
%   then the text file has 2N lines, with the odd lines giving the (x, y)
%   coordinates of the start of the line segments and the even lines giving
%   the end points. This function reads the file and puts the start and end
%   points into two matrices: start_points and end_points. If the
%   transposed_coordinates option is true, then the x and y coordinates are
%   swapped, as is usually the case in axisymmetric codes.
%   Reference: http://basilisk.fr/src/fractions.h#interface-output

    % Reads the matrix from the file
    interface_points = dlmread(interface_filename);
    
    % If transpose_coordinates is true, then we swap the x and y
    % coordinates in the matrix
    if exist('transpose_coordinates', 'var')
        if transpose_coordinates
            old_interface_points = interface_points;
            interface_points(:, 1) = old_interface_points(:, 2);
            interface_points(:, 2) = old_interface_points(:, 1);
        end
    end
    
    % Indexes the line segments    
    indices = 1 : 2 : length(interface_points);
    
    % Matrices of the start and end points of the line segments
    start_points = interface_points(indices, :);
    end_points = interface_points(indices + 1, :);
end