function [new_start_points, new_end_points] ...
    = coarsen_interface(start_points, end_points, first_seg, first_m, tol, ...
    remove_rest)
%coarsen_interface.m Coarsens the line segment plot of the interface by
%removing half of the segments
%   Function takes as input start_points and end_points, which
%   are Nx2 matrices of the starting and ending points of the N line
%   segments that are wanted to plotted. This function removes half of them
%   by, starting at first_seg, looks ahead at the line segment which
%   connects to its first closest neighbour, connects to it and then
%   removes the closest neighbour. The choice of first_seghere is
%   important, as the function will only remove the subsequent points, and
%   not the ones before it.

function value = point(seg, m)
%point Utility function for extracting the coordinates of a point, given
%the segment number "seg" and m, where m=1 is a start point and m=-1 is an
%end point.
    if m == 1
        value = start_points(seg, :);
    elseif m == -1
        value = end_points(seg, :);
    else
        disp("Invalid value of m");
        value = -1;
    end
end

function [seg_no, min_m, min_dist] = closest_point(curr_point, seg_nos)
    
    no_segs = length(seg_nos);

    start_dists = (curr_point - start_points(seg_nos, :)).^2;
    end_dists = (curr_point - end_points(seg_nos, :)).^2  ;
    
    dists_mat = zeros(2 * no_segs, 3);
    
    dists_mat(1 : no_segs, :) ...
        = [sqrt(start_dists(:, 1) + start_dists(:, 2)), ...
        seg_nos', ones(no_segs, 1)];
    
    dists_mat(no_segs + 1 : 2 * no_segs, :) ...
        = [sqrt(end_dists(:, 1) + end_dists(:, 2)), ...
        seg_nos', -ones(no_segs, 1)];
    
    [min_dist, min_dists_idx] = min(dists_mat(:, 1));
    
    seg_no = dists_mat(min_dists_idx, 2);
    min_m = dists_mat(min_dists_idx, 3);

end

curr_seg = first_seg; % Index of the starting segment
curr_m = first_m; % Indicates if the first point is a start_point or end_point

% Indices of the neighbouring segment and the segment after that
neighbour_seg = 0; 
next_seg = 0;

% Array to store the checked points
seg_nos = 1:length(start_points);
check_pts = ones(1, length(seg_nos));

% Array to store the points which are being kept
keep_pts = zeros(1, length(seg_nos));


% Loops until we can't find a neighbour or a next segment
while (neighbour_seg ~= -1) && (next_seg ~= -1)
    
    check_pts(curr_seg) = 0; % Removes the current segment from the search
    keep_pts(curr_seg) = 1; % Records that the current point is being kept

    % Point on the opposite side of the current segment, of which we are
    % trying to find another segment which is close to it
    curr_opp_pt = point(curr_seg, -curr_m);
    
    %% Find the neighbouring segment
    
    % Find the closest point to the segment
    [neighbour_seg, neighbour_m, min_dist] ...
        = closest_point(curr_opp_pt, seg_nos(check_pts == 1));

    if min_dist > tol
        % Neighbour not found, so break the loop
%         disp("Neighbour not found")
        break
    else
        % Saves the opposite point of the neighbour
        neighbour_opp_pt = point(neighbour_seg, -neighbour_m);
        if neighbour_opp_pt == -1
            break;
        end
        check_pts(neighbour_seg) = 0;
    end

    % Coarsen by setting the opposite point of the current segment to be
    % the opposite point of the neighbour segment
    if curr_m == 1
        end_points(curr_seg, :) = neighbour_opp_pt;
    else
        start_points(curr_seg, :) = neighbour_opp_pt;
    end
    
    % Remove the neighbour point from the start and end points by setting
    % them to the dummy value
    start_points(neighbour_seg, :) = [-1, -1];
    end_points(neighbour_seg, :) = [-1, -1];
    
    %% Find the next point
    % Find the closest point to the neighbour
    [next_seg, next_m, min_dist] ...
        = closest_point(neighbour_opp_pt, seg_nos(check_pts == 1));
    
    if min_dist > tol
        % Next segment not found, so break the loop
%         disp("Next not found with current segment ending at");
        curr_end_point = end_points(curr_seg, :);
        break
    elseif isempty(next_seg)
%         disp("Empty next seg");
        break;
    else
        % Sets the next current point to be this segment
        curr_seg = next_seg;
        curr_m = next_m;
    end
    
    
end

%% Remove the dummy points from the arrays
if remove_rest
    start_points = start_points(keep_pts == 1, :);
    end_points = end_points(keep_pts == 1, :);
end

new_start_points = start_points(start_points(:, 1) ~= -1, :);
new_end_points = end_points(end_points(:, 1) ~= -1, :);

%% Remove rest
% If the remove_rest option is set to be true, then all points that have
% not been found as a result of the coarsening are removed (this should in
% theory remove the trapped bubbles)


end


