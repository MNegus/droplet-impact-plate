%% pressure_comparison.m
% Script to compare the pressures on the plate for different plate widths,
% plus comparing to the analytical solution
%



% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions
% Here we specify the locations where the plate output files are stored. We
% expect a different directory for each simulation result.

% Parent directory where all the data is held
parent_directory = '/mnt/newarre/cantilever_paper_data/validation/box_width_validation/';


% Directory to save the figure(s)
results_directory = "Analysis";
results_directory = strcat(parent_directory, results_directory);

% Defines array with both directories in
data_directories = ["box_width_3", "box_width_6", "box_width_12"];

% Concatenates arrays to include parent directory
for k = 1 : length(data_directories)
    data_directories(k) = strcat(parent_directory, data_directories(k)); 
end

legend_entries = ["Box width = 3", "Box width = 6", "Box width = 12"];

%% Parameters
% Value of epsilon
eps = 1;

% Plate parameters
alpha = 1;
beta = 0;
gamma = 0;

% Initial drop height
initial_drop_height = 0.125;

% Impact time
impact_time = initial_drop_height;

% Maximum time
t_max = 0.5;

%% Wagner solution

% Plate displacement solution
[wagner_t, s_vals, sdot_vals, sddot_vals] = ...
    s_solution(t_max - impact_time, alpha, beta, gamma, eps);

% Interpolated functions for s and its derivatives
s = @(t) interp1(wagner_t, s_vals, t);
sdot = @(t) interp1(wagner_t, sdot_vals, t);
sddot = @(t) interp1(wagner_t, sddot_vals, t);

%% Pressure along plate
% Creates an animation of the pressure along the plate in time.This data is 
% stored in files "output_n.txt". The first column is the radial 
% coordinate, r. The second coordinate is z, vertical position. The third 
% is pressure, the fourth is the vertical velocity u_z and the fifth is the
% radial velocity u_r. 


% Reads the "times.txt" file from the first data directory. In theory this
% should be identical in all the cases
times = dlmread(sprintf('%s/cleaned_data/plate_outputs/times.txt', data_directories(3)));

% Position to start video at
start_pos =  125;
end_pos = 500;

output_range = start_pos : end_pos;
no_frames = end_pos - start_pos;


% Sets up figure
close all;
figure(1);
hold on;
grid on;
xlabel("$r$", "Interpreter", "latex", 'Fontsize',30);
ylabel("Pressure, $p$", "Interpreter", "latex", 'Fontsize', 30);
ax = gca;
ax.FontSize = 16;
set(gca,'TickLabelInterpreter','latex');
x_limits = [0, 2];
xlim(x_limits);
ylim([-1 1]);

% Creates animated line for the Wagner pressure
wagner_line = animatedline('color', [0 0 0], 'Linewidth', 1.5);

% Wagner turnover line
wagner_turnover_line = animatedline('color', 'red', 'Linestyle', '--', ...
    'LineWidth', 1.5);


% MAKE THIS SMARTER
% Creates animated lines for the numerical results
line1 = animatedline('Color', [0, 0.4470, 0.7410], 'Linewidth', 1.5);
line2 = animatedline('Color', [0.8500, 0.3250, 0.0980], 'Linewidth', 1.5);
line3 = animatedline('Color', [0.9290, 0.6940, 0.1250], 'Linewidth', 1.5);
% line4 = animatedline('Color', [0.4940, 0.1840, 0.5560], 'Linewidth', 1.5);
% line5 = animatedline('Color', [0.4660, 0.6740, 0.1880], 'Linewidth', 1.5);
animlines = [line1, line2, line3];

% Sets up the legend
L = legend(["Composite", "Wagner turnover point",legend_entries]);
set(L, 'Interpreter', 'latex');
set(L, 'FontSize', 15);

% Sets pixel size of the figure
width=800;
height=800;
set(gcf,'position',[10,10,width,height])

% Create the video writer with 5 fps
writerObj ...
    = VideoWriter(sprintf('%s/pressure_overlay.avi', results_directory));
writerObj.FrameRate = 5;
open(writerObj);

% Iterates over time
for m = start_pos : start_pos + no_frames -1

    t = times(m, 2); % Time

    for k = 1 : length(data_directories)
        if k == 1
            m = 10 * m;
        end
        
        % Loads in data from the text file
        output_matrix = ...
            dlmread(...
                sprintf('%s/cleaned_data/plate_outputs/output_%d.txt', ...
                    data_directories(k), m));
        if k == 1
            m = m / 10;
        end
        
        % Sorts in increasing order of r
        [~, sorted_idxs] = sort(output_matrix(:, 1));
        sorted_mat = output_matrix(sorted_idxs, :);

        % Saves values of r and pressure
        rs = sorted_mat(:, 1);
        ps = sorted_mat(:, 4);
        
        % Adds the pressure line
        clearpoints(animlines(k))
        addpoints(animlines(k), rs, ps);

    end
    % Wagner line
    if t > impact_time
        
        % Values of s
        s_val = s(t - impact_time);
        sdot_val = sdot(t - impact_time);
        sddot_val = sddot(t - impact_time);
        
        % Saves turnover point
        [d, ~, ~, ~] = s_dependents(t - impact_time, s_val, sdot_val, ...
            sddot_val);
        
        % Determines outer and composite solution
        [outer_rs, outer_ps, comp_rs, comp_ps] ...
            = outer_and_comp_pressure(t - impact_time, s_val, sdot_val, ...
                sddot_val, 0, 1.25 * d, 1);
        
        
        % Adds composite pressure to graph
%         clearpoints(wagner_line);
%         addpoints(wagner_line, comp_rs, comp_ps);

        % Draws a vertical line where Wagner theory turnover point is
        clearpoints(wagner_turnover_line);
        addpoints(wagner_turnover_line, d * [1 1], [0 100]);

    end

    title(sprintf("Pressure: Accelerating frame. $t$ = %.4f, $t_0$ = %.4f",...
        times(m, 2), impact_time), ...
        "Interpreter", "latex", 'Fontsize', 15);
    xlim(x_limits);
%     ylim([-2 15]);
    drawnow;
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
end

close(writerObj);
