%% data_analysis.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyses the data resulting from the Plate_validation simulations. It is
% assumed the data has been cleaned using the data_clean.sh utility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Definitions
% Here we specify the location of the data and other parameters required to
% conduct the analysis

% Adds the interface analysis to path
addpath("~/repos/plate-impact/data_analysis/interface_analysis");

% Master directory where all the data is stored
master_directory = '/scratch/Uniform_plate/cleaned_data/';

% Names of the individual directories where the data is stored
data_directories = ["uniform_plate"];

% Edits so the full directory address is given in data_directories
for k = 1:length(data_directories)
    data_directories(k) = strcat(master_directory, data_directories(k));
end

% Readable names to label the plots for each of the data directories
legend_entries = ["Plate"];

%% Volume conservation
% Plots the volume of the droplet volume fraction as a function of time for
% all of the simulation to check how well volume (or equivalently mass) is
% conserved

initial_vol = 4 * pi / 3; % Initial volume of the droplet

% Sets up figure
close(figure(1));
figure(1);
hold on;
grid on;
xlabel("$t$", "Interpreter", "latex", 'Fontsize', 15);
ylabel("Volume / Initial volume", "Interpreter", "latex", 'Fontsize', 15);
ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title("Volume conservation of droplet phase", "Interpreter", "latex", ...
    'Fontsize', 15);

% Iterates through the data directories. The "volumes.txt" files contain
% two columns, the first is time and the second is volume
for k = 1:length(data_directories)
    data_matrix = dlmread(strcat(data_directories(k), '/volumes.txt'));
    plot(data_matrix(:, 1), data_matrix(:, 2) / initial_vol);
end
% Sets up the legend
L = legend(legend_entries);
set(L, 'Interpreter', 'latex');
set(L, 'FontSize', 11);
set(L, 'Location', 'northwest');
print(gcf,'Figures/volume.png','-dpng','-r300');

%% Pressure along plate
% Creates an animation of the pressure along the plate in time. We can vary
% h from 0 to 2, where h is the number of cells above the plate the
% pressure is measured. This data is stored in files "output_n.txt". The
% first column is the radial coordinate, r. The second coordinate is h, the
% number of cells above the plate. The third is pressure, the fourth is the
% vertical velocity u_z and the fifth is the radial velocity u_r. 

% Turnover point
d = @(t) sqrt(3 * t);

% Wagner pressure
wagner_p = @(t, r) 3 ./ (pi * sqrt(3 * t - r.^2));

% Impact time
impact_time = 0.203;

h = 0; % Number of cells above plate

% Sets up figure
close(figure(2));
figure(2);
hold on;
grid on;
xlabel("$r$", "Interpreter", "latex", 'Fontsize',30);
ylabel("Pressure, $p$", "Interpreter", "latex", 'Fontsize', 30);
ax = gca;
ax.FontSize = 16;
set(gca,'TickLabelInterpreter','latex');
x_limits = [0, 1];
xlim(x_limits);
ylim([0 15]);

% Reads the "times.txt" file from the first data directory. In theory this
% should be identical in all the cases
times = dlmread(strcat(data_directories(1), '/plate_outputs/times.txt'));

% Create animated lines and place them in a matrix
line1 = animatedline('Color', [0    0.4470    0.7410]);
line2 = animatedline('Color', [0.8500    0.3250    0.0980]);
line3 = animatedline('Color', [0.9290    0.6940    0.1250]);
line4 = animatedline('Color', [0.4940    0.1840    0.5560]);
animlines = [line1, line2, line3, line4];

wagner_line = animatedline('Color', [0.4660    0.6740    0.1880]);


 % Sets up the legend
dir_arr = [1];
L = legend([legend_entries(dir_arr), 'Wagner']);
set(L, 'Interpreter', 'latex');
set(L, 'FontSize', 15);

width=800;
height=800;
set(gcf,'position',[10,10,width,height])
% create the video writer with 1 fps
writerObj = VideoWriter(sprintf('Videos/pressure_%d.avi', h));
writerObj.FrameRate = 5;
open(writerObj);
% Iterates over all times
for m = 205:750
    % Choses which data directories to show
    
    t = 0.001 * m; % Time
    
    for k = dir_arr
        % Loads in data from the text file
        output_matrix = dlmread(strcat(data_directories(k), ...
            '/plate_outputs/output_', num2str(m), '.txt'));
    
        % Isolates data with matching value of h
        h_column = output_matrix(:, 2); 
        refined_mat = output_matrix(h_column == h, :);

        % Sorts in increasing order of r
        [~, sorted_idxs] = sort(refined_mat(:, 1));
        sorted_mat = refined_mat(sorted_idxs, :);

        % Creates the animated line
        clearpoints(animlines(k))
        addpoints(animlines(k), sorted_mat(:, 1), sorted_mat(:, 3));
    end
   
    % Wagner line
    clearpoints(wagner_line);
    if t > impact_time
        wagner_rs = linspace(x_limits(1), 0.99 * d(t - impact_time), 1e3);
        wagner_vals = wagner_p(t - impact_time, wagner_rs);
        addpoints(wagner_line, wagner_rs, wagner_p(t - impact_time, ...
            wagner_rs));
    end
    ylim([0 ,max([max(wagner_vals), 10])]);
    xlim([0, max([1, 1.5 * d(t - impact_time)])]);
    
    title(sprintf("Pressure: $h$ = %d, $t$ = %g", h, times(m, 2)), ...
        "Interpreter", "latex", 'Fontsize', 15);
    drawnow;
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
end

close(writerObj);

%%
%