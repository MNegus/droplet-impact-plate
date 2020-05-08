%% pressure_analysis.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyses the data resulting from the uniform_plate simulations. It is
% assumed the data has been cleaned using the data cleaning utilities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;


%% Data definitions
% Here we specify the locations where the plate output files are stored. We
% expect a different directory for each simulation result.

% Directory where the resulting videos are to be stored
results_directory = "/mnt/newarre/analysis_tests";

% Master directory where all of the data is stored under (e.g. external
% hard drive location)
master_directory = "/mnt/newarre/analysis_tests/";

% Individual directory names under the master directory. We assume that the
% plate output files are stored under
% master_directory/data_directory(k)/cleaned_data
data_directories = ["wall_test", "stationary_prescribed", "coupled_test"];
no_dirs = length(data_directories); % Number of entries

% Adds the parent directory to the start of 
for k = 1 : length(data_directories)
    data_directories(k) = strcat(master_directory, data_directories(k)); 
end

% Readable names to label the plots for each of the data directories
legend_entries = ...
    ["Wall", "Embedded", "$\alpha = 0.01, \gamma = 1$"];

%% Parameters
% Parameters common to the simulations to aid in visualisation
initial_drop_height = 0.125; 

%% Pressure along plate
% Creates an animation of the pressure along the plate in time.This data is 
% stored in files "output_n.txt". The first column is the radial 
% coordinate, r. The second coordinate is z, vertical position. The third 
% is pressure, the fourth is the vertical velocity u_z and the fifth is the
% radial velocity u_r. 

% Stationary plate turnover point
d = @(t) sqrt(3 * t);

% Theoretical time of impact for a stationary plate
impact_time = initial_drop_height;

% Reads the "times.txt" file from the first data directory. In theory this
% should be identical in all the cases
times = dlmread(sprintf('%s/cleaned_data/plate_outputs/times.txt', data_directories(1)));

% Position to start video at
start_pos =  115;
end_pos = 400;
no_frames = end_pos - start_pos;

% Maximum pressure at each timestep
pmax = zeros(no_frames, length(data_directories) + 2);

% Sets up figure
figure(1);
hold on;
grid on;
xlabel("$r$", "Interpreter", "latex", 'Fontsize',30);
ylabel("Pressure, $p$", "Interpreter", "latex", 'Fontsize', 30);
ax = gca;
ax.FontSize = 16;
set(gca,'TickLabelInterpreter','latex');
x_limits = [0, 1];
xlim(x_limits);
ylim([-0.1 15]);

% Creates animated line for the Wagner pressure
wagner_line = animatedline('Color', [0    0.4470    0.7410], ...
    'Linestyle', '--', 'Linewidth', 1.5);

% MAKE THIS SMARTER
% Creates animated lines for the numerical results
line1 = animatedline('Color', [0.8500    0.3250    0.0980], 'Linewidth', 1.5);
line2 = animatedline('Color', [0.9290    0.6940    0.1250], 'Linewidth', 1.5);
line3 = animatedline('Color', [0.4940    0.1840    0.5560], 'Linewidth', 1.5);
animlines = [line1, line2, line3];

% Sets up the legend
L = legend(['Stationary Wagner pressure', legend_entries]);
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
        % Loads in data from the text file
        output_matrix = ...
            dlmread(...
                sprintf('%s/cleaned_data/plate_outputs/output_%d.txt', ...
                    data_directories(k), m));
        
        % Sorts in increasing order of r
        [~, sorted_idxs] = sort(output_matrix(:, 1));
        sorted_mat = output_matrix(sorted_idxs, :);

        % Saves values of r and pressure
        rs = sorted_mat(:, 1);
        ps = sorted_mat(:, 3);

        % Adds the pressure line
        clearpoints(animlines(k))
        addpoints(animlines(k), rs, ps);

        % Saves max value of pressure
        pmax(m - start_pos + 1, 1) = t; % Time
        pmax(m - start_pos + 1, 1 + k) = max(ps); % Computational pressure

    end
    % Wagner line
    if t > impact_time
        sigmas =  10.^linspace(-10, 5, 1e4);
        [wagner_rs, wagner_ps, outer_wagner_ps, wagner_pmax] ...
            = wagner_pressure(sigmas, t - impact_time, 0, 1);

        clearpoints(wagner_line);
        addpoints(wagner_line, wagner_rs, wagner_ps);

        % Draws a vertical line where the turnover point is
%         delete(wagner_turnover_line);
%         wagner_turnover_line = xline(sqrt(3 * (t - impact_time) * (1 - plate_velocity)), '--r' );
%         L = legend([legend_entry, 'Wagner pressure', 'Wagner turnover', 'Computational turnover']);

        % Records maximum Wagner pressure
        [pmax(m - start_pos + 1, length(data_directories) + 2), idx] = max(wagner_ps);
        wagner_rs(idx)
    end

    title(sprintf("Pressure, $t$ = %.3f, Theoretical impact time, $t$ = %.3f\n",...
        times(m, 2), impact_time), ...
        "Interpreter", "latex", 'Fontsize', 15);
    drawnow;
    frame = getframe(gcf);
    writeVideo(writerObj, frame);

end

close(writerObj);


% Creates plot of maximum pressure
figure(3);
hold on;
plot(pmax(:, 1), pmax(:, length(data_directories) + 2), '--'); % Wagner pressure
for k = 1 : length(data_directories)
    plot(pmax(:, 1), pmax(:, 1 + k)); % Computational pressure
end


ylim([0, max(1.5 * pmax(:, 1 + length(data_directories)))]);
grid on;
xlabel("$t$", "Interpreter", "latex", 'Fontsize',30);
ylabel("Max pressure, $p$", "Interpreter", "latex", 'Fontsize', 30);
ax = gca;
ax.FontSize = 10;
set(gca,'TickLabelInterpreter','latex');
title(sprintf("Combined max pressure, Wagner impact $t$ = %.3f", ...
    impact_time), "Interpreter", "latex", 'Fontsize', 15);

% xticks(linspace(min(pmax(:, 1)), max(pmax(:, 1)), 12))
% xtickformat('%.2f')
L = legend(["Wagner", legend_entries]);
set(L, 'Interpreter', 'latex');
set(L, 'FontSize', 11);
set(L, 'Location', 'northeast');

print(gcf, sprintf('%s/pmax.png', results_directory),'-dpng','-r300');