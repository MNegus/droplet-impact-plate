%% data_analysis.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyses the data resulting from the uniform_plate simulations. It is
% assumed the data has been cleaned using the data cleaning utilities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

% Adds the interface analysis to path
addpath("~/repos/plate-impact/data_analysis/interface_analysis");

%% Definitions
% Here we specify the location of the data and other parameters required to
% conduct the analysis

% Velocity of the plate
plate_velocity = 0;

% Time between plate outputs
plate_output_timestep = 1e-3;

% Parent directory where all the data is stored
parent_directory = '/mnt/newarre';

% Directory names inside the master directory
dir_names = ["force_test"];

% Adds the parent directory to the start of 
for k = 1 : length(dir_names)
    dir_names(k) = sprintf('%s/%s', parent_directory, dir_names(k));
end

% Directory to the graphs and videos 
analysis_directory = sprintf('%s/data_analysis', dir_names(1));

% Readable names to label the plots for each of the data directories
legend_entries = ["Solid wall"];


%% Pressure along plate
% Creates an animation of the pressure along the plate in time.This data is 
% stored in files "output_n.txt". The first column is the radial 
% coordinate, r. The second coordinate is z, vertical position. The third 
% is pressure, the fourth is the vertical velocity u_z and the fifth is the
% radial velocity u_r. 

% Stationary plate turnover point
d = @(t) sqrt(3 * t);

% Set true to adjust the Wagner impact time
adjusted = false;

% Impact time
if adjusted
    impact_time = 0.142;
else
    impact_time = (2.125 - 1 - 1) / (1 + plate_velocity);
end

% Reads the "times.txt" file from the first data directory. In theory this
% should be identical in all the cases
times = dlmread(sprintf('%s/cleaned_data/plate_outputs/times.txt', dir_names(1)));

% Position to start video at
start_pos =  floor(0.9 * impact_time / plate_output_timestep);
% start_pos = 1;

% Number of video frames
no_frames = 300; 

% Maximum pressure at each timestep
pmax = zeros(no_frames, length(dir_names) + 2);

% Force at each timestep
integrated_force = zeros(no_frames, length(dir_names) + 1);

% Sets up figure
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

% Creates animated line for the Wagner pressure
wagner_line = animatedline('Color', [0    0.4470    0.7410], 'Linestyle', '--');
outer_wagner_line = animatedline('Color', 'red', 'Linestyle', '--');

% Creates animated lines for the numerical results
line1 = animatedline('Color', [0.8500    0.3250    0.0980]);
line2 = animatedline('Color', [0.9290    0.6940    0.1250]);
line3 = animatedline('Color', [0.4940    0.1840    0.5560]);
animlines = [line1, line2, line3];


% Sets up the legend
L = legend(['Wagner pressure', 'Outer Wagner pressure', legend_entries]);
set(L, 'Interpreter', 'latex');
set(L, 'FontSize', 15);

% Sets pixel size of the figure
width=800;
height=800;
set(gcf,'position',[10,10,width,height])

% Create the video writer with 5 fps
if adjusted
    writerObj = VideoWriter(sprintf('%s/Videos/pressure_adjusted.avi', analysis_directory));
else
    writerObj = VideoWriter(sprintf('%s/Videos/pressure_unadjusted', analysis_directory), 'Uncompressed AVI');
end
writerObj.FrameRate = 5;
open(writerObj);


% Iterates over time
for m = start_pos: start_pos + no_frames -1

    t = times(m, 2); % Time

    for k = 1 : length(dir_names)
        % Loads in data from the text file
        output_matrix = ...
            dlmread(sprintf('%s/cleaned_data/plate_outputs/output_%d.txt', ...
                dir_names(k), m));
        
        % Sorts in increasing order of r
        [~, sorted_idxs] = sort(output_matrix(:, 1));
        sorted_mat = output_matrix(sorted_idxs, :);

        % Saves values of r and pressure
        rs = sorted_mat(:, 1);
        rs(2) - rs(1)
        ps = sorted_mat(:, 3);

        % Creates the animated line

        % Adds the pressure line
        clearpoints(animlines(k))
        addpoints(animlines(k), rs, ps);


        % Saves max value of pressure
        pmax(m - start_pos + 1, 1) = t; % Time
        pmax(m - start_pos + 1, 1 + k) = max(ps); % Computational pressure
        
        % Calculates the force using trapezoidal rule
        integrated_force(m - start_pos + 1, 1) = t;
        integrated_force(m - start_pos + 1, 1 + k) = ...
            trapz(rs, 2 * pi * rs .* ps);
    end
    % Wagner line
    if t > impact_time
        sigmas =  10.^linspace(-10, 5, 1e4);
        [wagner_rs, wagner_ps, outer_wagner_ps, wagner_pmax] = wagner_pressure(sigmas, t - impact_time, plate_velocity, 1);

        clearpoints(wagner_line);
        addpoints(wagner_line, wagner_rs, wagner_ps);
        
        clearpoints(outer_wagner_line);
        addpoints(outer_wagner_line, wagner_rs, outer_wagner_ps);

        % Draws a vertical line where the turnover point is
%         delete(wagner_turnover_line);
%         wagner_turnover_line = xline(sqrt(3 * (t - impact_time) * (1 - plate_velocity)), '--r' );
%         L = legend([legend_entry, 'Wagner pressure', 'Wagner turnover', 'Computational turnover']);

        % Records maximum Wagner pressure
        [pmax(m - start_pos + 1, length(dir_names) + 2), idx] = max(wagner_ps);
        wagner_rs(idx)
    end

    % Resets axes limits
%     ylim([0 , 20]);
%     xlim([0, max([1, 1.5 * d(t - impact_time)])]);

    title(sprintf("Combined pressures, $t$ = %.4f, Plate vel = %g, Wagner impact, $t$ = %.4f\n",...
        times(m, 2), plate_velocity, impact_time), ...
        "Interpreter", "latex", 'Fontsize', 15);
    drawnow;
    frame = getframe(gcf);
    writeVideo(writerObj, frame);

end

close(writerObj);


% Creates plot of maximum pressure
figure(3);
hold on;
plot(pmax(:, 1), pmax(:, length(dir_names) + 2), '--'); % Wagner pressure
for k = 1 : length(dir_names)
    plot(pmax(:, 1), pmax(:, 1 + k)); % Computational pressure
end


ylim([0, max(1.5 * pmax(:, 1 + length(dir_names)))]);
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

if adjusted
    print(gcf, sprintf('%s/Figures/pmax_adjusted.png', analysis_directory),'-dpng','-r300');
else
    print(gcf, sprintf('%s/Figures/pmax_unadjusted.png', analysis_directory),'-dpng','-r300');
end



%% Force on plate
% Plots the force on the plate as read from the cleaned log file
log_matrix = dlmread(sprintf('%s/cleaned_data/volumes.txt', dir_names(1)));
ts = log_matrix(:, 1);
Fs = log_matrix(:, 3);

% Force according to stationary plate Wagner theory
wagner_F = @(t) 6 * sqrt(3 * t);
wagner_ts = ts(ts > impact_time);

% Plot of force
figure(4);
hold on
% plot(integrated_force(:, 1), integrated_force(:, 2));

plot(wagner_ts, wagner_F(wagner_ts - impact_time));
plot(ts, Fs);


xlabel("t");
ylabel("Force");
title("Force on plate");
legend(["Wagner force", "Computational force"]);
print(gcf, sprintf('%s/Figures/force_comparison.png', analysis_directory),'-dpng','-r300');


%%
%