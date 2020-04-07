%% data_analysis.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyses the data resulting from the uniform_plate simulations. It is
% assumed the data has been cleaned using the data cleaning utilities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

%% Definitions
% Here we specify the location of the data and other parameters required to
% conduct the analysis

% Adds the interface analysis to path
addpath("~/repos/plate-impact/data_analysis/interface_analysis");

% Name of run
run_name = "Level 12";

% Master directory where all the data is stored
master_directory = '/mnt/newarre/level_12';

% Directory to store the outputs of the data analysis
analysis_directory = sprintf("%s/data_analysis", master_directory);

% Names of the individual directories where the data is stored
plate_velocity = -0.1;
% if plate_velocity < 0.1
%     data_directory = sprintf("%s/plate_vel_%.2f", master_directory, plate_velocity)
% else
%     data_directory = sprintf("%s/plate_vel_%.1f", master_directory, plate_velocity)
% end
cleaned_data_directory = sprintf("%s/cleaned_data", master_directory);

% Readable names to label the plots for each of the data directories
legend_entry = sprintf("Plate velocity = %.2f", plate_velocity);

%% Force on plate
% Plots the force on the plate as read from the cleaned log file
log_matrix = dlmread(sprintf('%s/volumes.txt', cleaned_data_directory));
ts = log_matrix(:, 1);
Fs = log_matrix(:, 3);
plot(ts, Fs);

%% Pressure along plate
% Creates an animation of the pressure along the plate in time.This data is 
% stored in files "output_n.txt". The first column is the radial 
% coordinate, r. The second coordinate is z, vertical position. The third 
% is pressure, the fourth is the vertical velocity u_z and the fifth is the
% radial velocity u_r. 

% Stationary plate turnover point
d = @(t) sqrt(3 * t);

% Impact time
impact_time = (2.125 - 1 - 1) / (1 + plate_velocity);

% Measured turnover point
comp_turnover_pts = dlmread(sprintf('%s/turnover_points.txt', analysis_directory));

% Reads the "times.txt" file from the first data directory. In theory this
% should be identical in all the cases
times = dlmread(strcat(cleaned_data_directory, '/plate_outputs/times.txt'));

% Position to start video at
start_pos =  floor(0.9 * impact_time * 1000)

no_frames = 100; % Number of video frames

% Maximum pressure at each timestep
pmax = zeros(no_frames, 3);

% Turn to false to not show graph
show_graph = true;


if show_graph == true 
    % Sets up figure
    close(figure(2));
    figure(2);
    
    % Vertical lines
    comp_turnover_line = xline(0);
    wagner_turnover_line = xline(0);
    
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

    % Create animated lines and place them in a matrix
    animline = animatedline('Color', [0    0.4470    0.7410]);

    wagner_line = animatedline('Color', [0.8500    0.3250    0.0980]);

     % Sets up the legend
    L = legend([legend_entry, 'Wagner pressure', 'Wagner turnover', 'Computational turnover']);
    set(L, 'Interpreter', 'latex');
    set(L, 'FontSize', 15);

    width=800;
    height=800;
    set(gcf,'position',[10,10,width,height])
    % create the video writer with 1 fps
    writerObj = VideoWriter(sprintf('%s/Videos/pressure_vel_%.2f.avi', analysis_directory, plate_velocity));
    writerObj.FrameRate = 5;
    open(writerObj);
end

% Iterates over time
for m = start_pos: start_pos + no_frames
    % Choses which data directories to show

    t = 0.001 * m; % Time

    % Loads in data from the text file
    output_matrix = dlmread(strcat(cleaned_data_directory, ...
        '/plate_outputs/output_', num2str(m), '.txt'));

    % Sorts in increasing order of r
    [~, sorted_idxs] = sort(output_matrix(:, 1));
    sorted_mat = output_matrix(sorted_idxs, :);

    % Saves values of r and pressure
    rs = sorted_mat(:, 1);
    ps = sorted_mat(:, 3);

    % Creates the animated line
    if show_graph == true
        % Adds the pressure line
        clearpoints(animline)
        addpoints(animline, rs, ps);
        
        % Adds the turnover point line
        delete(comp_turnover_line);
        comp_turnover_line = xline(comp_turnover_pts(m, 2), '--b');
    end

    % Saves max value of pressure
    pmax(m - start_pos + 1, 1) = t; % Time
    pmax(m - start_pos + 1, 2) = max(ps); % Computational pressure

    % Wagner line
    if show_graph == true
        if t > impact_time
            sigmas =  10.^linspace(-10, 5, 1e4);
            [wagner_rs, wagner_ps, wagner_pmax] = wagner_pressure(sigmas, t - impact_time, plate_velocity, 1);

            clearpoints(wagner_line);
            addpoints(wagner_line, wagner_rs, wagner_ps);

            % Draws a vertical line where the turnover point is
            delete(wagner_turnover_line);
            wagner_turnover_line = xline(sqrt(3 * (t - impact_time) * (1 - plate_velocity)), '--r' );
            L = legend([legend_entry, 'Wagner pressure', 'Wagner turnover', 'Computational turnover']);
            
            % Records maximum Wagner pressure
            [pmax(m - start_pos + 1, 3), idx] = max(wagner_ps);
            wagner_rs(idx)
        end
    end

    if show_graph == true
        ylim([0 , 20]);
        xlim([0, max([1, 1.5 * d(t - impact_time)])]);

        title(sprintf("%s, $t$ = %.3f, Plate vel = %g, Wagner impact = %.3f\n",...
            run_name, times(m, 2), plate_velocity, impact_time), ...
            "Interpreter", "latex", 'Fontsize', 15);
        drawnow;
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
    end
end
if show_graph == true
    close(writerObj);
end

% Creates plot of maximum pressure
figure(3);
hold on;
plot(pmax(:, 1), pmax(:, 2)); % Computational pressure
plot(pmax(:, 1), pmax(:, 3)); % Wagner pressure
ylim([0, max(1.5 * pmax(:, 2))]);
grid on;
xlabel("$t$", "Interpreter", "latex", 'Fontsize',30);
ylabel("Max pressure, $p$", "Interpreter", "latex", 'Fontsize', 30);
ax = gca;
ax.FontSize = 10;
set(gca,'TickLabelInterpreter','latex');
title(sprintf("%s. Max pressure for plate velocity = %g", run_name, ...
    plate_velocity), "Interpreter", "latex", 'Fontsize', 15);

% xticks(linspace(min(pmax(:, 1)), max(pmax(:, 1)), 12))
% xtickformat('%.2f')
L = legend(["Computational", "Wagner"]);
set(L, 'Interpreter', 'latex');
set(L, 'FontSize', 11);
set(L, 'Location', 'northeast');

print(gcf, sprintf('%s/Figures/pmax_vel_%g.png', analysis_directory, plate_velocity),'-dpng','-r300');
%     close(figure(3));








%%
%