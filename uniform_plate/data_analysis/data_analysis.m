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
plate_velocity= 0.1;
data_directory = sprintf("/scratch/Uniform_plate/cleaned_data/plate_vel_%.1f", plate_velocity)

% Readable names to label the plots for each of the data directories
legend_entry = sprintf("Plate velocity = %.1f", plate_velocity);

%% Volume conservation
% Plots the volume of the droplet volume fraction as a function of time for
% all of the simulation to check how well volume (or equivalently mass) is
% conserved
% 
% initial_vol = 4 * pi / 3; % Initial volume of the droplet
% 
% % Sets up figure
% close(figure(1));
% figure(1);
% hold on;
% grid on;
% xlabel("$t$", "Interpreter", "latex", 'Fontsize', 15);
% ylabel("Volume / Initial volume", "Interpreter", "latex", 'Fontsize', 15);
% ax = gca;
% ax.FontSize = 12;
% set(gca,'TickLabelInterpreter','latex');
% title("Volume conservation of droplet phase", "Interpreter", "latex", ...
%     'Fontsize', 15);
% 
% % Iterates through the data directories. The "volumes.txt" files contain
% % two columns, the first is time and the second is volume
% for k = 1:length(data_directories)
%     data_matrix = dlmread(strcat(data_directories(k), '/volumes.txt'));
%     plot(data_matrix(:, 1), data_matrix(:, 2) / initial_vol);
% end
% % Sets up the legend
% L = legend(legend_entries);
% set(L, 'Interpreter', 'latex');
% set(L, 'FontSize', 11);
% set(L, 'Location', 'northwest');
% print(gcf,'Figures/volume.png','-dpng','-r300');

%% Pressure along plate
% Creates an animation of the pressure along the plate in time. We can vary
% h from 0 to 2, where h is the number of cells above the plate the
% pressure is measured. This data is stored in files "output_n.txt". The
% first column is the radial coordinate, r. The second coordinate is h, the
% number of cells above the plate. The third is pressure, the fourth is the
% vertical velocity u_z and the fifth is the radial velocity u_r. 

% Turnover point
d = @(t) sqrt(3 * t);

% Impact time
impact_time = (2.125 - 1 - 1) / (1 - plate_velocity)


for h = [0, 1, 2]
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
    times = dlmread(strcat(data_directory, '/plate_outputs/times.txt'));

    % Create animated lines and place them in a matrix
    animline = animatedline('Color', [0    0.4470    0.7410]);

    wagner_line = animatedline('Color', [0.8500    0.3250    0.0980]);


     % Sets up the legend
    L = legend([legend_entry, 'Wagner']);
    set(L, 'Interpreter', 'latex');
    set(L, 'FontSize', 15);

    width=800;
    height=800;
    set(gcf,'position',[10,10,width,height])
    % create the video writer with 1 fps
    writerObj = VideoWriter(sprintf('Videos/pressure_h_%d_vel_%.1f.avi', h, plate_velocity));
    writerObj.FrameRate = 5;
    open(writerObj);

    % Iterates over all times
    start_pos = floor(0.9 * impact_time * 1000);
    for m = start_pos: start_pos + 100
        % Choses which data directories to show

        t = 0.001 * m; % Time

        % Loads in data from the text file
        output_matrix = dlmread(strcat(data_directory, ...
            '/plate_outputs/output_', num2str(m), '.txt'));

        % Isolates data with matching value of h
        h_column = output_matrix(:, 2); 
        refined_mat = output_matrix(h_column == h, :);

        % Sorts in increasing order of r
        [~, sorted_idxs] = sort(refined_mat(:, 1));
        sorted_mat = refined_mat(sorted_idxs, :);

        % Creates the animated line
        clearpoints(animline)
        addpoints(animline, sorted_mat(:, 1), sorted_mat(:, 3));

        % Wagner line
        clearpoints(wagner_line);
        if t > impact_time
            sigmas =  10.^linspace(-10, 5, 1e4);
            [wagner_rs, wagner_ps] = wagner_pressure(sigmas, t - impact_time, plate_velocity, 1);
            addpoints(wagner_line, wagner_rs, wagner_ps);
        end


        ylim([0 , 20]);
        xlim([0, max([1, 1.5 * d(t - impact_time)])]);

        title(sprintf("Pressure: $h$ = %d, $t$ = %g, Plate vel = %g\n", h, times(m, 2), plate_velocity), ...
            "Interpreter", "latex", 'Fontsize', 15);
        drawnow;
        frame = getframe(gcf);
        writeVideo(writerObj, frame);
    end

    close(writerObj);
    
end






%%
%