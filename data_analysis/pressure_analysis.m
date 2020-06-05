%% pressure_analysis.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyses the data resulting from the uniform_plate simulations. It is
% assumed the data has been cleaned using the data cleaning utilities.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

addpath('pressures');

%% Data definitions
% Here we specify the locations where the plate output files are stored. We
% expect a different directory for each simulation result.

% Parent directory where all of the data is stored under (e.g. external
% hard drive location)
parent_directory = "/mnt/newarre/a_1_test/";

% Directory where the resulting videos are to be stored
results_directory = sprintf("%s/Analysis", parent_directory);

% Individual directory names under the master directory. We assume that the
% plate output files are stored under
% master_directory/data_directory(k)/cleaned_data
% data_directories = ["a_0", "a_0.25", "a_0.5", "a_0.75", "a_1"];
data_directories = ["level_10", "level_11", "level_12", "level_13"];
% data_directories = ["bubble_attempt"];
no_dirs = length(data_directories); % Number of entries


% Adds the parent directory to the start of the data directories
for k = 1 : length(data_directories)
    data_directories(k) = strcat(parent_directory, data_directories(k)); 
end

% Readable names to label the plots for each of the data directories
% legend_entries = ["$a$ = 0", "$a$ = 0.25", "$a$ = 0.5", "$a$ = 0.75", ...
%     "$a$ = 1"];
legend_entries = ["Level 10", "Level 11", "Level 12", "Level 13"];

%% Parameters
% Physical parameters
rho_w = 998;
R0 = 1e-3;
U0 = 10;
T0 = R0 / U0;
Patm = 10^5;

p_dim = @(p) Patm + rho_w * U0^2 * p;
r_millimetre = @(r) R0 * r * 1000;
t_millisecond = @(t) T0 * t * 1000;

% (Constant) acceleration of the plate
a = 1.0;

% Displacement of s
s = @(t) 0.5 * a * t.^2;
sdot = @(t) a * t;
sddot = @(t) a;

% Parameters common to the simulations to aid in visualisation
initial_drop_height = 0.125; 

% If true, then pressures are cutoff at r = r_cutoff
cutoff = false;
% r_cutoff = @(t) 1.1 * d(t);

%% Pressure along plate
% Creates an animation of the pressure along the plate in time.This data is 
% stored in files "output_n.txt". The first column is the radial 
% coordinate, r. The second coordinate is z, vertical position. The third 
% is pressure, the fourth is the vertical velocity u_z and the fifth is the
% radial velocity u_r. 

% Theoretical time of impact for a stationary plate
impact_time = initial_drop_height;

% Reads the "times.txt" file from the first data directory. In theory this
% should be identical in all the cases
times = dlmread(sprintf('%s/cleaned_data/plate_outputs/times.txt', data_directories(1)));

% Position to start video at
start_pos =  115;
end_pos = 600;

output_range = start_pos : end_pos;
no_frames = end_pos - start_pos;

% Finds computational turnover points
read_comp_ds = dlmread(sprintf('%s/cleaned_data/turnover_points.txt', ...
    data_directories(end)));
comp_ds = read_comp_ds(:, 2:3);

% start_turnover_idx = find(read_comp_ds(:, 1) == start_pos);
% end_turnover_idx = find(read_comp_ds(:, 1) == end_pos);
% comp_ds = read_comp_ds(start_turnover_idx : end_turnover_idx, 2:3);
% comp_ds

% Maximum pressure at each timestep
pmax = zeros(no_frames, length(data_directories) + 2);

% Sets up figure
figure(1);
hold on;
grid on;
xlabel("$r$ / mm", "Interpreter", "latex", 'Fontsize',30);
ylabel("Pressure, $p$ / Pa", "Interpreter", "latex", 'Fontsize', 30);
ax = gca;
ax.FontSize = 16;
set(gca,'TickLabelInterpreter','latex');
x_limits = [0, 1.25];
xlim(r_millimetre(x_limits));
ylim(p_dim([-1 15]));

% Creates animated line for the Wagner pressure
outer_wagner_line = animatedline('color', 0.5 * [1 1 1], ...
    'Linestyle', '--', 'Linewidth', 1.5);
comp_wagner_line = animatedline('color', [0 0 0], 'Linewidth', 1.5);

% Wagner turnover line
wagner_turnover_line = animatedline('color', 'red', 'Linestyle', '--', ...
    'LineWidth', 1.5);

% Computational turnover line
comp_turnover_line = animatedline('color', 'blue', 'Linestyle', '--', ...
    'LineWidth', 1.5);

% MAKE THIS SMARTER
% Creates animated lines for the numerical results
line1 = animatedline('Color', [0, 0.4470, 0.7410], 'Linewidth', 1.5);
line2 = animatedline('Color', [0.8500, 0.3250, 0.0980], 'Linewidth', 1.5);
line3 = animatedline('Color', [0.9290, 0.6940, 0.1250], 'Linewidth', 1.5);
line4 = animatedline('Color', [0.4940, 0.1840, 0.5560], 'Linewidth', 1.5);
% line5 = animatedline('Color', [0.4660, 0.6740, 0.1880], 'Linewidth', 1.5);
animlines = [line1, line2, line3, line4];

% Sets up the legend
L = legend(["Outer", "Composite", "Wagner turnover point", ...
    "Computational turnover point",legend_entries]);
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

%%
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

        if (t >= impact_time) && cutoff
            % Configures cutoff r
            if comp_ds(m, 1) < 1e-3
                r_cutoff = d(t);
            else
                r_cutoff = comp_ds(m, 1);
            end
            ps = ps(rs <= r_cutoff);
            rs = rs(rs <= r_cutoff);
        end
        
        % Adds the pressure line
        clearpoints(animlines(k))
        addpoints(animlines(k), r_millimetre(rs), p_dim(ps));

        % Saves max value of pressure
        pmax(m - start_pos + 1, 1) = t_millisecond(t); % Time
        pmax(m - start_pos + 1, 1 + k) = p_dim(max(ps)); % Computational pressure

    end
    % Wagner line
    if t >= impact_time
        
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
        
        % Adds outer pressure to graph
        clearpoints(outer_wagner_line);
        addpoints(outer_wagner_line, r_millimetre(outer_rs), ...
            p_dim(outer_ps));
        
        % Cuts off composite pressure
        if cutoff
            comp_ps = comp_ps(comp_rs <= r_cutoff);
            comp_rs = comp_rs(comp_rs <= r_cutoff);
        end
        
        % Adds composite pressure to graph
        clearpoints(comp_wagner_line);
        addpoints(comp_wagner_line, r_millimetre(comp_rs), p_dim(comp_ps));

        % Draws a vertical line where Wagner theory turnover point is
        clearpoints(wagner_turnover_line);
        addpoints(wagner_turnover_line, ...
            r_millimetre(d * [1 1]), ...
            p_dim([0 100]));
        
        % Plots the computational turnover point
        clearpoints(comp_turnover_line);
        addpoints(comp_turnover_line, ...
            r_millimetre(comp_ds(m, 1)) * [1 1], p_dim([0 100]));
        
        % Records maximum Wagner pressure
        [pmax(m - start_pos + 1, length(data_directories) + 2), idx] ...
            = max(p_dim(comp_ps));
    end

    title(sprintf("Pressure: Accelerating frame. $t$ = %.4f ms, $t_0$ = %.4f ms",...
        t_millisecond(times(m, 2)), t_millisecond(impact_time)), ...
        "Interpreter", "latex", 'Fontsize', 15);
    xlim(r_millimetre(x_limits));
    ylim(p_dim([-2 15]));
    drawnow;
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
end

close(writerObj);


% % Creates plot of maximum pressure
% figure(3);
% hold on;
% plot(pmax(:, 1), pmax(:, length(data_directories) + 2), '--'); % Wagner pressure
% for k = 1 : length(data_directories)
%     plot(pmax(:, 1), pmax(:, 1 + k)); % Computational pressure
% end
% 
% 
% ylim([0, max(1.5 * pmax(:, 1 + length(data_directories)))]);
% grid on;
% xlabel("$t$", "Interpreter", "latex", 'Fontsize',30);
% ylabel("Max pressure, $p$", "Interpreter", "latex", 'Fontsize', 30);
% ax = gca;
% ax.FontSize = 10;
% set(gca,'TickLabelInterpreter','latex');
% title(sprintf("Combined max pressure, Wagner impact $t$ = %.3f", ...
%     impact_time), "Interpreter", "latex", 'Fontsize', 15);
% 
% % xticks(linspace(min(pmax(:, 1)), max(pmax(:, 1)), 12))
% % xtickformat('%.2f')
% L = legend(["Outer", "Composite", "Turnover point", legend_entries]);
% set(L, 'Interpreter', 'latex');
% set(L, 'FontSize', 11);
% set(L, 'Location', 'northeast');
% 
% print(gcf, sprintf('%s/pmax.png', results_directory),'-dpng','-r300');