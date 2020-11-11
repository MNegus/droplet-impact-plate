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

% Parent directories where all of the data is stored under (e.g. external
% hard drive location)
stationary_directory = "/media/michael/newarre/cantilever_paper_data/stationary_plate";
moving_directory = "/media/michael/newarre/cantilever_paper_data/gamma_varying/gamma_500";
% stationary_directory = "/scratch/negus/cant_paper_vid_comparisons/stationary";
% moving_directory = "/scratch/negus/cant_paper_vid_comparisons/moving";

% Directory where the resulting videos are to be stored
results_directory = "/media/michael/newarre/presentation_data";
% results_directory = "/scratch/negus/cant_paper_vid_comparisons/results";

% Readable names to label the plots for each of the data directories
legend_entries = ["Numerical stationary", "Numerical moving"];

%% Parameters

% Plate parameters
alpha = 2;
beta = 0;
gamma = 500;
eps = 1;

% Times
tmax = 0.8;

% Parameters common to the simulations to aid in visualisation
initial_drop_height = 0.125; 

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
times = dlmread(sprintf('%s/cleaned_data/plate_outputs/times.txt', stationary_directory));

% Analytical time variables
analytical_tvals = times(:, 2) - impact_time;
m0 = find(~analytical_tvals);
analytical_tvals = analytical_tvals(analytical_tvals > 0);

% Solves for the plate displacement
[~, s, sdot, sddot] = s_solution_alt(analytical_tvals, alpha, beta, gamma, eps);

% Position to start video at
start_pos =  1;
end_pos = 800;

output_range = start_pos : end_pos;
no_frames = end_pos - start_pos;

% Sets up figure
figure(1);
hold on;
grid on;
xlabel("$r$", "Interpreter", "latex", 'Fontsize',30);
ylabel("$p(r, -s(t), t)$", "Interpreter", "latex", 'Fontsize', 30);
ax = gca;
ax.FontSize = 16;
set(gca,'TickLabelInterpreter','latex');
x_limits = [-3, 3];
abs_pmax = 50;
y_limits = [-1, abs_pmax];
xlim(x_limits);
ylim(y_limits);

% Creates animated line for the Wagner pressure
stationary_wagner_line = animatedline('color', [0 0 0], 'Linewidth', 1.5, 'Linestyle', '--');
moving_wagner_line = animatedline('color', 0.5 * [1 1 1], 'Linewidth', 1.5, 'Linestyle', '--');

% % Wagner turnover line
% wagner_turnover_line = animatedline('color', 'red', 'Linestyle', '--', ...
%     'LineWidth', 1.5);

% MAKE THIS SMARTER
% Creates animated lines for the numerical results
% stationary_line = animatedline('Color', [0, 0.4470, 0.7410], 'Linewidth', 1.5);
% moving_line = animatedline('Color', [0.8500, 0.3250, 0.0980], 'Linewidth', 1.5);
stationary_line = animatedline('Color', [0 0 0], 'Linewidth', 1.5);
moving_line = animatedline('Color', 0.5 * [1 1 1], 'Linewidth', 1.5);


% Sets up the legend
L = legend(["Wagner stationary", "Wagner moving", legend_entries]);
set(L, 'Interpreter', 'latex');
set(L, 'FontSize', 15);

% Sets pixel size of the figure
width=2048;
height=257; % Either 257 or 512
set(gcf,'position',[10,10,width,height])
% set(gca,'LooseInset',get(gca,'TightInset'))

% Create the video writer with 5 fps
writerObj ...
    = VideoWriter(sprintf('%s/pressure_overlay.avi', results_directory));
writerObj.FrameRate = 25;
open(writerObj);

%%
% Iterates over timec

dim = [.15 .6 .3 .3];
str = "";
h = annotation('textbox',dim,'String',str,'FitBoxToText','on', "Interpreter", "latex", 'Fontsize', 15);

for m = start_pos : start_pos + no_frames -1

    t = times(m, 2); % Time

    d = sqrt(3 * t);
    
    %% Plots stationary computational
    output_matrix = dlmread(...
                sprintf('%s/cleaned_data/plate_outputs/output_%d.txt', ...
                    stationary_directory, m));
    
    % Sorts in increasing order of r
    [~, sorted_idxs] = sort(output_matrix(:, 1));
    sorted_mat = output_matrix(sorted_idxs, :);

    % Saves values of r and pressure
    rs = sorted_mat(:, 1);
    ps = sorted_mat(:, 3);
    
    % Adds the pressure line
    clearpoints(stationary_line)
    addpoints(stationary_line, -rs, ps);
    
    %% Plots moving computational
    output_matrix = dlmread(...
                sprintf('%s/cleaned_data/plate_outputs/output_%d.txt', ...
                    moving_directory, m));
    
    % Sorts in increasing order of r
    [~, sorted_idxs] = sort(output_matrix(:, 1));
    sorted_mat = output_matrix(sorted_idxs, :);

    % Saves values of r and pressure
    rs = sorted_mat(:, 1);
    ps = sorted_mat(:, 3);
    
    % Adds the pressure line
    clearpoints(moving_line)
    addpoints(moving_line, rs, ps);
    
  
    %% Wagner lines
    if t > impact_time
        
        %% Plots stationary line
        [~, ~, comp_rs, comp_ps] ...
            = outer_and_comp_pressure(t - impact_time, 0, 0, ...
                0, 0, 2 * d, 1);
            
        % Adds composite pressure to graph
        clearpoints(stationary_wagner_line);
        addpoints(stationary_wagner_line, -comp_rs, comp_ps);
        
        %% Plots moving line
        % Values of s
        s_val = s(m - m0);
        sdot_val = sdot(m - m0);
        sddot_val = sddot(m - m0);
        
        % Saves turnover point
        [d, ~, ~, ~] = s_dependents(t - impact_time, s_val, sdot_val, ...
            sddot_val);
        
        % Determines outer and composite solution
        [~, ~, comp_rs, comp_ps] ...
            = outer_and_comp_pressure(t - impact_time, s_val, sdot_val, ...
                sddot_val, 0, 2 * d, 1);
        
        % Adds composite pressure to graph
        clearpoints(moving_wagner_line);
        addpoints(moving_wagner_line, comp_rs, comp_ps);

    end

%     title(sprintf("$t$ = %.3f",... 
%         times(m, 2) - impact_time),  ...
%         "Interpreter", "latex", 'Fontsize', 15);
    
    str = sprintf("$t$ = %.3f", times(m, 2) - impact_time);
    set(h, 'String',str);
    xlim(x_limits);
    
    wagner_pmax = 3 / (8 * (t - impact_time));
    if (t > impact_time)
        if (wagner_pmax > 1.3 * abs_pmax)
            y_limits = [-0.1, abs_pmax];
        elseif (1.3 * wagner_pmax > 3)
            y_limits = [-0.1 * wagner_pmax, 1.3 *  wagner_pmax];
        else
            y_limits = [-0.1 * 3, 3];
        end
    end
    ylim(y_limits);
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