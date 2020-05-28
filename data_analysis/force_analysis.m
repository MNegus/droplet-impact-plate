%% force_analysis.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyses the force outputted by the simulations and compares them to the
% prediction given by Wagner theory. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('pressures');

%% Data definitions
% Here we specify the locations where the plate output files are stored. We
% expect a different directory for each simulation result.

% Parent directory where all of the data is stored under (e.g. external
% hard drive location)
% parent_directory = "/mnt/newarre/wall_force_test/";
% parent_directory = "/mnt/newarre/embed_force_test/";
parent_directory = "/mnt/newarre/acc_test/";

% Directory where the resulting videos are to be stored
results_directory = sprintf("%s/Analysis", parent_directory);

% Individual directory names under the master directory. We assume that the
% plate output files are stored under
% master_directory/data_directory(k)/cleaned_data
% data_directories = ["wall_test_10", "wall_test_11", "wall_test_12", "wall_test_13"];
% data_directories = ["embed_test_10"];
% data_directories = ["embed_test_10", "embed_test_11", "embed_test_12"]; 
% data_directories = ["wall_test_10"];
data_directories = ["constant_acc"];
no_dirs = length(data_directories); % Number of entries

% Adds the parent directory to the start of the data directories
for k = 1 : length(data_directories)
    data_directories(k) = strcat(parent_directory, data_directories(k)); 
end
% Readable names to label the plots for each of the data directories
% legend_entries = ["Level 10", "Level 11", "Level 12", "Level 13"];
% legend_entries = ["Level 10", "Level 11", "Level 12"];
% legend_entries = ["Level 10"];
legend_entries = ["Level 12"];

%% Parameters
% Physical parameters
rho_w = 998;
R0 = 1e-3;
U0 = 10;
T0 = R0 / U0;
Patm = 10^5;

p_dim = @(p) Patm + rho_w * U0^2 * p;
F_dim = @(F) rho_w * U0^2 * R0^2 * F;

% Parameters common to the simulations to aid in visualisation
initial_drop_height = 0.125; 

% (Constant) acceleration of the plate
a = 0.5;

% Position to start plots at 
start_pos =  5;
end_pos = 600;
no_frames = end_pos - start_pos;
plot_range = start_pos : end_pos;

% If true, then the jet is neglected in the force calculation
cutoff = false;

% Finds computational turnover points
if cutoff
    read_comp_ds = dlmread(sprintf('%s/cleaned_data/turnover_points.txt', ...
        data_directories(end)));
    comp_ds = read_comp_ds(:, 2:3);
end


%% Force comparison
% Loops over the outputs from the simulation and compares the force to the
% Wagner theory prediction

% Theoretical time of impact for a stationary plate
impact_time = initial_drop_height;

% Wagner theory force for stationary plate
outer_force = @(t) 6 * sqrt(3) * sqrt(t - impact_time);

% Array for calculated forces and pressures
pressure_force = zeros(no_frames, no_dirs); % Force contribution from pressure
stress_force = zeros(no_frames, no_dirs); % Force contribution from stress

% Cutoff pressure force
cutoff_pressure_force = zeros(no_frames, no_dirs); % Force contribution from pressure


% Times file
times = dlmread(sprintf('%s/cleaned_data/plate_outputs/times.txt',...
                        data_directories(1)));


% Integrated Wagner pressure
integrated_outer_pressure = zeros(no_frames, 1);
integrated_composite_pressure = zeros(no_frames, 1);
integrated_cutoff_pressure = zeros(no_frames, 1);
      
%% Iterates over time
for m = start_pos : end_pos
    
    % Saves time
    t = times(m, 2);
    
    % Integrated Wagner pressure
    if t >= impact_time
        % Saves turnover point
        [d, ~, ~, ~] = s_dependents(t - impact_time, 0, 0, 0);
        
        % Values of s for constant acceleration
        s = 0.5 * a * (t - impact_time)^2;
        sdot = a * (t - impact_time);
        sddot = a;
        
        [outer_rs, outer_ps, comp_rs, comp_ps] ...
            = outer_and_comp_pressure(t - impact_time, s, sdot, sddot, ...
                0, 1.25 * d, 1);
            
        integrated_outer_pressure(m - start_pos + 1) ...
            = trapz(outer_rs, 2 * pi * outer_rs .* outer_ps);
            
        integrated_composite_pressure(m - start_pos + 1) ...
            = trapz(comp_rs, 2 * pi * comp_rs .* comp_ps);
        
        if cutoff
            % Configures cutoff r
            if comp_ds(m, 1) < 1e-3
                r_cutoff = d;
            else
                r_cutoff = comp_ds(m, 1);
            end
            
            comp_ps = comp_ps(comp_rs <= r_cutoff);
            comp_rs = comp_rs(comp_rs <= r_cutoff);
            
            integrated_cutoff_pressure(m - start_pos + 1) ...
                = trapz(comp_rs, 2 * pi * comp_rs .* comp_ps);
        end
        
    end
    
    for k = 1 : length(data_directories)
        % Loads in data from the text file
        [k, m]
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
%         stress = sorted_mat(:, 4);
       
        % Calculate the forces using trapezoidal rule
        pressure_force(m - start_pos + 1, k) ...
            = trapz(rs, 2 * pi * rs .* ps);
%         stress_force(m - start_pos + 1,  k) ...
%             = trapz(rs, 2 * pi * rs .* stress);
        
        if (cutoff) && (t >= impact_time)
            % If jet is neglected, then remove any r values greater than a
            % specified cutoff
%             cutoff = sqrt(3 * (t - impact_time))
            ps = ps(rs <= r_cutoff);
            rs = rs(rs <= r_cutoff);
            
            if length(rs) == 1
                cutoff_pressure_force(m - start_pos + 1, k) = 0;
            else
                cutoff_pressure_force(m - start_pos + 1, k) ...
                    = trapz(rs, 2 * pi * rs .* ps);
            end
        end
        
    end
    
    
end



%% Plotting 
output_ts = times(plot_range, 2)
dimen_ts = 1000 * T0 * output_ts;

% Plots the outer, composite and output forces
close all;
figure(1);
hold on;
% plot(dimen_ts, F_dim(outer_force(output_ts)), 'color', 'black', ...
%     'Linestyle', '--', 'Linewidth', 1);
plot(dimen_ts, F_dim(integrated_outer_pressure), 'color', 'black', ...
    'Linestyle', '--', 'Linewidth', 1);
s = @(t) 0.5 * a * t.^2;
sdot = @(t) a * t;
sddot = @(t) a;
outer_F = @(t) 4 * sqrt(3) * (-sddot(t) .* (t - s(t)).^1.5 ...
    + (1 - sdot(t)).^2 .* sqrt(t - s(t)));
plot(dimen_ts, F_dim(outer_F(output_ts - impact_time)));
plot(dimen_ts, F_dim(integrated_composite_pressure), 'color', 'black', ...
    'Linestyle', '-.', 'Linewidth', 1);
for k = 1 : no_dirs
   output_mat = dlmread(sprintf('%s/cleaned_data/output.txt', ...
        data_directories(k)));
    plot(dimen_ts, F_dim(output_mat(plot_range, 3)));
end
grid on;
legend(["Integrated outer", "Analytical outer", "Composite", legend_entries], "Interpreter", "latex", ...
    'Location', 'northwest', 'fontsize', 12);
xlabel("$t$ / ms", "Interpreter", "latex", "Fontsize", 15);
ylabel("Force / N", "Interpreter", "latex", "Fontsize", 15);
ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title("Force on plate: Moving frame. $a$ = 0.5", "Interpreter", "latex", ...
    'Fontsize', 14);
print(gcf, sprintf('%s/force_comparison.png', results_directory), ...
 '-dpng', '-r300');

% % Plots the cutoff forces
% figure(2);
% hold on;
% plot(dimen_ts, F_dim(outer_force(output_ts)), 'color', 'black', ...
%     'Linestyle', '--', 'Linewidth', 1);
% plot(dimen_ts, F_dim(integrated_cutoff_pressure), 'color', 'black', ...
%     'Linestyle', '-.', 'Linewidth', 1);
% for k = 1 : no_dirs
%     plot(dimen_ts, F_dim(cutoff_pressure_force(:, k)));
% end
% grid on;
% legend(["Outer", "Composite", legend_entries], "Interpreter", "latex", ...
%     'Location', 'northwest', 'fontsize', 12);
% xlabel("$t$ / ms", "Interpreter", "latex", "Fontsize", 15);
% ylabel("Force / N", "Interpreter", "latex", "Fontsize", 15);
% ax = gca;
% ax.FontSize = 12;
% set(gca,'TickLabelInterpreter','latex');
% title("Force on plate: Solid wall. Cutoff at turnover.", "Interpreter", "latex", ...
%     'Fontsize', 14);
% print(gcf, sprintf('%s/cutoff_force_comparison.png', results_directory), ...
%  '-dpng', '-r300');

    
    
    
    
    
    
    

