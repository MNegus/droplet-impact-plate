%% force_analysis.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyses the force outputted by the simulations and compares them to the
% prediction given by Wagner theory. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('pressures');
addpath('forces');

%% Data definitions
% Here we specify the locations where the plate output files are stored. We
% expect a different directory for each simulation result.

% Parent directory where all of the data is stored under (e.g. external
% hard drive location)
parent_directory = "/mnt/newarre/const_acc_f_output/";

% Directory where the resulting videos are to be stored
results_directory = sprintf("%s/Analysis", parent_directory);

% Individual directory names under the master directory. We assume that the
% plate output files are stored under
% master_directory/data_directory(k)/cleaned_data
data_directory = "level_12";

% Adds the parent directory to the start of the data directory
data_directory = strcat(parent_directory, data_directory); 

% Readable names to label the plots for each of the data directories
legend_entry = "Level 12";


%% Parameters
% Physical parameters
rho_w = 998;
R0 = 1e-3;
U0 = 5;
T0 = R0 / U0;
Patm = 10^5;

eps = 1;

p_dim = @(p) Patm + rho_w * U0^2 * p;
F_dim = @(F) rho_w * U0^2 * R0^2 * F;

% Parameters common to the simulations to aid in visualisation
initial_drop_height = 0.125; 

% Position to start plots at 
start_pos =  5;
end_pos = 600;
no_frames = end_pos - start_pos;
plot_range = start_pos : end_pos;

% (Constant) acceleration of the plate
a = 0.5;

% Displacement of s
s = @(t) 0.5 * a * t.^2;
sdot = @(t) a * t;
sddot = @(t) a;

% Finds computational turnover points
read_comp_ds = dlmread(sprintf('%s/cleaned_data/turnover_points.txt', ...
    data_directory));
comp_ds = read_comp_ds(:, 2:3);


%% Force comparison
% Loops over the outputs from the simulation and compares the force to the
% Wagner theory prediction

% Theoretical time of impact for a stationary plate
impact_time = initial_drop_height;

% Wagner theory force for moving plate
outer_force = @(t) 2 * sqrt(3 * (t - s(t))) ...
    .* (3 * (1 - sdot(t)).^2 - 2 * sddot(t) .* (t - s(t)));

% Arrays for calculated forces 
% Total force contribution from pressure
total_pressure_force = zeros(no_frames, 1); 

% Total force contribution from stress
total_stress_force = zeros(no_frames, 1); 

% Force contribution from pressure in the jet
jet_pressure_force = zeros(no_frames, 1); 

% Force contribution from entrapped bubbles
bubble_pressure_force = zeros(no_frames, 1); 

% Force contribution from bulk of the liquid
bulk_pressure_force = zeros(no_frames, 1);

% Times file
times = dlmread(sprintf('%s/cleaned_data/plate_outputs/times.txt',...
                        data_directory));

% % Wagner s-solution
% t_max = times(end_pos, 2)
% [t_vals, s_vals, sdot_vals, sddot_vals] = s_solution(t_max, alpha, beta, gamma, eps);
% t_vals = t_vals(1 : end - 1);
% wagner_force = composite_force(t_vals, s_vals, sdot_vals, sddot_vals, eps);

% Integrated Wagner pressure
total_composite_force = zeros(no_frames, 1);
jet_composite_force = zeros(no_frames, 1);
bulk_composite_force = zeros(no_frames, 1);
wagner_ds = zeros(no_frames, 1);
      
%% Iterates over time
for m = start_pos : end_pos
    
    % Saves time
    t = times(m, 2);
    
    % Integrated Wagner pressure
    if t >= impact_time
        s_val = s(t - impact_time);
        sdot_val = sdot(t - impact_time);
        sddot_val = sddot(t - impact_time);
        
        % Saves turnover point
        [d, ~, ~, ~] = s_dependents(t - impact_time, s_val, sdot_val, ...
            sddot_val);
        wagner_ds(m) = d;
        
        [outer_rs, outer_ps, comp_rs, comp_ps] ...
            = outer_and_comp_pressure(t - impact_time, s_val, sdot_val, ...
            sddot_val, 0, 1.25 * d, 1);
            
        % Integrates all of the composite solution
        total_composite_force(m - start_pos + 1) ...
            = trapz(comp_rs, 2 * pi * comp_rs .* comp_ps);
        
        % Configures cutoff r
        r_cutoff = d;
    
        % Values of r and p in the bulk solution
        bulk_ps = comp_ps(comp_rs <= r_cutoff);
        bulk_rs = comp_rs(comp_rs <= r_cutoff);

        bulk_composite_force(m - start_pos + 1) ...
            = trapz(bulk_rs, 2 * pi * bulk_rs .* bulk_ps);
        
        % Values of r and p in the jet solution
        jet_ps = comp_ps(comp_rs > r_cutoff);
        jet_rs = comp_rs(comp_rs > r_cutoff);
        
        jet_composite_force(m - start_pos + 1) ...
            = trapz(jet_rs, 2 * pi * jet_rs .* jet_ps);
        
    end
    
    
    % Loads in data from the text file
    m
    output_matrix = ...
        dlmread(...
            sprintf('%s/cleaned_data/plate_outputs/output_%d.txt', ...
                data_directory, m));

    % Sorts in increasing order of r
    [~, sorted_idxs] = sort(output_matrix(:, 1));
    sorted_mat = output_matrix(sorted_idxs, :);

    % Saves values of r and pressure
    rs = sorted_mat(:, 1);
    ps = sorted_mat(:, 3);
    stress = sorted_mat(:, 4);

    % Calculate the forces using trapezoidal rule
    total_pressure_force(m - start_pos + 1) ...
        = trapz(rs, 2 * pi * rs .* ps);
    total_stress_force(m - start_pos + 1) ...
        = trapz(rs, 2 * pi * rs .* stress);

    % Calculates component forces after the impact time
    if (t >= impact_time)
        if (comp_ds(m, 1) == 0) 
            % If turnover point not defined, set the forces components to
            % be zero
            jet_pressure_force(m - start_pos + 1) = 0;
            bubble_pressure_force(m - start_pos + 1) = 0;
            bulk_pressure_force(m - start_pos + 1) = 0;
        else
            % Set the cutoff point to be the turnover point
            r_cutoff = comp_ds(m, 1);
            
            % Calculate force due to pressure in jet region
            jet_ps = ps(rs > r_cutoff);
            jet_rs = rs(rs > r_cutoff);
            jet_pressure_force(m - start_pos + 1) ...
                = trapz(jet_rs, 2 * pi * jet_rs .* jet_ps);
            
            % Values of the volume fraction on surface (1 for liquid, 0 for
            % air)
            fs = sorted_mat(:, 5); 
            
            % Calculate force due to bubbles within the turnover curve
            bubble_rs = rs((rs <= r_cutoff) & (fs == 0));
            bubble_ps = ps((rs <= r_cutoff) & (fs == 0));
            if length(bubble_rs) >= 2
                bubble_pressure_force(m - start_pos + 1) ...
                    = trapz(bubble_rs, 2 * pi * bubble_rs .* bubble_ps);
            else
                bubble_pressure_force(m - start_pos + 1) = 0;
            end
            
            % Calculate force due to liquid within the turnover curve
            bulk_rs = rs((rs <= r_cutoff) & (fs == 1));
            bulk_ps = ps((rs <= r_cutoff) & (fs == 1));
            bulk_pressure_force(m - start_pos + 1) ...
                = trapz(bulk_rs, 2 * pi * bulk_rs .* bulk_ps);
        end
        
        
    end
        

    
    
end



%% Plotting 
output_ts = times(plot_range, 2);
dimen_ts = 1000 * T0 * output_ts;

% Plots the outer, composite and output forces
close all;
figure(1);
hold on;
% plot(dimen_ts, F_dim(outer_force(output_ts - impact_time)), 'color', 'black', ...
%     'Linestyle', '--', 'Linewidth', 1);
plot(dimen_ts, F_dim(total_composite_force),'color', 'black', ...
    'Linestyle', '-', 'Linewidth', 1); 
plot(dimen_ts, F_dim(bulk_composite_force), 'color', 'black', ...
    'Linestyle', '--', 'Linewidth', 1); 
plot(dimen_ts, F_dim(jet_composite_force), 'color', 'black', ...
    'Linestyle', '-.', 'Linewidth', 1); 

plot(dimen_ts, F_dim(total_pressure_force), 'Linewidth', 1);
% plot(dimen_ts, F_dim(total_stress_force));
plot(dimen_ts, F_dim(bulk_pressure_force), 'Linewidth', 1);
plot(dimen_ts, F_dim(jet_pressure_force), 'Linewidth', 1);
plot(dimen_ts, F_dim(bubble_pressure_force), 'Linewidth', 1);

grid on;
legend(["Composite", "Bulk composite", "Jet composite", ...
    "Total pressure force", "Bulk pressure force", ...
    "Jet pressure force", "Bubble pressure force"], ...
    "Interpreter", "latex", ...
    'Location', 'northwest', 'fontsize', 10);

xlabel("$t$ / ms", "Interpreter", "latex", "Fontsize", 15);
ylabel("Force / N", "Interpreter", "latex", "Fontsize", 15);
ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title("Force on plate: Accelerating frame. Level 12. $a \equiv 1250$ m/s$^2$", "Interpreter", "latex", ...
    'Fontsize', 14);
ylim([-0.01 0.1])
print(gcf, sprintf('%s/component_force_comparison.png', results_directory), ...
 '-dpng', '-r300');

%% PLOT COMPARISON OF TURNOVER POINTS
figure(2);
hold on
plot(dimen_ts, wagner_ds(plot_range, 1), 'color', 'black', 'Linewidth', 1);
plot(dimen_ts, comp_ds(plot_range, 1), 'linewidth', 1);

legend(["Wagner", "Computational"], ...
    "Interpreter", "latex", ...
    'Location', 'northwest', 'fontsize', 12);


xlabel("$t$ / ms", "Interpreter", "latex", "Fontsize", 15);
ylabel("$d(t)$ / mm", "Interpreter", "latex", "Fontsize", 15);
grid on;
ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title("Turnover points: Accelerating frame. Level 12. $a \equiv 1250$ m/s$^2$", "Interpreter", "latex", ...
    'Fontsize', 14);    
print(gcf, sprintf('%s/turnover_points.png', results_directory), ...
 '-dpng', '-r300');
    
    
    
    
    

