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
parent_directory = "/mnt/newarre/cantilever_paper_data/stationary_plate";

% Directory where the resulting videos are to be stored
results_directory = sprintf("%s/Analysis", parent_directory);

% Individual directory names under the master directory. We assume that the
% plate output files are stored under
% master_directory/data_directory(k)/cleaned_data
data_directories = [""];
no_dirs = length(data_directories); % Number of entries

% Adds the parent directory to the start of the data directories
for k = 1 : length(data_directories)
    data_directories(k) = strcat(parent_directory, data_directories(k)); 
end
% Readable names to label the plots for each of the data directories
legend_entries = ["Stationary plate"];


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
start_pos =  100;
end_pos = 200;
no_frames = end_pos - start_pos;
plot_range = start_pos : end_pos;

% (Constant) acceleration of the plate
a = 0;

% Displacement of s
s = @(t) 0.5 * a * t.^2;
sdot = @(t) a * t;
sddot = @(t) a;

% % Plate parameters
% alpha = 1;
% beta = 0;
% gamma = 0;

% If true, then we separately calculate the force only up until the
% turnover point
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


% Array for calculated forces and pressures
full_pressure_force = zeros(no_frames, no_dirs); % Force contribution from pressure
full_stress_force = zeros(no_frames, no_dirs); % Force contribution from stress

% Cutoff pressure force
cutoff_pressure_force = zeros(no_frames, no_dirs); % Force contribution from pressure


% Times file
times = dlmread(sprintf('%s/cleaned_data/plate_outputs/times.txt',...
                        data_directories(1)));

% Wagner s-solution
wagner_t = times(times > impact_time) - impact_time;
s_vals = s(wagner_t);
sdot_vals = sdot(wagner_t);
sddot_vals = sddot(wagner_t);
wagner_outer_F = outer_force(wagner_t, s_vals, sdot_vals, sddot_vals, eps);
wagner_composite_F = composite_force(wagner_t, s_vals, sdot_vals, sddot_vals, eps);

% % Integrated Wagner pressure
integrated_composite_pressure = zeros(no_frames, 1);
integrated_cutoff_pressure = zeros(no_frames, 1);
      
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
        
        [outer_rs, outer_ps, comp_rs, comp_ps] ...
            = outer_and_comp_pressure(t - impact_time, s_val, sdot_val, ...
            sddot_val, 0, 1.25 * d, 1);
            
        integrated_composite_pressure(m - start_pos + 1) ...
            = trapz(comp_rs, 2 * pi * comp_rs .* comp_ps);
        
        if cutoff
            % Configures cutoff r to be theoretical turnover point
            r_cutoff = d;
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
        stress = sorted_mat(:, 4);
       
        % Calculate the forces using trapezoidal rule
        full_pressure_force(m - start_pos + 1, k) ...
            = trapz(rs, 2 * pi * rs .* ps);
        full_stress_force(m - start_pos + 1,  k) ...
            = trapz(rs, 2 * pi * rs .* stress);
        
        if (cutoff) && (t >= impact_time)
            % If jet is neglected, then remove any r values greater than a
            % specified cutoff
            
            % Determine cutoff
            if (comp_ds(m, 1) == 0) 
                cutoff_pressure_force(m - start_pos + 1, k) = 0;
            else
            
                ps = ps(rs <= comp_ds(m, 1));
                rs = rs(rs <= comp_ds(m, 1));

                if length(rs) == 1
                    cutoff_pressure_force(m - start_pos + 1, k) = 0;
                else
                    cutoff_pressure_force(m - start_pos + 1, k) ...
                        = trapz(rs, 2 * pi * rs .* ps);
                end
            end
        end
        
    end
    
    
end



%% Plotting 
output_ts = times(plot_range, 2)
dimen_ts = 1000 * T0 * output_ts;

% % Plots the outer, composite and output forces
% close all;
% figure(1);
% hold on;
% plot(1000 * T0 * t_vals, F_dim(outer_force), 'color', 'black', ...
%     'Linestyle', '--', 'Linewidth', 1);
% plot(dimen_ts(1: end - 1), F_dim(integrated_composite_pressure), 'color', 'black', ...
%     'Linestyle', '-.', 'Linewidth', 1);
% for k = 1:no_dirs
%     plot(dimen_ts, F_dim(full_pressure_force(:, k)));
% end
% 
% grid on;
% legend(["Outer force", "Composite force", legend_entries], ...
%     "Interpreter", "latex", "location", "southeast", "Fontsize", 12);
% xlabel("$t$ / ms", "Interpreter", "latex", "Fontsize", 15);
% ylabel("Force / N", "Interpreter", "latex", "Fontsize", 15);
% ax = gca;
% ax.FontSize = 12;
% set(gca,'TickLabelInterpreter','latex');
% title("Force on plate: Accelerating frame. $a \equiv 1250$ m/s$^2$", "Interpreter", "latex", ...
%     'Fontsize', 14);
% ylim([-0.01 0.1])
% print(gcf, sprintf('%s/full_force_comparison.png', results_directory), ...
%  '-dpng', '-r300');
% 

output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(1)));
simulation_ts = output_mat(:, 1);
simulation_Fs = output_mat(:, 2);

figure(2);
hold on;
plot(simulation_ts, simulation_Fs);
plot(output_ts, full_pressure_force);

    
    
    
    
    
    

