%% plate_analysis.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyses the plate position outputted by the simulations and compares
% them to the Wagner theory prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters

% Physical parameters
rho_w = 998;
R0 = 1e-3;
U0 = 5;
T0 = R0 / U0;
Patm = 10^5;

% Plate ODE coefficients
alpha = 1;
beta = 0;
gamma = 40;
eps = 1;

% Initial drop height
initial_drop_height = 0.125;

impact_time = initial_drop_height;

t_max = 0.6;

% Solid wall Wagner solution
solid_t = linspace(1e-6, t_max - impact_time, 1e3)';
wagner_s = zeros(length(solid_t), 1);
wagner_sdot = zeros(length(solid_t), 1);
wagner_sddot = zeros(length(solid_t), 1);
solid_F = composite_force(solid_t, wagner_s, wagner_sdot, wagner_sddot, eps);

% Moving plate Wagner solution
[wagner_t, wagner_s, wagner_sdot, wagner_sddot] = s_solution(t_max + impact_time, alpha, ...
    beta, gamma, eps);
wagner_F = composite_force(wagner_t, wagner_s, wagner_sdot, wagner_sddot, eps);


%% Data definitions
% Here we specify the locations where the plate output files are stored. We
% expect a different directory for each simulation result.

% Parent directory where all of the data is stored under (e.g. external
% hard drive location)
parent_directory = "/mnt/newarre/supervision_runs_15_july/";
% parent_directory = "/mnt/newarre/gamma_0_runs/";

% Directory where the resulting videos are to be stored
results_directory = sprintf("%s/Analysis", parent_directory);

% Individual directory names under the master directory. We assume that the
% plate output files are stored under
% master_directory/data_directory(k)/cleaned_data
data_directories = ["alpha_1_gamma_40", "solid_wall_peak_removal"];
no_dirs = length(data_directories); % Number of entries

% Adds the parent directory to the start of the data directories
for k = 1 : length(data_directories)
    data_directories(k) = strcat(parent_directory, data_directories(k)); 
end
% Readable names to label the plots for each of the data directories
legend_entries = ["$\alpha = 1, \gamma = 40$", "Solid wall"];


%% Force comparison
close all;
figure(1);
hold on;

for k = 1 : length(data_directories)
    output_mat = dlmread(sprintf('%s/cleaned_data/output.txt',...
                            data_directories(k)));
    comp_ts = output_mat(:, 1);
    comp_Fs = output_mat(:, 3);
    plot(comp_ts, comp_Fs);
end

plot(wagner_t + impact_time, wagner_F, 'color', [0 0 0]);
plot(solid_t + impact_time, solid_F, 'color', [0 0 0], 'Linestyle', '--');

grid on;
legend([legend_entries, "Composite", "Solid wall composite"], ...
    "Interpreter", "latex", "location", "northwest", "Fontsize", 10);
xlabel("$t$", "Interpreter", "latex", "Fontsize", 15);
ylabel("Force", "Interpreter", "latex", "Fontsize", 15);
ax = gca;
ax.FontSize = 12;
xlim([0 t_max]);
% ylim([0 3.5]);
set(gca,'TickLabelInterpreter','latex');
titlestr = ['Force on plate. $\alpha$ = ', num2str(alpha), ...
    ', $\beta = $', num2str(beta), ', $\gamma = $', num2str(gamma)];
% titlestr = ['Force on plate. $\beta = $', num2str(beta), ', $\gamma = $', num2str(gamma)];


title(titlestr, "Interpreter", "latex", ...
    'Fontsize', 14);
print(gcf, sprintf('%s/force_comparison.png', results_directory), ...
 '-dpng', '-r300');


% %% Peak explanation
% close all;
% start_pos = 1;
% peak_thresh = 4;
% influence = 0.1;
% output_mat = dlmread(sprintf('%s/cleaned_data/output.txt',...
%                             data_directories(1)));
% ts = output_mat(start_pos : end, 1);
% raw_force = output_mat(start_pos : end, 2);
% used_force = output_mat(start_pos : end, 3);
% avg = output_mat(start_pos : end, 4);
% std = output_mat(start_pos : end, 5);
% 
% figure(2);
% hold on;
% plot(ts, raw_force);
% plot(ts, used_force, 'color', [0 0 0]);
% 
% upper_bound = avg + peak_thresh * std;
% lower_bound = avg - peak_thresh * std;
% % plot(ts, lower_bound, 'color', [0 0 0]);
% % plot(ts, upper_bound, 'color', [0 0 0]);
% 
% 
% shade(ts, upper_bound, ts, lower_bound,'FillType',[1 2]);
% xlabel("t");
% ylabel("F");
% legend(["Raw", "Used", "Upper bound", "Lower bound", "Valid region"]);
% 
