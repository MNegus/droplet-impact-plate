%% force_comparison.m
% Script to compare the forces on the plate for where the initial height
% has been varied. The times are shifted so the theoretical time of impact
% is always the same
%

% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions


% Parent directory where all the data is held
parent_directory = '/mnt/newarre/cantilever_paper_data/gamma_varying/';

% Directory to save the figure(s)
analysis_directory = "Analysis";
analysis_directory = strcat(parent_directory, analysis_directory);

% Range of alphas
gammas = [0, 10, 100, 500, 1000];

% Defines arrays for all the values of alpha
data_directories = string(length(gammas));
legend_entries = string(length(gammas));

for k = 1 : length(gammas)
    gamma = gammas(k);
    
    data_directories(k) = [parent_directory, '/gamma_', num2str(gamma)];
    
    legend_entries(k) = ['$\gamma =$ ', num2str(gamma)] ;
end

% Stationary plate data
stationary_plate_directory = '/mnt/newarre/cantilever_paper_data/stationary_plate';


%% Parameters

% Value of epsilon
eps = 1;

% Plate parameters
alpha = 2;
beta = 0;

% Initial drop height 
initial_drop_heights = 0.125;

% Impact time 
impact_time = initial_drop_heights;

% Maximum time
t_max = 0.8;

%% Wagner solution
% Stationary plate Wagner solution
wagner_t = linspace(1e-7, t_max - impact_time, 1e3);
s = zeros(size(wagner_t));
sdot = zeros(size(wagner_t));
sddot = zeros(size(wagner_t));

% Composite force solution
wagner_force = composite_force(wagner_t, s, sdot, sddot, eps);

%% Plotting force comparisons
close all;

fig = figure(1);
hold on;

N = length(data_directories);

start_color = 0.75;
end_color = 0.3;

m = (end_color - start_color) / (N - 1);
c = 0.5 * (start_color + end_color - (N + 1) * m);

for k = 1 : length(data_directories)
   output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(k)));
   
   line_color = (m * k + c) * [1 1 1];
   
   ts = output_mat(:, 1);
   Fs = output_mat(:, 3);
   
   % Rescale ts with the impact time
   ts = ts - impact_time;
   
   figure(1);
   plot(ts, Fs, 'Displayname', legend_entries(k), 'Linewidth', 2, ...
       'color', line_color);
end

% Stationary plate solution
output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", ...
    stationary_plate_directory));
ts = output_mat(:, 1);
ts = ts - impact_time;
Fs = output_mat(:, 3);
plot(ts, Fs, 'Linewidth', 3, 'color', 'black', ...
    'Displayname', ['Stationary', newline, 'plate']);

% Analytical solution
plot(wagner_t, wagner_force, ...
    'color', 'black', 'Linewidth', 3, 'Linestyle', '--', ...
    'Displayname', ['Analytical', newline, 'solution']);


% x limits
xlim([-impact_time t_max - impact_time]);


% Arrow for increasing gamma
X = [0.45 0.3];
Y = [0.3 0.65];
annotation('arrow', X, Y);

% Arrow label
txt = '$\gamma$';
text(0.05, 3.15, txt, "Interpreter", "Latex", "Fontsize", 30);

% legend("Interpreter", "latex", "location", "northwest", "Fontsize", 14);

grid on;
xlabel("$t$", "Interpreter", "latex");
ylabel("$F(t)$", 'Interpreter', 'latex');
set(gca, 'XTick', 0 : 0.2 : t_max);
set(gca, 'YTick', 0 : 1 : 5 );

ax = gca;
ax.FontSize = 30;
set(gca,'TickLabelInterpreter','latex');
% title(['Force on plate: $\alpha =$ ', ...
%     num2str(alpha), ', $\beta =$ ' num2str(beta)], "Interpreter", "latex", ...
%     'Fontsize', 14);

set(gcf, 'Position',  [0, 0, 500, 700]);
ax = gca;

plot_name = sprintf("%s/gamma_force_comparison.png", analysis_directory);
exportgraphics(ax, plot_name, 'resolution', 300);

% print(gcf, sprintf("%s/force_comparison.png", analysis_directory), ...
%     '-dpng', '-r300');