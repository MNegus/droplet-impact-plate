%% force_comparison.m
% Script to compare the forces on the plate for where the initial height
% has been varied. The times are shifted so the theoretical time of impact
% is always the same
%

% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions


% Parent directory where all the data is held
parent_directory = '/mnt/newarre/cantilever_paper_data/alpha_varying/';

% Directory to save the figure(s)
analysis_directory = "Analysis";
analysis_directory = strcat(parent_directory, analysis_directory);

% Range of alphas
alphas = [1, 2, 5, 10, 20, 100];

% Defines arrays for all the values of alpha
data_directories = string(length(alphas));
legend_entries = string(length(alphas));

for k = 1 : length(alphas)
    alpha = alphas(k);
    
    data_directories(k) = [parent_directory, '/alpha_', num2str(alpha)];
    
    legend_entries(k) = ['$\alpha =$ ', num2str(alpha)] ;
end

% Stationary plate data
stationary_plate_directory = '/mnt/newarre/cantilever_paper_data/stationary_plate';


%% Parameters

% Value of epsilon
eps = 1;

% Plate parameters
beta = 0;
gamma = 0;

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

for k = 1 : length(data_directories)
   output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(k)));
   
   ts = output_mat(:, 1);
   Fs = output_mat(:, 3);
   
   % Rescale ts with the impact time
   ts = ts - impact_time;
   
   figure(1);
   plot(ts, Fs, 'Displayname', legend_entries(k), 'Linewidth', 2);
end



% Stationary plate solution
output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", ...
    stationary_plate_directory));
ts = output_mat(:, 1);
ts = ts - impact_time;
Fs = output_mat(:, 3);
plot(ts, Fs, 'Linewidth', 3, 'color', 0.3 * [1 1 1], ...
    'Displayname', ['Stationary', newline, 'plate']);

% Analytical solution
plot(wagner_t, wagner_force, ...
    'color', 'black', 'Linewidth', 3, 'Linestyle', '--', ...
    'Displayname', ['Analytical', newline, 'solution']);


% x limits
xlim([-impact_time  t_max - impact_time]);


% Arrow for increasing alpha
X = [0.5 0.6];
Y = [0.25 0.70];
annotation('arrow', X, Y);

% Arrow label
txt = '$\alpha$';
text(0.375, 3.4, txt, "Interpreter", "Latex", "Fontsize", 30);

% legend("Interpreter", "latex", "location", "northwest", "Fontsize", 12);

grid on;
xlabel("$t$", "Interpreter", "latex");
ylabel("$F(t)$", 'Interpreter', 'latex');
set(gca, 'XTick', 0 : 0.2 : t_max);
set(gca, 'YTick', 0 : 1 : 5 );
ax = gca;
ax.FontSize = 30;
set(gca,'TickLabelInterpreter','latex');
% title(['Force on plate: $\beta =$ ', ...
%     num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
%     'Fontsize', 14);
set(gcf, 'Position',  [0, 0, 500, 700]);
ax = gca;

plot_name = sprintf("%s/alpha_force_comparison.png", analysis_directory);
pause(0.5);
exportgraphics(ax, plot_name, 'resolution', 300);

