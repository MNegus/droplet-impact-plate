%% force_comparison.m
% Script to compare the displacements s. 
% The times are shifted so the theoretical time of impact
% is always the same
%

% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions


% Parent directory where all the data is held
parent_directory = '/mnt/newarre/cantilever_paper_data/';

% Directory to save the figure(s)
analysis_directory = "special_case/Analysis";
analysis_directory = strcat(parent_directory, analysis_directory);

% Defines array with both directories in
data_directories = ["gamma_varying/gamma_500"];

% Concatenates arrays to include parent directory
for k = 1 : length(data_directories)
    data_directories(k) = strcat(parent_directory, data_directories(k)); 
end

legend_entries = ["Computational"];


%% Parameters

% Value of epsilon
eps = 1;

% Plate parameters
alpha = 2;
beta = 0;
gamma = 500;

% Initial drop height 
initial_drop_heights = 0.125;

% Impact time 
impact_time = initial_drop_heights;

% Maximum computational time
t_max = 0.8;

%% Wagner solution

% Plate displacement solution
[wagner_t, s, sdot, sddot] = s_solution(t_max - impact_time, alpha, beta, gamma, eps);

% Composite force solution
wagner_force = composite_force(wagner_t, s, sdot, sddot, eps);

%% Plotting
close all;

% % Force comparison
% figure(1);
% hold on;

% s displacement comparison
figure(2);
hold on;
% 
% % sdot comparison
% figure(3);
% hold on;
% 
% % sddot comparison
% figure(4);
% hold on;


output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(1)));

ts = output_mat(:, 1);
Fs = output_mat(:, 3);
ss = output_mat(:, 6);
sdots = output_mat(:, 7);
sddots = output_mat(:, 8);

% Shift time so t = 0 happens at impact time
ts = ts - impact_time;

% figure(1);
% h(1) = plot(ts, Fs, 'Linewidth', 2, 'color', 0.25 * [1 1 1]);
% h(2) = plot(ts, gamma * ss, 'Linewidth', 2, 'color', 0.75 * [1 1 1]);
   
figure(2);
plot(ts, ss, 'color', 0.25 * [1 1 1], 'Linewidth', 2);
%    
% figure(3);
% plot(ts, sdots, 'color', 0.5 * [1 1 1], 'Linewidth', 2);
% 
% figure(4);
% plot(ts, sddots, 'color', 0.5 * [1 1 1], 'Linewidth', 2);



% % Forces comparisons
% figure(1);
% h(3) = plot(wagner_t, composite_force(wagner_t, s, sdot, sddot, eps), ...
%     'Linewidth', 2, 'Linestyle', '--', 'color', 'black');
% h(4) = plot(wagner_t, gamma * s, ...
%     'Linewidth', 2, 'Linestyle', '-.', 'color', 'black');
% legend(h([1, 3, 2, 4]), ["$F(t)$: Analytical", "$F(t)$: Numerical", ...
%     "$\gamma s(t)$: Analytical", "$\gamma s(t)$: Numerical"], ...
%     "Interpreter", "latex", "location", "northwest", "Fontsize", 12);
% xlim([-impact_time t_max - impact_time]);
% grid on;
% xlabel("$t$", "Interpreter", "latex");
% ylabel("Force", 'Interpreter', 'latex');
% ax = gca;
% ax.FontSize = 12;
% % ylim([-0.8 0.8]);
% % set(gca, 'XTick', 0 : 0.2 : t_max);
% % set(gca, 'YTick', -0.8 : 0.4 : 0.8);
% set(gca,'TickLabelInterpreter','latex');
% % title(['Plate acceleration: $\alpha =$ ', num2str(alpha), ', $\beta =$ ', ...
% %     num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
% %     'Fontsize', 14);
% % print(gcf, sprintf("%s/plate_acceleration_comparison.png", analysis_directory), ...
% %     '-dpng', '-r300');
% set(gcf, 'Position',  [0, 0, 700, 250]);
% ax = gca;
% plot_name = sprintf("%s/force_comparison.png", analysis_directory);
% pause(0.1);
% exportgraphics(ax, plot_name, 'resolution', 300);

% s comparisons
figure(2);
plot(wagner_t, s, 'color', 'black', 'Linewidth', 2, 'Linestyle', ':');
% Add scatter for desired labelled points
t_labeled = [0.13, 0.37, 0.62];
s_labeled = interp1(ts, ss, t_labeled);
% scatter(t_labeled, s_labeled, '*');
% legend([legend_entries, "Analytical"], ...
%     "Interpreter", "latex", "location", "northwest", "Fontsize", 12);
xlim([-impact_time t_max - impact_time]);
ylim([0 1.1 * max(s)]);
set(gca, 'XTick', -impact_time : impact_time : t_max - impact_time);
grid on;
xlabel("$t$", "Interpreter", "latex");
ylabel("$s(t)$", 'Interpreter', 'latex');
ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');

% title(['Plate displacement: $\alpha =$ ', num2str(alpha), ', $\beta =$ ', ...
%     num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
%     'Fontsize', 14);
% print(gcf, sprintf("%s/plate_displacement_comparison.png", analysis_directory), ...
%     '-dpng', '-r300');
set(gcf, 'Position',  [0, 0, 700, 200]);
plot_name = sprintf("%s/plate_displacement_comparison.png", analysis_directory);
pause(0.1);
exportgraphics(ax, plot_name, 'resolution', 300);
% 
% % sdot comparisons
% figure(3);
% plot(wagner_t, sdot, ...
%     'color', 'black', 'Linewidth', 2, 'Linestyle', '--');
% % legend([legend_entries, "Analytical"], ...
% %     "Interpreter", "latex", "location", "northeast", "Fontsize", 12);
% 
% 
% xlim([-impact_time t_max - impact_time]);
% ylim([-0.08 0.08]);
% set(gca, 'XTick', 0 : 0.2 : t_max - impact_time);
% set(gca, 'YTick', -0.08 : 0.04 : 0.08);
% grid on;
% xlabel("$t$", "Interpreter", "latex");
% ylabel("$\dot{s}(t)$", 'Interpreter', 'latex');
% 
% ax = gca;
% ax.FontSize = 20;
% set(gca,'TickLabelInterpreter','latex');
% % title(['Plate velocity: $\alpha =$ ', num2str(alpha), ', $\beta =$ ', ...
% %     num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
% %     'Fontsize', 14);
% % print(gcf, sprintf("%s/plate_velocity_comparison.png", analysis_directory), ...
% %     '-dpng', '-r300');
% set(gcf, 'Position',  [0, 0, 500, 300]);
% ax = gca;
% plot_name = sprintf("%s/plate_velocity_comparison.png", analysis_directory);
% pause(0.1);
% exportgraphics(ax, plot_name, 'resolution', 300);
% 
% % sddot comparisons
% figure(4);
% plot(wagner_t, sddot, ...
%     'color', 'black', 'Linewidth', 2, 'Linestyle', '--');
% % legend([legend_entries, "Analytical"], ...
% %     "Interpreter", "latex", "location", "northeast", "Fontsize", 12);
% xlim([-impact_time t_max - impact_time]);
% grid on;
% xlabel("$t$", "Interpreter", "latex");
% ylabel("$\ddot{s}(t)$", 'Interpreter', 'latex');
% ax = gca;
% ax.FontSize = 20;
% ylim([-0.8 0.8]);
% set(gca, 'XTick', 0 : 0.2 : t_max);
% set(gca, 'YTick', -0.8 : 0.4 : 0.8);
% set(gca,'TickLabelInterpreter','latex');
% % title(['Plate acceleration: $\alpha =$ ', num2str(alpha), ', $\beta =$ ', ...
% %     num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
% %     'Fontsize', 14);
% % print(gcf, sprintf("%s/plate_acceleration_comparison.png", analysis_directory), ...
% %     '-dpng', '-r300');
% set(gcf, 'Position',  [0, 0, 500, 300]);
% ax = gca;
% plot_name = sprintf("%s/plate_acceleration_comparison.png", analysis_directory);
% pause(0.1);
% exportgraphics(ax, plot_name, 'resolution', 300);
% 
% 
