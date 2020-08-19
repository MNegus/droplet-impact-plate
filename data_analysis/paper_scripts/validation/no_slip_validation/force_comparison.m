%% force_comparison.m
% Script to compare the forces on the plate for where the initial height
% has been varied. The times are shifted so the theoretical time of impact
% is always the same
%

% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions


% Parent directory where all the data is held
parent_directory = '/mnt/newarre/cantilever_paper_data/';

% Directory to save the figure(s)
analysis_directory = "validation/no_slip_validation/Analysis";
analysis_directory = strcat(parent_directory, analysis_directory);

% Defines array with both directories in
data_directories = ["alpha_varying/alpha_100", ...
    "validation/no_slip_validation/alpha_100_beta_0_gamma_0"];

% Concatenates arrays to include parent directory
for k = 1 : length(data_directories)
    data_directories(k) = strcat(parent_directory, data_directories(k)); 
end

legend_entries = ["No-slip", "Free-slip"];

%% Parameters

% Value of epsilon
eps = 1;

% Plate parameters
alpha = 100;
beta = 0;
gamma = 0;

% Initial drop height 
initial_drop_heights = 0.125;

% Impact time 
impact_time = initial_drop_heights;

% Maximum time
t_max = 0.5;

%% Wagner solution

% Plate displacement solution
[wagner_t, s, sdot, sddot] = s_solution(t_max - impact_time, alpha, beta, gamma, eps);

% Composite force solution
wagner_force = composite_force(wagner_t, s, sdot, sddot, eps);

%% Plotting
close all;

% Force comparison
figure(1);
hold on;

% % s displacement comparison
% figure(2);
% hold on;

% % sdot comparison
% figure(3);
% hold on;
% 
% % sddot comparison
% figure(4);
% hold on;

for k = 1 : length(data_directories)
   output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(k)));
   
   ts = output_mat(:, 1);
   Fs = output_mat(:, 3);
   ss = output_mat(:, 6);
   sdots = output_mat(:, 7);
   sddots = output_mat(:, 8);
   
   figure(1);
   plot(ts, Fs);
   
%    figure(2);
%    plot(ts, ss);
   
%    figure(3);
%    plot(ts, sdots);
%    
%    figure(4);
%    plot(ts, sddots);
end

% Force comparison
figure(1);
plot(wagner_t + impact_time, wagner_force, ...
    'color', 'black', 'Linewidth', 2, 'Linestyle', '--');
legend([legend_entries, "Analytical"], ...
    "Interpreter", "latex", "location", "northwest", "Fontsize", 12);


xlim([0 t_max]);
grid on;
xlabel("$t$", "Interpreter", "latex");
ylabel("$F(t)$ ", 'Interpreter', 'latex');

ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title(['Force on plate: $\alpha =$ ', num2str(alpha), ', $\beta =$ ', ...
    num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
    'Fontsize', 14);
print(gcf, sprintf("%s/force_comparison.png", analysis_directory), ...
    '-dpng', '-r300');

% s comparisons
% figure(2);
% plot(wagner_t + impact_time, s, ...
%     'color', 'black', 'Linewidth', 2, 'Linestyle', '--');
% % legend([legend_entries, "Analytical"], ...
% %     "Interpreter", "latex", "location", "northwest", "Fontsize", 12);
% xlim([0 t_max]);
% ylim([0 1.1 * max(s)]);
% grid on;
% xlabel("$t$", "Interpreter", "latex");
% ylabel("$s(t)$", 'Interpreter', 'latex');
% ax = gca;
% ax.FontSize = 12;
% set(gca,'TickLabelInterpreter','latex');
% % title(['Plate displacement: $\alpha =$ ', num2str(alpha), ', $\beta =$ ', ...
% %     num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
% %     'Fontsize', 14);
% % print(gcf, sprintf("%s/plate_displacement_comparison.png", analysis_directory), ...
% %     '-dpng', '-r300');
% set(gcf, 'Position',  [0, 0, 700, 150]);
% plot_name = sprintf("%s/plate_displacement_comparison.png", analysis_directory);
% pause(0.5);
% exportgraphics(ax, plot_name, 'resolution', 300);

% % sdot comparisons
% figure(3);
% plot(wagner_t + impact_time, sdot, ...
%     'color', 'black', 'Linewidth', 2, 'Linestyle', '--');
% legend([legend_entries, "Analytical"], ...
%     "Interpreter", "latex", "location", "northeast", "Fontsize", 12);
% 
% 
% xlim([0 t_max]);
% grid on;
% xlabel("$t$", "Interpreter", "latex");
% ylabel("$\dot{s}(t)$", 'Interpreter', 'latex');
% 
% ax = gca;
% ax.FontSize = 12;
% set(gca,'TickLabelInterpreter','latex');
% % title(['Plate velocity: $\alpha =$ ', num2str(alpha), ', $\beta =$ ', ...
% %     num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
% %     'Fontsize', 14);
% print(gcf, sprintf("%s/plate_velocity_comparison.png", analysis_directory), ...
%     '-dpng', '-r300');
% 
% % sddot comparisons
% figure(4);
% plot(wagner_t + impact_time, sddot, ...
%     'color', 'black', 'Linewidth', 2, 'Linestyle', '--');
% legend([legend_entries, "Analytical"], ...
%     "Interpreter", "latex", "location", "northeast", "Fontsize", 12);
% 
% 
% xlim([0 t_max]);
% grid on;
% xlabel("$t$", "Interpreter", "latex");
% ylabel("$\ddot{s}(t)$", 'Interpreter', 'latex');
% 
% ax = gca;
% ax.FontSize = 12;
% set(gca,'TickLabelInterpreter','latex');
% % title(['Plate acceleration: $\alpha =$ ', num2str(alpha), ', $\beta =$ ', ...
% %     num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
% %     'Fontsize', 14);
% print(gcf, sprintf("%s/plate_acceleration_comparison.png", analysis_directory), ...
%     '-dpng', '-r300');