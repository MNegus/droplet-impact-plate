%% force_comparison.m
% Script to compare the forces on the plate for different amounts of
% damping, ranging from no damping, under-damped, critically damped and
% over-damped 
%

% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Parameters

% Physical parameters
rho_w = 998;
R0 = 1e-3;
U0 = 5;
T0 = R0 / U0;
Patm = 10^5;

% Dimensional functions
F_newton = @(F) rho_w * U0^2 * R0^2 * F;
t_milli = @(t) T0 * 1000 * t;
r_milli = @(r) R0 * 1000 * r;

% Value of epsilon
eps = 1;

% Plate parameters
alpha = 1;
betas = [0, 7.07, 14.14, 28.28, 56.57];
gamma = 200;

% Initial drop height
initial_drop_height = 0.125;

% Impact time
impact_time = initial_drop_height;

% Maximum time
t_max = 0.5;


%% Data definitions

% Parent directory where all the data is held
parent_directory = '/mnt/newarre/cantilever_paper_data/damping_gamma_200/';

% Directory to save the figure(s)
analysis_directory = "Analysis";
analysis_directory = strcat(parent_directory, analysis_directory);

% Defines arrays for all the values of beta
data_directories = string(size(betas));
legend_entries = string(size(betas));
for k = 1 : length(betas)
    beta = betas(k);
    
    data_directories(k) = [parent_directory, '/beta_', num2str(beta)]
%     
%     data_directories(k) = sprintf("%s/alpha_%d_beta_%d_gamma_%d", ...
%         parent_directory, alpha, beta, gamma);
    
    legend_entries(k) = ['$\beta =$ ', num2str(beta)] ;
end

%% Plotting
close all;

% % Force comparison
% figure(1);
% hold on;
% 
% % s displacement comparison
% figure(2);
% hold on;
% 
% % sdot comparison
% figure(3);
% hold on;
% 
% % sddot comparison
% figure(4);
% hold on;

for k = 1 : length(betas)
    beta = betas(k);
    
    % Compute Wagner solution
    % Plate displacement solution
    [wagner_t, s, sdot, sddot] = s_solution(t_max - impact_time, alpha, beta, gamma, eps);

    % Composite force solution
    wagner_force = composite_force(wagner_t, s, sdot, sddot, eps);
    
    output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(k)));

    ts = output_mat(:, 1);
    Fs = output_mat(:, 3);
    ss = output_mat(:, 6);
    sdots = output_mat(:, 7);
    sddots = output_mat(:, 8);

    % Force comparison
    figure(1);
    plot(t_milli(ts), F_newton(Fs), 'Displayname', "Computational");
    hold on;
    plot(t_milli(wagner_t + impact_time), F_newton(wagner_force), ...
        'Displayname', "Analytical", ...
        'Linestyle', '--', 'Linewidth', 1, 'Color', 'black');
    hold off;
    legend("Interpreter", "latex", "location", "northwest", "Fontsize", 12);
    xlim([0 t_milli(t_max)]);
    grid on;
    xlabel("$t^*$ / ms", "Interpreter", "latex");
    ylabel("$F^*(t^*)$ / N ", 'Interpreter', 'latex');
    ax = gca;
    ax.FontSize = 12;
    set(gca,'TickLabelInterpreter','latex');
    title(['Force on plate: $\alpha =$ ', ...
        num2str(alpha), '$, \beta =$', num2str(beta), ...
        ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
        'Fontsize', 14);
    print(gcf, sprintf("%s/force_comparisons/force_comparison_beta_%.2f.png", analysis_directory, beta), ...
        '-dpng', '-r300');

    % s comparison
    figure(2);
    plot(t_milli(ts), r_milli(ss), 'Displayname', "Computational");
    hold on;
    plot(t_milli(wagner_t + impact_time), r_milli(s), ...
        'Displayname', "Analytical", ...
        'Linestyle', '--', 'Linewidth', 1, 'Color', 'black');
    hold off;
    legend("Interpreter", "latex", "location", "northwest", "Fontsize", 12);
    xlim([0 t_milli(t_max)]);
    grid on;
    xlabel("$t^*$ / ms", "Interpreter", "latex");
    ylabel("$s^*(t^*)$ / mm", 'Interpreter', 'latex');
    ax = gca;
    ax.FontSize = 12;
    set(gca,'TickLabelInterpreter','latex');
    title(['Plate displacement: $\alpha =$ ', ...
        num2str(alpha), '$, \beta =$ ', num2str(beta), ...
        ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
        'Fontsize', 14);
    print(gcf, sprintf("%s/s_comparisons/s_comparison_beta_%.2f.png", analysis_directory, beta), ...
        '-dpng', '-r300');
    
    % sdot comparison
    figure(3);
    plot(t_milli(ts), sdots * U0, 'Displayname', "Computational");
    hold on;
    plot(t_milli(wagner_t + impact_time), sdot * U0, ...
        'Displayname', "Analytical", ...
        'Linestyle', '--', 'Linewidth', 1, 'Color', 'black');
    hold off;
    legend("Interpreter", "latex", "location", "northwest", "Fontsize", 12);
    xlim([0 t_milli(t_max)]);
    grid on;
    xlabel("$t^*$ / ms", "Interpreter", "latex");
    ylabel("$\dot{s}^*(t^*)$ / m/s", 'Interpreter', 'latex');
    ax = gca;
    ax.FontSize = 12;
    set(gca,'TickLabelInterpreter','latex');
    title(['Plate velocity: $\alpha =$ ', ...
        num2str(alpha), '$, \beta =$ ', num2str(beta), ...
        ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
        'Fontsize', 14);
    print(gcf, sprintf("%s/sdot_comparisons/sdot_comparison_beta_%.2f.png", analysis_directory, beta), ...
        '-dpng', '-r300');
    
    % sddot comparison
    figure(3);
    plot(t_milli(ts), sddots * U0 / T0, 'Displayname', "Computational");
    hold on;
    plot(t_milli(wagner_t + impact_time), sddot * U0 / T0, ...
        'Displayname', "Analytical", ...
        'Linestyle', '--', 'Linewidth', 1, 'Color', 'black');
    hold off;
    legend("Interpreter", "latex", "location", "northeast", "Fontsize", 12);
    xlim([0 t_milli(t_max)]);
    grid on;
    xlabel("$t^*$ / ms", "Interpreter", "latex");
    ylabel("$\ddot{s}^*(t^*)$ / m/s$^2$", 'Interpreter', 'latex');
    ax = gca;
    ax.FontSize = 12;
    set(gca,'TickLabelInterpreter','latex');
    title(['Plate acceleration: $\alpha =$ ', ...
        num2str(alpha), '$, \beta =$ ', num2str(beta), ...
        ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
        'Fontsize', 14);
    print(gcf, sprintf("%s/sddot_comparisons/sddot_comparison_beta_%.2f.png", analysis_directory, beta), ...
        '-dpng', '-r300');
    
    
    
%     figure(2);
%     plot(t_milli(ts), r_milli(ss), 's', 'Displayname', legend_entries(k));
% 
%     figure(3);
%     plot(t_milli(ts), sdots * U0, 's', 'Displayname', legend_entries(k));
% 
%     figure(4);
%     plot(t_milli(ts), sddots * U0 / T0, 's', 'Displayname', legend_entries(k));
end
% 
% % Force comparison
% figure(1);
% legend("Interpreter", "latex", "location", "northwest", "Fontsize", 12);
% xlim([t_milli(-0.1) t_milli(t_max)]);
% grid on;
% xlabel("$t^*$ / ms", "Interpreter", "latex");
% ylabel("$F^*(t^*)$ / N ", 'Interpreter', 'latex');
% ax = gca;
% ax.FontSize = 12;
% set(gca,'TickLabelInterpreter','latex');
% title(['Force on plate: $\alpha =$ ', ...
%     num2str(alpha), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
%     'Fontsize', 14);
% print(gcf, sprintf("%s/force_comparison.png", analysis_directory), ...
%     '-dpng', '-r300');
% 
% % s comparisons
% figure(2);
% legend([legend_entries], ...
%     "Interpreter", "latex", "location", "northwest", "Fontsize", 12);
% xlim([0 t_milli(t_max)]);
% grid on;
% xlabel("$t^*$ / ms", "Interpreter", "latex");
% ylabel("$s^*(t^*)$ / mm", 'Interpreter', 'latex');
% ax = gca;
% ax.FontSize = 12;
% set(gca,'TickLabelInterpreter','latex');
% title(['Plate displacement: $\alpha =$ ', ...
%     num2str(alpha), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
%     'Fontsize', 14);
% print(gcf, sprintf("%s/plate_displacement_comparison.png", analysis_directory), ...
%     '-dpng', '-r300');
% 
% 
% % sdot comparisons
% figure(3);
% legend([legend_entries], ...
%     "Interpreter", "latex", "location", "northwest", "Fontsize", 12);
% xlim([0 t_milli(t_max)]);
% grid on;
% xlabel("$t^*$ / ms", "Interpreter", "latex");
% ylabel("$\dot{s}^*(t^*)$ / m/s", 'Interpreter', 'latex');
% ax = gca;
% ax.FontSize = 12;
% set(gca,'TickLabelInterpreter','latex');
% title(['Plate velocity: $\alpha =$ ', ...
%     num2str(alpha), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
%     'Fontsize', 14);
% print(gcf, sprintf("%s/plate_velocity_comparison.png", analysis_directory), ...
%     '-dpng', '-r300');
% 
% % sddot comparisons
% figure(4);
% legend([legend_entries], ...
%     "Interpreter", "latex", "location", "northwest", "Fontsize", 12);
% xlim([t_milli(-0.05) t_milli(t_max)]);
% grid on;
% xlabel("$t^*$ / ms", "Interpreter", "latex");
% ylabel("$\ddot{s}^*(t^*)$ / m/s$^2$", 'Interpreter', 'latex');
% ax = gca;
% ax.FontSize = 12;
% set(gca,'TickLabelInterpreter','latex');
% title(['Plate acceleration: $\alpha =$ ', ...
%     num2str(alpha), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
%     'Fontsize', 14);
% print(gcf, sprintf("%s/plate_acceleration_comparison.png", analysis_directory), ...
%     '-dpng', '-r300');