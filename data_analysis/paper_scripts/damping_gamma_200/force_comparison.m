%% force_comparison.m
% Script to compare the forces on the plate for different amounts of
% damping, ranging from no damping, under-damped, critically damped and
% over-damped 
%

% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Parameters

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

% Force comparison
figure(1);
hold on;

% s displacement comparison
% figure(2);
% hold on;

for k = 1 : length(betas)
    beta = betas(k);
    
%     % Compute Wagner solution
%     % Plate displacement solution
%     [wagner_t, s, sdot, sddot] = s_solution(t_max - impact_time, alpha, beta, gamma, eps);
% 
%     % Composite force solution
%     wagner_force = composite_force(wagner_t, s, sdot, sddot, eps);
    
    output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(k)));

    ts = output_mat(:, 1);
    Fs = output_mat(:, 3);
    ss = output_mat(:, 6);

    figure(1);
    plot(ts, Fs);
%     plot(impact_time + wagner_t, wagner_force, 'color', 'black');
% 
%     figure(2);
%     plot(ts, ss);
end

% % Stationary plate solution
% stationary_plate_dir = '/mnt/newarre/cantilever_paper_data/stationary_plate';
% output_mat = dlmread(sprintf('%s/cleaned_data/output.txt', stationary_plate_dir));
% ts = output_mat(:, 1);
% Fs = output_mat(:, 3);
% figure(1);
% plot(ts, Fs, 'color', 'black', 'Linewidth', 1);

figure(1);
% plot(wagner_t + impact_time, wagner_force, ...
%     'color', 0.5 * [1 1 1], 'Linewidth', 2, 'Linestyle', '--');

% legend([legend_entries, "Stationary plate"],...
%     "Interpreter", "latex", "location", "northwest", "Fontsize", 12);


xlim([0 t_max]);
ylim([0 4.5]);
grid on;
xlabel("$t$", "Interpreter", "latex");
ylabel("$F(t)$", 'Interpreter', 'latex');

ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title(['Force on plate: $\alpha =$ ', num2str(alpha), ...
    ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
    'Fontsize', 14);
print(gcf, sprintf("%s/force_comparison.png", analysis_directory), ...
    '-dpng', '-r300');


% figure(2);
% plot(wagner_t + impact_time, s);