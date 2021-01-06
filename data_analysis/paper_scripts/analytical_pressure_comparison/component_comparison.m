%% analytical_pressure_comparison


% Adds analytical pressures to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));


%% Parameters
close all;
eps = 1;
t0 = 1e-9;
tvals = linspace(t0, t_max + t0, 1e4);
impact_time = 0.125;
t_max = 0.8 - impact_time;
% output_range = 1 : 800;
% origin_num_ts = 1e-3 * output_range - impact_time;

% Plate parameters, a set of triples for alpha, beta, gamma
params = [[2, 0, 0]; [2, 0, 100]; [2, 2 * sqrt(2 * 100), 100]];

% Line styles
line_styles = ["--", ":", "-."];

% Analysis directory
analysis_directory ...
    = sprintf("/home/michael/Documents/supplementary_material/composite_pressure_comparison/overall_comparison");

% Defines data directories for all the values
no_params = 3;
% data_directories = string(no_params);
% data_directories(1) = "/media/michael/newarre/cantilever_paper_data/alpha_varying/alpha_2";
% data_directories(2) = "/media/michael/newarre/cantilever_paper_data/beta_varying/beta_0";
% data_directories(3) = "/media/michael/newarre/cantilever_paper_data/beta_varying/beta_28.28";


% % Stationary plate data
% stationary_plate_directory = '/media/michael/newarre/cantilever_paper_data/stationary_plate';
% output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", ...
%     stationary_plate_directory));
% ts = output_mat(:, 1);
% stationary_ts = ts - impact_time;
% stationary_ss = output_mat(:, 6);
% turnover_mat = dlmread(sprintf("%s/cleaned_data/turnover_points.txt", ...
%     stationary_plate_directory));
% stationary_turnover_ts = 1e-3 * turnover_mat(:, 1) - impact_time; 
% stationary_ds = turnover_mat(:, 2);
% num_p0s = dlmread(sprintf("%s/cleaned_data/origin_pressure.txt", ...
%         stationary_plate_directory));
%     
%% Creates figures
% Plate displacement comparison
% figure(1); 
subplot(1, 3, 2);
hold on;

% Origin pressure comparison
% figure(2);
subplot(1, 3, 1);
hold on;
[ds, ddots, dddots, Js] = s_dependents(tvals, 0, 0, 0);
% p0_stationary = sqrt(3 ./ tvals) / (pi * eps);
p0_stationary = (4 * ds / (3 * pi)) .* (2 * ddots.^2 + ds .* dddots) / eps;
plot(tvals, p0_stationary, 'color', 0 * [1 1 1], 'Linewidth', 2); % Plot stationary value
% plot(origin_num_ts, num_p0s, 'color', 0 * [1 1 1]);


% Turnover point comparison
% figure(3);
subplot(1, 3, 3);
hold on;
plot(tvals, sqrt(3 * tvals).^2, 'color', 0 * [1 1 1], ...
    'Linewidth', 2); % Plot stationary value
% plot(stationary_turnover_ts, stationary_ds.^2, 'color', 0 * [1 1 1]);


%% Loops over plate parameters
for idx = 1 : 3
    
%     line_color = (0.2 + 0.1 * idx) * [1 1 1];
    line_color = 'black';
    
    alpha = params(idx, 1); beta = params(idx, 2); gamma = params(idx, 3);
    
%     output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", ...
%         data_directories(idx)));
%     ts = output_mat(:, 1);
%     num_ts = ts - impact_time;
%     num_ss = output_mat(:, 6);
    
%     turnover_mat = dlmread(sprintf("%s/cleaned_data/turnover_points.txt", ...
%         data_directories(idx)));
%     turnover_ts = 1e-3 * turnover_mat(:, 1) - impact_time; 
%     num_ds = turnover_mat(:, 2);
%     
%     num_p0s = dlmread(sprintf("%s/cleaned_data/origin_pressure.txt", ...
%         data_directories(idx)));
%     

    %% Solve for displacement and dependents
    [tvals, ss, sdots, sddots] = s_solution_alt(tvals, alpha, beta, gamma, eps);
    [ds, ddots, dddots, Js] = s_dependents(tvals, ss, sdots, sddots);
    
    %% Plot plate displacement
%     figure(1);
    subplot(1, 3, 2);
    plot(tvals, ss, 'linestyle', line_styles(idx), 'color', line_color, ...
        'Linewidth', 2);
%     plot(num_ts, num_ss, 'color', line_color);
    
    %% Plot origin pressure
%     figure(2);
    subplot(1, 3, 1);
    p0_moving = (4 * ds / (3 * pi)) .* (2 * ddots.^2 + ds .* dddots) / eps;
    plot(tvals, p0_moving, 'Linestyle', line_styles(idx), ...
        'color', line_color, 'Linewidth', 2);
%     plot(origin_num_ts, num_p0s, 'color', line_color);
    
    %% Plot turnover point
%     figure(3);
    subplot(1, 3, 3);
    plot(tvals, ds.^2, 'color', line_color, ...
        'Linestyle', line_styles(idx), 'Linewidth', 2);
%     plot(turnover_ts, num_ds.^2, 'color', line_color);
    
    
end

%% Plate displacement
% figure(1);
subplot(1, 3, 2);
xlabel("$t$", "Interpreter", "latex");
ylabel("$s(t)$", "Interpreter", "latex");
grid on;
plot_name = sprintf("%s/plate_displacement.png", analysis_directory);
ax = gca;
ax.FontSize = 30;
set(gca,'TickLabelInterpreter','latex');
set(gca, 'XTick', 0 : 0.2 : t_max);
% ax = gca;
% exportgraphics(ax, plot_name, 'resolution', 300);

%% Origin pressure
% figure(2);
subplot(1, 3, 1);
ylim([-1, 4]);
ylabel("$p(0, -s(t), t)$", "Interpreter", "latex");
xlabel("$t$", "Interpreter", "latex");
grid on;
plot_name = sprintf("%s/origin_pressure.png", analysis_directory);
ax = gca;
ax.FontSize = 30;
set(gca,'TickLabelInterpreter','latex');
set(gca, 'XTick', 0 : 0.2 : t_max);
% ax = gca;
% exportgraphics(ax, plot_name, 'resolution', 300);

%% Turnover point
% figure(3);
subplot(1, 3, 3);
ylabel("$d(t)^2$", "Interpreter", "latex");
xlabel("$t$", "Interpreter", "latex");
grid on;
plot_name = sprintf("%s/turnover_point.png", analysis_directory);
ax = gca;
ax.FontSize = 30;
set(gca,'TickLabelInterpreter','latex');
% set(gcf, 'Position',  [0, 0, 500, 700]);
set(gca, 'XTick', 0 : 0.2 : t_max);
% L = legend(["Stationary plate", ...
%     "Moving plate: $\alpha = 2, \beta = 0, \gamma = 0$", ...
%     "Moving plate: $\alpha = 2, \beta = 0, \gamma = 100$", ...
%     "Moving plate: $\alpha = 2, \beta = 28.28, \gamma = 100$"], ...
%     'location', 'southoutside', 'interpreter', 'latex', ...
%     'Numcolumns', 2);
% newPosition = [0.5 0.5 0 0];
% newUnits = 'normalized';
% set(L,'Position', newPosition,'Units', 'normalized');
% ax = gca;
% exportgraphics(ax, plot_name, 'resolution', 300);

%% Save triple figure
set(gcf, 'Position',  [0, 0, 1500, 700]);
plot_name = sprintf("%s/analytical_result_comparison.png", analysis_directory);
f = gcf;
exportgraphics(f, plot_name, 'resolution', 300);
