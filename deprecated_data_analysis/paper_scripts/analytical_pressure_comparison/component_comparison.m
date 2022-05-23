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

fontsize = 25;

% Plate parameters, a set of triples for alpha, beta, gamma
alpha = 0.1;
gamma = 25;
beta = 2 * sqrt(alpha * gamma);
params = [[alpha, 0, 0]; [alpha, 0, gamma]; [alpha, beta, gamma]];
% % 
% alpha = 2;
% gamma = 100;
% beta = 2 * sqrt(alpha * gamma)
% params = [[alpha, 0, 0]; [alpha, 0, gamma]; [alpha, beta, gamma]];

% Line styles
line_styles = ["--", ":", "-."];

% Analysis directory
analysis_directory ...
    = sprintf("/home/michael/Documents/supplementary_material/composite_pressure_comparison/overall_comparison_alpha_%g", alpha);

% Defines data directories for all the values
no_params = 3;

%% Creates figures
% Origin pressure comparison
subplot(1, 4, 1);
hold on;
[ds, ddots, dddots, Js] = s_dependents(tvals, 0, 0, 0);
p0_stationary = (4 * ds / (3 * pi)) .* (2 * ddots.^2 + ds .* dddots) / eps;
plot(tvals, p0_stationary, 'color', 0 * [1 1 1], 'Linewidth', 2); % Plot stationary value

% Maximum pressure comparison
subplot(1, 4, 2);
hold on;
plot(tvals, ddots.^2 / (2 * eps^2), 'color', 0 * [1 1 1], ...
    'Linewidth', 2);

% Plate displacement comparison
subplot(1, 4, 3);
hold on;

% Turnover point comparison
subplot(1, 4, 4);
hold on;
plot(tvals, sqrt(3 * tvals), 'color', 0 * [1 1 1], ...
    'Linewidth', 2); % Plot stationary value


%% Loops over plate parameters
for idx = 1 : 3
    
    line_color = 'black';
    
    alpha = params(idx, 1); beta = params(idx, 2); gamma = params(idx, 3);
    
    %% Solve for displacement and dependents
    [tvals, ss, sdots, sddots] = s_solution_alt(tvals, alpha, beta, gamma, eps);
    [ds, ddots, dddots, Js] = s_dependents(tvals, ss, sdots, sddots);
    
    %% Plot origin pressure
    subplot(1, 4, 1);
    p0_moving = (4 * ds / (3 * pi)) .* (2 * ddots.^2 + ds .* dddots) / eps;
    plot(tvals, p0_moving, 'Linestyle', line_styles(idx), ...
        'color', line_color, 'Linewidth', 2);
    
    %% Plot the maximum pressure
    subplot(1, 4, 2);
    plot(tvals, ddots.^2 / (2 * eps^2), 'color', line_color, ...
        'Linestyle', line_styles(idx), 'Linewidth', 2);
    
    %% Plot plate displacement
    subplot(1, 4, 3);
    plot(tvals, ss, 'linestyle', line_styles(idx), 'color', line_color, ...
        'Linewidth', 2);
    
    %% Plot turnover point
    subplot(1, 4, 4);
    plot(tvals, ds, 'color', line_color, ...
        'Linestyle', line_styles(idx), 'Linewidth', 2);
    
end


%% Origin pressure
subplot(1, 4, 1);
ylim([-1, 4]);
ylabel("$p(0, -s(t), t)$", "Interpreter", "latex");
xlabel("$t$", "Interpreter", "latex");
grid on;
ax = gca;
ax.FontSize = fontsize;
set(gca,'TickLabelInterpreter','latex');
set(gca, 'XTick', 0 : 0.2 : t_max);

if alpha == 0.1
    subplot(1, 4, 1);
    rectangle('position',[0 -1 0.1 5], 'Curvature',0.3, ...
       'Edgecolor', 0.5 * [1 1 1], 'Linewidth', 2);
end

%% Maximum pressure
subplot(1, 4, 2);
ylim([-1, 7]);
ylabel("$p_{max}(t)$", "Interpreter", "latex");
xlabel("$t$", "Interpreter", "latex");
grid on;
ax = gca;
ax.FontSize = fontsize;
set(gca,'TickLabelInterpreter','latex');
set(gca, 'XTick', 0 : 0.2 : t_max);
if alpha == 0.1
    subplot(1, 4, 2);
    height = 8
    rectangle('position',[0 -1 0.1 height], 'Curvature',0.3, ...
       'Edgecolor', 0.5 * [1 1 1], 'Linewidth', 2);
end



%% Plate displacement
subplot(1, 4, 3);
xlabel("$t$", "Interpreter", "latex");
ylabel("$s(t)$", "Interpreter", "latex");
ylim([0, 0.475]);
grid on;
ax = gca;
ax.FontSize = fontsize;
set(gca,'TickLabelInterpreter','latex');
set(gca, 'XTick', 0 : 0.2 : t_max);

%% Turnover point
% figure(3);
subplot(1, 4, 4);
ylabel("$d(t)$", "Interpreter", "latex");
xlabel("$t$", "Interpreter", "latex");
grid on;
ax = gca;
ax.FontSize = fontsize;
% set(gca, 'Yscale', 'log');
set(gca,'TickLabelInterpreter','latex');
set(gca, 'XTick', 0 : 0.2 : t_max);


%% Save triple figure
set(gcf, 'Position',  [0, 0, 1800, 400]);
plot_name = sprintf("%s/component_comparison_alpha_%g.png", analysis_directory, alpha);
f = gcf;
exportgraphics(f, plot_name, 'resolution', 300);

%% Legend
close(figure(2));
figure(2);
hold on;
plot([0 1], [0, 1], 'color', line_color);
for q = 1 : 3
    plot([0 1], [0, 1], 'color', line_color, ...
        'Linestyle', line_styles(q), 'Linewidth', 1);
end

if alpha == 2
    L = legend(["Stationary plate", ...
        "Moving plate: $\beta = 0, \gamma = 0$", ...
        "Moving plate: $\beta = 0, \gamma = 100$", ...
        "Moving plate: $\beta = 28.28, \gamma = 100$"], ...
        'location', 'northoutside', 'interpreter', 'latex', ...
        'Numcolumns', 2);
else
    L = legend(["Stationary plate", ...
    "Moving plate: $\beta = 0, \gamma = 0$", ...
    "Moving plate: $\beta = 0, \gamma = 25$", ...
    "Moving plate: $\beta = 3.16, \gamma = 25$"], ...
    'location', 'northoutside', 'interpreter', 'latex', ...
    'Numcolumns', 2);
end
plot_name = sprintf("%s/legend_full.png", analysis_directory);
ax = gca;
exportgraphics(ax, plot_name, 'resolution', 500);
