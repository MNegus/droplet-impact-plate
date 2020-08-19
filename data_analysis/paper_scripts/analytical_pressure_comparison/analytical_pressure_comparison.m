%% analytical_pressure_comparison


% Adds analytical pressures to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

% Analysis directory
analysis_directory ...
    = "/home/michael/Documents/cantilever_paper_data/Data figures";


%% Parameters
eps = 0.1;
for t = [1, 10]
    s = 0;
    sdot = 0;
    sddot = 0;

    % s depdendents
    [d, ddot, dddot, J] = s_dependents(t, 0, 0, 0);

    r_max = 1.25 * eps * d;


    %% Outer pressure
    rhats = linspace(0, d, 1e3);
    outer_rs = eps * rhats;
    outer_ps = outer_pressure(rhats, d, ddot, dddot, eps);

    %% Produces etas for inner and composite
    indexat = @(expr, index) expr(index);
    inner_r = @(eta) indexat(inner_pressure(eta, J, d, ddot, eps), 1);

    % Determines minimum and maximum value of eta
    r_min = 0;
    options = optimoptions('fsolve', 'OptimalityTolerance', 1e-10);
    min_fun = @(eta) inner_r(eta) - r_min;
    min_eta = fsolve(min_fun, -10, options);

    max_fun = @(eta) inner_r(eta) - r_max;
    max_eta = fsolve(max_fun, 1e4, options);

    % Creates array for the etas (bunched near min_eta)
    etas = max_eta ...
        + (1 - exp(-linspace(100, 0, 1e4))) ...
            * (min_eta - max_eta);

    %% Inner pressure
    [inner_rs, inner_ps] = inner_pressure(etas, J, d, ddot, eps);

    %% Composite pressure
    [composite_rs, composite_ps] ...
        = composite_pressure(etas, d, ddot, dddot, J, eps);

    %% Plotting
    close all;
    figure(1);
    hold on;
    h(3) = plot(composite_rs, composite_ps, 'Linewidth', 3, ...
        'color', 0.5 * [1 1 1]);
    h(1) = plot(outer_rs, outer_ps, 'Linewidth', 1.5, 'Linestyle', '--', ...
        'color', 'black');
    h(2) = plot(inner_rs, inner_ps, 'Linewidth', 1.5, 'Linestyle', ':', ...
        'color', 'black');
    h(4) = xline(eps * d, 'Linestyle', '--');

    legend(h([1 : 3]), ["Outer solution", "Inner solution", "Composite solution"], ...
        "Interpreter", "latex", "location", "northwest", "Fontsize", 12);

    ylim([0 1.75 * max(composite_ps)]);
    xlim([0 1.25 * eps * d]);

    grid on;
    xlabel("$r$", "Interpreter", "latex");
    ylabel("$p$", 'Interpreter', 'latex');
    ax = gca;
    ax.FontSize = 12;
    set(gca,'TickLabelInterpreter','latex');
    set(gcf, 'Position',  [0, 0, 400, 300]);
    ax = gca;

    plot_name = sprintf("%s/analytical_pressure_comparison_%g.png", analysis_directory, t);
    pause(0.1);
    exportgraphics(ax, plot_name, 'resolution', 300);
end
