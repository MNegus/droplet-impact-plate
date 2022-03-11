%% analytical_pressure_comparison


% Adds analytical pressures to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));


%% Parameters
close all;
eps = 1;
impact_time = 0.125
t_max = 0.8 - impact_time;
t0 = 1e-9;
tvals = linspace(t0, t_max + t0, 1e4);

% Plate parameters, a set of triples for alpha, beta, gamma
params = [[2, 0, 0]; [2, 0, 100]; [2, 2 * sqrt(2 * 100), 100]];

% Line styles
line_styles = ["--", ":", "-."];

% Analysis directory
analysis_directory ...
    = sprintf("/home/michael/Documents/supplementary_material/composite_pressure_comparison/overall_comparison");


%% Time definitions
timestep_range = 1 : floor(length(tvals) / 6)
squared_range = timestep_range.^2
squared_range = squared_range(squared_range < length(timestep_range))
tvals(squared_range(4 : 6 : end))


%% Loop over parameters
for idx = 1 : 3
    % Saves parameters
    alpha = params(idx, 1); beta = params(idx, 2); gamma = params(idx, 3);
    
    % Solve for s and dependents
    [tvals, ss, sdots, sddots] = s_solution_alt(tvals, alpha, beta, gamma, eps);
%     [ds, ddots, dddots, Js] = s_dependents(tvals, ss, sdots, sddots);

    %% Loops over time
    for m = squared_range(4 : 6 : end)
        t = tvals(m)
        %% Loops between stationary and moving plate
        for k = [-1, 1]
            if k == -1
                % Stationary plate case
                s = 0; sdot = 0; sddot = 0;
                linestyle = "-";
                linewidth = 0.5;
            else
                % Moving plate case
                s = ss(m); sdot = sdots(m); sddot = sddots(m);
                linestyle = line_styles(idx);
                linewidth = 1;
            end
            
            % s depdendents
            [d, ddot, dddot, J] = s_dependents(t, s, sdot, sddot);
            
            % Maximum value of r
            r_max = 1.25 * eps * d;
            
            % Solves for the composite pressure
            [~, ~, comp_rs, comp_ps] = outer_and_comp_pressure(t, s, ...
                sdot, sddot, 0, r_max, eps);
            
            % Plots the pressure and turnover point
            subplot(3, 1, idx);
            hold on;
            plot(k * comp_rs, comp_ps, 'Linewidth', 0.5, ...
                'color', 'black', 'Linestyle', linestyle, 'linewidth', linewidth);
            xline(k * eps * d, 'Linestyle', '--', 'Color', 0.5 * [1 1 1]);
            drawnow;
            
        end
    end
    
    subplot(3, 1, idx);
    ylim([-2 52]);
    xlim([-0.66 0.66]);
    grid on;
    ylabel("$p(r, -s(t), t)$", 'Interpreter', 'latex');
    if idx == 3
        set(gca,'xtickMode', 'auto')
        xlabel("$r$", "Interpreter", "latex");
    end
    
    if idx == 1
        txt = '$t = 0.001$';
        text(-0.001, 48, txt, "Interpreter", "Latex", "Fontsize", 9);
        
        txt = '$t = 0.108$';
        text(0.48, 15, txt, "Interpreter", "Latex", "Fontsize", 9);
    end
    
    ax = gca;
    ax.FontSize = 9;
    set(gca,'TickLabelInterpreter','latex');
end
%% Overall fig options
subplot(3, 1, 1);
ax = gca;
ax.Position = [0.1 0.74 0.8 0.2];

subplot(3, 1, 2);
ax = gca;
ax.Position = [0.1 0.53 0.8 0.2];

subplot(3, 1, 3);
ax = gca;
ax.Position = [0.1 0.32 0.8 0.2];

% for shift = [0, 0.3, 0.6]
for shift = 0.6
    % Arrow for increasing time
    X = [0.545 0.82];
    Y = [0.31 + shift, 0.17 + shift];
    annotation('arrow', X, Y, 'headwidth', 5, 'headlength', 5);

%     X = [0.47 0.16];
%     Y = [0.38 + shift, 0.145 + shift];
%     annotation('arrow', X, Y);
end


% set(gcf, 'Position',  [0, 0, 650, 650]);
plot_name = sprintf("%s/mirrored_pressure.png", analysis_directory);
f = gcf;
exportgraphics(f, plot_name, 'resolution', 300);


