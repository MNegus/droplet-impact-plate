%% analytical_pressure_comparison


% Adds analytical pressures to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));


%% Parameters
close all;
eps = 1;
% t_max = 0.8 - impact_time;
t_max = 0.1;
t0 = 1e-9;
tvals = linspace(t0, t_max + t0, 1e3);
tvals(10)

% Plate parameters, a set of triples for alpha, beta, gamma
alpha = 0.1;
gamma = 25;
beta = 2 * sqrt(alpha * gamma)
params = [[alpha, 0, 0]; [alpha, 0, gamma]; [alpha, beta, gamma]];

% Line styles
line_styles = ["--", ":", "-."];

% Analysis directory
analysis_directory ...
    = sprintf("/home/michael/Documents/supplementary_material/composite_pressure_comparison/overall_comparison_alpha_%g", alpha);


%% Time definitions

% timestep_range = 1 : 1000
% % squared_range = timestep_range.^2
% squared_range = squared_range(squared_range < length(timestep_range))
% tvals(squared_range(4 : 6 : end))


%% Loop over parameters
for idx = 1 : 3
    % Saves parameters
    alpha = params(idx, 1); beta = params(idx, 2); gamma = params(idx, 3);
    
    % Solve for s and dependents
    [tvals, ss, sdots, sddots] = s_solution_alt(tvals, alpha, beta, gamma, eps);
%     [ds, ddots, dddots, Js] = s_dependents(tvals, ss, sdots, sddots);
    
    figure(idx);
    hold on;
    
    %% Loops over time
%     for m = squared_range(4 : 6 : end)
    for m = 10 : 100 : length(tvals)
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
            figure(idx);
            plot(k * comp_rs, comp_ps, 'Linewidth', 0.5, ...
                'color', 'black', 'Linestyle', linestyle, 'linewidth', linewidth);
            xline(k * eps * d, 'Linestyle', '--', 'Color', 0.5 * [1 1 1]);
            drawnow;
            
        end
    end
    
    figure(idx);
    ylim([-2 20]);
%     xlim([-0.66 0.66]);
    xlim([-0.6, 0.6]);
    grid on;
    ylabel("$p(r, -s(t), t)$", 'Interpreter', 'latex');
    
    set(gca,'TickLabelInterpreter','latex');
    if idx == 3
        set(gca,'xtickMode', 'auto')
        xlabel("$r$", "Interpreter", "latex");
    else
        ticklabels = get(gca,'XTickLabel');
        ticklabels_new = cell(size(ticklabels));
        for i = 1:length(ticklabels)
            ticklabels_new{i} = ['\color{red} ' ticklabels{i}];
        end
        % set the tick labels
        set(gca, 'XTickLabel', ticklabels_new);
    end
    
    if idx == 1
        % Arrow for increasing time
        X = [0.65 0.8];
        Y = [0.8, 0.32];
        annotation('arrow', X, Y, 'headwidth', 5, 'headlength', 5); 
        
        txt = '$t$';
        text(0.3, 13, txt, "Interpreter", "Latex");
    end
    
%     if idx == 1
%         txt = '$t = 0.001$';
%         text(-0.001, 48, txt, "Interpreter", "Latex", "Fontsize", 9);
%         
%         txt = '$t = 0.108$';
%         text(0.48, 15, txt, "Interpreter", "Latex", "Fontsize", 9);
%     end
    
    ax = gca;
    ax.FontSize = 9;
    
    set(gcf, 'Position',  [0, 0, 650, 125]);
    plot_name = sprintf("%s/mirrored_pressure_%d.png", analysis_directory, idx);
    f = gcf;
    exportgraphics(f, plot_name, 'resolution', 300);

end
%% Overall fig options
% subplot(3, 1, 1);
% ax = gca;
% ax.Position = [0.1 0.74 0.8 0.2];
% 
% subplot(3, 1, 2);
% ax = gca;
% ax.Position = [0.1 0.53 0.8 0.2];
% 
% subplot(3, 1, 3);
% ax = gca;
% ax.Position = [0.1 0.32 0.8 0.2];
% 
% % for shift = [0, 0.3, 0.6]
% for shift = 0.6
%     % Arrow for increasing time
%     X = [0.545 0.82];
%     Y = [0.31 + shift, 0.17 + shift];
%     annotation('arrow', X, Y, 'headwidth', 5, 'headlength', 5);
% 
% %     X = [0.47 0.16];
% %     Y = [0.38 + shift, 0.145 + shift];
% %     annotation('arrow', X, Y);
% end


% set(gcf, 'Position',  [0, 0, 650, 650]);
% plot_name = sprintf("%s/mirrored_pressure.png", analysis_directory);
% f = gcf;
% exportgraphics(f, plot_name, 'resolution', 300);


