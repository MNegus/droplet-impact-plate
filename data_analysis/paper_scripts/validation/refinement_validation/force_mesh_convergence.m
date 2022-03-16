% force_comparison.m
% Script to compare the forces on the plate for where the initial height
% has been varied. The times are shifted so the theoretical time of impact
% is always the same
%

% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions
close all;
fig_no = 1;


close all;
for validation_case = ["stationary", "moving"]

    % Parent directory where all the data is held
    parent_directory = sprintf('/media/michael/newarre/supplementary_material/validation/%s_validation/', ...
        validation_case);

    % Directory to save the figure(s)
    analysis_directory = "/home/michael/Documents/supplementary_material/validation";
    % analysis_directory = strcat(parent_directory, analysis_directory);

    % Defines array with both directories in
    data_directories = ["level_10", "level_11", "level_12", "level_13", "level_14"];

    % Concatenates arrays to include parent directory
    for k = 1 : length(data_directories)
        data_directories(k) = strcat(parent_directory, data_directories(k)); 
    end

%     legend_entries = ["Max level = 10", "Max level = 11", "Max level = 12", ...
%         "Max level = 13", "Max level = 14"];
    legend_entries = ["$m = 10$", "$m = 11$", "$m = 12$", "$m = 13$", "$m = 14$"];
    % Line colors
    start_color = 0.2;
    end_color = 0.6;
    no_dirs = length(data_directories);
    m = (end_color - start_color) / (no_dirs - 1);
    c = 0.5 * (start_color + end_color - (no_dirs + 1) * m);
    line_colors = (m * (1 : no_dirs) + c)' .* [1 1 1];

    %% Parameters

    % Value of epsilon
    eps = 1;

    % Plate parameters alpha = 2, beta = 7.07, gamma = 100
    alpha = 2;
    beta = 7.07;
    gamma = 100;

    % Impact time 
    impact_time = 0.125;

    % Maximum time
    t_max = 0.8 - impact_time;
    
    %% Set up convergence testing
    % Saves the highest refinement data
    output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(end)));
    ts = output_mat(:, 1);
    Fs_refined = output_mat(:, 3);

    % Computational time values
    tvals = ts - impact_time;
    
    % Number of timesteps
    N = length(tvals)
    
    mean_errors = zeros(length(data_directories) - 1, 1);

    %% Plotting
    % Force comparison
    figure(fig_no);
    hold on;

    p(no_dirs) = plot(ts - impact_time, Fs_refined, 'color', line_colors(no_dirs, :), ...
        'Linewidth', 1.5);
    
    for k = length(data_directories) - 1 : -1 : 1
        output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(k)));

        ts = output_mat(:, 1);
        Fs = output_mat(:, 3);

        p(k) = plot(ts - impact_time, Fs, 'color', line_colors(k, :), 'Linewidth', 1.5);

        mean_errors(k) = norm(Fs - Fs_refined) / sqrt(N);
    end
    
    ax1 = gca;
    legend(p(1:no_dirs), legend_entries, ...
        "Interpreter", "latex", "location", "northwest", ...
        "Fontsize", 14);


    xlim([-impact_time - 0.21 t_max]);
    grid on;
    xlabel("$t$", "Interpreter", "latex");
    ylabel("$F_m(t)$", 'Interpreter', 'latex');
    
    font_size = 16;
    ax1.FontSize = font_size;
    set(gca,'TickLabelInterpreter','latex');
    % title(['Force on plate: $\alpha =$ ', num2str(alpha), ', $\beta =$ ', ...
    %     num2str(beta), ', $\gamma =$ ' num2str(gamma)], "Interpreter", "latex", ...
    %     'Fontsize', 14);


    ax_width = 0.28;
    ax2 = axes('Position',[.625 .25 ax_width 1.3 * ax_width]);
    
    box on;
    hold on;
    plot(10 : 13, mean_errors, 'color', 'black');
    sz = 75;
    scatter(10 : 13, mean_errors, sz, line_colors(1 : 4, :), 'filled');
    xlim([9.5, 13.5]);
    ylim([0 0.08]);
%     set(gca, 'Yscale', 'log');
    set(gca,'TickLabelInterpreter','latex');
    xticks([10, 11, 12, 13]);
    xlabel("$m$", "Interpreter", "latex");
    ylabel("$||F_m - F_{14}||_2$", "Interpreter", "latex");
    ax2.FontSize = font_size;
    grid on;
    
    set(gcf, 'Position',  [0, 0, 450, 500]);
    
    plot_filename = sprintf("%s/%s_force_validation.png", analysis_directory, validation_case);
    pause(0.5);
    print(gcf, plot_filename, ...
        '-dpng', '-r300');
    
    
    fig_no = fig_no + 1;

end
