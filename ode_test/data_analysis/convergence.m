%% convergence.m
% Determines the convergence of the numerical scheme solving the following
% ode:
% s''(t) + 2 nu s'(t) + omega0^2 s(t) = omega0^2 * q
% 

%% Parameters
raw_data_dir = "~/repos/plate-impact/ode_test/raw_data";
analysis_dir = "~/repos/plate-impact/ode_test/data_analysis";

dt0 = 0.1; % First timestep
no_runs = 13; % Number of runs

nu = 1.; 
omega0 = 10.; 
q = 1.; 
omega = sqrt(omega0^2 - nu^2);
MAX_TIME = 5;

dts = zeros(no_runs, 1);
max_errors = zeros(no_runs, 1);

% Analytic solution
t_analytic = linspace(0, MAX_TIME, 1e3);
s_analytic = @(t) ...
    q * (1 - exp(-nu * t) .* (cos(omega * t) + (nu / omega) * sin(omega * t)));

%% Finds the max errors
for run_no = 0 : no_runs - 1
% for run_no = 8
    dts(run_no + 1) = dt0 / 2^run_no;
    
    output_matrix = dlmread(sprintf('%s/output_%d.txt', raw_data_dir, run_no));
    ts = output_matrix(:, 1);
    numeric = output_matrix(:, 2);
    analytic = output_matrix(:, 3);
    
    errors = abs(numeric - analytic);
    max_errors(run_no + 1) = max(errors);
end

p = polyfit(log(dts), log(max_errors), 1);
fit = @(dt) exp(p(2)) * dt.^p(1);

figure;
loglog(dts, fit(dts), dts, max_errors, 's');
grid on;

legend(sprintf("$y = %.2f x + %.2f$", p(1), p(2)), "Numerical", ...
    "Interpreter", "Latex", "Fontsize", 12, 'location', 'northwest'); 
xlabel("$\Delta t$", "Interpreter", "latex", 'Fontsize', 15);
ylabel("Error", "Interpreter", "latex", 'Fontsize', 15);
title("Convergence of numerical scheme for ODE", "Interpreter", "latex", 'Fontsize', 15);
set(gca,'TickLabelInterpreter', 'latex', 'Fontsize', 13);

%% Create animation of solution
% Creates video object
writerObj = VideoWriter(sprintf('%s/solution.avi', analysis_dir));
writerObj.FrameRate = 1;
open(writerObj);


figure(1);

for run_no = 0 : no_runs - 1
    % Reads values from output files
    output_matrix = dlmread(sprintf('%s/output_%d.txt', raw_data_dir, run_no));
    ts = output_matrix(:, 1);
    s_numeric = output_matrix(:, 2);
    
    % Plots analytic solution
    figure(1);
    plot(t_analytic, s_analytic(t_analytic), 'Linewidth', 5, 'Color', [0.75 0.75 0.75]);
    hold on;
    
    % Plots numerical solution
    plot(ts, s_numeric, 'color', 'red');
    hold off
    
    % Sets up axes
    grid on;
    legend(["Analytic", "Numeric"], 'Interpreter', 'latex');
    title(sprintf('Solution for timestep = $0.1 / 2^{%d}$', run_no), ...
        'Interpreter', 'latex', 'Fontsize', 14);
    xlabel("$t$", "Interpreter", "latex", 'Fontsize', 15)
    ylabel("$s(t)$", "Interpreter", "latex", 'Fontsize', 15);
    set(gca,'TickLabelInterpreter', 'latex', 'Fontsize', 13);
    
    % Draws to video writer
    drawnow;
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
end
close(writerObj);

%% Creates animation of error
% Creates video object
writerObj = VideoWriter(sprintf('%s/error.avi', analysis_dir));
writerObj.FrameRate = 1;
open(writerObj);

figure(2);

for run_no = 0 : no_runs - 1
    % Reads values from output files
    output_matrix = dlmread(sprintf('%s/output_%d.txt', raw_data_dir, run_no));
    ts = output_matrix(:, 1);
    s_numeric = output_matrix(:, 2);
    
    % Plots error
    plot(ts, s_numeric - s_analytic(ts));
    hold off
    
    % Sets up axes
    grid on;
    title(sprintf('Error for timestep size = $0.1 / 2^{%d}$', run_no), ...
        'Interpreter', 'latex', 'Fontsize', 14);
    xlabel("$t$", "Interpreter", "latex", 'Fontsize', 15)
    ylabel("Numerical - Analytical", "Interpreter", "latex", 'Fontsize', 15);
    set(gca,'TickLabelInterpreter', 'latex', 'Fontsize', 13);
    
    % Draws to video writer
    drawnow;
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
end
close(writerObj);
