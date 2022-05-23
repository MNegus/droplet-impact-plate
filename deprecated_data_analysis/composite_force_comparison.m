%% composite_force_comparison.m
% This script compares two different ways of determining the composite
% force. First is to get the composite pressure equations and integrate 
% numerically using MATLAB. The second is integrating analytically. The
% analytical integration still requires some numerical solving, given that
% the lower limit of the integral in terms of eta has to be determined, but
% this is a much cheaper computation than the trapeze integration

% Parameters
eps = 1e-2; % Value of epsilon
tvals = linspace(1e-6, 10, 1e3)'; % Time interval for plotting
no_points = 1e4; % Number of spatial points to integrate over


% Forms of d(t), its derivatives ddot(t) and dddot(t) and J(t)
d = @(t) sqrt(3 * t);
ddot = @(t) sqrt(3) ./ (2 * sqrt(t));
dddot = @(t) - sqrt(3) ./(4 * t.^(3/2));
J = @(t) 2 * d(t).^3 / (9 * pi);


% Arrays to store the values of the integrated forces
outer_force_integral = zeros(length(tvals), 1);
inner_force_integral = zeros(length(tvals), 1);
composite_force_integral = zeros(length(tvals), 1);

% Solution for eta_0 
eta_0_full_fun = @(eta, t) eps * d(t) ...
    - (eps^3 * J(t) / pi) * (exp(2 * eta) + 4 * exp(eta) + 2 * eta + 1);
eta_0_vals = zeros(length(tvals), 1);

% Analytical pressure functions
p_outer = @(rhat, t) real((1/eps) * (4 * d(t) * dddot(t) * sqrt(d(t)^2 - rhat.^2) / (3 * pi) ...
    - 4 * (rhat.^2 - 2 * d(t)^2) * ddot(t)^2 ./ (3 * pi * sqrt(d(t)^2 - rhat.^2))));
p_inner = @(eta, t) (1 / eps^2) * 2 * ddot(t)^2 * exp(eta) ./ (1 + exp(eta)).^2;
tilde_r = @(eta, t) -(J(t)/pi) * (exp(2 * eta) + 4 * exp(eta) + 2 * eta + 1);
inner_r = @(eta, t) eps * d(t) + eps^3 * tilde_r(eta, t);
p_overlap = @(r, t) 2 * sqrt(2) * d(t)^(3/2) * ddot(t)^2 ...
        ./ (3 * pi * sqrt(eps) * sqrt(eps * d(t) - r));

%%
% Loops over all times
for k = 1 : length(tvals)
    
    % Value of t
    t = tvals(k);
    
    % Values of r
    rs = eps * d(t) + eps^3 * tilde_r(etas, t);

    % Outer force
    
    outer_integrand = @(r) 2 * pi * r .* p_outer(r / eps, t);
    outer_force_integral(k) = quadgk(outer_integrand, 0, eps * d(t));
    
    
    % Inner force. Done by creating an interpolate function
    
    % Solving for eta_0
    eta_0_fun = @(eta) eta_0_full_fun(eta, t);
    eta_0_vals(k) = fsolve(eta_0_fun, 10);
    
    % Array for etas
    eta_min = -100;
    etas = linspace(eta_0_vals(k), eta_min, no_points);
    
    % Interpolated inner
    interp_p_inner = @(eta_query) interp1(etas, p_inner(etas, t), eta_query);
    
    % Calculation of force
    inner_integrand = @(eta) 2 * pi * inner_r(eta, t) .* interp_p_inner(eta) * eps^3;
    inner_force_integral(k) = quadgk(inner_integrand, eta_min, eta_0_vals(k)); 

    % Composite array
    p_composite = p_outer(rs / eps, t) + p_inner(etas, t) - p_overlap(rs, t);
    
    % Forces
%     inner_force_integral(k) = trapz(rs, 2 * pi * rs .* p_inner(etas, t));
    composite_force_integral(k) = trapz(rs, 2 * pi * rs .* p_composite);
end


%% Analytical
outer_force = eps * (8/9) * d(tvals).^3 .* (4 * ddot(tvals).^2 + dddot(tvals) .* d(tvals));
inner_force = ((8 * eps^4 * ddot(tvals).^2 .* J(tvals).^2 .* exp(eta_0_vals)) / pi) ...
        .* (pi * d(tvals) ./ (eps^2 * J(tvals)) + 1 - exp(2 * eta_0_vals) / 3 ...
            - 2 * exp(eta_0_vals) - 2 * eta_0_vals);
overlap_force = eps * 16 * sqrt(2) * d(tvals).^3 .* ddot(tvals).^2 / 9;
composite_force = outer_force + inner_force -overlap_force;


%% Force comparisons
close all;
analysis_dir = '/home/michael/Desktop/Wednesday 22nd July/Quadgk force integral';

figure(1);
hold on;
plot(tvals, outer_force, 'color', 0.75 * [1 1 1], 'Linewidth', 3);
plot(tvals, outer_force_integral);
grid on;
xlabel("$t$", "Interpreter", "latex");
ylabel("$F(t)$", 'Interpreter', 'latex');
legend(["Analytical", "Trapezoidal"], ...
    "Interpreter", "latex", "location", "northwest", "Fontsize", 12);
ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title(['Outer region force. $\epsilon = $', num2str(eps)], "Interpreter", "latex", ...
    'Fontsize', 14);
print(gcf, sprintf("%s/outer_force_comparison.png", analysis_dir), ...
    '-dpng', '-r300');

figure(2);
hold on;
plot(tvals, inner_force, 'color', 0.75 * [1 1 1], 'Linewidth', 3);
plot(tvals, inner_force_integral);
grid on;
xlabel("$t$", "Interpreter", "latex");
ylabel("$F(t)$", 'Interpreter', 'latex');
legend(["Analytical", "Trapezoidal"], ...
    "Interpreter", "latex", "location", "northwest", "Fontsize", 12);
ax = gca;
ax.FontSize = 12;
set(gca,'TickLabelInterpreter','latex');
title(['Inner region force. $\epsilon = $', num2str(eps)], "Interpreter", "latex", ...
    'Fontsize', 14);
print(gcf, sprintf("%s/inner_force_comparison.png", analysis_dir), ...
    '-dpng', '-r300');

% figure(3);
% hold on;
% plot(tvals, composite_force, ...
%     'color', 0.75 * [1 1 1], 'Linewidth', 3);
% plot(tvals, composite_force_integral);
% grid on;
% xlabel("$t$", "Interpreter", "latex");
% ylabel("$F(t)$", 'Interpreter', 'latex');
% legend(["Analytical", "Trapezoidal"], ...
%     "Interpreter", "latex", "location", "northwest", "Fontsize", 12);
% ax = gca;
% ax.FontSize = 12;
% set(gca,'TickLabelInterpreter','latex');
% title(['Composite force. $\epsilon = $', num2str(eps)], "Interpreter", "latex", ...
%     'Fontsize', 14);
% print(gcf, sprintf("%s/composite_force_comparison.png", analysis_dir), ...
%     '-dpng', '-r300');
% 
% %% Computational comparison
% impact_time = 0.125;
% output_filename = '/mnt/newarre/supervision_runs_15_july/solid_wall_peak_removal/cleaned_data/output.txt';
% output_mat = dlmread(output_filename);
% simulation_ts = output_mat(:, 1);
% simulation_Fs = output_mat(:, 3);
% 
% figure(4);
% hold on;
% plot(simulation_ts, simulation_Fs);
% plot(tvals + impact_time, outer_force, 'Linestyle', '--', ...
%     'Linewidth', 1, 'color', [0 0 0]);
% plot(tvals + impact_time, composite_force, ...
%     'Linewidth', 1, 'color', [0 0 0] );
% grid on;
% xlabel("$t$", "Interpreter", "latex");
% ylabel("$F(t)$", 'Interpreter', 'latex');
% legend(["Computational", "Outer", "Composite"], ...
%     "Interpreter", "latex", "location", "northwest", "Fontsize", 12);
% ax = gca;
% ax.FontSize = 12;
% set(gca,'TickLabelInterpreter','latex');
% title(['Force on plate: Stationary plate'], "Interpreter", "latex", ...
%     'Fontsize', 14);
% print(gcf, sprintf("%s/computational_comparison.png", analysis_dir), ...
%     '-dpng', '-r300');
