%% comparison.m
% Script to compare the output of the C code to the analytical solution of
% the ODE:
% alpha s''(t) + beta s'(t) + gamma s(t) = F(t)

%% Load in the numerical solution
log_name = "/home/michael/repos/plate-impact/ode_test/code/ode/log"

log_matrix = dlmread(log_name)

times = log_matrix(:, 1);
s_numeric = log_matrix(:, 2);

%% Parameters
alpha = 1;
beta = 1;
gamma = 1;
omega = 1.2;
force = @(t) exp(sin(omega * t.^1.5));

%% ode45 solution
s_arr0 = [0, 0]; % Initial condition
dsdt = @(t, s) ...
    [s(2); ...
    (1/alpha) * (force(t) - beta * s(2) - gamma * s(1))];
tspan = [min(times), max(times)];
[t_ode45, s_ode45] = ode45(dsdt, times, s_arr0);


% %% Analytical solution: 
% % alpha = beta = gamma = 1, F(t) = q, 
% q = 1;
% s_analytic = @(t) ...
%     q * (1 - exp(-t/2) ...
%     .* (cos(sqrt(3) * t / 2) + (1/sqrt(3)) * sin(sqrt(3) * t / 2)));
% 
% %% Analytical solution: 
% % alpha = 1, beta = nu, gamma = omega0^2, F = omega0^2 F0 sin(omega t)
% % (Damped, forced harmonic oscillator)
% nu = 0.1;
% F0 = 1;
% omega0 = 1;
% omega = 1.1;
% 
% 
% % Analytic solution is s(t) = R sin(omega t - phi)
% phip = atan(nu * omega / (omega0^2 - omega^2))
% Rp = omega0^2 * F0 / sqrt((omega0^2 - omega^2)^2 + nu^2 * omega^2)
% if omega > omega0
%     Rp = - Rp;
% end
% 
% phih = atan(omega0 * tan(phip) / omega)
% if phip == 0
%     Rh = -Rp
% else
%     Rh = -Rp * sin(phip) / sin(phih)
% end
% 
% s_analytic = @(t) ...
%     Rh * sin(omega0 * t - phih) + Rp * sin(omega * t - phip);
% 
% %% nu = 0
% nu = 0;
% F0 = 1;
% omega0 = 1;
% omega = 1.1;
% 
% s_analytic = @(t) (omega0 * F0) / (omega0^2 - omega^2) ...
%     * (-omega * sin(omega0 * t) + omega0 * sin(omega * t));
%% Comparison plot
figure(1);
hold on;
plot(times, s_numeric);
% times = linspace(0, 100, 1e4);
% plot(times, s_analytic(times));
plot(t_ode45, s_ode45(:, 1));
xlabel("t");
ylabel("s(t)");
legend(["Numeric", "ode45"]);
% title("alpha, beta, gamma, q = 1");

%% Error plot
figure(2);
% plot(times, s_numeric - s_analytic(times));
plot(times, s_numeric - s_ode45(:, 1));
xlabel("t");
ylabel("Error");
title("Error between numeric and analytic");
