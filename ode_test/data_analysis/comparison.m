%% comparison.m
% Script to compare the output of the C code to the analytical solution of
% the ODE:
% alpha s''(t) + beta s'(t) + gamma s(t) = F(t)
% where F(t) = const.

%% Load in the numerical solution
log_name = "/home/michael/repos/plate-impact/ode_test/code/ode/log"

log_matrix = dlmread(log_name)

times = log_matrix(:, 1);
s_numeric = log_matrix(:, 2);

%% Parameters
% We say that alpha = 1, beta = 2 nu, gamma = omega0^2 and F = omega0^2 q,
% so we only consider nu, omega0 and q
nu = 1;
omega0 = 10;
q = 1;


%% ode45 solution
beta = 2 * nu;
gamma = omega0^2;
force = @(t) omega0^2 * q;

s_arr0 = [0, 0]; % Initial condition
dsdt = @(t, s) ...
    [s(2); ...
    (1/alpha) * (force(t) - beta * s(2) - gamma * s(1))];
times = linspace(0, 10, 1e3);
tspan = [min(times), max(times)];
[t_ode45, s_ode45] = ode45(dsdt, times, s_arr0);


%% Analytical solution
% Some extra constants needed
omega = sqrt(omega0^2 - nu^2);
s_analytic = @(t) q * (1 - exp(-nu * t) ...
    .* (cos(omega * t) + (nu / omega) * sin(omega * t)));

%%
figure;
hold on;
plot(t_ode45, s_ode45(:, 1));
plot(times, s_analytic(times));

%% Comparison plot
figure(1);
hold on;
plot(times, s_numeric);
% times = linspace(0, 100, 1e4);
% plot(times, s_analytic(times));
plot(t_ode45, s_ode45(:, 1));
xlabel("t");
ylabel("s(t)");
legend(["Basilisk", "ode45"]);
title("$\alpha, \beta, \gamma = 1$,  $F = \exp(\sin(\omega t)), \omega = 1.2$", "Interpreter", "latex");
print(gcf, '/mnt/newarre/ode_test/comparison.png','-dpng','-r300');


%% Error plot
figure(2);
% plot(times, s_numeric - s_analytic(times));
plot(times, s_numeric - s_ode45(:, 1));
xlabel("t");
ylabel("Error");
title("Error between Basilisk and ode45");
print(gcf, '/mnt/newarre/ode_test/error.png','-dpng','-r300');
