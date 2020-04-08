

q = 1;

s_analytic = @(t) ...
    q * (1 - exp(-t/2) ...
    .* (cos(sqrt(3) * t / 2) + (1/sqrt(3)) * sin(sqrt(3) * t / 2)));

log_name = "/home/michael/repos/plate-impact/ode_test/code/ode/log"

log_matrix = dlmread(log_name)

times = log_matrix(:, 1);
s_numeric = log_matrix(:, 2);

figure(1);
hold on;
plot(times, s_numeric);
plot(times, s_analytic(times));
xlabel("t");
ylabel("s(t)");
legend(["Numeric", "Analytic"]);
title("alpha, beta, gamma, q = 1");

figure(2);
plot(times, s_numeric - s_analytic(times));
xlabel("t");
ylabel("Error");
title("Error between numeric and analytic");
