addpath("pressures");

% Parameters
eps = 1;
t = 0.1;
s  = 0;
sdot = 0;
sddot = 0;

[d, ddot, dddot, J] = s_dependents(t, s, sdot, sddot);
% d =  sqrt(3 * t);
% ddot = sqrt(3) ./ (2 * sqrt(t));
% dddot =  - sqrt(3) ./(4 * t.^(3/2));
J = 2 * d^3 / (9 * pi);
% J = 2 * t^(3/2) / (sqrt(3) * pi);


% Inner pressure
eta_max = 10;
eta_min = -log((pi * d / J )/ eps^2);
etas = linspace(eta_min, eta_max, 1e3);
[rs, inner_p] = inner_pressure(etas, J, d, ddot, eps);

% Outer pressure
outer_p = outer_pressure(rs / eps, d, ddot, dddot, eps);

% Overlap pressure
overlap_p = overlap_pressure(rs, d, ddot, eps);

% Composite pressure
comp_p = inner_p + outer_p - overlap_p;

figure(1);
hold on;
plot(rs, outer_p);
plot(rs, inner_p);
plot(rs, comp_p);
legend("Outer", "Inner", "Composite");