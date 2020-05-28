function [outer_rs, outer_ps, comp_rs, comp_ps] ...
    = outer_and_comp_pressure(t, s, sdot, sddot, r_min, r_max, eps)
%%OUTER_AND_COMP_PRESSURE Returns the outer and composite pressure

% Loads in s dependent terms
[d, ddot, dddot, J] = s_dependents(t, s, sdot, sddot);

%% Determines outer region pressure
tol = 1e-8;
delta_max = log(d / tol);
idxs = linspace(0, delta_max, 1e3);
rhats = d * (1 - exp(-idxs));
outer_rs = eps * rhats;
outer_ps = outer_pressure(rhats, d, ddot, dddot, eps);

%% Determines composite pressure

% Anonymous function for r in terms of eta
indexat = @(expr, index) expr(index);
inner_r = @(eta) indexat(inner_pressure(eta, J, d, ddot, eps), 1);

% Determines minimum and maximum value of eta
options = optimoptions('fsolve', 'OptimalityTolerance', 1e-10);
min_fun = @(eta) inner_r(eta) - r_min;
min_composite_eta = fsolve(min_fun, -10, options);
    
max_fun = @(eta) inner_r(eta) - r_max;
max_composite_eta = fsolve(max_fun, 1e4, options);

% Creates array for the composite etas (bunched near min_eta)
composite_etas = max_composite_eta ...
    + (1 - exp(-linspace(100, 0, 1e4))) ...
        * (min_composite_eta - max_composite_eta);
    
% Saves composite pressure
[comp_rs, comp_ps] ...
    = composite_pressure(composite_etas, d, ddot, dddot, J, eps);

end