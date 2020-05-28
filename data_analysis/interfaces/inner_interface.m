function [rs, hs] = inner_interface(r_min, r_max, no_points, d, s, J, eps)
%INNER_INTERFACE Returns the interface location in the inner region

    % Inner variables for r and h, parameterised by eta 
    inner_r = @(eta) (J / pi) * (exp(eta) - eta - 1);
    inner_h = @(eta) J * (1 + 4 * exp(eta / 2) / pi);
    
    % Function for the non-dimensionalised r and h, parameterised by eta
    r = @(eta) eps * d + eps^3 * inner_r(eta);
    h = @(eta) - eps^2 * s + eps^3 * inner_h(eta);
    
    % Finds min and max value of eta
    options = optimoptions('fsolve', 'OptimalityTolerance', 1e-8);
    min_eta_fun = @(eta) r(eta) - r_min;
    min_eta = fsolve(min_eta_fun, -10)
    
    max_eta_fun = @(eta) r(eta) - r_max;
    max_eta = fsolve(max_eta_fun, 10)
    
    % Creates etas matrix
    etas = max_eta + (1 - exp(-linspace(100, 0, no_points))) ...
        * (min_eta - max_eta);
    
    % Returns rs and hs
    rs = r(etas);
    hs = h(etas);

end