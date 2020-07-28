function F = inner_force(t_vals, s_vals, sdot_vals, sddot_vals, eps)
%INNER_FORCE Returns the value of force from the inner-region 
    % Function uses the analytical expression for the force in the
    % inner-region given values of time (t_vals) and the position of s and
    % its derivatives.
    
    
    % Finds d, its derivatives and the values of J as functions of t and s
    [d_vals, ddot_vals, ~, J_vals] ...
        = s_dependents(t_vals, s_vals, sdot_vals, sddot_vals);
    
    %% Finds the value of eta_0, the upper limit of the integral
    % We find eta_0 by using a different parameterisation defined by sigma
    % = exp(2 * eta), solve for sigma_0 and then convert back to get eta_0.
    
    % Solver options
    options = optimoptions('fsolve', 'TolFun', 1e-10, 'TolX', 1e-10);
    
    % For eps -> 0, we have exp(2 eta) dominating
    sigma_0_guess = pi * d_vals ./ (eps^2 * J_vals);
    
    % Function to solve
    sigma_0_fun = @(sigmas) sigmas + 4 * sqrt(sigmas) + log(sigmas) ...
        - pi * d_vals ./ (eps^2 * J_vals);
    
    % Use fsolve to find sigma_0
    sigma_0_vals = fsolve(sigma_0_fun, sigma_0_guess, options);
    
    % Convert sigma_0 into eta_0
    eta_0_vals = 0.5 * log(sigma_0_vals);
    
    %% Determine F from the analytical expression
    F = (8 * eps^4 / pi) * ddot_vals.^2 .* J_vals.^2 .* exp(eta_0_vals) ...
        .* (pi * d_vals ./ (eps^2 * J_vals) + 1 ...
        - exp(2 * eta_0_vals) / 3 - 2 * exp(eta_0_vals) - 2 * eta_0_vals);
        
end