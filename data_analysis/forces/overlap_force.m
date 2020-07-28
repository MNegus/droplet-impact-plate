function F = overlap_force(t_vals, s_vals, sdot_vals, sddot_vals, eps)
%OVERLAP_FORCE Returns the value of force from the overlap region
    % Function uses the analytical expression for the force in the
    % overlap-region given values of time (t_vals) and the position of s and
    % its derivatives.
    
    % Finds d and its derivatives as functions of t and s
    [d_vals, ddot_vals, ~, ~] ...
        = s_dependents(t_vals, s_vals, sdot_vals, sddot_vals); 
    
    % Determine F from the analytical expression
    F = (16 * sqrt(2) * eps / 9) * d_vals.^3 .* ddot_vals.^2;
end