function F = outer_force(t_vals, s_vals, sdot_vals, sddot_vals, eps)
%OUTER_FORCE Returns the value of force from the outer-region 
    % Function uses the analytical expression for the force in the
    % outer-region given values of time (t_vals) and the position of s and
    % its derivatives.
    
    % Finds d and its derivatives as functions of t and s
    [d_vals, ddot_vals, dddot_vals, ~] ...
        = s_dependents(t_vals, s_vals, sdot_vals, sddot_vals);
    
    
    % Determine F from the analytical expression
    F = (8 * eps / 9) * d_vals.^3 ...
        .* (4 * ddot_vals.^2 + dddot_vals .* d_vals);
end