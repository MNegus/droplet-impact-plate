function F = composite_force(t_vals, s_vals, sdot_vals, sddot_vals, eps)
%COMPOSITE_FORCE Returns the value of the composite force
    % Function uses the analytical expression for the composite expansion
    % force given the values of the outer, inner and overlap force
    
    F = outer_force(t_vals, s_vals, sdot_vals, sddot_vals, eps) ...
        + inner_force(t_vals, s_vals, sdot_vals, sddot_vals, eps) ...
        - overlap_force(t_vals, s_vals, sdot_vals, sddot_vals, eps);
end