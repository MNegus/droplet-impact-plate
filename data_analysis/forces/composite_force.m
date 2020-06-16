function F = composite_force(t_vals, s_vals, sdot_vals, sddot_vals, eps)
%COMPOSITE_FORCE Returns the force from the integrated composite pressure
    
    F = zeros(length(t_vals));
    
    for k = 1 : length(t_vals)
        t = t_vals(k);
        s = s_vals(k);
        sdot = sdot_vals(k);
        sddot = sddot_vals(k);
        

        % Finds d to work out maximum r value
        [d, ~, ~, ~] = s_dependents(t, s, sdot, sddot);

        % Minimum and maximum values of r to integrate over
        r_min = 0;
        r_max = 1.25 * eps * d;

        % Calculate composite rs and ps
        [~, ~, comp_rs, comp_ps] ...
            = outer_and_comp_pressure(t, s, sdot, sddot, r_min, r_max, eps);

        % Integrates the composite expansion
        F(k) = trapz(comp_rs, 2 * pi * comp_rs .* comp_ps);
    end

end