function [d, ddot, dddot, J] = s_dependents(t, s, sdot, sddot)
%S_DEPENDENTS Returns values of the functions that are dependent on s
%   Values for the turnover point (d), its first derivative (ddot), second
%   derivative (dddot), and the jet-thickness (J), are returned as
%   functions of the time (t), cantilever displacement (s), its first
%   derivative (sdot) and second derivative (sddot).

% Turnover point 
d = sqrt(3 * (t - s));

% First derivative of turnover point
ddot = (sqrt(3)/2) * (1 - sdot) ./ sqrt(t - s);

% Second derivative of turnover point
dddot = - (sqrt(3)/4) * ((1 - sdot).^2 + 2 * (t - s) .* sddot)...
    ./ (t - s).^(3/2);

% Jet thickness
J = 2 * (t - s).^(3/2) / (sqrt(3) * pi);

end