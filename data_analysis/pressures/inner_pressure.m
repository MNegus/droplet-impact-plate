function [r, p] = inner_pressure(eta, J, d, ddot, eps)
%INNER_PRESSURE Returns the inner pressure
%   Returns the pressure in the inner region as a function of the parameter
%   eta, which can go from -infinity to +infinity. Requires input of the
%   jet-thickness J, the turnover point d, the time derivative of the
%   turnover point ddot and eps, the value of epsilon. 

% Inner variable for r
tilde_r = (J / pi) * (eta - 1 - exp(-eta) - 4 * exp(-eta/2));

% Value for r
r = eps * d + eps^3 * tilde_r;

% Value for pressure
p = (1 / eps^2) * 2 * ddot^2 * exp(-eta / 2) ./ (1 + exp(-eta / 2)).^2;

end