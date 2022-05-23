function p = outer_pressure(rhat, d, ddot, dddot, eps)
%OUTER_PRESSURE Returns the outer pressure
%   Returns the pressure in the outer region as a function of the outer
%   variable rhat. Requires input of the turnover point, d, and its first
%   and second derivatives ddot and dddot. eps is the value of epsilon.

p = (1/eps) * (4 * (2 * d^2 - rhat.^2) * ddot^2 ...
        ./ (3 * pi * sqrt(d^2 - rhat.^2)) ...
    + 4 * d * dddot * sqrt(d^2 - rhat.^2) / (3 * pi));
end