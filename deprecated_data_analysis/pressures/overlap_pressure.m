function p = overlap_pressure(r, d, ddot, eps)
%OVERLAP_PRESSURE Returns the overlap pressure between the inner and outer
%   Returns the overlap function required to build a composite solution of
%   the pressure between the outer and inner regions. Requires the input of
%   r, the turnover point d and its derivative ddot, and the value of
%   epsilon, eps. 

p = 2 * sqrt(2) * d^(3/2) * ddot^2 ...
    ./ (3 * pi * sqrt(eps) * sqrt(eps * d - r));
end