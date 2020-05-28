function [r, p] = composite_pressure(eta, d, ddot, dddot, J, eps)
%COMPOSITE_PRESSURE Returns the composite pressure
%   Returns the composite pressure, which is the sum of the outer and inner
%   pressure, minus the overlap function. Requires the input of the
%   parameter eta, the turnover point and its first and second derivatives
%   (d, ddot and dddot), the jet-thickness J and the value of epsilon, eps.

% Finds r as a function of eta, and the inner pressure
[r, inner_p] = inner_pressure(eta, J, d, ddot, eps);

% Produces the composite
p = real(outer_pressure(r / eps, d, ddot, dddot, eps)) + inner_p ...
    - real(overlap_pressure(r, d, ddot, eps));

end