function [rs, hs] = outer_interface(rhats, t, d, s, eps)
%OUTER_INTERFACE Returns the location of the interface in the outer region
%   Returns the location, hs, of the interface at rs. Input is rhats (the
%   outer region coordinate), time t, turnover point d, plate dispacement s
%   and epsilon, eps. 

% Calculates rs
rs = eps * rhats;

% Calculates hs
hs = eps^2 * (0.5 * rhats.^2 - t ...
    + (2/pi) * (asin(d ./ rhats) .* (t - s - 0.5 * rhats.^2) ...
        + 0.5 * d * sqrt(rhats.^2 - d^2)));
end