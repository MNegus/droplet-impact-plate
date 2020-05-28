function [rs, hs] = jet_interface(taus, ds, ddots, dddots, Js, t, s, eps)
%JET_INTERFACE Returns the location of the interface of the jet

% Calculates rs
rs = eps * (2 * ddots .* (t - taus) + ds);

% Calculates hs
hs = - eps^2 * s ...
    + eps^3 * (ddots .* Js) ./ (ddots - 2 * dddots .* (t - taus));
end