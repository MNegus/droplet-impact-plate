%% interface_analysis

addpath("interfaces");

% Wagner parameters
eps = 1;
% tvals = linspace(1e-2, 1, 1e2);
tvals = 0.5;
s = @(t) 0;
sdot = @(t) 0;
sddot = @(t) 0;

% Plotting parameters
outer_line = animatedline('color', 0.5 * [1 1 1], ...
    'Linestyle', '--', 'Linewidth', 1.5);
jet_line = animatedline('color', 0.5 * [1 1 1], ...
    'Linestyle', ':', 'Linewidth', 1.5);
inner_line = animatedline('Linewidth', 1.5);

writerObj = VideoWriter('interface_video.avi');
writerObj.FrameRate = 5;
open(writerObj);

figure(1);
for t = tvals
    % s dependent terms
    [d, ddot, dddot, J] = s_dependents(t, s(t), sdot(t), sddot(t));
    
    % Axes limits
    r_max = 2 * eps * d;
    
    figure(1);
    
    % Outer interface plotting
    rhats = linspace(d, r_max / eps, 1e2);
    [outer_rs, outer_hs] = outer_interface(rhats, t, d, s(t), eps);

    clearpoints(outer_line);
    addpoints(outer_line, outer_rs, outer_hs);
    
    % Jet interface plotting
    tau_min_fun = @(tau) tau_min_full_fun(tau, s(tau), sdot(tau), ...
        sddot(tau), t, r_max / eps);
    tau_min = fsolve(tau_min_fun, 1e-2);
    
    taus = linspace(tau_min, t, 1e3);
    [ds, ddots, dddots, Js] ...
        = s_dependents(taus, s(taus), sdot(taus), sddot(taus));
    [jet_rs, jet_hs] ...
        = jet_interface(taus, ds, ddots, dddots, Js, t, s(t), eps);

    clearpoints(jet_line);
    addpoints(jet_line, jet_rs, jet_hs);
    
    % Inner interface plotting
    [inner_rs, inner_hs] = inner_interface(0.8 * eps * d, r_max, 1000, ...
        d, s(t), J, eps);
    
    clearpoints(inner_line);
    addpoints(inner_line, inner_rs, inner_hs);
    
    
    
    % Plotting limits
    xlim([0 1.5 * max(outer_rs)]);
%     ylim([0 1.5 * max(outer_hs)]);
    
    legend(["Outer", "Jet"]);
    drawnow;
    
end

function residue = tau_min_full_fun(tau, s, sdot, sddot, t, r_bar_max)
    [d, ddot, ~, ~] = s_dependents(tau, s, sdot, sddot);
    
    residue = 2 * ddot * (t - tau) + d - r_bar_max;

end