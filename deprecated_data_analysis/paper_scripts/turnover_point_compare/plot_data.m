%% Data definitions

% Parent directory where all the data is held
parent_directory = '/mnt/newarre/cantilever_paper_data/alpha_varying/';

% Directory to save the figure(s)
turnover_directory = "turnover_points";
turnover_directory = strcat(parent_directory, turnover_directory);

% Range of alphas
alphas = [2, 5, 10, 20, 100];

% Timestep value
dt = 1e-3;

t_max = 0.6;


figure(1);
hold on;
for k = 1:length(alphas)
    alpha = alphas(k);
    output_dir = dlmread(sprintf("%s/alpha_%d.txt", ...
        turnover_directory, alpha));
    
    ts = dt * output_dir(:, 1);
    ds = output_dir(:, 2);
    
    % Shortens according to t_max
    ds = ds(ts < t_max);
    ts = ts(ts < t_max);
    
    plot(ts, ds);
    
end