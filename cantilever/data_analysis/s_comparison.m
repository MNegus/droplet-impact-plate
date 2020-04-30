%% s_comparison.m
% Compares the output s values from the simulation to the theoretical 
% values given by Wagner theory 

%% Parameters
% Parameters the specific simulation was run with
alpha = 0.01;
beta = 0;
gamma = 0;
impact_time = 0.115; % GUESS

% Location of data
parent_directory = "/mnt/newarre/low_alpha";

%% Loads in data
output_matrix = ...
    dlmread(sprintf("%s/cleaned_data/volumes.txt", parent_directory));
times = output_matrix(:, 1);
s_num = output_matrix(:, 4);
sdot_num = output_matrix(:, 5);
sddot_num = output_matrix(:, 5);

%% Wagner solution (assuming beta = gamma = 0)
t_wag = linspace(impact_time, max(times), 1e3);
zero_fun = @(s) alpha * s - (8 * sqrt(3) / 5) * ((t_wag - impact_time) - s).^(5/2);
s_wag = fsolve(zero_fun, zeros(size(t_wag)));


%% Plots s
figure(1);
hold on;
plot(times, s_num);
plot(t_wag, s_wag);
ylabel("s(t)");
xlabel("t");
legend(["Numerical", "Wagner"]);
title("Plotting $s(t)$", "Interpreter", "latex");
print(gcf, sprintf('%s/data_analysis/Figures/s.png', parent_directory),'-dpng','-r300');


%% Plots sdot
figure(2);
hold on;
plot(times, sdot_num);
plot(t_wag(2:end), diff(s_wag) ./ diff(t_wag - impact_time) );
ylabel("s'(t)");
xlabel("t");
legend(["Numerical", "Wagner"]);
title("Plotting $s'(t)$", "Interpreter", "latex");
print(gcf, sprintf('%s/data_analysis/Figures/s_dot.png', parent_directory),'-dpng','-r300');

%% Plots sddot
sdot = diff(s_wag) ./ diff(t_wag - impact_time);
sddot = diff(sdot) ./ diff(t_wag(2:end) - impact_time);
figure(3);
plot(times, sddot_num);
hold on
% plot(t_wag(3:end), sddot);
ylabel("s''(t)");
xlabel("t");
% legend(["Numerical", "Wagner"]);
title("Plotting $s''(t)$", "Interpreter", "latex");
print(gcf, sprintf('%s/data_analysis/Figures/s_ddot.png', parent_directory),'-dpng','-r300');
