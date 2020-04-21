%% comparison.m
% Script to compare the output of the C code which integrates the following
% pressure about a box of width L:
% p(r, t) = exp(-(r-t))

% Parent directory
parent_directory = "~/repos/plate-impact/force_test";

% Reads the file
log_matrix = dlmread(sprintf("%s/code/force/force_file.txt", parent_directory));
ts = log_matrix(:, 1);
F_analytic = log_matrix(:, 3);
F_numeric = log_matrix(:, 2);

% Plots the results
figure(1);
hold on;
plot(ts, F_analytic);
plot(ts, F_numeric);
xlabel("t");
ylabel("Force");
legend(["Analytic", "Numeric"]);
title("Comparison of force on plate");
print(gcf, '/home/michael/repos/plate-impact/force_test/data_analysis/comparison.png','-dpng','-r300');

figure(2);
hold on;
plot(ts, F_analytic - F_numeric);
xlabel("t");
ylabel("Error");
title("Error between analytic and numeric force");
print(gcf, '/home/michael/repos/plate-impact/force_test/data_analysis/error.png','-dpng','-r300');