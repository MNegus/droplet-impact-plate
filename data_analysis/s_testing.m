alpha = 1;
gamma = 1;
t_max = 0.8;

figure(1);
hold on;
% [t, s, sdot, sddot] = s_solution(t_max, alpha, beta, gamma, eps);
% plot(t, s);
gamma = 1000;
[t, s, sdot, sddot] = s_solution(t_max, alpha, beta, gamma, eps)
plot(t, s);
