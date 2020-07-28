

addpath("forces");

eps = 1;
t_max = 0.8;
alpha = 1;
beta = 0;
gamma = 200;

beta_critical = 2 * sqrt(alpha * gamma);

close all;

figure(1);
hold on;

figure(2);
hold on;

figure(3);
hold on;

figure(4);
hold on;
for beta = [0, 0.25 * beta_critical, 0.5 * beta_critical, beta_critical, 2.0 * beta_critical]
% for beta = [0, 0.5 * beta_critical, beta_critical, 2.0 * beta_critical]
   [t, s, sdot, sddot] = s_solution(t_max, alpha, beta, gamma, eps);
   Fs = composite_force(t, s, sdot, sddot, eps);
   figure(1);
   plot(t, Fs);
   
   figure(2);
   plot(t, s);
   
   figure(3);
   plot(t, sdot);
   
   figure(4);
   plot(t, sddot);
end

% Stationary
s = zeros(size(t));
sdot = zeros(size(t));
sddot = zeros(size(t));
Fs = composite_force(t, s, sdot, sddot, eps);
figure(1);
plot(t, Fs);


figure(1);
legend("0", "0.1", "0.5", "1", "2", "Stationary");

figure(2);
legend("0", "0.1", "0.5", "1", "2", "Stationary");

figure(3);
legend("0", "0.1", "0.5", "1", "2", "Stationary");

figure(4);
legend("0", "0.1", "0.5", "1", "2", "Stationary");

% %%
% 
% addpath("forces");
% 
% eps = 1;
% t_max = 1.0;
% alpha = 1;
% beta = 0;
% gamma = 50;
% 
% [t, s, sdot, sddot] = s_solution(t_max, alpha, beta, gamma, eps);
% Fs = composite_force(t, s, sdot, sddot, eps);
% 
% %%
% figure(1);
% plot(t, Fs);
% 
% figure(2);
% plot(t, s);
% 
% figure(3);
% plot(t, sdot);
% 
% figure(4);
% plot(t, sddot);