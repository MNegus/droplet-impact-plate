%% force_comparison.m
% Script to compare the displacements s. 
% The times are shifted so the theoretical time of impact
% is always the same
%

% Adds analytical forces to the path
addpath(genpath("~/repos/plate-impact/data_analysis"));

%% Data definitions


% Parent directory where all the data is held
parent_directory = '/media/michael/newarre/cantilever_paper_data/';

% Directory to save the figure(s)
results_directory = "/home/michael/Documents/supplementary_material/plate_displacement_animation";

% Defines array with both directories in
data_directories = ["gamma_varying/gamma_500"];

% Concatenates arrays to include parent directory
for k = 1 : length(data_directories)
    data_directories(k) = strcat(parent_directory, data_directories(k)); 
end

legend_entries = ["Computational"];


%% Parameters

% Value of epsilon
eps = 1;

% Plate parameters
alpha = 2;
beta = 0;
gamma = 500;

% Initial drop height 
initial_drop_heights = 0.125;

% Impact time 
impact_time = initial_drop_heights;

% Maximum computational time
t_max = 0.8;

%% Numerical solution 
output_mat = dlmread(sprintf("%s/cleaned_data/output.txt", data_directories(1)));
num_times = output_mat(:, 1);
num_times = num_times - impact_time;
m0 = find(~num_times)
num_s = output_mat(:, 6);

%% Wagner solution
analytical_tvals = num_times(num_times > 0);

% Plate displacement solution
[wagner_t, s, sdot, sddot] = s_solution_alt(analytical_tvals, alpha, beta, gamma, eps);

wagner_s = [zeros(m0 + 1, 1); s];

%% Animated graph
% Create the video writer with 25 fps
writerObj ...
    = VideoWriter(sprintf('%s/original_animation.avi', results_directory));
writerObj.FrameRate = 25;
open(writerObj);

close(figure(2));
figure(2);
width=2048;
height=260; % Either 257 or 512
set(gcf,'position',[10,10,width,height])

analytical_line = animatedline('Color', 0.5 * [1 1 1], 'Linewidth', 1.5, ...
        'Linestyle', '--');
numerical_line = animatedline('Color', 0.5 * [1 1 1], 'Linewidth', 1.5);
scatter_line = animatedline('Marker','o');

L = legend(["Analytical", "Numerical"]);
set(L, 'Interpreter', 'latex');
set(L, 'FontSize', 20);
set(L, 'Location', 'Northwest');
set(L, 'Numcolumns', 2);

set(gca, 'XTick', -impact_time : impact_time : t_max - impact_time);
grid on;
xlabel("$t$", "Interpreter", "latex", 'Fontsize',30);
ylabel("$s(t)$", 'Interpreter', 'latex', 'Fontsize',30);
ax = gca;
ax.FontSize = 20;
set(gca,'TickLabelInterpreter','latex');
set(gca, 'XTick', -impact_time : impact_time : t_max - impact_time);
set(gcf,'color','w');


for m = 1 : 10 : length(num_times)
    tvals = num_times(1 : m);
    wagner_svals = wagner_s(1 : m);
    num_svals = num_s(1 : m);

    clearpoints(analytical_line);
    addpoints(analytical_line, tvals, wagner_svals);

    clearpoints(numerical_line);
    addpoints(numerical_line, tvals, num_svals);
    
    clearpoints(scatter_line);
    addpoints(scatter_line, tvals(m), num_svals(m));

    
    xlim([min(num_times) max(num_times)]);
    ylim([0 1.1 * max(s)]);
    
    drawnow;
    frame = getframe(gcf);
    writeVideo(writerObj, frame);

end

close(writerObj);


