%% pressure_comparison
% Compares the outputted pressure with the Wagner prediction

% Location of data
parent_directory = "/mnt/newarre/low_alpha";

impact_time = 0.13;

%% Load in times and s values
output_matrix = ...
    dlmread(sprintf("%s/cleaned_data/volumes.txt", parent_directory));
times = output_matrix(:, 1);
s_num = output_matrix(:, 4);
sdot_num = output_matrix(:, 5);
sddot_num = output_matrix(:, 5);


figure(1);

width=800;
height=800;
set(gcf,'position',[10,10,width,height])
% create the video writer with 1 fps
writerObj = VideoWriter(sprintf('%s/data_analysis/Videos/pressure.avi', parent_directory));
writerObj.FrameRate = 5;
open(writerObj);

for m = 100 : 200
    
    t = times(m); % Time value
    
    % Loads in pressure file
    pressure_matrix = ...
        dlmread(sprintf("%s/cleaned_data/plate_outputs/output_%d.txt", ...
            parent_directory, m));
        
    % Sorts in increasing order of r
    [~, sorted_idxs] = sort(pressure_matrix(:, 1));
    sorted_mat = pressure_matrix(sorted_idxs, :);

    % Saves values of r and pressure
    rs = sorted_mat(:, 1);
    ps = sorted_mat(:, 3);
    
    if t > impact_time
       sigmas =  10.^linspace(-10, 5, 1e4);
       [wagner_rs, wagner_ps, wagner_pmax] = wagner_pressure(sigmas, ...
           t - impact_time, s_num(m), sdot_num(m), sddot_num(m), 1);
    end
    
    figure(1);
    plot(rs, ps);
    hold on;
    if t > impact_time
        plot(wagner_rs, wagner_ps);
    end
    hold off;
    xlim([0, 1]);
    ylim([0, 10]);
    ylabel("p");
    xlabel("r");
    legend(["Numerical", "Semi-Wagner"]);
    title(sprintf("Semi-Wagner. t = %.3f", t));
    drawnow;
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    
end
close(writerObj);


