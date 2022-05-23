function origin_pressure(parent_dir, output_range)
    
    p0 = zeros(length(output_range), 1);
    for k = output_range
        output_mat = dlmread(sprintf("%s/cleaned_data/plate_outputs/output_%d.txt", ...
           parent_dir, k));
        rs = output_mat(:, 1);
        ps = output_mat(:, 3);
        [min_r, min_idx] = min(rs);
        [k, min_r]
        p0(k) = ps(min_idx);
    end
    
    dlmwrite(sprintf('%s/cleaned_data/origin_pressure.txt', parent_dir), ...
        p0);

    plot(output_range, p0);
end