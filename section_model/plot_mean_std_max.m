function [H1, H2, H3, H4, H5] = plot_mean_std_max(n, mean_vec, std_vec, max_vec)
H1 = plot(1:100:n, mean_vec(1:100:end), 'LineWidth', 2, 'Color', 'k');
H2 = plot(
    1:100:n,
    [max(mean_vec(1:100:end) - 0.5*std_vec(1:100:end), 0); mean_vec(1:100:end) + 0.5*std_vec(1:100:end)],
    'LineWidth', 2, 'Color', 'b'
);
H3 = plot(
    1:100:n,
    [max(mean_vec(1:100:end) - std_vec(1:100:end), 0); mean_vec(1:100:end) + std_vec(1:100:end)],
    'LineWidth', 2, 'Color', 'm'
);
H4 = plot(1:100:n, max_vec(1:100:end), 'LineWidth', 2, 'Color', 'g');
H5 = plot(1:100:n, 1/3*ones(size(1:100:n)), 'LineWidth', 2, 'k--');
end
