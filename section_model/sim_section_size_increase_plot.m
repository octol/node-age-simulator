clear all
load dat/sim_section_size_increase_net_age_12.dat

network_size = number_of_sections(1) * start_section_size(1);
assert(all((number_of_sections.*start_section_size - network_size) == 0))

for ii = 1:length(number_of_sections)
    section_size = sum(nodes{ii}.active, 2)';

    malicious_fraction_mean(ii) = mean(stats{ii}.malicious_nodes_fraction);
    malicious_fraction_std(ii) = std(stats{ii}.malicious_nodes_fraction);
    malicious_fraction_max(ii) = max(stats{ii}.malicious_nodes_fraction);

    malicious_elders_fraction_mean(ii) = mean(stats{ii}.malicious_elders_fraction);
    malicious_elders_fraction_std(ii) = std(stats{ii}.malicious_elders_fraction);
    malicious_elders_fraction_max(ii) = max(stats{ii}.malicious_elders_fraction);
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vary section size at end time
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on
plot(min_section_size, malicious_fraction_mean, 'LineWidth', 2);
plot(min_section_size, malicious_fraction_mean + malicious_fraction_std, 'LineWidth', 2);
plot(min_section_size, malicious_fraction_max, 'LineWidth', 2);
plot(min_section_size, 1/3*ones(size(min_section_size)),'k--', 'LineWidth', 2);
hold off
xlabel('Section size')
ylabel('Malicious nodes / section')
legend({'\mu', '\sigma', 'max'})
title([
    'iter: ',num2str(network_iterations), ', ',...
    'a_n: ', num2str(initial_network_age), ', ',...
    'nodes: ', num2str(network_size), ', ',...
    'adversary: ', num2str(fraction_of_new_nodes_are_malicious)
])
filename = [...
    'section_model_malicious_nodes_vs_section_size',...
    '_age_', num2str(initial_network_age), ...
    '_adversary_', num2str(fraction_of_new_nodes_are_malicious), ...
    '_n_', num2str(network_iterations), ...
    '_size_', num2str(network_size), ...
    '.png']
print(filename, '-dpng');

figure(2); clf;
hold on
plot(min_section_size, malicious_elders_fraction_mean, 'LineWidth', 2);
plot(min_section_size, malicious_elders_fraction_mean + malicious_elders_fraction_std, 'LineWidth', 2);
plot(min_section_size, malicious_elders_fraction_max, 'LineWidth', 2);
plot(min_section_size, 1/3*ones(size(min_section_size)),'k--', 'LineWidth', 2);
hold off
xlabel('Section size')
ylabel('Malicious elders / section')
legend({'\mu', '\sigma', 'max'})
title([
    'iter: ',num2str(network_iterations), ', ',...
    'a_n: ', num2str(initial_network_age), ', ',...
    'nodes: ', num2str(network_size), ', ',...
    'adversary: ', num2str(fraction_of_new_nodes_are_malicious)
])
filename = [...
    'section_model_malicious_elders_vs_section_size',...
    '_age_', num2str(initial_network_age), ...
    '_adversary_', num2str(fraction_of_new_nodes_are_malicious), ...
    '_n_', num2str(network_iterations), ...
    '_size_', num2str(network_size), ...
    '.png']
print(filename, '-dpng');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot against iterations
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii = 1:length(number_of_sections)

    section_stats = stats{ii};
    n = network_iterations;

    figure(2 + ii); clf;

    subplot(2,2,1)
    hold on
    [H1,H2,H3,H4,H5] = plot_mean_std_max(n,
        section_stats.malicious_nodes_fraction_mean,
        section_stats.malicious_nodes_fraction_std,
        section_stats.malicious_nodes_fraction_max
    );
    hold off
    xlabel("Iterations");
    ylabel("Mal. nodes / sect");
    axis([0 n 0 0.7])
    title([
        'Sections (#/min): ', num2str(number_of_sections(ii)),'/', num2str(min_section_size(ii)),...
        ', adversary: ', num2str(fraction_of_new_nodes_are_malicious)]
    )
    drawnow;

    subplot(2,2,2)
    hold on
    [H1,H2,H3,H4,H5] = plot_mean_std_max(n,
        section_stats.malicious_node_age_fraction_mean,
        section_stats.malicious_node_age_fraction_std,
        section_stats.malicious_node_age_fraction_max
    );
    hold off
    xlabel("Iterations");
    ylabel("Mal. node age / sect");
    axis([0 n 0 0.7])
    legend([H1, H2(1), H3(1), H4(1)], '\mu', '0.5\sigma', '\sigma', 'max', 'Location', 'NorthEastOutside');
    title([
        'a_n=', num2str(initial_network_age),...
        ', elders=', num2str(num_of_elders(ii))]
    )
    drawnow;

    subplot(2,2,3)
    hold on
    [H1,H2,H3,H4,H5] = plot_mean_std_max(n,
        section_stats.malicious_elders_fraction_mean,
        section_stats.malicious_elders_fraction_std,
        section_stats.malicious_elders_fraction_max
    );
    hold off
    xlabel("Iterations");
    ylabel("Mal. elders / sect");
    axis([0 n 0 0.7])
    drawnow;

    subplot(2,2,4)
    hold on
    [H1,H2,H3,H4,H5] = plot_mean_std_max(n,
        section_stats.malicious_age_fraction_mean,
        section_stats.malicious_age_fraction_std,
        section_stats.malicious_age_fraction_max
    );
    hold off
    xlabel("Iterations");
    ylabel("Mal. elder age / sect");
    axis([0 n 0 0.7])
    drawnow;

    filename = [...
        'section_model_malicious_per_section',...
        '_age_', num2str(initial_network_age), ...
        '_adversary_', num2str(fraction_of_new_nodes_are_malicious), ...
        '_section_size_', num2str(min_section_size(ii)),...
        '_no_sections_', num2str(number_of_sections(ii)),...
        '_elders_', num2str(num_of_elders(ii)), ...
        '.png']
    print(filename, '-dpng');
end

