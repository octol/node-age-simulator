clear all
load dat/sim_network_growth_net_age_12_adv_0.2_section_100_elders_100.dat

for ii = 1:length(network_growth_rate)

    section_stats = stats{ii};
    current_number_of_sections = size(nodes{ii}.active, 1);

    figure(ii); clf;

    subplot(2,2,1)
    hold on
    [H1,H2,H3,H4,H5] = plot_mean_std_max(network_iterations,
        section_stats.malicious_nodes_fraction_mean,
        section_stats.malicious_nodes_fraction_std,
        section_stats.malicious_nodes_fraction_max
    );
    hold off
    xlabel("Iterations");
    ylabel("Mal. nodes / sect");
    axis([0 network_iterations 0 0.7])
    title([
        'Sections (#/min): ', num2str(number_of_sections),'/', num2str(min_section_size),...
        ', adversary: ', num2str(fraction_of_new_nodes_are_malicious),
        'Current sections: ', num2str(current_number_of_sections),', ',...]
        'Growth rate: ', num2str(network_growth_rate(ii))]
    )
    drawnow;

    subplot(2,2,2)
    hold on
    [H1,H2,H3,H4,H5] = plot_mean_std_max(network_iterations,
        section_stats.malicious_node_age_fraction_mean,
        section_stats.malicious_node_age_fraction_std,
        section_stats.malicious_node_age_fraction_max
    );
    hold off
    xlabel("Iterations");
    ylabel("Mal. node age / sect");
    axis([0 network_iterations 0 0.7])
    legend([H1, H2(1), H3(1), H4(1)], '\mu', '0.5\sigma', '\sigma', 'max', 'Location', 'NorthEastOutside');
    title([
        'a_n=', num2str(initial_network_age),...
        ', elders=', num2str(num_of_elders)]
    )
    drawnow;

    subplot(2,2,3)
    hold on
    [H1,H2,H3,H4,H5] = plot_mean_std_max(network_iterations,
        section_stats.malicious_elders_fraction_mean,
        section_stats.malicious_elders_fraction_std,
        section_stats.malicious_elders_fraction_max
    );
    hold off
    xlabel("Iterations");
    ylabel("Mal. elders / sect");
    axis([0 network_iterations 0 0.7])
    drawnow;

    subplot(2,2,4)
    hold on
    [H1,H2,H3,H4,H5] = plot_mean_std_max(network_iterations,
        section_stats.malicious_age_fraction_mean,
        section_stats.malicious_age_fraction_std,
        section_stats.malicious_age_fraction_max
    );
    hold off
    xlabel("Iterations");
    ylabel("Mal. elder age / sect");
    axis([0 network_iterations 0 0.7])
    drawnow;

    filename = [...
        'section_model_with_growth_malicious_per_section',...
        '_net_growth_', num2str(network_growth_rate(ii)), ...
        '_age_', num2str(initial_network_age), ...
        '_adversary_', num2str(fraction_of_new_nodes_are_malicious), ...
        '_section_size_', num2str(min_section_size),...
        '_no_sections_', num2str(number_of_sections),...
        '_elders_', num2str(num_of_elders), ...
        '.png']
    print(filename, '-dpng');
end

