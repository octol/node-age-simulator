clear all
%load dat/sim_section_size_increase_net_age_12.dat
load sim_section_size_increase_large_dataset_fix_elders_10_net_age_12_adv_0.3.dat

network_size = number_of_sections(1) * start_section_size(1);
%assert(all((number_of_sections.*start_section_size - network_size) == 0))
%max(abs(number_of_sections.*start_section_size - network_size))

N = 1:network_iterations;

% Sample points
n = [1 2 5 10 20]*network_iterations/20;

for ii = 1:length(number_of_sections)
    section_size = sum(nodes{ii}.active, 2)';

    malicious_fraction_max(ii, 1) = max(stats{ii}.malicious_nodes_fraction_max(1:n(1)));
    malicious_fraction_max(ii, 2) = max(stats{ii}.malicious_nodes_fraction_max(1:n(2)));
    malicious_fraction_max(ii, 3) = max(stats{ii}.malicious_nodes_fraction_max(1:n(3)));
    malicious_fraction_max(ii, 4) = max(stats{ii}.malicious_nodes_fraction_max(1:n(4)));
    malicious_fraction_max(ii, 5) = max(stats{ii}.malicious_nodes_fraction_max(1:n(5)));

    malicious_elders_fraction_max(ii, 1) = max(stats{ii}.malicious_elders_fraction_max(1:n(1)));
    malicious_elders_fraction_max(ii, 2) = max(stats{ii}.malicious_elders_fraction_max(1:n(2)));
    malicious_elders_fraction_max(ii, 3) = max(stats{ii}.malicious_elders_fraction_max(1:n(3)));
    malicious_elders_fraction_max(ii, 4) = max(stats{ii}.malicious_elders_fraction_max(1:n(4)));
    malicious_elders_fraction_max(ii, 5) = max(stats{ii}.malicious_elders_fraction_max(1:n(5)));
end

% Find threshold iterations
first_iteration_above_threshold = network_iterations*ones(size(number_of_sections));
for ii = 1:length(number_of_sections)
    first_above_threshold = min(find(stats{ii}.malicious_nodes_fraction_max > 1/3));
    if length(first_above_threshold) > 0
        first_iteration_above_threshold(ii) = first_above_threshold;
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Max values when varying section size for different times
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); clf;
hold on
plot(min_section_size, malicious_fraction_max(:,1), 'LineWidth', 2);
plot(min_section_size, malicious_fraction_max(:,2), 'LineWidth', 2);
plot(min_section_size, malicious_fraction_max(:,3), 'LineWidth', 2);
plot(min_section_size, malicious_fraction_max(:,4), 'LineWidth', 2);
plot(min_section_size, malicious_fraction_max(:,5), 'LineWidth', 2);
plot(min_section_size, 1/3*ones(size(min_section_size)),'k--', 'LineWidth', 2);
hold off
xlabel('Section size')
ylabel('Malicious nodes / section')
legend({'1k','2k','5k','10k','20k'})
%legend({'\mu', '\sigma', 'max'})
axis([0 200 0 0.7])
title([
    'iter: ',num2str(network_iterations), ', ',...
    'a_n: ', num2str(initial_network_age), ', ',...
    'nodes: ', num2str(network_size), ', ',...
    'fix elders: ', num2str(num_of_elders(1)), ', ',...
    'adversary: ', num2str(fraction_of_new_nodes_are_malicious)
])
filename = [...
    'section_model_max_malicious_nodes_vs_section_size',...
    '_fix_elders_', num2str(num_of_elders(1)), ...
    '_age_', num2str(initial_network_age), ...
    '_adversary_', num2str(fraction_of_new_nodes_are_malicious), ...
    '_n_', num2str(network_iterations), ...
    '_size_', num2str(network_size), ...
    '.png'
]
print(filename, '-dpng');

figure(2); clf;
hold on
plot(min_section_size, malicious_elders_fraction_max(:,1), 'LineWidth', 2);
plot(min_section_size, malicious_elders_fraction_max(:,2), 'LineWidth', 2);
plot(min_section_size, malicious_elders_fraction_max(:,3), 'LineWidth', 2);
plot(min_section_size, malicious_elders_fraction_max(:,4), 'LineWidth', 2);
plot(min_section_size, malicious_elders_fraction_max(:,5), 'LineWidth', 2);
plot(min_section_size, 1/3*ones(size(min_section_size)),'k--', 'LineWidth', 2);
hold off
xlabel('Section size')
ylabel('Malicious elders / section')
legend({'1k','2k','5k','10k','20k'})
axis([0 200 0 0.7])
title([
    'iter: ',num2str(network_iterations), ', ',...
    'a_n: ', num2str(initial_network_age), ', ',...
    'nodes: ', num2str(network_size), ', ',...
    'fix elders: ', num2str(num_of_elders(1)), ', ',...
    'adversary: ', num2str(fraction_of_new_nodes_are_malicious)
])
filename = [...
    'section_model_max_malicious_elders_vs_section_size',...
    '_fix_elders_', num2str(num_of_elders(1)), ...
    '_age_', num2str(initial_network_age), ...
    '_adversary_', num2str(fraction_of_new_nodes_are_malicious), ...
    '_n_', num2str(network_iterations), ...
    '_size_', num2str(network_size), ...
    '.png'
]
print(filename, '-dpng');

figure(3); clf;
plot(min_section_size, first_iteration_above_threshold, 'LineWidth', 2);
xlabel('Section size')
ylabel('Iterations until first stalled section')
title([
    'a_n: ', num2str(initial_network_age), ', ',...
    'nodes: ', num2str(network_size), ', ',...
    'fix elders: ', num2str(num_of_elders(1)), ', ',...
    'adversary: ', num2str(fraction_of_new_nodes_are_malicious)
])
filename = [...
    'section_model_first_stalled_iteration',...
    '_fix_elders_', num2str(num_of_elders(1)), ...
    '_age_', num2str(initial_network_age), ...
    '_adversary_', num2str(fraction_of_new_nodes_are_malicious), ...
    '_n_', num2str(network_iterations), ...
    '_size_', num2str(network_size), ...
    '.png'
]
%print(filename, '-dpng');
