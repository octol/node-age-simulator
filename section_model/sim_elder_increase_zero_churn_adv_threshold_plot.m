clear all
%load dat/sim_section_size_increase_net_age_12.dat
load sim_elder_increase_zero_churn_adv_min_section_100_net_age_12_adv_0.05.dat

first_iteration_above_threshold = network_iterations*ones(size(num_of_elders));
for ii = 1:length(num_of_elders)
    first_above_threshold = min(find(stats{ii}.malicious_elders_fraction_max > 1/10));
    if length(first_above_threshold) > 0
        first_iteration_above_threshold(ii) = first_above_threshold;
    end
end

figure(1); clf;
plot(num_of_elders, first_iteration_above_threshold./2^12, 'LineWidth', 2);
xlabel('Elder size')
%ylabel('Iterations until first stalled section (elders)')
ylabel('Time until first stalled section / initial network age')
%legend({'0.1 adv','0.2 adv','0.25 adv','0.3 adv'})
title([
    'a_n: ', num2str(initial_network_age), ', ',...
    'nodes: ', num2str(network_size), ', ',...
    'section size: ', num2str(num_of_elders(1))
])
filename = [...
    'section_model_first_stalled_iteration_fix_elders_',...
    '_min_section_', num2str(min_section_size), ...
    '_age_', num2str(initial_network_age), ...
    '_n_', num2str(network_iterations), ...
    '_size_', num2str(network_size), ...
    '.png'
]
%print(filename, '-dpng');
