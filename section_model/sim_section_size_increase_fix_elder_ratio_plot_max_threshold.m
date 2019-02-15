clear all
%load dat/sim_section_size_increase_net_age_12.dat
load sim_section_size_increase_fix_elder_ratio_10_net_age_12_adv_0.1.dat

network_size = number_of_sections(1) * start_section_size(1);
%assert(all((number_of_sections.*start_section_size - network_size) == 0))
%max(abs(number_of_sections.*start_section_size - network_size))

% Find threshold iterations
jj = 1;
first_iteration_above_threshold(jj,:) = network_iterations*ones(size(number_of_sections));
for ii = 1:length(number_of_sections)
    first_above_threshold = min(find(stats{ii}.malicious_elders_fraction_max > 1/3));
    if length(first_above_threshold) > 0
        first_iteration_above_threshold(jj, ii) = first_above_threshold;
    end
end

load sim_section_size_increase_fix_elder_ratio_10_net_age_12_adv_0.2.dat
jj = 2;
first_iteration_above_threshold(jj,:) = network_iterations*ones(size(number_of_sections));
for ii = 1:length(number_of_sections)
    first_above_threshold = min(find(stats{ii}.malicious_elders_fraction_max > 1/3));
    if length(first_above_threshold) > 0
        first_iteration_above_threshold(jj, ii) = first_above_threshold;
    end
end

load sim_section_size_increase_fix_elder_ratio_10_net_age_12_adv_0.25.dat
jj = 3;
first_iteration_above_threshold(jj,:) = network_iterations*ones(size(number_of_sections));
for ii = 1:length(number_of_sections)
    first_above_threshold = min(find(stats{ii}.malicious_elders_fraction_max > 1/3));
    if length(first_above_threshold) > 0
        first_iteration_above_threshold(jj, ii) = first_above_threshold;
    end
end

load sim_section_size_increase_fix_elder_ratio_10_net_age_12_adv_0.3.dat
jj = 4;
first_iteration_above_threshold(jj,:) = network_iterations*ones(size(number_of_sections));
for ii = 1:length(number_of_sections)
    first_above_threshold = min(find(stats{ii}.malicious_elders_fraction_max > 1/3));
    if length(first_above_threshold) > 0
        first_iteration_above_threshold(jj, ii) = first_above_threshold;
    end
end



figure(1); clf;
plot(min_section_size, first_iteration_above_threshold./2^12, 'LineWidth', 2);
xlabel('Section size')
%ylabel('Iterations until first stalled section (elders)')
ylabel('Time until first stalled section / initial network age')
legend({'0.1 adv','0.2 adv','0.25 adv','0.3 adv'})
title([
    'a_n: ', num2str(initial_network_age), ', ',...
    'nodes: ', num2str(network_size), ', ',...
    'fix elder ratio: ', num2str(num_of_elders(1)/min_section_size(1))
])
filename = [...
    'section_model_first_stalled_iteration_combined',...
    '_fix_elder_ratio_', num2str(num_of_elders(1)/min_section_size(1)), ...
    '_age_', num2str(initial_network_age), ...
    '_n_', num2str(network_iterations), ...
    '_size_', num2str(network_size), ...
    '.png'
]
print(filename, '-dpng');
