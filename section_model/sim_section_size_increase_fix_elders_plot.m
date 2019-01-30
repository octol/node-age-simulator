clear all
load dat/sim_section_size_increase_fix_elders.dat

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
    'section_model_malicious_nodes_vs_section_size_fix_elders',...
    '_age_', num2str(initial_network_age), ...
    '_adversary_', num2str(fraction_of_new_nodes_are_malicious), ...
    '_n_', num2str(network_iterations), ...
    '_size_', num2str(network_size), ...
    '.png'];
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
    'section_model_malicious_elders_vs_section_size_fix_elders',...
    '_age_', num2str(initial_network_age), ...
    '_adversary_', num2str(fraction_of_new_nodes_are_malicious), ...
    '_n_', num2str(network_iterations), ...
    '_size_', num2str(network_size), ...
    '.png'];
print(filename, '-dpng');

%figure(2); clf;
%hist(malicious_per_section{end},100);
