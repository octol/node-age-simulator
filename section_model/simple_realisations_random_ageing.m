clear all

network_iterations = 20000;
init_iterations = 0;
initial_network_age = 12;

network_size = 100000;
min_section_size = 100;
start_section_size = 2*min_section_size
max_section_size = 5*min_section_size
number_of_sections = ceil(network_size ./ start_section_size)
num_of_elders = 10*ones(size(number_of_sections))
fraction_of_new_nodes_are_malicious = 0.20;

[nodes,stats] = run_section_model(
    number_of_sections,
    start_section_size,
    min_section_size,
    max_section_size,
    initial_network_age,
    num_of_elders,
    network_iterations,
    init_iterations,
    fraction_of_new_nodes_are_malicious
);

figure(1)
filename = [...
    'section_model_random_ageing',...
    '_age_', num2str(initial_network_age), ...
    '_adversary_', num2str(fraction_of_new_nodes_are_malicious), ...
    '_section_size_', num2str(min_section_size(1)),...
    '_no_sections_', num2str(number_of_sections(1)),...
    '_elders_', num2str(num_of_elders(1)), ...
    '.png']
%print(filename, '-dpng');
