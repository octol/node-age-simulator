clear all

network_iterations = 20000;
init_iterations = 0;
initial_network_age = 12;

network_size = 100000;
min_section_size = 100
start_section_size = 2*min_section_size
max_section_size = 5*min_section_size
number_of_sections = network_size ./ start_section_size
num_of_elders = min_section_size;
fraction_of_new_nodes_are_malicious = 0.20
zero_churn_adversary = false;
network_growth_rate = 2e-4;

[nodes,stats] = run_section_model(
    number_of_sections,
    start_section_size,
    min_section_size,
    max_section_size,
    initial_network_age,
    num_of_elders,
    network_iterations,
    init_iterations,
    fraction_of_new_nodes_are_malicious,
    zero_churn_adversary,
    network_growth_rate
);
