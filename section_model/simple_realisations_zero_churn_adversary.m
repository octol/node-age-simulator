clear all

network_iterations = 20000;
init_iterations = 0;
initial_network_age = 16;

% Scenario A
number_of_sections = 5000;
start_section_size = 20;
max_section_size = 50;
min_section_size = 10;
num_of_elders = min_section_size;
fraction_of_new_nodes_are_malicious = 0.20;
zero_churn_adversary = true;

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
    zero_churn_adversary
);
