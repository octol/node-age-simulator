clear all
figure(1); clf
figure(2); clf

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

% Scenario B (double section size)
%number_of_sections = 2500;
%start_section_size = 40;
%max_section_size = 100;
%min_section_size = 20;
%num_of_elders = min_section_size;
%fraction_of_new_nodes_are_malicious = 0.20;

% Scenario C (double network size)
%number_of_sections = 10000;
%start_section_size = 20;
%max_section_size = 50;
%min_section_size = 10;
%num_of_elders = min_section_size;
%fraction_of_new_nodes_are_malicious = 0.20;

% Scenario D (further increase section size)
%number_of_sections = 1000;
%start_section_size = 100;
%max_section_size = 250;
%min_section_size = 50;
%num_of_elders = min_section_size;
%fraction_of_new_nodes_are_malicious = 0.20;

% Scenario E (further increase section size)
%number_of_sections = 500;
%start_section_size = 200;
%max_section_size = 500;
%min_section_size = 100;
%num_of_elders = min_section_size;
%fraction_of_new_nodes_are_malicious = 0.20;

% Scenario F (0.5 attacker)
%number_of_sections = 500;
%start_section_size = 200;
%max_section_size = 500;
%min_section_size = 100;
%num_of_elders = min_section_size;
%fraction_of_new_nodes_are_malicious = 0.50;

% Scenario V
%number_of_sections = 5000;
%start_section_size = 200;
%max_section_size = 500;
%min_section_size = 100;
%num_of_elders = min_section_size;
%fraction_of_new_nodes_are_malicious = 0.20;

nodes = run_model_with_sections(
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
