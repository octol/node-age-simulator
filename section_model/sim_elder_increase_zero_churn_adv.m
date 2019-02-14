clear all

%network_iterations = 20000;
network_iterations = 8000;
init_iterations = 0;
initial_network_age = 12;

% Test increasing section size
network_size = 100000;
min_section_size = 100
start_section_size = 2*min_section_size
max_section_size = 5*min_section_size
number_of_sections = ceil(network_size ./ start_section_size)
num_of_elders = [10 20 30 40 50 60 70 80 90]*ones(size(number_of_sections))
fraction_of_new_nodes_are_malicious = 0.05
zero_churn_adversary = true;

for ii = 1:length(num_of_elders)
    fprintf('Running with num_of_elders: %d\n', num_of_elders(ii));
    [nodes{ii}, stats{ii}] = run_section_model(
        number_of_sections,
        start_section_size,
        min_section_size,
        max_section_size,
        initial_network_age,
        num_of_elders(ii),
        network_iterations,
        init_iterations,
        fraction_of_new_nodes_are_malicious,
        zero_churn_adversary
    );
end

filename = [
    'sim_elder_increase_zero_churn_adv',...
    '_min_section_', num2str(min_section_size),...
    '_net_age_', num2str(initial_network_age),...
    '_adv_',num2str(fraction_of_new_nodes_are_malicious),...
    '.dat'
]
save(filename)
