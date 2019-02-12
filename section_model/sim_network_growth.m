clear all

network_iterations = 20000;
init_iterations = 0;
initial_network_age = 12;

% Test increasing section size
network_size = 100000;
min_section_size = [100]
%min_section_size = [10 20 50 100]
%min_section_size = [8 10 13 17 20 30 40 50 60 70 80 90 100 120 140 160 180 200]
start_section_size = 2*min_section_size
max_section_size = 5*min_section_size
number_of_sections = network_size ./ start_section_size
num_of_elders = min_section_size;
fraction_of_new_nodes_are_malicious = 0.20
zero_churn_adversary = false;
network_growth_rate = [0 1e-5 1e-4];

for ii = 1:length(network_growth_rate)
    fprintf('Running with network_growth_rate: %d\n', network_growth_rate(ii));
    [nodes{ii}, stats{ii}] = run_section_model(
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
        network_growth_rate(ii)
    );
end

filename = [
    'sim_network_growth',...
    '_net_age_',num2str(initial_network_age),...
    '_adv_',num2str(fraction_of_new_nodes_are_malicious),...
    '_section_',num2str(min_section_size),...
    '_elders_',num2str(num_of_elders),...
    '.dat']
save(filename)
