clear all

network_iterations = 20000;
init_iterations = 0;
initial_network_age = 12;

% Test increasing section size
network_size = 100000;
%min_section_size = [10 20 50 100]
min_section_size = [100 200 300 400 500]
start_section_size = 2*min_section_size
max_section_size = 5*min_section_size
number_of_sections = ceil(network_size ./ start_section_size)
%num_of_elders = 10*ones(size(number_of_sections))
num_of_elders = min_section_size ./ 10
fraction_of_new_nodes_are_malicious = 0.3;

for ii = 1:length(number_of_sections)
    fprintf('Running with number_of_sections: %d\n', number_of_sections(ii));
    [nodes{ii}, stats{ii}] = run_section_model(
        number_of_sections(ii),
        start_section_size(ii),
        min_section_size(ii),
        max_section_size(ii),
        initial_network_age,
        num_of_elders(ii),
        network_iterations,
        init_iterations,
        fraction_of_new_nodes_are_malicious
    );
end

filename = [
    'sim_section_size_increase_fix_elder_ratio_', num2str(num_of_elders(1)),...
    '_net_age_',num2str(initial_network_age),...
    '_adv_',num2str(fraction_of_new_nodes_are_malicious),...
    '.dat'
]
save(filename)
