clear all

network_iterations = 1000;
init_iterations = 0;
initial_network_age = 16;

% Test increasing section size
number_of_sections = [500 250 100 50];
start_section_size = [20 40 100 200];
max_section_size = [50 100 250 500];
min_section_size = [10 20 50 100];
num_of_elders = min_section_size;
fraction_of_new_nodes_are_malicious = 0.20;

for ii = 1:length(number_of_sections)
    fprintf('Running with number_of_sections: %d\n', number_of_sections(ii));
    figure(1); clf
    figure(2); clf
    nodes{ii} = run_model_with_sections(
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

size(sum(nodes{1}.malicious, 2));
