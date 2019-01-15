clear all

network_iterations = 200;
init_iterations = 100;

number_of_sections = 10000;
start_section_size = 10;
max_section_size = 30;
part_of_empty_section_slots_filled_per_iteration = 8; % means: 1/8
fraction_of_new_nodes_are_malicious = 0.10;
max_durability = 9;

section_stalled_threshold = 1/3;
section_compromised_threshold = 2/3;

% Initialise network with min_section_size everywhere and only honest nodes
% with zero age
nodes.work = zeros(number_of_sections, max_section_size);
nodes.age = zeros(number_of_sections, max_section_size);
nodes.durability = randi([1,max_durability], number_of_sections, max_section_size);
nodes.malicious = zeros(number_of_sections, max_section_size);
nodes.active = logical(zeros(number_of_sections, max_section_size));
nodes.active(:,1:start_section_size) = ones(number_of_sections, start_section_size);

figure(1); clf
figure(2); clf

% Evolve network before starting
for n = 1:network_iterations
    % Part of the network does work
    nodes_that_does_work = logical(randi([0,1], number_of_sections, max_section_size));
    nodes_that_does_work = and(nodes_that_does_work, nodes.active);
    nodes.work(nodes_that_does_work) += 1;

    % Increase age and relocate if work == 2^n
    nodes_to_age = and(rem(log2(nodes.work), 1) == 0, nodes_that_does_work);
    nodes.age(nodes_to_age) += 1;
    nodes_to_age_indices = find(nodes_to_age);

    % Find target slots
    % Note: this means we relocate towards emptier sections
    node_slots_available = find(nodes.active == false);
    assert(length(nodes_to_age_indices) < length(node_slots_available));
    I = randperm(length(node_slots_available));
    I = I(1:length(nodes_to_age_indices));

    % Relocate
    nodes.work(node_slots_available(I)) = nodes.work(nodes_to_age_indices);
    nodes.age(node_slots_available(I)) = nodes.age(nodes_to_age_indices);
    nodes.durability(node_slots_available(I)) = nodes.durability(nodes_to_age_indices);
    nodes.malicious(node_slots_available(I)) = nodes.malicious(nodes_to_age_indices);
    nodes.active(node_slots_available(I)) = nodes.active(nodes_to_age_indices);
    nodes.work(nodes_to_age_indices) = 0;
    nodes.age(nodes_to_age_indices) = 0;
    nodes.durability(nodes_to_age_indices) = 0;
    nodes.malicious(nodes_to_age_indices) = 0;
    nodes.active(nodes_to_age_indices) = 0;

    % Join new nodes
    node_slots_available = find(nodes.active == false);
    I = randperm(length(node_slots_available));
    I = I(1:end/part_of_empty_section_slots_filled_per_iteration);
    nodes.work(node_slots_available(I)) = 0;
    nodes.age(node_slots_available(I)) = 0;
    nodes.durability(node_slots_available(I)) = randi([1,max_durability], length(I), 1);
    if n > init_iterations
        nodes.malicious(node_slots_available(I)) = logical(rand(length(I), 1) < fraction_of_new_nodes_are_malicious);
    end
    nodes.active(node_slots_available(I)) = true;

    % Randomly drop nodes depending on durability
    nodes_to_drop = (rand(number_of_sections, max_section_size).*nodes.active) < 1./(nodes.durability + 1);
    nodes.work(nodes_to_drop) = 0;
    nodes.age(nodes_to_drop) = 0;
    nodes.durability(nodes_to_drop) = 0;
    nodes.malicious(nodes_to_drop) = 0;
    nodes.active(nodes_to_drop) = false;

    % Only start plotting once we've established a steady-state (honest) network
    if n > init_iterations
        m = n - init_iterations;
        section_size = sum(nodes.active, 2)';
        section_size_mean(m) = mean(section_size);
        section_size_std(m) = std(section_size);
        network_size(m) = sum(section_size);

        section_size_malicious = sum(nodes.malicious, 2)';
        section_size_malicious_mean(m) = mean(section_size_malicious);
        section_size_malicious_std(m) = std(section_size_malicious);
        netsize_size_malicious(m) = sum(section_size_malicious);

        section_work = sum(nodes.work, 2)';
        section_work_mean(m) = mean(section_work);
        section_work_std(m) = std(section_work);
        network_work(m) = sum(section_work);

        section_work_malicious = sum(nodes.work.*nodes.malicious, 2)';
        section_work_malicious_mean(m) = mean(section_work_malicious);
        section_work_malicious_std(m) = std(section_work_malicious);
        network_work_malicious(m) = sum(section_work_malicious);
        network_work_malicious_fraction(m) = network_work_malicious(m) / network_work(m);

        section_load = section_size_malicious ./ section_size;
        section_load_mean(m) = mean(section_load);
        section_load_std(m) = std(section_load);

        stalled_sections(m) = sum(section_load > section_stalled_threshold) / number_of_sections;
        compromised_sections(m) = sum(section_load > section_compromised_threshold) / number_of_sections;

        fprintf('section size (mean/std): %d / %d \n', section_size_mean(m), section_size_std(m));
        fprintf('  section work (mean/std): %d / %d \n', section_work_mean(m), section_work_std(m));
        figure(1)
        hold on
        plot(network_work_malicious_fraction, 'b-', 'linewidth',2);
        title('\Sigma(w_{malicious}) / \Sigma(w)');
        xlabel('Iterations')
        drawnow
        hold off

        figure(2)
        hold on
        plot(stalled_sections, 'b-', 'linewidth', 2);
        plot(compromised_sections, 'b-', 'linewidth', 2);
        title('Fraction of stalled and compromised sections');
        xlabel('Iterations')
        drawnow
        hold off
    end
end

% figure(1)
% print -dpng fraction_malicious_work.png
% figure(2)
% print -dpng fraction_stalled_sections.png
