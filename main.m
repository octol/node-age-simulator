clear all

function s = get_sizes(section)
    for ii = 1:length(section)
        s(ii) = length(section(ii).node);
    end
end

function index = random_index(indices)
    index = indices(randi([1, length(indices)]));
end

network_iterations = 1000;
number_of_sections = 10000;
min_section_size = 8;
start_section_size = 8;
max_section_size = 20;
total_start_nodes = start_section_size * number_of_sections;
total_nodes = total_start_nodes;
total_bad_nodes = 0;

% Initialise network with min_section_size everywhere and only honest nodes
% with zero age
nodes.work = zeros(number_of_sections, max_section_size);
nodes.age = zeros(number_of_sections, max_section_size);
nodes.durability = randi([1,3], number_of_sections, max_section_size);
nodes.evil = zeros(number_of_sections, max_section_size);
nodes.active = logical(zeros(number_of_sections, max_section_size));
nodes.active(:,1:start_section_size) = ones(number_of_sections, start_section_size);

% Evolve network before starting
for n = 1:network_iterations
    % Part of the network does work
    nodes_that_does_work = logical(randi([0,1], number_of_sections, max_section_size));
    nodes.work(and(nodes_that_does_work, nodes.active)) += 1;

    % Increase age and relocate if work == 2^n
    nodes_to_age = rem(log2(nodes.work), 1) == 0;
    nodes.age(nodes_to_age) += 1;
    nodes_to_age_indices = find(nodes_to_age);

    % Find target slots
    % Note: this means we relocate towards emptier section
    node_slots_available = find(nodes.active == false);
    assert(length(nodes_to_age_indices) < length(node_slots_available));
    I = randperm(length(node_slots_available));
    I = I(1:length(nodes_to_age_indices));

    % Relocate
    nodes.work(node_slots_available(I)) = nodes.work(nodes_to_age_indices);
    nodes.age(node_slots_available(I)) = nodes.age(nodes_to_age_indices);
    nodes.durability(node_slots_available(I)) = nodes.durability(nodes_to_age_indices);
    nodes.evil(node_slots_available(I)) = nodes.evil(nodes_to_age_indices);
    nodes.active(node_slots_available(I)) = nodes.active(nodes_to_age_indices);
    nodes.work(nodes_to_age_indices) = 0;
    nodes.age(nodes_to_age_indices) = 0;
    nodes.durability(nodes_to_age_indices) = 0;
    nodes.evil(nodes_to_age_indices) = 0;
    nodes.active(nodes_to_age_indices) = 0;

    % Join new nodes
    node_slots_available = find(nodes.active == false);
    I = randperm(length(node_slots_available));
    I = I(1:end/8);
    nodes.work(node_slots_available(I)) = 0;
    nodes.age(node_slots_available(I)) = 0;
    nodes.durability(node_slots_available(I)) = randi([1,3], length(I), 1);
    if n > 100
        nodes.evil(node_slots_available(I)) = logical(rand(length(I), 1) < 0.1);
    end
    nodes.active(node_slots_available(I)) = true;

    % Randomly drop nodes depending on durability
    nodes_to_drop = ((rand(number_of_sections, max_section_size).*nodes.active) ./ nodes.durability) > 0.3;
    nodes.work(nodes_to_drop) = 0;
    nodes.age(nodes_to_drop) = 0;
    nodes.durability(nodes_to_drop) = 0;
    nodes.evil(nodes_to_drop) = 0;
    nodes.active(nodes_to_drop) = false;

    section_spread = sum(nodes.active, 2)';
    section_load = sum(nodes.evil, 2)';

    fprintf('total_nodes: %d, total_bad_nodes: %d, max(section_load): %d\n', sum(section_spread), sum(section_load), max(section_load));
end
