clear all

function s = get_sizes(section)
    for ii = 1:length(section)
        s(ii) = length(section(ii).node);
    end
end

function index = random_index(indices)
    index = indices(randi([1, length(indices)]));
end

number_of_sections = 100;
min_section_size = 8;
start_section_size = 12;
total_start_nodes = start_section_size * number_of_sections;
total_nodes = total_start_nodes;
total_bad_nodes = 0;

% Initialise network with min_section_size everywhere and only honest nodes
% with zero age
for section_ii = 1:number_of_sections
    for node_ii = 1:start_section_size
        section(section_ii).node(node_ii).work = 0;
        section(section_ii).node(node_ii).age = 0;
        section(section_ii).node(node_ii).durability = randi([1, 3]);
        section(section_ii).node(node_ii).evil = 0;
    end
end

% Evolve network before starting
N = 1000;
for n = 1:N
    % Spread out total_nodes/3 amount of work
    % Generate total_nodes number of index tuples (section_ii, node_ii)
    total_work = total_nodes/3;
    for ii = 1:total_work
        section_ii = randi([1, number_of_sections]);
        node_ii = randi([1, length(section(section_ii).node)]);

        section(section_ii).node(node_ii).work += 1;

        % Increase age and relocate if work == 2^N
        log2_value = log2(section(section_ii).node(node_ii).work);
        if abs(rem(log2_value, 1)) < 1e3*eps
            section(section_ii).node(node_ii).age += 1;

            % Target index should be a among the smallest sections
            section_sizes = get_sizes(section);

            % Alt 1: random among section below min_section_size
            small_section_indices = find(section_sizes < min_section_size);
            % Alt 2: random among the smallest sections
            %minimal_section_size = min(section_sizes);
            %small_section_indices = find(section_sizes == minimal_section_size);

            if length(small_section_indices) > 0
                target_ii = random_index(small_section_indices);
            else
                target_ii = randi([1, number_of_sections]);
            end
            if target_ii != section_ii
                section(target_ii).node = [section(target_ii).node section(section_ii).node(node_ii)];
                section(section_ii).node(node_ii) = [];
            end
        end
    end

    % Randomly join new nodes
    if total_nodes < total_start_nodes
        new_nodes = total_start_nodes - total_nodes;
        for ii = 1:new_nodes
            section_ii = randi([1, number_of_sections]);
            node_ii = length(section(section_ii).node) + 1;

            evil = (n > 20) * (rand() < 0.1);
            total_bad_nodes += evil;

            section(section_ii).node(node_ii).work = 0;
            section(section_ii).node(node_ii).age = 0;
            section(section_ii).node(node_ii).durability = randi([1, 3]);
            section(section_ii).node(node_ii).evil = evil;
        end
        total_nodes += new_nodes;
    end

    % Randomly drop nodes depending on durability
    for ii = 1:total_nodes/10
        section_ii = randi([1, number_of_sections]);
        node_ii = randi([1, length(section(section_ii).node)]);

        chance_to_drop = rand()/section(section_ii).node(node_ii).durability;
        if chance_to_drop > 0.5
            total_bad_nodes -= section(section_ii).node(node_ii).evil;
            total_nodes -= 1;
            section(section_ii).node(node_ii) = [];
            assert(length(section(section_ii).node) > 0, 'section size dropped to zero!')
        end
    end

    % Max evil section load
    for ii = 1:length(section)
        section_size = length(section(ii).node);
        s = sum([section(ii).node.evil]);
        section_load(ii) = s/section_size;
    end

    fprintf('total_nodes: %d, total_bad_nodes: %d, max(section_load): %d\n', total_nodes, total_bad_nodes, max(section_load));
end
