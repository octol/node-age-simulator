clear all

network_iterations = 30000;
init_iterations = 0;

number_of_sections = 1000;
initial_network_age = 16;

% Scenario A
start_section_size = 16;
max_section_size = 40;
min_section_size = 8;
fraction_of_new_nodes_are_malicious = 0.20;

% Scenario B
%start_section_size = 32;
%max_section_size = 80;
%min_section_size = 16;
%fraction_of_new_nodes_are_malicious = 0.10;

% Scenario C
%start_section_size = 16;
%max_section_size = 40;
%min_section_size = 8;
%fraction_of_new_nodes_are_malicious = 0.20;

% Scenario D
%start_section_size = 32;
%max_section_size = 80;
%min_section_size = 16;
%fraction_of_new_nodes_are_malicious = 0.20;

% Scenario E
%init_iterations = 0;
%start_section_size = 128;
%max_section_size = 320;
%min_section_size = 64;
%fraction_of_new_nodes_are_malicious = 0.25;

section_stalled_threshold = 1/3;
section_compromised_threshold = 2/3;

% Initialise network with min_section_size everywhere and only honest nodes
% with zero age
nodes.work = zeros(number_of_sections, max_section_size);
nodes.age = zeros(number_of_sections, max_section_size);
nodes.malicious = zeros(number_of_sections, max_section_size);
nodes.active = logical(zeros(number_of_sections, max_section_size));
nodes.active(:,1:start_section_size) = ones(number_of_sections, start_section_size);

nodes_active_indices = find(nodes.active);
num_of_age_buckets = 16;
age_bucket_size = numel(nodes_active_indices) / num_of_age_buckets;

% Flat spread of ages
% Start nodes at age 4 (they're adults, not infants)
nodes.work(nodes_active_indices) = round(2.^((initial_network_age-4)*rand(size(nodes_active_indices)) + 4));
nodes.age(nodes.active) = floor(log2(nodes.work(nodes.active)));

% Evolve network before starting
for n = 1:network_iterations
    % All nodes does one unit of work
    nodes.work(nodes.active) += 1;

    % Increase age and relocate if work == 2^n
    nodes_to_age = rem(log2(nodes.work), 1) == 0;
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
    nodes.malicious(node_slots_available(I)) = nodes.malicious(nodes_to_age_indices);
    nodes.active(node_slots_available(I)) = nodes.active(nodes_to_age_indices);
    nodes.work(nodes_to_age_indices) = 0;
    nodes.age(nodes_to_age_indices) = 0;
    nodes.malicious(nodes_to_age_indices) = false;
    nodes.active(nodes_to_age_indices) = false;

    % Randomly drop nodes according to 1/w
    network_size_before_drop = sum(sum(nodes.active));
    nodes_to_drop = and(rand(number_of_sections, max_section_size) < 1./nodes.work, nodes.active);
    nodes.work(nodes_to_drop) = 0;
    nodes.age(nodes_to_drop) = 0;
    nodes.malicious(nodes_to_drop) = false;
    nodes.active(nodes_to_drop) = false;

    % Join new nodes
    nodes_to_add = network_size_before_drop - sum(sum(nodes.active));
    %fraction_of_network_dropped = nodes_to_add / network_size_before_drop

    % Node joins should flow towards the smallest sections, but to be on the
    % safe side add nodes to the smallest sections first.
    % Quite an ugly workaround, but this all of course comes from performance
    % considerations.
    % Duplication galore!
    small_sections = sum(nodes.active, 2) < min_section_size;
    number_of_small_sections = sum(small_sections);
    while sum(small_sections) > 0
        i = find(small_sections, 1);
        j = find(nodes.active(i,:) == false, 1);
        nodes.work(i,j) = 2^4;
        nodes.age(i,j) = round(log2(nodes.work(i,j)));
        if n > init_iterations
            nodes.malicious(i,j) = logical(rand() < fraction_of_new_nodes_are_malicious);
        else
            nodes.malicious(i,j) = false;
        end
        nodes.active(i,j) = true;
        small_sections = sum(nodes.active, 2) < min_section_size;
        nodes_to_add -= 1;
    end

    % Add the rest
    node_slots_available = find(nodes.active == false);
    I = randperm(length(node_slots_available));
    assert(length(I) > nodes_to_add);
    I = I(1:nodes_to_add);
    nodes.work(node_slots_available(I)) = 2^4;
    nodes.age(node_slots_available(I)) = round(log2(2^4));
    if n > init_iterations
        nodes.malicious(node_slots_available(I)) = logical(rand(length(I), 1) < fraction_of_new_nodes_are_malicious);
    else
        nodes.malicious(node_slots_available(I)) = false;
    end
    nodes.active(node_slots_available(I)) = true;
    assert(floor(log2(nodes.work(nodes.active))) == nodes.age(nodes.active));

    % Collect statistics
    section_size = sum(nodes.active, 2)';
    section_size_mean(n) = mean(section_size);
    section_size_std(n) = std(section_size);
    network_size(n) = sum(section_size);

    section_size_malicious = sum(nodes.malicious, 2)';
    section_size_malicious_mean(n) = mean(section_size_malicious);
    section_size_malicious_std(n) = std(section_size_malicious);
    netsize_size_malicious(n) = sum(section_size_malicious);

    section_work = sum(nodes.work, 2)';
    section_work_mean(n) = mean(section_work);
    section_work_std(n) = std(section_work);
    network_work(n) = sum(section_work);

    section_work_malicious = sum(nodes.work.*nodes.malicious, 2)';
    section_work_malicious_mean(n) = mean(section_work_malicious);
    section_work_malicious_std(n) = std(section_work_malicious);
    network_work_malicious(n) = sum(section_work_malicious);
    network_work_malicious_fraction(n) = network_work_malicious(n) / network_work(n);

    section_load = section_size_malicious ./ section_size;
    section_load_mean(n) = mean(section_load);
    section_load_std(n) = std(section_load);

    section_work_load = section_work_malicious ./ section_work;
    section_work_load_mean(n) = mean(section_work_load);
    section_work_load_std(n) = std(section_work_load);

    stalled_sections(n) = sum(section_load > section_stalled_threshold) / number_of_sections;
    compromised_sections(n) = sum(section_load > section_compromised_threshold) / number_of_sections;

    stalled_sections_work(n) = sum(section_work_load > section_stalled_threshold) / number_of_sections;
    compromised_sections_work(n) = sum(section_work_load > section_compromised_threshold) / number_of_sections;

    if mod(n, 100) == 0
        fprintf('section size (mean/std): %d / %d \n', section_size_mean(n), section_size_std(n));
        fprintf('  section work (mean/std): %d / %d \n', section_work_mean(n), section_work_std(n));

        figure(1)
        subplot(2,2,1)
        hist(nodes.age(nodes.active), 50);
        xlabel('Age');
        title(['Sections (number/mean/std/min): ',...
            num2str(number_of_sections),'/',...
            num2str(section_size_mean(n)),'/',...
            num2str(section_size_std(n),2),'/',...
            num2str(min_section_size)]
        )
        drawnow;

        subplot(2,2,2)
        hold on
        H1 = plot(1:n, section_load_mean, 'LineWidth', 2, 'Color', 'k');
        H2 = plot(
            1:n, [max(section_load_mean - 0.5*section_load_std, 0); section_load_mean + 0.5*section_load_std],
            'LineWidth', 2, 'Color', 'b'
        );
        H3 = plot(
            1:n, [max(section_load_mean - section_load_std, 0); section_load_mean + section_load_std],
            'LineWidth', 2, 'Color', 'm'
        );
        hold off
        legend([H1, H2(1), H3(1)], '\mu', '0.5\sigma', '\sigma', 'Location', 'Northwest');
        title('Malicious nodes per section')
        drawnow;

        subplot(2,2,3)
        hold on
        H1 = plot(1:n, section_work_load_mean, 'LineWidth', 2, 'Color', 'k');
        H2 = plot(
            1:n, [max(section_work_load_mean - 0.5*section_work_load_std,0); section_work_load_mean + 0.5*section_work_load_std],
            'LineWidth', 2, 'Color', 'b'
        );
        H3 = plot(
            1:n, [max(section_work_load_mean - section_work_load_std,0); section_work_load_mean + section_work_load_std],
            'LineWidth', 2, 'Color', 'm'
        );
        hold off
        legend([H1, H2(1), H3(1)], '\mu', '0.5\sigma', '\sigma', 'Location', 'Northwest');
        title('Malicious work per section')
        drawnow;

        subplot(2,2,4)
        plot(1:n, stalled_sections, 'LineWidth', 2, 1:n, stalled_sections_work, 'LineWidth', 2);
        legend('Stallable section (nodes)', 'Stallable sections (work)','Location','Northwest');
        title(['Adversary: ', num2str(fraction_of_new_nodes_are_malicious)]);
        drawnow
    end
end

%figure(1)
%print(["section_model_age_",num2str(initial_network_age), "_adversary_",num2str(fraction_of_new_nodes_are_malicious), "_section_size_", num2str(min_section_size),".png"],'-dpng');
