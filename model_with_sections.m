clear all
figure(1); clf
figure(2); clf

function nodes = initialise_network(number_of_sections, max_section_size, start_section_size)
    nodes.work = zeros(number_of_sections, max_section_size);
    nodes.age = zeros(number_of_sections, max_section_size);
    nodes.malicious = zeros(number_of_sections, max_section_size);
    nodes.active = logical(zeros(number_of_sections, max_section_size));
    nodes.active(:,1:start_section_size) = ones(number_of_sections, start_section_size);
    nodes.elder = logical(zeros(number_of_sections, max_section_size));
end

function nodes = initialise_nodes(nodes, initial_network_age, num_of_elders)
    % Flat spread of ages
    % Start nodes at age 4 (they're adults, not infants)
    nodes_active_indices = find(nodes.active);
    nodes.work(nodes_active_indices) = round(2.^((initial_network_age-4)*rand(size(nodes_active_indices)) + 4));
    nodes.age(nodes.active) = floor(log2(nodes.work(nodes.active)));

    for s = 1:size(nodes.work,1)
        [sorted_work,I] = sort(nodes.work(s,:),'descend');
        nodes.elder(s,I(1:num_of_elders)) = true;
        if numel(I) > num_of_elders
            nodes.elder(s,I(num_of_elders+1:end)) = false;
        end
    end
end

function [nodes, I] = increase_age(nodes)
    nodes_to_age = rem(log2(nodes.work), 1) == 0;
    nodes.age(nodes_to_age) += 1;
    I = find(nodes_to_age);
end

function nodes = relocate(nodes, nodes_to_age_indices)
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
end

function nodes = churn(nodes)
    nodes_to_drop = and(rand(size(nodes.active)) < 1./nodes.work, nodes.active);
    nodes.work(nodes_to_drop) = 0;
    nodes.age(nodes_to_drop) = 0;
    nodes.malicious(nodes_to_drop) = false;
    nodes.active(nodes_to_drop) = false;
end

function nodes = join_new(nodes, nodes_to_add, min_section_size, add_malicious_nodes, fraction_of_new_nodes_are_malicious)
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
        if add_malicious_nodes
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
    if add_malicious_nodes
        nodes.malicious(node_slots_available(I)) = logical(rand(length(I), 1) < fraction_of_new_nodes_are_malicious);
    else
        nodes.malicious(node_slots_available(I)) = false;
    end
    nodes.active(node_slots_available(I)) = true;
end

function nodes = assign_elder_status(nodes, num_of_elders)
    [sorted_work,I] = sort(nodes.work, 2, 'descend');
    nodes.elder = logical(zeros(size(nodes.active)));
    ind = sub2ind (size(nodes.elder), repmat(1:rows(nodes.elder), num_of_elders, 1), I(:, 1:num_of_elders)');
    nodes.elder(ind) = true;
    assert(all(sum(nodes.elder,2) == num_of_elders))
end

function section_stats = collect_section_statistics(n, nodes, section_stalled_threshold, num_of_elders)
    number_of_sections = size(nodes.active, 1);

    % Section size statistics
    section_stats.size = sum(nodes.active, 2)';
    section_stats.size_mean(n) = mean(section_stats.size);
    section_stats.size_std(n) = std(section_stats.size);
    section_stats.size_malicious = sum(nodes.malicious.*nodes.active, 2)';
    section_stats.size_load = section_stats.size_malicious ./ section_stats.size;
    section_stats.size_load_max(n) = max(section_stats.size_load);
    section_stats.size_load_mean(n) = mean(section_stats.size_load);
    section_stats.size_load_std(n) = std(section_stats.size_load);
    section_stats.stalled_size(n) = sum(section_stats.size_load > section_stalled_threshold) / number_of_sections;

    % Section elders statistics
    section_stats.elders_malicious = sum(nodes.elder.*nodes.malicious.*nodes.active, 2)';
    section_stats.elders_load = section_stats.elders_malicious ./ num_of_elders;
    section_stats.elders_load_max(n) = max(section_stats.elders_load);
    section_stats.elders_load_mean(n) = mean(section_stats.elders_load);
    section_stats.elders_load_std(n) = std(section_stats.elders_load);
    section_stats.stalled_elders(n) = sum(section_stats.elders_load > section_stalled_threshold) / number_of_sections;

    % Section work statistics
    section_stats.work = sum(nodes.work.*nodes.elder.*nodes.active, 2)';
    section_stats.work_malicious = sum(nodes.work.*nodes.elder.*nodes.malicious.*nodes.active, 2)';
    section_stats.work_load = section_stats.work_malicious ./ section_stats.work;
    section_stats.work_load_max(n) = max(section_stats.work_load);
    section_stats.work_load_mean(n) = mean(section_stats.work_load);
    section_stats.work_load_std(n) = std(section_stats.work_load);
    section_stats.stalled_work(n) = sum(section_stats.work_load > section_stalled_threshold) / number_of_sections;

    % Section age statistics
    section_stats.age = sum(nodes.age.*nodes.elder.*nodes.active, 2)';
    section_stats.age_malicious = sum(nodes.age.*nodes.elder.*nodes.malicious.*nodes.active, 2)';
    section_stats.age_load = section_stats.age_malicious ./ section_stats.age;
    section_stats.age_load_max(n) = max(section_stats.age_load);
    section_stats.age_load_mean(n) = mean(section_stats.age_load);
    section_stats.age_load_std(n) = std(section_stats.age_load);
    section_stats.stalled_age(n) = sum(section_stats.age_load > section_stalled_threshold) / number_of_sections;
end

function plot_statistics(n, section_stats, min_section_size, fraction_of_new_nodes_are_malicious, initial_network_age, num_of_elders)
    number_of_sections = size(section_stats.size, 2);
    assert(size(section_stats.size, 1) == 1)

    figure(1)

    subplot(2,2,1)
    hold on
    H1 = plot(1:100:n, section_stats.size_load_mean(1:100:end), 'LineWidth', 2, 'Color', 'k');
    H2 = plot(
        1:100:n,
        [max(section_stats.size_load_mean(1:100:end) - 0.5*section_stats.size_load_std(1:100:end), 0);...
            section_stats.size_load_mean(1:100:end) + 0.5*section_stats.size_load_std(1:100:end)],
        'LineWidth', 2, 'Color', 'b'
    );
    H3 = plot(
        1:100:n,
        [max(section_stats.size_load_mean(1:100:end) - section_stats.size_load_std(1:100:end), 0); ...
            section_stats.size_load_mean(1:100:end) + section_stats.size_load_std(1:100:end)],
        'LineWidth', 2, 'Color', 'm'
    );
    H4 = plot(1:100:n, section_stats.size_load_max(1:100:end), 'LineWidth', 2, 'Color', 'g');
    hold off
    xlabel("Iterations");
    ylabel("Fraction");
    legend([H1, H2(1), H3(1), H4(1)], '\mu', '0.5\sigma', '\sigma', 'max', 'Location', 'Northwest');
    title('Malicious nodes per section')
    drawnow;

    subplot(2,2,2)
    hold on
    H1 = plot(1:100:n, section_stats.elders_load_mean(1:100:end), 'LineWidth', 2, 'Color', 'k');
    H2 = plot(
        1:100:n,
        [max(section_stats.elders_load_mean(1:100:end) - 0.5*section_stats.elders_load_std(1:100:end),0);
            section_stats.elders_load_mean(1:100:end) + 0.5*section_stats.elders_load_std(1:100:end)],
        'LineWidth', 2, 'Color', 'b'
    );
    H3 = plot(
        1:100:n,
        [max(section_stats.elders_load_mean(1:100:end) - section_stats.elders_load_std(1:100:end),0);
            section_stats.elders_load_mean(1:100:end) + section_stats.elders_load_std(1:100:end)],
        'LineWidth', 2, 'Color', 'm'
    );
    H4 = plot(1:100:n, section_stats.elders_load_max(1:100:end), 'LineWidth', 2, 'Color', 'g');
    hold off
    xlabel("Iterations");
    ylabel("Fraction");
    legend([H1, H2(1), H3(1), H4(1)], '\mu', '0.5\sigma', '\sigma', 'max', 'Location', 'Northwest');
    title('Malicious elders per section')
    drawnow;

    subplot(2,2,3)
    hold on
    H1 = plot(1:100:n, section_stats.age_load_mean(1:100:end), 'LineWidth', 2, 'Color', 'k');
    H2 = plot(
        1:100:n,
        [max(section_stats.age_load_mean(1:100:end) - 0.5*section_stats.age_load_std(1:100:end),0);...
            section_stats.age_load_mean(1:100:end) + 0.5*section_stats.age_load_std(1:100:end)],
        'LineWidth', 2, 'Color', 'b'
    );
    H3 = plot(
        1:100:n,
        [max(section_stats.age_load_mean(1:100:end) - section_stats.age_load_std(1:100:end),0);...
            section_stats.age_load_mean(1:100:end) + section_stats.age_load_std(1:100:end)],
        'LineWidth', 2, 'Color', 'm'
    );
    H4 = plot(1:100:n, section_stats.age_load_max(1:100:end), 'LineWidth', 2, 'Color', 'g');
    hold off
    xlabel("Iterations");
    ylabel("Fraction");
    legend([H1, H2(1), H3(1), H4(1)], '\mu', '0.5\sigma', '\sigma', 'max', 'Location', 'Northwest');
    title('Malicious age per section')
    drawnow;

    subplot(2,2,4)
    hold on
    H1 = plot(1:100:n, section_stats.work_load_mean(1:100:end), 'LineWidth', 2, 'Color', 'k');
    H2 = plot(
        1:100:n,
        [max(section_stats.work_load_mean(1:100:end) - 0.5*section_stats.work_load_std(1:100:end),0);...
            section_stats.work_load_mean(1:100:end) + 0.5*section_stats.work_load_std(1:100:end)],
        'LineWidth', 2, 'Color', 'b'
    );
    H3 = plot(
        1:100:n,
        [max(section_stats.work_load_mean(1:100:end) - section_stats.work_load_std(1:100:end),0);...
            section_stats.work_load_mean(1:100:end) + section_stats.work_load_std(1:100:end)],
        'LineWidth', 2, 'Color', 'm'
    );
    H4 = plot(1:100:n, section_stats.work_load_max(1:100:end), 'LineWidth', 2, 'Color', 'g');
    hold off
    xlabel("Iterations");
    ylabel("Fraction");
    legend([H1, H2(1), H3(1), H4(1)], '\mu', '0.5\sigma', '\sigma', 'max', 'Location', 'Northwest');
    title('Malicious work per section')
    drawnow;

    %subplot(2,2,4)
    figure(2)
    plot(
        1:100:n, section_stats.stalled_size(1:100:end), 'LineWidth', 2,
        1:100:n, section_stats.stalled_elders(1:100:end), 'LineWidth', 2,
        1:100:n, section_stats.stalled_age(1:100:end), 'LineWidth', 2,
        1:100:n, section_stats.stalled_work(1:100:end), 'LineWidth', 2
    );
    legend(
        'Stallable sections (nodes)',
        'Stallable sections (elders)',
        'Stallable sections (age)',
        'Stallable sections (work)',
        'Location','Northwest'
    );
    title(['Sections (#/\mu/\sigma/min): ',...
        num2str(number_of_sections),'/',...
        num2str(section_stats.size_mean(n)),'/',...
        num2str(section_stats.size_std(n),2),'/',...
        num2str(min_section_size),...
        ', adversary: ', num2str(fraction_of_new_nodes_are_malicious),...
        ', a_n=', num2str(initial_network_age),...
        ', elders=', num2str(num_of_elders)]
    )
    xlabel("Iterations");
    ylabel("Fraction");
    drawnow
end

network_iterations = 20000;
init_iterations = 0;
initial_network_age = 16;

section_stalled_threshold = 1/3;
fraction_of_existing_malicious_nodes_less_than_fraction_of_new_ones = true;

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

% Initialise network with min_section_size everywhere and only honest nodes
% with zero age
nodes = initialise_network(number_of_sections, max_section_size, start_section_size);
nodes = initialise_nodes(nodes, initial_network_age, num_of_elders);
assert(size(nodes.work, 1) == number_of_sections);

% Evolve network before starting
for n = 1:network_iterations
    tic
    % All nodes does one unit of work
    nodes.work(nodes.active) += 1;

    % Increase age and relocate if work == 2^n
    [nodes, nodes_to_age_indices] = increase_age(nodes);
    nodes = relocate(nodes, nodes_to_age_indices);

    % Randomly drop nodes according to 1/w
    network_size_before_drop = sum(sum(nodes.active));
    nodes = churn(nodes);

    % Join new nodes
    nodes_to_add = network_size_before_drop - sum(sum(nodes.active));
    %fraction_of_network_dropped = nodes_to_add / network_size_before_drop

    add_malicious_nodes = fraction_of_existing_malicious_nodes_less_than_fraction_of_new_ones && n > init_iterations;
    nodes = join_new(nodes, nodes_to_add, min_section_size, add_malicious_nodes, fraction_of_new_nodes_are_malicious);
    assert(floor(log2(nodes.work(nodes.active))) == nodes.age(nodes.active));

    % Assign elder status
    nodes = assign_elder_status(nodes, num_of_elders);

    % If we should stop adding malicious nodes
    fraction_of_existing_malicious_nodes_less_than_fraction_of_new_ones =  ...
        sum(sum(nodes.malicious))/sum(sum(nodes.active)) < fraction_of_new_nodes_are_malicious;

    section_stats = collect_section_statistics(n, nodes, section_stalled_threshold, num_of_elders);

    if mod(n, 100) == 1
        fprintf('section size (mean/std): %d / %d \n', section_stats.size_mean(n), section_stats.size_std(n));

        plot_statistics(n, section_stats, min_section_size, fraction_of_new_nodes_are_malicious, initial_network_age, num_of_elders);
    end
    toc
end

%figure(1)
%print(["section_model_malicious_per_section_age_",num2str(initial_network_age), "_adversary_",num2str(fraction_of_new_nodes_are_malicious), "_section_size_", num2str(min_section_size),"_no_sections_",num2str(number_of_sections),"_elders_",num2str(num_of_elders),".png"],'-dpng');
%figure(2)
%print(["section_model_stallable_sections_age_",num2str(initial_network_age), "_adversary_",num2str(fraction_of_new_nodes_are_malicious), "_section_size_", num2str(min_section_size),"_no_sections_",num2str(number_of_sections),"_elders_",num2str(num_of_elders),".png"],'-dpng');

