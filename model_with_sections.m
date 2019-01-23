clear all
figure(1); clf
figure(2); clf

network_iterations = 20000;
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
%fraction_of_new_nodes_are_malicious = 0.20;

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

% Initialise network with min_section_size everywhere and only honest nodes
% with zero age
nodes.work = zeros(number_of_sections, max_section_size);
nodes.age = zeros(number_of_sections, max_section_size);
nodes.malicious = zeros(number_of_sections, max_section_size);
nodes.active = logical(zeros(number_of_sections, max_section_size));
nodes.active(:,1:start_section_size) = ones(number_of_sections, start_section_size);
nodes.elder = logical(zeros(number_of_sections, max_section_size));

nodes_active_indices = find(nodes.active);
num_of_age_buckets = 16;
age_bucket_size = numel(nodes_active_indices) / num_of_age_buckets;

% Flat spread of ages
% Start nodes at age 4 (they're adults, not infants)
nodes.work(nodes_active_indices) = round(2.^((initial_network_age-4)*rand(size(nodes_active_indices)) + 4));
nodes.age(nodes.active) = floor(log2(nodes.work(nodes.active)));

assert(size(nodes.work,1) == number_of_sections);
for s = 1:size(nodes.work,1)
    [sorted_work,I] = sort(nodes.work(s,:),'descend');
    nodes.elder(s,I(1:min_section_size)) = true;
    if numel(I) > min_section_size
        nodes.elder(s,I(min_section_size+1:end)) = false;
    end
end

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

    % Assign elder status
    [sorted_work,I] = sort(nodes.work, 2, 'descend');
    nodes.elder = logical(zeros(number_of_sections, max_section_size));
    ind = sub2ind (size(nodes.elder), repmat(1:rows(nodes.elder),min_section_size,1), I(:,1:min_section_size)');
    nodes.elder(ind) = true;
    assert(all(sum(nodes.elder,2) == min_section_size))

    % Section size statistics
    section_size = sum(nodes.active, 2)';
    section_size_mean(n) = mean(section_size);
    section_size_std(n) = std(section_size);
    section_size_malicious = sum(nodes.malicious.*nodes.active, 2)';
    section_size_load = section_size_malicious ./ section_size;
    section_size_load_mean(n) = mean(section_size_load);
    section_size_load_std(n) = std(section_size_load);
    stalled_sections_size(n) = sum(section_size_load > section_stalled_threshold) / number_of_sections;

    % Section elders statistics
    section_elders_malicious = sum(nodes.elder.*nodes.malicious.*nodes.active, 2)';
    section_elders_load = section_elders_malicious ./ min_section_size;
    section_elders_load_mean(n) = mean(section_elders_load);
    section_elders_load_std(n) = std(section_elders_load);
    stalled_sections_elders(n) = sum(section_elders_load > section_stalled_threshold) / number_of_sections;

    % Section work statistics
    section_work = sum(nodes.work.*nodes.active, 2)';
    section_work_malicious = sum(nodes.work.*nodes.malicious.*nodes.active, 2)';
    section_work_load = section_work_malicious ./ section_work;
    section_work_load_mean(n) = mean(section_work_load);
    section_work_load_std(n) = std(section_work_load);
    stalled_sections_work(n) = sum(section_work_load > section_stalled_threshold) / number_of_sections;

    % Section age statistics
    section_age = sum(nodes.age.*nodes.active, 2)';
    section_age_malicious = sum(nodes.age.*nodes.malicious.*nodes.active, 2)';
    section_age_load = section_age_malicious ./ section_age;
    section_age_load_mean(n) = mean(section_age_load);
    section_age_load_std(n) = std(section_age_load);
    stalled_sections_age(n) = sum(section_age_load > section_stalled_threshold) / number_of_sections;

    if mod(n, 100) == 0
        fprintf('section size (mean/std): %d / %d \n', section_size_mean(n), section_size_std(n));

        figure(1)

        subplot(2,2,1)
        hold on
        H1 = plot(1:100:n, section_size_load_mean(1:100:end), 'LineWidth', 2, 'Color', 'k');
        H2 = plot(
            1:100:n,
            [max(section_size_load_mean(1:100:end) - 0.5*section_size_load_std(1:100:end), 0);...
             section_size_load_mean(1:100:end) + 0.5*section_size_load_std(1:100:end)],
            'LineWidth', 2, 'Color', 'b'
        );
        H3 = plot(
            1:100:n,
            [max(section_size_load_mean(1:100:end) - section_size_load_std(1:100:end), 0); ...
             section_size_load_mean(1:100:end) + section_size_load_std(1:100:end)],
            'LineWidth', 2, 'Color', 'm'
        );
        hold off
        xlabel("Iterations");
        ylabel("Fraction");
        legend([H1, H2(1), H3(1)], '\mu', '0.5\sigma', '\sigma', 'Location', 'Northwest');
        title('Malicious nodes per section')
        drawnow;

        subplot(2,2,2)
        hold on
        H1 = plot(1:100:n, section_elders_load_mean(1:100:end), 'LineWidth', 2, 'Color', 'k');
        H2 = plot(
            1:100:n,
            [max(section_elders_load_mean(1:100:end) - 0.5*section_elders_load_std(1:100:end),0);
             section_elders_load_mean(1:100:end) + 0.5*section_elders_load_std(1:100:end)],
            'LineWidth', 2, 'Color', 'b'
        );
        H3 = plot(
            1:100:n,
            [max(section_elders_load_mean(1:100:end) - section_elders_load_std(1:100:end),0);
             section_elders_load_mean(1:100:end) + section_elders_load_std(1:100:end)],
            'LineWidth', 2, 'Color', 'm'
        );
        hold off
        xlabel("Iterations");
        ylabel("Fraction");
        legend([H1, H2(1), H3(1)], '\mu', '0.5\sigma', '\sigma', 'Location', 'Northwest');
        title('Malicious elders per section')
        drawnow;

        subplot(2,2,3)
        hold on
        H1 = plot(1:100:n, section_age_load_mean(1:100:end), 'LineWidth', 2, 'Color', 'k');
        H2 = plot(
            1:100:n,
            [max(section_age_load_mean(1:100:end) - 0.5*section_age_load_std(1:100:end),0);...
             section_age_load_mean(1:100:end) + 0.5*section_age_load_std(1:100:end)],
            'LineWidth', 2, 'Color', 'b'
        );
        H3 = plot(
            1:100:n,
            [max(section_age_load_mean(1:100:end) - section_age_load_std(1:100:end),0);...
             section_age_load_mean(1:100:end) + section_age_load_std(1:100:end)],
            'LineWidth', 2, 'Color', 'm'
        );
        hold off
        xlabel("Iterations");
        ylabel("Fraction");
        legend([H1, H2(1), H3(1)], '\mu', '0.5\sigma', '\sigma', 'Location', 'Northwest');
        title('Malicious age per section')
        drawnow;

        subplot(2,2,4)
        hold on
        H1 = plot(1:100:n, section_work_load_mean(1:100:end), 'LineWidth', 2, 'Color', 'k');
        H2 = plot(
            1:100:n,
            [max(section_work_load_mean(1:100:end) - 0.5*section_work_load_std(1:100:end),0);...
             section_work_load_mean(1:100:end) + 0.5*section_work_load_std(1:100:end)],
            'LineWidth', 2, 'Color', 'b'
        );
        H3 = plot(
            1:100:n,
            [max(section_work_load_mean(1:100:end) - section_work_load_std(1:100:end),0);...
             section_work_load_mean(1:100:end) + section_work_load_std(1:100:end)],
            'LineWidth', 2, 'Color', 'm'
        );
        hold off
        xlabel("Iterations");
        ylabel("Fraction");
        legend([H1, H2(1), H3(1)], '\mu', '0.5\sigma', '\sigma', 'Location', 'Northwest');
        title('Malicious work per section')
        drawnow;

        %subplot(2,2,4)
        figure(2)
        plot(
            1:100:n, stalled_sections_size(1:100:end), 'LineWidth', 2,
            1:100:n, stalled_sections_elders(1:100:end), 'LineWidth', 2,
            1:100:n, stalled_sections_age(1:100:end), 'LineWidth', 2,
            1:100:n, stalled_sections_work(1:100:end), 'LineWidth', 2
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
            num2str(section_size_mean(n)),'/',...
            num2str(section_size_std(n),2),'/',...
            num2str(min_section_size),...
            ', adversary: ', num2str(fraction_of_new_nodes_are_malicious),...
            ', a_n=', num2str(initial_network_age)]
        )
        xlabel("Iterations");
        ylabel("Fraction");
        drawnow
    end
end

%figure(1)
%print(["section_model_malicious_per_section_age_",num2str(initial_network_age), "_adversary_",num2str(fraction_of_new_nodes_are_malicious), "_section_size_", num2str(min_section_size),".png"],'-dpng');
%figure(2)
%print(["section_model_stallable_sections_age_",num2str(initial_network_age), "_adversary_",num2str(fraction_of_new_nodes_are_malicious), "_section_size_", num2str(min_section_size),".png"],'-dpng');
