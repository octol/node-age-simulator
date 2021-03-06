clear all

%init_network_iterations = 10000;
init_network_iterations = 0;
network_iterations = 20000;
network_size = 100000;
init_network_age = 16;
fraction_of_new_nodes_are_malicious = 0.2;

% Initial setup
nodes.work = round(2.^((init_network_age-4)*rand(1,network_size) + 4));
nodes.age = floor(log2(nodes.work));
nodes.malicious = logical(zeros(size(nodes.work)));

% Run without adversaries first for awhile to check the work distribution assumption
for n = 1:init_network_iterations
    % All nodes does 1 unit of work w
    nodes.work += 1;
    nodes.age = floor(log2(nodes.work));

    % Reset nodes according to 1/w
    nodes_resetting = 1./nodes.work > rand(1, network_size);
    nodes.work(nodes_resetting) = 16;
    nodes.age(nodes_resetting) = log2(nodes.work(nodes_resetting));

    fraction_of_network_resetting(n) = sum(nodes_resetting) / network_size;
    fraction_of_work_resetting(n) = sum(nodes.work(nodes_resetting)) / sum(nodes.work);

    if mod(n, 100) == 0
        figure(1)
        hist(nodes.work, 50);
        xlabel('work');
        title(['Nodes: ',num2str(network_size),', Initial network age: ', num2str(init_network_age),', Iterations: ', num2str(n)])
        drawnow

        figure(2)
        hist(log2(nodes.work), 50);
        xlabel('log2(work)');
        title(['Nodes: ',num2str(network_size),', Initial network age: ', num2str(init_network_age),', Iterations: ', num2str(n)])
        drawnow

        figure(3)
        hist(nodes.age, 50);
        xlabel('age = floor(log2(work))')
        title(['Nodes: ',num2str(network_size),', Initial network age: ', num2str(init_network_age),', Iterations: ', num2str(n)])
        drawnow

        figure(4)
        semilogy(1:n, fraction_of_network_resetting, 1:n, fraction_of_work_resetting)
        legend("Network reset rate", "Work reset rate");
        grid on
        title(['Nodes: ',num2str(network_size),', Initial network age: ', num2str(init_network_age),', Iterations: ', num2str(n)])
        drawnow
    end
end

%figure(1)
%print -dpng simple_model_work_dist.png
%figure(2)
%print -dpng simple_model_log_work_dist.png
%figure(3)
%print -dpng simple_model_age_dist.png
%figure(4)
%print -dpng simple_model_network_reset_rate.png

fprintf("-- Adding adversary --\n");

% Keep track of how many nodes are lost
nodes.initial = logical(ones(size(nodes.work)));

% Add adversary
for n = 1:network_iterations
    % All nodes does 1 unit of work w
    nodes.work += 1;
    nodes.age = floor(log2(nodes.work));

    % Reset nodes according to 1/w
    nodes_resetting = 1./nodes.work > rand(1, network_size);
    nodes.work(nodes_resetting) = 16;
    nodes.age(nodes_resetting) = log2(nodes.work(nodes_resetting));
    nodes.malicious(nodes_resetting) = logical(rand(numel(find(nodes_resetting)),1) < fraction_of_new_nodes_are_malicious);
    nodes.initial(nodes_resetting) = false;

    % Collect network work stats
    network_work = sum(nodes.work);
    network_age = sum(nodes.age);
    malicious_work = sum(nodes.work(nodes.malicious));
    malicious_age = sum(nodes.age(nodes.malicious));
    frac_malicious_work(init_network_iterations + n) = malicious_work / network_work;
    frac_malicious_age(init_network_iterations + n) = malicious_age / network_age;

    fraction_of_nodes_are_malicious(init_network_iterations + n) = sum(nodes.malicious) / length(nodes.malicious);

    % Collect elder stats
    % Let's assume ~50% of adults are elders
    [sorted_age,I] = sort(nodes.age);
    elder_work = sum(nodes.work(I(end/2:end)));
    malicious_nodes_work = nodes.work.*nodes.malicious;
    elder_work_malicious = sum(malicious_nodes_work(I(end/2:end)));
    frac_malicious_elder_work(init_network_iterations + n) = elder_work_malicious / elder_work;

    fraction_of_network_resetting(init_network_iterations + n) = sum(nodes_resetting) / network_size;
    fraction_of_work_resetting(init_network_iterations + n) = sum(nodes.work(nodes_resetting)) / sum(nodes.work);

    fraction_of_initial_nodes_lost(init_network_iterations + n) = 1 - sum(nodes.initial)/length(nodes.initial);

    if mod(n, 100) == 0
        figure(1)
        hist(nodes.work, 50);
        xlabel('work');
        drawnow

        figure(2)
        hist(log2(nodes.work), 50);
        xlabel('log2(work)');
        drawnow

        figure(3)
        hist(nodes.age, 50);
        xlabel('age = floor(log2(work))')
        drawnow

        figure(4)
        semilogy(
            1:length(fraction_of_network_resetting), fraction_of_network_resetting,
            1:length(fraction_of_work_resetting), fraction_of_work_resetting
        )
        legend({'Network reset rate', 'Work reset rate'}, 'Location', 'West');
        xlabel('Iteration')
        grid on
        title(['Nodes: ',num2str(network_size),', Initial network age: ', num2str(init_network_age)])
        drawnow

        figure(5)
        N = 1:length(fraction_of_network_resetting);
        plot(
            N, fraction_of_nodes_are_malicious, 'LineWidth', 2,
            N, frac_malicious_age, 'LineWidth', 2,
            N, frac_malicious_work, 'LineWidth', 2
        );
        xlabel('Iteration')
        ylabel('Fraction')
        legend({'Malicious nodes', 'Malicious age', 'Malicious work'}, 'Location', 'NorthWest');
        title([
            'Nodes: ',num2str(network_size),...
            ', Initial network age: ', num2str(init_network_age),...
            ' , Adversary: ', num2str(fraction_of_new_nodes_are_malicious)
        ])
        drawnow
    end
end

%figure(4)
%print(['simple_model_network_reset_rate_with_attack_network_age_',num2str(init_network_age),'_adversary_',num2str(fraction_of_new_nodes_are_malicious),'.png'],'-dpng')
%figure(5)
%print(['simple_model_malicious_fractions_network_age_',num2str(init_network_age),'_adversary_',num2str(fraction_of_new_nodes_are_malicious),'.png'],'-dpng')
