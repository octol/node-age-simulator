clear all

init_network_iterations = 10000;
network_iterations = 10000;
network_size = 100000;

% Initial setup
nodes.work = 2.^linspace(3,8,network_size);
nodes.age = floor(log2(nodes.work));
nodes.malicious = logical(zeros(size(nodes.work)));

% Fully mature network
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
        semilogy(1:n, fraction_of_network_resetting, 1:n, fraction_of_work_resetting)
        legend("Network reset rate", "Work reset rate");
        grid on
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
    nodes.malicious(nodes_resetting) = logical(rand(numel(find(nodes_resetting)),1) < 0.1);
    nodes.initial(nodes_resetting) = false;

    % Collect network work stats
    network_work = sum(nodes.work);
    malicious_work = sum(nodes.work(nodes.malicious));
    frac_malicious_work(init_network_iterations + n) = malicious_work / network_work;

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
        figure(4)
        semilogy(1:length(fraction_of_network_resetting), fraction_of_network_resetting, 1:length(fraction_of_work_resetting), fraction_of_work_resetting)
        legend("Network reset rate", "Work reset rate");
        grid on
        drawnow

        figure(5)
        N = 1:length(fraction_of_network_resetting);
        plot(N, fraction_of_nodes_are_malicious, 'LineWidth', 2, N, frac_malicious_work, 'LineWidth', 2, N, frac_malicious_elder_work, 'LineWidth', 2);
        xlabel("Iteration")
        legend({"Malicious nodes", "Malicious work", "Malicious elder work"},"Location", "NorthWest");
        drawnow
    end
end

%figure(4)
%print -dpng simple_model_network_reset_rate_with_attack.png
%figure(5)
%print -dpng simple_model_malicious_fractions.png
