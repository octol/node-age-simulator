clear all

init_network_iterations = 1000;
network_iterations = 10000;
%network_size = 10;
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

    if mod(n, 10) == 0
        figure(1)
        %hist(nodes.work,100);
        hist(log2(nodes.work),100);
        drawnow
    end
end

fprintf("Adding adversary\n");

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

    % Collect network work stats
    network_work = sum(nodes.work);
    malicious_work = sum(nodes.work(nodes.malicious));
    frac_malicious_work(n) = malicious_work / network_work;

    % Collect elder stats
    % Let's assume ~50% of adults are elders
    [sorted_age,I] = sort(nodes.age);
    elder_work = sum(nodes.work(I(end/2:end)));
    malicious_nodes_work = nodes.work.*nodes.malicious;
    elder_work_malicious = sum(malicious_nodes_work(I(end/2:end)));
    frac_malicious_elder_work(n) = elder_work_malicious / elder_work;

    if mod(n, 10) == 0
        fraction_of_network_resetting = sum(nodes_resetting) / network_size
        fraction_of_work_resetting = sum(nodes.work(nodes_resetting)) / sum(nodes.work)

        figure(1)
        %hist(nodes.work,100);
        hist(log2(nodes.work),100);
        xlabel("Node age")
        drawnow

        figure(2)
        plot(frac_malicious_work, 'b-','LineWidth',2);
        xlabel("Iteration")
        title("Fraction of malicious work")
        drawnow

        figure(3)
        plot(frac_malicious_elder_work, 'b-','LineWidth',2);
        xlabel("Iteration")
        title("Fraction of malicious elder work")
        drawnow
    end
end
