clear all

network_iterations = 10000;
network_size = 1000000;

% Initial setup
%nodes.work = randi([4,32], 1, network_size);
nodes.work = 4*ones(1,network_size);

for n = 1:network_iterations
    % All nodes does 1 unit of work w
    nodes.work += 1;

    % Reset nodes according to 1/w
    nodes_resetting = 1./nodes.work > rand(1, network_size);
    nodes.work(nodes_resetting) = 4;

    figure(1)
    %hist(nodes.work,100);
    hist(log2(nodes.work),100);
    drawnow
end



