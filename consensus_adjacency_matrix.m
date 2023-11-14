function A = consensus_adjacency_matrix(rcomm,y)
%{
Let "n" be the number of agents.

--- Inputs ---
y: n by 2 array holding output data for all agents
   e.g., y(j,2) is the second output for agent-j

--- Outputs ---
A:  n by n adjacency matrix - see notes for details on construction
%}

% get the number of agents, as the number of rows of y
n_agents = size(y, 1);

% initialize A to be all zeros
A = zeros(n_agents, n_agents);

% OPINION algorithm
% calculate adjacency based on radius
for i = 1:n_agents
    for j = 1:n_agents
       if i ~= j
            agent1 = [y(i,1), y(i,2)];
            agent2 = [y(j,1), y(j,2)];
            vector = agent2 - agent1;
            magnitude = norm(vector);
            if magnitude <= rcomm
               A(i,j) = 1;

end