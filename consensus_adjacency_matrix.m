function A = consensus_adjacency_matrix(rcom,y)
% OPINION algorithm
% calculate adjacency based on radius

%Determines number of rows in the input matrix
n_agents = size(y, 1);

%Initializes a nxn matrix filled with zeros
A = zeros(n_agents);

%Fill the matrix y with a 1 if the radius is close enough
for i = 1:n_agents
    for j = 1:n_agents
       if i ~= j
            agent1 = [y(i,1), y(i,2)];
            agent2 = [y(j,1), y(j,2)];
            vector = agent2 - agent1;
            magnitude = norm(vector);
            if magnitude <= rcom
               A(i,j) = 1;
            end
       end
    end
end
