function A = opinion_adjacency_matrix(rComms, p0)
%{
Let "n" be the number of agents.

--- Inputs ---
rComms: n by 1 array holding radii of communication for each node

p0: n by 2 array holding X and Y position values
   e.g., y(j,2) is the second output for agent-j

--- Outputs ---
A:  n by n adjacency matrix - see notes for details on construction
%}

n = length(p0);      % number of agents

% communication params
K = 5;
SIGMA = 0.5;
BETA = 0.9;

A = NaN(n); 
for i=1:n
    currentRadiusComm = rComms(i);
    for j=1:n
        a2a = sqrt((p0(i,1) - p0(j,1))^2 + (p0(i,2) - p0(j,2))^2);
        if currentRadiusComm > 0
            % simple communication radius
            if a2a < currentRadiusComm
                A(i,j) = 1;
            else 
                A(i,j) = 0;
            end
        else
            % distance decay
            A(i,j) = K / (SIGMA^2 + a2a^2)^BETA;
        end
    end
end


end