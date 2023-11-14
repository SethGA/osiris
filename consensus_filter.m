function x1 = consensus_filter(x0, L, t, dt)
%{
Let "n" be the number of agents.

--- Inputs ---
x0: n by 2 array holding consensus state information for timestep i
    e.g., x0(j,1) is the first consensus state for agent-j

L:  n by n Laplacian matrix

t and dt: the current time and timestep, respectively. 

--- Outputs ---
x1: n by 2 array holding consensus state information for timestep i+1
    Same structure as x0.
%}


end