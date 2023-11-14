function [y1, E1] = update_agents(y0, ydot, x1, E0, dt)
%{
Let "n" be the number of agents.

--- Inputs ---
y0: n by 2 array holding all system outputs at timestep i
    e.g., y0(j,2) is the second output for agent-j

ydot: n by 2 array holding all system output time derivatives at timestep i
    e.g., ydot(j,2) is the second output time derivative for agent-j

x1: n by 2 array holding all consensus states for all agents at timestep i+1
    e.g., x1(j,2) is the second consensus state for agent-j

E0:  1 by n array holding energy data for all agents at timestep i
    e.g., E0(j) is the energy of agent-j at timestep i  

dt: timestep

--- Outputs ---
y1:  n by 2 array holding all system outputs at timestep i+1
     Same structure as y0.

E1: 1 by n array holding energy data for all agents at timestep i+1.
    Same structure as E0.
%}


end