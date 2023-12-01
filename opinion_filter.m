function x1 = opinion_filter(x0, L, t, dt)
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

n = length(x0);
tau = 0 * ones(n,1);      % n by 1 signal delay vector


x1 = x0 - L * x0 *dt;
for i = 1:n
    %%The stationary points where agents do not move.
    if all(x0(i,:) == [10 20]) || all(x0(i,:) == [20 19])|| all(x0(i,:) == [25 12])|| all(x0(i,:) == [33 10])
        x1(i,:) = x0(i,:);
    end
    if tau(i) > t
        x1(i,:) = x0(i,:);
    end
end

end
