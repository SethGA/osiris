function y1 = update_agents(y0, ydot, x1, dt, lNoise, rNoise)
%{
Let "n" be the number of agents.

--- Inputs ---
y0: n by 2 array holding all system outputs at timestep i
    e.g., y0(j,2) is the second output for agent-j

ydot: n by 2 array holding all system output time derivatives at timestep i
    e.g., ydot(j,2) is the second output time derivative for agent-j

x1: n by 2 array holding all consensus states for all agents at timestep i+1
    e.g., x1(j,2) is the second consensus state for agent-j

dt: timestep

lNoise: n by 1 array holding all left noise values. 

rNoise: n by 1 array holding all right noise values.

--- Outputs ---
y1:  n by 2 array holding all system outputs at timestep i+1
     Same structure as y0.

E1: 1 by n array holding energy data for all agents at timestep i+1.
    Same structure as E0.
%}

y1 = x1;

% Noise restriction example --------------------------------------------
% Here Noise will act as a source of skeptizism. 
% ie. If an agent is under right influence but subject to left Noise,
% Movement will be damped.
n = length(y0);     % number of agents
dvec = ydot * dt;
if (sum(lNoise) || sum(rNoise))
     lNoise = -(lNoise-1); rNoise = -(rNoise-1);
     for i=1:n
         if (dvec(i,1)>0), dvec(i,1) = dvec(i,1)*lNoise(i);
         else, dvec(i,1) = dvec(i,1)*rNoise(i); end
     end
end

% With velocity restriction --------------------------------------------
MAXD = 0.2;
if MAXD > 0
    n = length(y0);     % number of agents
    dmag = sqrt(sum(dvec.^2,2));
    for i=1:n
        r = dmag(i) / MAXD;
        if r > 1
            y1(i,:) = y0(i,:) + dvec(i,:) / r;
        end
    end
end

end