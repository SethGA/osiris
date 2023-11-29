%% HEADER 
%{
APSC 200 MODULE P2 - MTHE - CONSENSUS ALGORITHM TEMPLATE
Author: H. Wilson, B. Van Eden, Math & Eng., Queen's University
Release Date: 2023-06-01

Inputs:     opinion_input.csv = initial values for all agents
Outputs:    opinion_output.csv = time vs state data
            opinion_output.mat = time, state and network data 
                           (req'd for annimation)
            Time vs state and arena plots

Functions (user-defined):   opinion_adjacency_matrix.m
                            opinion_filter.m
                            update_agents.m

The I/O descriptions can be found in the preamble of each function.

animate_network.m generates a visualization of the network over the course
of this simulation.  It can only be run once this script has been
successfully executed.

--- Revision History ---
2023-06-01: Initial release

***************************************************************************
This script simulates the opinion algorithm outlined in the APSC 200 MTHE 
Course Manual. The script has been designed with transparency in mind.
(i.e., It allows you to see how all of the pieces fit together.) 

For this script to execute successfully, two tasks must be completed:
1: All of the functions listed above must be built.
2: The simulation parameters (marked by UPPERCASE letters) and the initial 
   values (specified in opinion_input.csv or opinion_input_2.csv) must be updated to match your
   application. 

Advanced users are welcome to alter this script to meet their needs; 
however, it will be done at their own risk. (i.e., Do not rely on the
instructor/TAs to help with large structual modifications to this script.)
***************************************************************************
%}

%% PARAMETERS - can edit w/o comprimising script execution

% simulation parameters
% Changes the time in which the simulation is ran.
TFINAL = 200;
% Not to sure set but I think it is how many steps each agent takes
NSTEPS = 10000;

N_DIMENSIONS = 2;
% 1 or 2 depending on Application

% plot toggles - set to 0 to suppress plot
SHOW_ARENA = 1;         % 2D plot showing agent paths
SHOW_OUTPUT = 1;        % output vs time plot
SHOW_OPINION = 1;     % consensus state vs time plot

%% SETUP

% importing initial values -----------------------------------------------
file = 'opinion_input';
if (N_DIMENSIONS == 2), file = strcat(file, '_2'); end

initval = readmatrix(file);  
nagents = length(initval);

rNoise = zeros(nagents, 1);
lNoise = zeros(nagents, 1);
if(N_DIMENSIONS == 1),lNoise = initval(:,3);rNoise = initval(:,4);end

% pre-allocation & initialization -----------------------------------------
P = zeros(NSTEPS, nagents, 2);% all output data 
for i = 1:N_DIMENSIONS
    P(:,:,i) = ones(NSTEPS,1) * initval(:,i+1)';
end
% P(i,j,k) is agent-j's k-th output at the i-th timestep; k = {1,2}
% If Dimension is 1, the P(i,j,2) will be zero

Q = zeros(NSTEPS, nagents, 2);         % all consensus state data
for i = 1:N_DIMENSIONS
    Q(:,:,i) = ones(NSTEPS,1) * initval(:,i+1)';
end
% same structure as P but for consensus state

Pdot = zeros(NSTEPS, nagents, 2);
% same structure as P but for time derivative of output

G = nan(nagents, nagents, NSTEPS-1);  % all network data
% G(:,:,i) is the Laplacian matrix at the i-th timestep


%% SIMULATION
t = linspace(0,TFINAL,NSTEPS)'; 
tstep = t(2);
rComms = initval(:,1);

for i = 2:NSTEPS-1
    % repackaging for easy use
    p0 = [P(i,:,1)' P(i,:,2)'];
    q0 = [Q(i,:,1)' Q(i,:,2)'];
    pdot = [Pdot(i,:,1)' Pdot(i,:,2)'];
    
    % ============================================================= WEEK 9
    % identifying which agents can communicate with one another
    G(:,:,i) = opinion_adjacency_matrix(rComms, p0); 
    
    D = diag(sum(G(:,:,i), 2));        % degree matrix
    L = D - G(:,:,i);               % Laplacian matrix

    % ============================================================= WEEK 10
    % updating consensus state of all agents
    q1 = opinion_filter(q0, L, t(i), tstep);
    
    % ============================================================= WEEK 11
    % updating output of all agents while 

    p1 = update_agents(p0, pdot, q1, tstep, lNoise, rNoise); %

    % saving new output and consensus state for plotting
    for k = 1:nagents
        for j = 1:2
            P(i+1,k,j) = p1(k,j)';
            Q(i+1,k,j) = q1(k,j)';
            Pdot(i+1,k,j) = (P(i,k,j) - P(i-1,k,j))/tstep;
        end
    end
end


%% PLOTTING
if (SHOW_ARENA)
    my_colours = ["#FF0000", "#00FF00", "#0000FF", "#00FFFF", ...
    "#FF00FF", "#FFFF00", "#0072BD", "#D95319", "#EDB120", ...
    "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"];

    figure 
    hold on
    
    % text offset for agent index
    tos = 0.01*[max(P(:,:,1),[],"all") max(P(:,:,1),[],"all")]; 

    for i=1:nagents
        % start = agent index
        clr = my_colours(mod(i,length(my_colours))+1);  
        text(P(1,i,1)+tos(1), P(1,i,2)+tos(2), ...
            num2str(i), 'Color', clr);     

        % agent path
        plot(P(:,i,1), P(:,i,2), "Color", clr);       

        % end = star
        plot(P(end,i,1), P(end,i,2), "Marker", "pentagram",... 
            "Color",clr,"MarkerFaceColor",clr);         
    end
    xlabel("p1");
    ylabel("p2");
    hold off
end

legend_labels = cell(nagents,1);
for i = 1:nagents
    legend_labels{i} = strcat("agent ", num2str(i));
end

if (SHOW_OUTPUT)    
    figure
    subplot(211)
    plot(t, P(:,:,1))
    title("Output")
    xlabel("t")
    ylabel("p1")
    legend(legend_labels);
    
    subplot(212)
    plot(t, P(:,:,2))
    xlabel("t")
    ylabel("p2")
end

if (SHOW_OPINION) 
    figure
    subplot(211)
    plot(t, Q(:,:,1))
    title("Consensus State")
    xlabel("t")
    ylabel("q1")
    legend(legend_labels);
    
    subplot(212)
    plot(t, Q(:,:,2))
    xlabel("t")
    ylabel("q2")
end

%% EXPORTING DATA

% --- Matlab - for arena animation ---
save("opinion_output.mat", "t", "P", "G");  

% --- Excel/Other - for further analysis ---
p1_ = P(:,:,1);
p2_ = P(:,:,2);
q1_ = Q(:,:,1);
q2_ = Q(:,:,2);
output_table = [array2table(t) array2table(p1_) array2table(p2_) ...
    array2table(q1_) array2table(q2_)];
writetable(output_table,"opinion_output.csv"); 

%{
***************************** CSV headers ******************************
You will notice the headers in the output csv are a awkward to read.
Lets take "q2_5" as an example:
  2 - indicates it is the second element (y dir) of the consensus state
  5 - indicates that this for agent 5
%}
