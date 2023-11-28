%% animate_network.m
%{
The purpose of this script is to allow the user to visualize how the
network evolves over the course of the simulation.  

                            *** NOTE *** 
Either run_lloyds.m, run_consensus.m, or run_opinion.m must have been successfully 
executed BEFORE animate_network.m is run.  This script uses the t, P and G 
arrays to produce the visualization and so the data must already exist in
the Workspace. 
%}
figure
%% PARAMETERS - can change w/o compromising script execution
MAX_EDGE_WEIGHT = 1;
UPDATE_LAG = 0.05;       % in seconds
R_COMM_EDGE_WEIGHT = 1;
%% SIMULATION
n=length(t)-1;
for i = 2:n 
    Gt = digraph(G(:,:,i),'omitselfloops');
    wts = MAX_EDGE_WEIGHT * Gt.Edges.Weight/max(Gt.Edges.Weight);
    p = plot(Gt,'LineWidth',wts);
    p.XData = P(i,:,1);
    p.YData = P(i,:,2);
    %plot the radius of communication for each agent
    for j = 1:nagents
        centers = [p.XData(j) p.YData(j)];
        viscircles(centers,rComms(j),'Color','black', 'LineWidth', R_COMM_EDGE_WEIGHT);
        text(centers(1,1)+0.75*rComms(j),centers(1,2)+0.75*rComms(j), string(j)); %label radius of communication
    end
    title(strcat('t = ',num2str(t(i),2)));
    xlabel('x');
    ylabel('y');
    drawnow;
    pause(UPDATE_LAG);
end