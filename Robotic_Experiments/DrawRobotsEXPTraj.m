%% Visualize the trajectory of robotic experiments for BOC-based interaction and random-based interaction.
clear;clc;
addpath("../Data/Robots Data/")
addpath("../Utility/")
load("robotic_experiments_data_BOC.mat")
load("robotic_experiments_data_random.mat")
episode = 1;
DrawTraj(robots_traj_BOC{episode}{1,1});
DrawTraj(robots_traj_rand{episode}{1,1});
function DrawTraj(actors)
    G = struct;
    G.actor = actors;
    G.simStep = length(G.actor{1}.memory);
    G.maxID = 50;
    dt = 0.5;
    for t = 1:G.simStep
        posDir = [];
        p_posDir = [];
        for i = 1:G.maxID
            heading(t,i) = vel2heading_deg(G.actor{i}.memory(t,[3,4]));
            posDir(i,[1,2,3,4]) = G.actor{i}.memory(t,[1,2,3,4]);
            if t>=2
                p_posDir(i,[1,2,3,4]) = G.actor{i}.memory(t-1,[1,2,3,4]);
                speed(t,i) = norm(posDir(i,[1,2])-p_posDir(i,[1,2]))/dt;
                rotRate(t,i) = acosd(dot(posDir(i,[3,4]),p_posDir(i,[3,4])))/dt;
            end
        end
    end
    figure;
    r = 30;          
    arrow_scale = 30;   
    for i = 1:G.maxID
        pos = G.actor{i}.pose;
        vel = G.actor{i}.vel;
        tailTraj = G.actor{i}.memory(:,[1,2]);
        if ~isnan(pos)
            quiver(pos(1),pos(2),arrow_scale*vel(1),arrow_scale*vel(2),0,'k','linewidth',1); hold on;
            line(tailTraj(2:end,1),tailTraj(2:end,2),'linestyle','-','linewidth',0.5); hold on;
            rectangle('Position', [pos(1)-r, pos(2)-r, r*2, r*2], 'Curvature', [1 1]); hold on;
        end
    end
    % stable tracking regions
    rectangle('Position', [-2790,-2890,5380, 5680], 'Curvature', [0 0],'linewidth',2); hold on;
    % experimental arena border
    rectangle('Position', [-2100,-2350,4200, 4700], 'Curvature', [0 0],'linewidth',3,'edgecolor',[0,0,1]); hold on;
    if isfield(G,'r_sense')
        pos = G.actor{1}.pose;
        r = G.r_sense;
        rectangle('Position', [pos(1)-r, pos(2)-r, r*2, r*2], 'Curvature', [1 1],'edgecolor',[0,1,0]); hold on;
    end
    box on; grid on; axis equal; 
end