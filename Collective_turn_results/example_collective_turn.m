%% replay the example BOC 
clear;clc
addpath("../Utility/")
N = 100;
activateTime = 50;
turning_angle = 180; % 90
folder_name = "../Data/Simulation Data/collective turn/BOC_example_" + num2str(turning_angle);
[acc, r_area] = replay_collective_turn_snapshots(folder_name, 1, activateTime, turning_angle);
turning_angle = 90;
folder_name = "../Data/Simulation Data/collective turn/BOC_example_" + num2str(turning_angle);
[acc, r_area] = replay_collective_turn_snapshots(folder_name, 1, activateTime, turning_angle);
%% %% replay the example rand
clear;clc
addpath("../Utility/")
N = 100;
activateTime = 50;
turning_angle = 180; % 90
folder_name = "../Data/Simulation Data/collective turn/Random_example_" + num2str(turning_angle);
[acc, r_area] = replay_collective_turn_snapshots(folder_name, 1, activateTime, turning_angle);
turning_angle = 90;
folder_name = "../Data/Simulation Data/collective turn/Random_example_" + num2str(turning_angle);
[acc, r_area] = replay_collective_turn_snapshots(folder_name, 1, activateTime, turning_angle);
