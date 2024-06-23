%% replay BOC-based interaction in collective spin, Fig2(a-c) in main text
clear;clc;
addpath("../Utility")
data_folder = "../Data/Simulation Data/collective spin/BOC_example/";
activateTime = 25;
[heading, op, transfer_Ang] = replay_collective_spin_snapshots(data_folder, 1, activateTime);
%% replay random-based interaction in collective spin, Fig2(e-f) in main text
clear;clc;
addpath("../Utility")
data_folder = "../Data/Simulation Data/collective spin/randam_example/";
activateTime = 25;
[heading, op, transfer_Ang] = replay_collective_spin_snapshots(data_folder, 1, activateTime);