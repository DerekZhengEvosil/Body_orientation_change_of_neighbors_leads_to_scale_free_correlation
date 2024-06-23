addpath("../Data/Robots Data/")
addpath("../Utility/")
%% BOC-based interaction
BOC = load("EXPData_BOC_example.mat");
analysisRobots_performance(BOC.G)