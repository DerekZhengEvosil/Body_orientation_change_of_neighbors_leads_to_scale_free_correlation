addpath("../Data/Robots Data/")
addpath("../Utility/")
%% Random-based interaction
rand = load("EXPData_random_example.mat");
analysisRobots_performance(rand.G)