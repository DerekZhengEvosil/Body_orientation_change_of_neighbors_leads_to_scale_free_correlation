%% mean acc with rand
clear;clc
addpath("../Utility/")
load("../Data/Robots Data/BOC_resp_info.mat")
load("../Data/Robots Data/random_resp_info.mat")
type_list_draw = ["BOC", "Random"];
Proj_acc = BOC_resp_info{1};
rand_acc = random_resp_info{1};
violin_data = [Proj_acc  rand_acc];
violin_label = [reshape(repmat(type_list_draw(1), 1, size(Proj_acc , 2)), 1, []),... 
                reshape(repmat(type_list_draw(2), 1, size(rand_acc, 2)), 1, [])];
figure
figSize_L = 6;
figSize_W = 8;
set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
vs = violinplot(violin_data, violin_label, 'ShowMean', true, 'ShowMedian', false, 'MarkerSize', 15);
ylabel("accuracy")
set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
% mean resp with rand
type_list_draw = ["BOC", "Random"];
Proj_area = BOC_resp_info{2};
rand_area = random_resp_info{2};
violin_data = [Proj_area, rand_area];
violin_label = [reshape(repmat(type_list_draw(1), 1, size(Proj_area , 2)), 1, []),... 
                reshape(repmat(type_list_draw(2), 1, size(rand_area, 2)), 1, [])];
figure
figSize_L = 6;
figSize_W = 8;
set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
vs = violinplot(violin_data, violin_label, 'ShowMean', true, 'ShowMedian', false, 'MarkerSize', 15);
ylabel("responsiveness")
set(gca, 'Fontname', 'helvetica', 'FontSize', 15)

type_list_draw = ["BOC", "Random"];
figure
figSize_L = 12;
figSize_W = 5;
set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
Maxstep = 299;
Seg = 5;
cyctime = 0.5;
time_x = [1:Seg:Maxstep];
proj_acc = BOC_resp_info{3};
rand_acc = random_resp_info{3};
shadedErrorBar(time_x.* cyctime, proj_acc(:,time_x), {@mean,@std}, 'lineProps', ...
                                                    {'-ob','Color', [11, 92, 175]./255, ...
                                                    'markerfacecolor',[11, 92, 175]./255}, ...
                                                    'patchSaturation',0.1)
hold on
shadedErrorBar(time_x.* cyctime, rand_acc(:,time_x), {@mean,@std}, 'lineprops', ...
                                                    {'-or','Color', [242, 41, 79]./255, ...
                                                    'markerfacecolor',[242, 41, 79]./255}, ...
                                                    'patchSaturation',0.1)
xlabel("time (s)")
ylabel("accuracy")
box on
legend(type_list_draw,'Box','off')
set(gca, 'Fontname', 'helvetica', 'FontSize', 15)