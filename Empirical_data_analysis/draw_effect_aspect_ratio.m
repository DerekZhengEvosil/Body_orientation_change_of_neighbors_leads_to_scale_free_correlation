%% combine the Nc = 8/10 except front preference in occ,corresponding to FigureS8 in the SI
clear;
clc;
addpath("../Data/")
addpath("../Data/Empirical Data/")
draw_alpha_list = [0, 1];
Nc_list = [8, 10];
b_list = [0.1, 0.2, 0.3, 0.4, 0.5];
spearman_cell = {};
dataType = "Body Orientation Change"; % dataType = ["Distance" or "Bearing change" or "Body Orientation Change"]
draw_alpha = 1;
for b = b_list
    b_idx = find(b == b_list);
    spearman_turning_all = [];
    for Nc = Nc_list
        if dataType == "Distance"
            load("../Data/Empirical Data/Nc=" + num2str(Nc) + "/coarse_data_seg_LF_DS_occ_relation_revised" + "_b_" + num2str(b) + ".mat")
        elseif dataType == "Bearing change"
            load("../Data/Empirical Data/Nc=" + num2str(Nc) + "/coarse_data_seg_LF_BC_occ_relation_revised" + "_b_" + num2str(b) + ".mat")
        elseif dataType == "Body Orientation Change"
            load("../Data/Empirical Data/Nc=" + num2str(Nc) + "/coarse_data_seg_LF_BOC_occ_relation_revised" + "_b_" + num2str(b) + ".mat")
        end
        alpha_idx = find(draw_alpha == draw_alpha_list);
        spearman_turning = LF_MS_relation.spearman_turning{alpha_idx};
        spearman_turning_all = [spearman_turning_all;spearman_turning];
    end
    spearman_cell{b_idx} = spearman_turning_all;
end
figure
figSize_L = 7;
figSize_W = 12;
set(gcf,'Units','centimeter','Position',[5 5 figSize_L figSize_W]);
group = [repmat({"0.1"}, length(spearman_cell{1}), 1); ...
            repmat({"0.2"}, length(spearman_cell{2}), 1); ...
            repmat({"0.3"}, length(spearman_cell{3}), 1);...
            repmat({"0.4"}, length(spearman_cell{4}), 1);...
            repmat({"0.5"}, length(spearman_cell{5}), 1);];
boxplot([spearman_cell{1};spearman_cell{2};spearman_cell{3};spearman_cell{4};spearman_cell{5}], group,'symbol','')
hold on
yline(0,'LineStyle','--','LineWidth',1)
xlabel("aspect ratio");
ylabel("$\rho$",'Interpreter','latex');
title("$\alpha$=" + num2str(draw_alpha),'Interpreter','latex')
set(gca, 'Fontname', 'helvetica', 'FontSize', 15)