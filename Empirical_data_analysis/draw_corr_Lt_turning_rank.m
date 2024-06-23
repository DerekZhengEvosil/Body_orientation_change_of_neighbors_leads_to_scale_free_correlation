%% Draw the correlation between the Lt and turning rank,corresponding to Figure S2 in the SI
Nc = 10;
load("../Data/Empirical Data/Nc=" + num2str(Nc) + "/coarse_data_seg_LF_turning_rank_relation.mat")
figure
figSize_L = 15;
figSize_W = 9;
set(gcf,'Units','centimeter','Position',[5 5 figSize_L figSize_W]);
myColor = [16 95 177;241 43 79;28 156 72];
for i = 1
    bar_spearman_mat_tmp = LF_MS_relation.bar_spearman_mat{i};
    bar_spearman_mat_tmp = bar_spearman_mat_tmp(~isnan(bar_spearman_mat_tmp));
    hist = histogram(bar_spearman_mat_tmp,'Normalization','probability');
    hist.EdgeColor = [1,1,1];
    hist.FaceColor = myColor(i,:)/255;
    hist.FaceAlpha = 0.15;
    hist.EdgeAlpha = 0.15;
    hold on
    plot(hist.BinEdges(1:length(hist.Values)) + hist.BinWidth/2, smooth(hist.Values,5),'Linewidth',3,'Color',myColor(i,:)/255)
end
hold on 
xline(0,"LineWidth",3,"Color",[0,0,0])
xlabel("\rho")
ylabel("P(\rho)")
set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
%% 
Nc = 8;
load("../Data/Empirical Data/Nc=" + num2str(Nc) + "/coarse_data_seg_LF_turning_rank_relation.mat")
figure
figSize_L = 15;
figSize_W = 9;
set(gcf,'Units','centimeter','Position',[5 5 figSize_L figSize_W]);
myColor = [16 95 177;241 43 79;28 156 72];
for i = 1
    bar_spearman_mat_tmp = LF_MS_relation.bar_spearman_mat{i};
    bar_spearman_mat_tmp = bar_spearman_mat_tmp(~isnan(bar_spearman_mat_tmp));
    hist = histogram(bar_spearman_mat_tmp,'Normalization','probability');
    hist.EdgeColor = [1,1,1];
    hist.FaceColor = myColor(i,:)/255;
    hist.FaceAlpha = 0.15;
    hist.EdgeAlpha = 0.15;
    hold on
    plot(hist.BinEdges(1:length(hist.Values)) + hist.BinWidth/2, smooth(hist.Values,5),'Linewidth',3,'Color',myColor(i,:)/255)
end
hold on 
xline(0,"LineWidth",3,"Color",[0,0,0])
xlabel("\rho")
ylabel("P(\rho)")
set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
