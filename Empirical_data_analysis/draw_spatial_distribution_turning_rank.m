%% Draw the spatial distribution of turning rank, corresponding to Figure S3-S4 in the main text
clear;clc;
addpath("../Utility/")
addpath("../Data/")
addpath("../Data/Empirical Data/")
load("../Data/Empirical Data/Nc=8/turning_rank_distribution.mat")
figure
figSize_L = 20;
figSize_W = 50;
for k = 1:numel(turning_rank_tmp)
    subplot(5,2,k)
    set(gcf,'Units','centimeter','Position',[5 5 figSize_L figSize_W]);
    scatter(all_pos(1,:)',all_pos(2,:)',10,'MarkerEdgeColor',[127 127 127]/255,...
                  'MarkerFaceColor',[127 127 127]/255)
    hold on
    box on
    dscatter(rank_pos_id{k}(1,:)', rank_pos_id{k}(2,:)')
    hold on 
    yline(0,"--k","LineWidth",2)
    colormap(jet)
    grid on
    xlabel("x(mm)")
    ylabel("y(mm)")
    legend(["others"],"Box","off")
    h = colorbar;
    h.Label.String = 'probability density estimate';
    title("turning rank = " + num2str(k))
    axis equal
    set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
end
%%
clear;clc;
addpath("../Utility/")
addpath("../Data/")
load("../Data/Empirical Data/Nc=10/turning_rank_distribution.mat")
figure
figSize_L = 20;
figSize_W = 50;
for k = 1:numel(turning_rank_tmp)
    subplot(5,2,k)
    set(gcf,'Units','centimeter','Position',[5 5 figSize_L figSize_W]);
    scatter(all_pos(1,:)',all_pos(2,:)',10,'MarkerEdgeColor',[127 127 127]/255,...
                  'MarkerFaceColor',[127 127 127]/255)
    hold on
    box on
    dscatter(rank_pos_id{k}(1,:)',rank_pos_id{k}(2,:)')
    hold on 
    yline(0,"--k","LineWidth",2)
    colormap(jet)
    grid on
    xlabel("x(mm)")
    ylabel("y(mm)")
    legend(["others"],"Box","off")
    h = colorbar;
    h.Label.String = 'probability density estimate';
    title("turning rank = " + num2str(k))
    axis equal
    set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
end