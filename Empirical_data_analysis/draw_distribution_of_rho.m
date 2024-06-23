%% Draw the distribution of rho for the Body Orientation Change,corresponding to Figure1-(e) in the main text
clear;
clc;
addpath("../Data/")
addpath("../Data/Empirical Data/")
draw_alpha_list = [0, 1, 10];
Nc_list = [8, 10];
figure("Name","Body Orientation Change")
figSize_L = 10;
figSize_W = 12;
set(gcf,'Units','centimeter','Position',[5 5 figSize_L figSize_W]);
cnt = 0;
color_list = [223 41 47; 64 168 39; 28 60 145]./255;
for draw_alpha = draw_alpha_list
    cnt = cnt + 1;
        spearman_turning_all = [];
        for Nc = Nc_list
            load("../Data/Empirical Data/Nc=" + num2str(Nc) + "/coarse_data_seg_LF_BOC_occ_relation_revised_b_0.1.mat")
            alpha1_list = draw_alpha_list;
            for i = 1:size(LF_MS_relation.spearman_turning,2)
                if alpha1_list(i) == draw_alpha
                    spearman_turning = LF_MS_relation.spearman_turning{i};
                    spearman_turning_all = [spearman_turning_all;spearman_turning];
                end
            end
        end
        subplot(3,1,cnt)
        hist(cnt) = histogram(spearman_turning_all,'Normalization','probability');
        hist(cnt).EdgeColor = [1,1,1];
        hist(cnt).FaceColor = color_list(cnt,:);
        hist(cnt).FaceAlpha = 0.15;
        hist(cnt).EdgeAlpha = 0.15;
        hold on
        H2(cnt) = plot(hist(cnt).BinEdges(1:length(hist(cnt).Values)) + hist(cnt).BinWidth/2, ...
            smooth(hist(cnt).Values,5), 'Linewidth', 3,'Color', color_list(cnt,:));
        xline(0,"LineWidth",3,"Color",[0,0,0])
        ylabel("P(\rho)")
        xlabel("\rho")
        set(get(get(hist(cnt),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        legend(H2(cnt), "\alpha=" + draw_alpha_list(cnt),'Box','off') 
        set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
end
%% Draw the distribution of rho for the Distance,corresponding to Figure1-(f) in the main text
clear;
clc;
draw_alpha_list = [0, 1, 10];
Nc_list = [8, 10];
addpath("../Data/")
figure("Name","Distance")
figSize_L = 10;
figSize_W = 12;
set(gcf,'Units','centimeter','Position',[5 5 figSize_L figSize_W]);
cnt = 0;
color_list = [223 41 47; 64 168 39; 28 60 145]./255;
for draw_alpha = draw_alpha_list
    cnt = cnt + 1;
        spearman_turning_all = [];
        for Nc = Nc_list
            load("../Data/Empirical Data/Nc=" + num2str(Nc) + "/coarse_data_seg_LF_DS_occ_relation_revised_b_0.1.mat")
            alpha1_list = draw_alpha_list;
            for i = 1:size(LF_MS_relation.spearman_turning,2)
                if alpha1_list(i) == draw_alpha
                    spearman_turning = LF_MS_relation.spearman_turning{i};
                    spearman_turning_all = [spearman_turning_all;spearman_turning];
                end
            end
        end
        subplot(3,1,cnt)
        hist(cnt) = histogram(spearman_turning_all,'Normalization','probability');
        hist(cnt).EdgeColor = [1,1,1];
        hist(cnt).FaceColor = color_list(cnt,:);
        hist(cnt).FaceAlpha = 0.15;
        hist(cnt).EdgeAlpha = 0.15;
        hold on
        H2(cnt) = plot(hist(cnt).BinEdges(1:length(hist(cnt).Values)) + hist(cnt).BinWidth/2, ...
            smooth(hist(cnt).Values,5), 'Linewidth', 3,'Color', color_list(cnt,:));
        xline(0,"LineWidth",3,"Color",[0,0,0])
        ylabel("P(\rho)")
        xlabel("\rho")
        set(get(get(hist(cnt),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        legend(H2(cnt), "\alpha=" + draw_alpha_list(cnt),'Box','off') 
        set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
end
%% Draw the distribution of rho for the Bearing Change,corresponding to Figure1-(g) in the main text
clear;
clc;
draw_alpha_list = [0, 1, 10];
Nc_list = [8, 10];
addpath("../Data/")
figure("Name","Bearing Change")
figSize_L = 10;
figSize_W = 12;
set(gcf,'Units','centimeter','Position',[5 5 figSize_L figSize_W]);
cnt = 0;
color_list = [223 41 47; 64 168 39; 28 60 145]./255;
for draw_alpha = draw_alpha_list
    cnt = cnt + 1;
        spearman_turning_all = [];
        for Nc = Nc_list
            load("../Data/Empirical Data/Nc=" + num2str(Nc) + "/coarse_data_seg_LF_BC_occ_relation_revised_b_0.1.mat")
            alpha1_list = draw_alpha_list;
            for i = 1:size(LF_MS_relation.spearman_turning,2)
                if alpha1_list(i) == draw_alpha
                    spearman_turning = LF_MS_relation.spearman_turning{i};
                    spearman_turning_all = [spearman_turning_all;spearman_turning];
                end
            end
        end
        subplot(3,1,cnt)
        hist(cnt) = histogram(spearman_turning_all,'Normalization','probability');
        hist(cnt).EdgeColor = [1,1,1];
        hist(cnt).FaceColor = color_list(cnt,:);
        hist(cnt).FaceAlpha = 0.15;
        hist(cnt).EdgeAlpha = 0.15;
        hold on
        H2(cnt) = plot(hist(cnt).BinEdges(1:length(hist(cnt).Values)) + hist(cnt).BinWidth/2, ...
            smooth(hist(cnt).Values,5), 'Linewidth', 3,'Color', color_list(cnt,:));
        xline(0,"LineWidth",3,"Color",[0,0,0])
        ylabel("P(\rho)")
        xlabel("\rho")
        set(get(get(hist(cnt),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        legend(H2(cnt), "\alpha=" + draw_alpha_list(cnt),'Box','off') 
        set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
end