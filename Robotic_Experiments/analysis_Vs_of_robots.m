%% The analysis of information transfer speed with the increasing group size
addpath("../Utility/")
load('../Data/Robots Data/ascend_delay_all.mat');
load('../Data/Robots Data/ascend_trans_dis_all.mat');
figure
figSize_L = 10;
figSize_W = 3.5;
set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
N = 50;
Type_list = ["BOC", "Rand"];
slope_cell = {};
slope = [];
for Type = Type_list
    Type_idx = find(Type == Type_list);
    for idx = 1:2
        ascend_delay = ascend_delay_all{idx, Type_idx};
        ascend_trans_dis = ascend_trans_dis_all{idx, Type_idx};
        x = nanmean(ascend_delay,1)';
        y = nanmean(ascend_trans_dis, 1)';
        std_x = nanstd(ascend_delay, 1)';
        std_y = nanstd(ascend_trans_dis, 1)';
        if Type == "BOC"
            linear_region = find(x>1&x<50);
            linear_start = linear_region(1);
            linear_end = linear_region(end);
        else
            linear_region = find(x>1&x<50);
            linear_start = linear_region(1);
            linear_end = linear_region(end);
        end
        % 定义线性模型
        linearModel = fit(x(linear_start:linear_end), y(linear_start:linear_end), 'poly1');  % 'poly1' 表示线性模型
        % 提取拟合的系数
        coefficients = coeffvalues(linearModel);
        % 提取斜率和截距
        slope(Type_idx, idx) = coefficients(1);
        intercept = coefficients(2);
        % 生成拟合的直线
        x_fit = linspace(min(x), max(x), 100);  % 生成用于绘制直线的 x 值
        y_fit = polyval(coefficients, x_fit);    % 计算对应的 y 值
        % 创建散点图
        if Type_idx == 1
            draw_x_idx = find(x>0&x<20);
            if idx == 1
                plot(x(draw_x_idx(1:4:end)),y(draw_x_idx(1:4:end)),'ko', 'LineWidth',2,'MarkerEdgeColor',[11,92,175]./255, ...
                                                            'MarkerFaceColor',[11,92,175]./255,'MarkerSize',5)
            else
                plot(x(draw_x_idx(1:4:end)),y(draw_x_idx(1:4:end)),'k^', 'LineWidth',2,'MarkerEdgeColor',[11,92,175]./255, ...
                                                            'MarkerFaceColor',[11,92,175]./255,'MarkerSize',5)
            end
        else
            draw_x_idx = find(x>0&x<20);
            if idx == 1
                plot(x(draw_x_idx),y(draw_x_idx),'ko', 'LineWidth',2,'MarkerEdgeColor',[205,62,21]./255, ...
                                                            'MarkerFaceColor',[205,62,21]./255,'MarkerSize',5)
            else
                plot(x(draw_x_idx),y(draw_x_idx),'k^', 'LineWidth',2,'MarkerEdgeColor',[205,62,21]./255, ...
                                                            'MarkerFaceColor',[205,62,21]./255,'MarkerSize',5)
            end
        end
        hold on;  % 保持图形以便后续添加直线
    end
end
xlabel("turning delay (s)")
ylabel("$d_{i}^{turn}$(mm)",'Interpreter','latex')
legend(["BOC turn 1","BOC turn 2", "random turn 1", "random turn 2"],'Box','off')
% legend_list_1 = "slope = " + round(slope(1,:)',1);
% legend_list_2 = "slope = " + round(slope(2,:)',1);
% legend([legend_list_1;legend_list_2],'Box','off')
set(gca, 'Fontname', 'helvetica', 'FontSize', 7)