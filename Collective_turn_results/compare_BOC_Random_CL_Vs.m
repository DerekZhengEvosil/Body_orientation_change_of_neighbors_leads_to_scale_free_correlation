%% Correlation Length 
clear;clc
addpath("../Utility/")
load( "../Data/Simulation Data/collective turn/BOC_corr_data_all.mat")
N_list = [10 20 30 40 50 60 70 80 90 100];
turning_angle_list = [90, 180];
mean_r_all = turn_all_data{4};
mean_C_all = turn_all_data{5};
mean_GroupSZ_all = turn_all_data{6};
mean_CL_all = turn_all_data{7};
for turning_angle = turning_angle_list
    for N = N_list
        N_idx = find(N == N_list);
        turning_angle_idx = find(turning_angle == turning_angle_list);
        mean_r = mean_r_all{N_idx, turning_angle_idx};
        mean_C = mean_C_all{N_idx, turning_angle_idx};
        mean_FS = mean_GroupSZ_all{N_idx, turning_angle_idx};
        CL{N_idx, turning_angle_idx} = mean_CL_all{N_idx, turning_angle_idx};
        FS{N_idx, turning_angle_idx} = mean_FS;
    end
end
N_list = [10 20 30 40 50 60 70 80 90 100];
draw_N = [20 30 40 50 60 70 80 90 100];
color_list = [[228,26,28]./255; [228,26,28]./255];
figure
figSize_L = 12;
figSize_W = 5;
set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
turning_angle_label = ["\pi/2$", "\pi$"];
for turning_angle = turning_angle_list
    turning_angle_idx = find(turning_angle == turning_angle_list);
    subplot(1,2,turning_angle_idx)
    for i = 1:length(draw_N)
        idx = find(draw_N(i) == N_list);
        FS_tmp = (FS{idx, turning_angle_idx});
        CL_tmp = (CL{idx, turning_angle_idx});
        SZ(i) = nanmean(FS_tmp(:));
        CorL(i) = nanmean(CL_tmp(:));
        yneg(i) = std(CL_tmp(:));
        ypos(i) = std(CL_tmp(:));
        xneg(i) = std(FS_tmp(:));
        xpos(i) = std(FS_tmp(:));
    end
    if turning_angle_idx == 1
        h_b(turning_angle_idx) = errorbar(SZ, CorL, yneg, ypos, xneg, xpos, "o", "color", color_list(turning_angle_idx,:), "MarkerSize", 6,...
         "MarkerEdgeColor", color_list(turning_angle_idx,:), "MarkerFaceColor", color_list(turning_angle_idx,:),'LineWidth',1);
    else
        h_b(turning_angle_idx) = errorbar(SZ, CorL, yneg, ypos, xneg, xpos, "^", "color", color_list(turning_angle_idx,:), "MarkerSize", 6,...
         "MarkerEdgeColor", color_list(turning_angle_idx,:), "MarkerFaceColor", color_list(turning_angle_idx,:),'LineWidth',1);
    end
%     h_e(turning_angle_idx) = errorbar(SZ, CorL, yneg, ypos, "o", "color", color_list(turning_angle_idx,:), "MarkerSize", 8,...
%      "MarkerEdgeColor", color_list(turning_angle_idx,:), "MarkerFaceColor", color_list(turning_angle_idx,:),'LineWidth',1.5);
    hold on
    [xData, yData] = prepareCurveData(SZ, CorL);
    % Set up fittype and options.
    ft = fittype( 'poly1' );
    % Fit model to data.
    [fitresult, gof] = fit(xData, yData, ft);
    x_fit = [xData(1)-5:xData(end)+300];
    y_fit = feval(fitresult, x_fit);
    plot(x_fit, y_fit,'color',color_list(turning_angle_idx,:), 'linewidth',1.5);
%     hold on 
%     yline(10,'k--','lineWidth',2)
    xlabel("Flock Size",'Interpreter','latex')
    ylabel("Correlation Length", 'Interpreter','latex')
%     legend(h_e(turning_angle_idx), "$\theta_{info}=" + turning_angle_label(turning_angle_idx),'box','off','Interpreter','latex')
%     set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
    set(gca, 'Fontname', 'helvetica', 'FontSize', 9)
end
hold on
load("../Data/Simulation Data/collective turn/Random_corr_data_all.mat")
N_list = [10 20 30 40 50 60 70 80 90 100];
turning_angle_list = [90, 180];
mean_r_all = turn_all_data{4};
mean_C_all = turn_all_data{5};
mean_GroupSZ_all = turn_all_data{6};
mean_CL_all = turn_all_data{7};
for turning_angle = turning_angle_list
    for N = N_list
        N_idx = find(N == N_list);
        turning_angle_idx = find(turning_angle == turning_angle_list);
        mean_r = mean_r_all{N_idx, turning_angle_idx};
        mean_C = mean_C_all{N_idx, turning_angle_idx};
        mean_FS = mean_GroupSZ_all{N_idx, turning_angle_idx};
        CL{N_idx, turning_angle_idx} = mean_CL_all{N_idx, turning_angle_idx};
        FS{N_idx, turning_angle_idx} = mean_FS;
    end
end
N_list = [10 20 30 40 50 60 70 80 90 100];
draw_N = [20 30 40 50 60 70 80 90 100];
color_list = [[55,126,184]./255; [55,126,184]./255];
turning_angle_label = ["\pi/2$", "\pi$"];
for turning_angle = turning_angle_list
    turning_angle_idx = find(turning_angle == turning_angle_list);
    subplot(1,2,turning_angle_idx)
    for i = 1:length(draw_N)
        idx = find(draw_N(i) == N_list);
        SZ(i) = nanmean(FS{idx, turning_angle_idx});
        CorL(i) = nanmean(CL{idx, turning_angle_idx});
        yneg(i) = std(CL{idx, turning_angle_idx});
        ypos(i) = std(CL{idx, turning_angle_idx});
        xneg(i) = std(FS{idx, turning_angle_idx});
        xpos(i) = std(FS{idx, turning_angle_idx});
    end
    CorL(isnan(CorL)) = 1;
    if turning_angle_idx == 1
        h_r(turning_angle_idx) = errorbar(SZ, CorL, yneg, ypos, xneg, xpos, "o",  "color", color_list(turning_angle_idx,:), "MarkerSize", 6,...
         "MarkerEdgeColor", color_list(turning_angle_idx,:), "MarkerFaceColor", color_list(turning_angle_idx,:),'LineWidth',1);
    else
        h_r(turning_angle_idx) = errorbar(SZ, CorL, yneg, ypos, xneg, xpos, "^",  "color", color_list(turning_angle_idx,:), "MarkerSize", 6,...
         "MarkerEdgeColor", color_list(turning_angle_idx,:), "MarkerFaceColor", color_list(turning_angle_idx,:),'LineWidth',1);
    end
%     h_e(turning_angle_idx) = errorbar(SZ, CorL, yneg, ypos, "o", "color", color_list(turning_angle_idx,:), "MarkerSize", 8,...
%      "MarkerEdgeColor", color_list(turning_angle_idx,:), "MarkerFaceColor", color_list(turning_angle_idx,:),'LineWidth',1.5);
    hold on 
    [xData, yData] = prepareCurveData(SZ, CorL);
    % Set up fittype and options.
    ft = fittype( 'smoothingspline' );
    % Fit model to data.
    [fitresult, gof] = fit(xData, yData, ft);
    x_fit = [xData(1)-5:xData(end)+10];
    y_fit = feval(fitresult, x_fit);
    plot(x_fit, y_fit,'color',color_list(turning_angle_idx,:), 'linewidth',1.5);
    xlim([1200, 3500])
    xlabel("Flock Size",'Interpreter','latex')
    ylabel("Correlation Length", 'Interpreter','latex')
    set(gca, 'Fontname', 'helvetica', 'FontSize', 9)
end
% legend([h_b(1), h_b(2), h_r(1), h_r(2)], ["BOC", "BOC", "Random", "Random"],'box','off','Interpreter','latex')
% set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
%% information transfer 
clear;clc
addpath("../Utility/")
load( "../Data/Simulation Data/collective turn/BOC_corr_data_all.mat")
N_list = [10 20 30 40 50 60 70 80 90 100];
% draw_N_list = [50 100 200];
draw_N_list = [10 20 30 40 50 60 70 80 90 100];
turning_angle_list = [90, 180];
slope = {};
transfer_Ang_all = turn_all_data{1};
ascend_delay_all = turn_all_data{2};
ascend_trans_dis_all = turn_all_data{3};
for turning_angle = turning_angle_list
    turning_angle_idx = find(turning_angle == turning_angle_list);
    N_cnt = 0;
    for N = draw_N_list
        N_cnt = N_cnt + 1;
        N_idx = find(N == N_list);
        ascend_delay = ascend_delay_all{N_idx, turning_angle_idx};
        ascend_trans_dis = ascend_trans_dis_all{N_idx, turning_angle_idx};
        right_idx = find(ascend_delay(1,:) >= 0);
        coefficients_list = [];
        cnt = 0;
        for r_idx = right_idx
            cnt = cnt + 1;
            x = ascend_delay(:,r_idx);
            y = ascend_trans_dis(:,r_idx);
            linear_start = floor(N * 0.3);
            linear_end = floor(N * 0.75);
            linearModel = fit(x(linear_start:linear_end), y(linear_start:linear_end), 'poly1');  % 'poly1' 表示线性模型
            coefficients = coeffvalues(linearModel);
            intercept = coefficients(2);
            coefficients_list(cnt) = coefficients(1);
        end
        slope{N_cnt, turning_angle_idx} = coefficients_list;
    end
end
N_list = [10 20 30 40 50 60 70 80 90 100];
color_list = [[228,26,28]./255; [228,26,28]./255];
turning_angle_label = ["\pi/2$", "\pi$"];
mean_slope = cellfun(@nanmean, slope) * 60;
err_slope = cellfun(@nanstd, slope) * 60;
figure
figSize_L = 12;
figSize_W = 3;
set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
for turning_angle = turning_angle_list
    turning_angle_idx = find(turning_angle == turning_angle_list);
    subplot(1,2,turning_angle_idx)
    if turning_angle_idx == 1
        h_b(turning_angle_idx) = errorbar(N_list, mean_slope(:, turning_angle_idx), err_slope(:, turning_angle_idx), "o", 'Color','k', 'MarkerEdgeColor','k','MarkerFaceColor', ...
                            color_list(turning_angle_idx,:),'MarkerSize', 8,'LineWidth',1);
    else
        h_b(turning_angle_idx) = errorbar(N_list, mean_slope(:, turning_angle_idx), err_slope(:, turning_angle_idx), "^", 'Color','k', 'MarkerEdgeColor','k','MarkerFaceColor', ...
                            color_list(turning_angle_idx,:),'MarkerSize', 8,'LineWidth',1);
    end
    hold on
    [xData, yData] = prepareCurveData(N_list, mean_slope(:,turning_angle_idx) );
    % Set up fittype and options.
    ft = fittype( 'poly1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    x_fit = [N_list(1)-10:N_list(end)+10];
    y_fit = feval(fitresult, x_fit);
    plot(x_fit, y_fit,'color',color_list(turning_angle_idx,:), 'linewidth',2.5);
    hold on
    xlabel("$N$",'interpreter','latex')
    ylabel("$V_s$",'Interpreter','latex')
    set(gca, 'Fontname', 'helvetica', 'FontSize', 9)
end

load("../Data/Simulation Data/collective turn/Random_corr_data_all.mat")
N_list = [10 20 30 40 50 60 70 80 90 100];
draw_N_list = [10 20 30 40 50 60 70 80 90 100];
turning_angle_list = [90, 180];
slope = {};
transfer_Ang_all = turn_all_data{1};
ascend_delay_all = turn_all_data{2};
ascend_trans_dis_all = turn_all_data{3};
for turning_angle = turning_angle_list
    turning_angle_idx = find(turning_angle == turning_angle_list);
    N_cnt = 0;
    for N = draw_N_list
        N_cnt = N_cnt + 1;
        N_idx = find(N == N_list);
        ascend_delay = ascend_delay_all{N_idx, turning_angle_idx};
        ascend_trans_dis = ascend_trans_dis_all{N_idx, turning_angle_idx};
        right_idx = find(ascend_delay(1,:) >= 0);
        coefficients_list = [];
        cnt = 0;
        for r_idx = right_idx
            cnt = cnt + 1;
            x = ascend_delay(:,r_idx);
            y = ascend_trans_dis(:,r_idx);
            x(isnan(x)) = 0;
            y(isnan(y)) = 0;
            linear_start = floor(N * 0.3);
            linear_end = floor(N * 0.75);
            linearModel = fit(x(linear_start:end), y(linear_start:end), 'poly1');  % 'poly1' 表示线性模型
            coefficients = coeffvalues(linearModel);
            intercept = coefficients(2);
            coefficients_list(cnt) = coefficients(1);
        end
        data_liers = coefficients_list(:); 
        slope{N_cnt, turning_angle_idx} = data_liers;
    end
end
figure
figSize_L = 12;
figSize_W = 3;
set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
N_list = [10 20 30 40 50 60 70 80 90 100];
color_list = [[55,126,184]./255; [55,126,184]./255];
turning_angle_label = ["\pi/2$", "\pi$"];
mean_slope = cellfun(@nanmean, slope) * 60;
err_slope = cellfun(@nanstd, slope) * 60;
for turning_angle = turning_angle_list
    turning_angle_idx = find(turning_angle == turning_angle_list);
    subplot(1,2,turning_angle_idx)
    if turning_angle_idx == 1
        h_r(turning_angle_idx) = errorbar(N_list, mean_slope(:, turning_angle_idx), err_slope(:, turning_angle_idx), "o", 'Color','k', 'MarkerEdgeColor','k','MarkerFaceColor', ...
                            color_list(turning_angle_idx,:),'MarkerSize', 8,'LineWidth',1);
    else
        h_r(turning_angle_idx) = errorbar(N_list, mean_slope(:, turning_angle_idx), err_slope(:, turning_angle_idx), "^", 'Color','k', 'MarkerEdgeColor','k','MarkerFaceColor', ...
                            color_list(turning_angle_idx,:),'MarkerSize', 8,'LineWidth',1);
    end
    hold on
    [xData, yData] = prepareCurveData(N_list, mean_slope(:,turning_angle_idx) );
    % Set up fittype and options.
    ft = fittype( 'poly1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );
    x_fit = [N_list(1)-10:N_list(end)+10];
    y_fit = feval(fitresult, x_fit);
    plot(x_fit, y_fit,'color',color_list(turning_angle_idx,:), 'linewidth',2.5);
    hold on
    xlabel("$N$",'interpreter','latex')
    ylabel("$V_s$",'Interpreter','latex')
    set(gca, 'Fontname', 'helvetica', 'FontSize', 9)
end