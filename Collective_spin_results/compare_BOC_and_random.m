clear;clc;
N_list = [30, 50, 75, 100, 150, 200];
rand_cs_data = "../Data/Simulation Data/collective spin/Random_cs_CL_Cr_data.mat";
boc_cs_data = "../Data/Simulation Data/collective spin/BOC_cs_CL_Cr_data.mat";
load(boc_cs_data)
CL_cell = collective_spin_CL_Cr_data{1};
GroupSZ_cell = collective_spin_CL_Cr_data{2};
figure
figSize_L = 10;
figSize_W = 4.5;
set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
for i = 1:length(N_list)
    SZ(i) = nanmean(GroupSZ_cell{i});
    CorL(i) = nanmean(CL_cell{i});
    yneg(i) = nanstd(CL_cell{i});
    ypos(i) = nanstd(CL_cell{i});
    xneg(i) = nanstd(GroupSZ_cell{i});
    xpos(i) = nanstd(GroupSZ_cell{i});
end

e_b = errorbar(SZ, CorL, yneg, ypos, xneg, xpos, "o", "color",[213,62,79]./255, "MarkerSize", 6,...
    "MarkerEdgeColor", [213,62,79]./255, "MarkerFaceColor", [213,62,79]./255,'LineWidth',1);
hold on
[xData, yData] = prepareCurveData( SZ, CorL );
% Set up fittype and options.
ft = fittype( 'poly1' );
% Fit model to data.
[fitresult, gof] = fit(xData, yData, ft);
x_fit = [xData(1)-5:xData(end)+10];
y_fit = feval(fitresult, x_fit);
h = plot(x_fit, y_fit,'color',[213,62,79]./255, 'linewidth',1.5);
hold on
load(rand_cs_data)
CL_cell = collective_spin_CL_Cr_data{1};
GroupSZ_cell = collective_spin_CL_Cr_data{2};
for i = 1:length(N_list)
    SZ(i) = nanmean(GroupSZ_cell{i});
    CorL(i) = nanmean(CL_cell{i});
    yneg(i) = std(CL_cell{i});
    ypos(i) = std(CL_cell{i});
    xneg(i) = std(GroupSZ_cell{i});
    xpos(i) = std(GroupSZ_cell{i});
end
hold on
e_r = errorbar(SZ, CorL, yneg, ypos, xneg, xpos, "o", "color",[50,136,189]./255, "MarkerSize", 6,...
    "MarkerEdgeColor", [50,136,189]./255, "MarkerFaceColor", [50,136,189]./255,'LineWidth',1.5);
hold on
[xData, yData] = prepareCurveData( SZ, CorL );
% Set up fittype and options.
ft = fittype( 'poly1' );
% Fit model to data.
[fitresult, gof] = fit(xData, yData, ft);
x_fit = [xData(1)-5:xData(end)+10];
y_fit = feval(fitresult, x_fit);
h = plot(x_fit, y_fit,'color',[50,136,189]./255, 'linewidth',1.5);
xlim([1500, 4800])
xlabel("Flock Size",'Interpreter','latex')
ylabel("Correlation Length", 'Interpreter','latex')
legend([e_b, e_r], ["BOC","Random"],'box','off')
set(gca, 'Fontname', 'helvetica', 'FontSize', 9)
%% 
clear;clc;
N_list = [30, 50, 75, 100, 150, 200];
rand_Vs_data = "../Data/Simulation Data/collective spin/Random_Vs_info.mat";
boc_Vs_data = "../Data/Simulation Data/collective spin/BOC_Vs_info.mat";
load(boc_Vs_data)
figure
figSize_L = 10;
figSize_W = 4.5;
set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
mean_slope = cellfun(@nanmean, Vs_info);
std_slope = cellfun(@nanstd, Vs_info);
e_b = errorbar(N_list, mean_slope, std_slope, "^", 'Color', 'k', 'MarkerEdgeColor', 'k','MarkerFaceColor', ...
                        [213,62,79]/255, 'MarkerSize', 10,'LineWidth',1);
hold on
[xData, yData] = prepareCurveData(N_list, mean_slope);
% Set up fittype and options.
ft = fittype( 'power2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [2.76284081529424 0.482165172368064 -0.358059417250643];
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
x_fit = [N_list(1)-10:N_list(end)+50];
y_fit = feval(fitresult, x_fit);
plot(x_fit, y_fit,'color',[213,62,79]./255, 'linewidth',1.5);
hold on
load(rand_Vs_data)
mean_slope = cellfun(@nanmean, Vs_info);
std_slope = cellfun(@nanstd, Vs_info);
e_r = errorbar(N_list, mean_slope, std_slope, "^", 'Color', 'k', 'MarkerEdgeColor', 'k','MarkerFaceColor', ...
                        [50,136,189]/255, 'MarkerSize', 10,'LineWidth',1);
hold on
[xData, yData] = prepareCurveData(N_list, mean_slope);
% Set up fittype and options.
ft = fittype( 'power2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [2.76284081529424 0.482165172368064 -0.358059417250643];
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
x_fit = [N_list(1)-10:N_list(end)+50];
y_fit = feval(fitresult, x_fit);
plot(x_fit, y_fit,'color',[50,136,189]./255, 'linewidth',1.5);
xlabel("$N$",'interpreter','latex')
ylabel("$V_s$",'Interpreter','latex')
legend([e_b, e_r],["BOC", "Random"],'Interpreter','latex','box','off') %  = 0.1864 * N^{0.7347} + 19.29
set(gca, 'Fontname', 'helvetica', 'FontSize', 9)


