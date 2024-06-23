%% Correlation function as a function of distance r
clear;clc;
N_list = [30, 50, 75, 100, 150, 200];
figure;
figSize_L = 5.5;
figSize_W = 5.5;
set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
my_color = [254,227,145
254,196,79
254,153,41
236,112,20
204,76,2
153,52,4]./255;
cs_data = "../Data/Simulation Data/collective spin/Random_cs_CL_Cr_data.mat";
load(cs_data)
r = collective_spin_CL_Cr_data{3};
Cr = collective_spin_CL_Cr_data{4};
for N = N_list
    N_idx = find(N == N_list);
    mean_r = nanmean(r{N_idx});
    mean_Cr = nanmean(Cr{N_idx});
    plot(mean_r, mean_Cr, 'o', 'Color',my_color(N_idx,:),'MarkerEdgeColor','none','MarkerFaceColor',my_color(N_idx,:),'MarkerSize',3)
    hold on
end
yline(0, '-','linewidth',1);
box on
xlabel("$r$",'Interpreter','latex')
ylabel("Correlation Function")
set(gca, 'Fontname', 'helvetica', 'FontSize', 9)
legend("N = " + N_list','box','off','location','best')
%% Correlation length increases with distance r
CL_cell = collective_spin_CL_Cr_data{1};
GroupSZ_cell = collective_spin_CL_Cr_data{2};
figure
figSize_L = 24;
figSize_W = 8;
set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
for i = 1:length(N_list)
    SZ(i) = nanmean(GroupSZ_cell{i});
    CorL(i) = nanmean(CL_cell{i});
    yneg(i) = nanstd(CL_cell{i});
    ypos(i) = nanstd(CL_cell{i});
    xneg(i) = nanstd(GroupSZ_cell{i});
    xpos(i) = nanstd(GroupSZ_cell{i});
end
errorbar(SZ, CorL, yneg, ypos, xneg, xpos, "o", "color",[238 0 9]./255, "MarkerSize", 10,...
    "MarkerEdgeColor", [238 0 9]./255, "MarkerFaceColor", [238 0 9]./255,'LineWidth',1.5)
hold on
[xData, yData] = prepareCurveData( SZ, CorL );
% Set up fittype and options.
ft = fittype( 'poly1' );
% Fit model to data.
[fitresult, gof] = fit(xData, yData, ft);
x_fit = [xData(1)-5:xData(end)+10];
y_fit = feval(fitresult, x_fit);
h = plot(x_fit, y_fit,'color','k', 'linewidth',2.5);
xlabel("Flock Size")
ylabel("Correlation Length")
set(gca, 'Fontname', 'helvetica', 'FontSize', 15)