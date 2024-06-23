%% The analysis of information transfer speed with the increasing group size
addpath("../Utility/")
load('../Data/Robots Data/information_trans_data.mat');
figure
figSize_L = 10;
figSize_W = 5;
set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
N = 50;
Type_list = ["BOC", "Rand"];
slope_cell = {};
slope = [];
mean_delay = information_trans_data{1};
mean_trans_dis =  information_trans_data{2};
std_delay = information_trans_data{3};
std_trans_dis = information_trans_data{4};
for Type = Type_list
    Type_idx = find(Type == Type_list);
    for idx = 1:2
        x = mean_delay{idx, Type_idx};
        y = mean_trans_dis{idx, Type_idx};
        std_x = std_delay{idx, Type_idx};
        std_y =std_trans_dis{idx, Type_idx};
        linear_region = find(x>5&x<15);
        linear_start = linear_region(1);
        linear_end = linear_region(end);
        linearModel = fit(x(linear_start:linear_end), y(linear_start:linear_end), 'poly1'); 
        coefficients = coeffvalues(linearModel);
        slope(Type_idx, idx) = coefficients(1);
        intercept = coefficients(2);
        x_fit = linspace(min(x), max(x), 100); 
        y_fit = polyval(coefficients, x_fit);    
        if Type_idx == 1
            if idx == 1
                plot(x,y,'ko', 'LineWidth',2,'MarkerEdgeColor',[11,92,175]./255, ...
                                                            'MarkerFaceColor',[11,92,175]./255,'MarkerSize',5)
            else
                plot(x,y,'k^', 'LineWidth',2,'MarkerEdgeColor',[11,92,175]./255, ...
                                                            'MarkerFaceColor',[11,92,175]./255,'MarkerSize',5)
            end
        else
            if idx == 1
                plot(x,y,'ko', 'LineWidth',2,'MarkerEdgeColor',[205,62,21]./255, ...
                                                            'MarkerFaceColor',[205,62,21]./255,'MarkerSize',5)
            else
                plot(x,y,'k^', 'LineWidth',2,'MarkerEdgeColor',[205,62,21]./255, ...
                                                            'MarkerFaceColor',[205,62,21]./255,'MarkerSize',5)
            end
        end
        hold on; 
    end
end
xlabel("turning delay (s)")
ylabel("$d_{trans}$(mm)",'Interpreter','latex')
legend_list_1 = "slope = " + round(slope(1,:)',1);
legend_list_2 = "slope = " + round(slope(2,:)',1);
legend([legend_list_1;legend_list_2],'Box','off')