%% acc on bar
clear;clc
type_list = ["projection"];
load("../Data/Simulation Data/collective turn/Random_acc_all.mat")
for type_tmp = type_list
    type_id = find(type_tmp == type_list);
    N_list = [10 20 30 40 50 60 70 80 90 100];
    turning_angle_list = ["\pi/2", "\pi"];
    drawResults = struct;
    cnt = 0;
    for turning_angle_tmp = turning_angle_list
        turning_idx = find(turning_angle_tmp == turning_angle_list);
        acc_N_tmp = acc_all_mt_type{turning_idx};
        for N_tmp = N_list
            N_idx = find(N_tmp == N_list);
            r_tmp = acc_N_tmp{N_idx}; 
            for i = 1:length(r_tmp) % 1:length(r_tmp)
                cnt = cnt + 1;
                drawResults.big_categ{cnt} = num2str(N_tmp);
                drawResults.result(cnt) = r_tmp(i);
                drawResults.colors{cnt} = num2str(turning_angle_tmp);
            end 
        end
    end
    g = gramm('x',drawResults.big_categ','y',drawResults.result', 'color', drawResults.colors');
    g.stat_summary('geom',{'bar','black_errorbar'},'setylim','true');
    % g.stat_violin('normalization','width');
    g.set_text_options('interpreter', 'tex','base_size',15)
    g.axe_property('box', 'on'); % 将'box'参数设置为'on'，表示显示框边线
    g.set_names('x', 'N', 'y', 'accuracy','color', "Turning angle",'label',{}) 
    g.set_title('');
    g.axe_property('YLim', [0, 1]);
    g.set_order_options('x',0);
    figure;
    figSize_L = 16;
    figSize_W = 6;
    set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
    g.draw();
end
%% responsivity on bar
clear;clc
type_list = ["projection"];
load("../Data/Simulation Data/collective turn/Random_resp_all.mat")
for type_tmp = type_list
    type_id = find(type_tmp == type_list);
    N_list = [10 20 30 40 50 60 70 80 90 100];
    turning_angle_list = ["\pi/2", "\pi"];
    drawResults = struct;
    cnt = 0;
    for turning_angle_tmp = turning_angle_list
        turning_idx = find(turning_angle_tmp == turning_angle_list);
        r_all_N_tmp = r_all_mt_type{turning_idx};
        for N_tmp = N_list
            N_idx = find(N_tmp == N_list);
            r_tmp = r_all_N_tmp{N_idx}; 
            for i = 1:length(r_tmp) %1:length(r_tmp)
                cnt = cnt + 1;
                drawResults.big_categ{cnt} = num2str(N_tmp);
                drawResults.result(cnt) = r_tmp(i);
                drawResults.colors{cnt} = num2str(turning_angle_tmp);
            end 
        end
    end
    g = gramm('x',drawResults.big_categ','y',drawResults.result', 'color', drawResults.colors');
    g.stat_summary('geom',{'bar','black_errorbar'},'setylim','true');
    g.axe_property('box', 'on'); % 将'box'参数设置为'on'，表示显示框边线
    % g.stat_violin('normalization','width');
    g.set_text_options('interpreter', 'tex','base_size',15)
    g.set_names('x', 'N', 'y', 'R','color', "Turning angle",'label',{}) 
    % g.set_title('');
    g.axe_property('YLim', [0, 2]); % 将min_value和max_value替换为所需的y轴范围值
    g.set_order_options('x',0);
    figure;
    figSize_L = 16;
    figSize_W = 6;
    set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
    g.draw();
end