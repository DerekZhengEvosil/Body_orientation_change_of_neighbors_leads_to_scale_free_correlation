%% The analysis of correlation length with the increasing robotic flock size
N_list = [20, 30, 50];
addpath("../Utility/")
load("../Data/Robots Data/SZ_cell.mat")
load("../Data/Robots Data/CL_cell.mat")
color_list = [[239 52 61]./255; [60 52 193]./255];
type_list = ["BOC","Random"];
cnt = 0;
drawResults = struct;
for type = type_list
    type_idx = find(type == type_list);
    for N = N_list
        N_idx = find(N == N_list);
        SZ = SZ_cell{N_idx, type_idx};
        CorL = CL_cell{N_idx, type_idx};
        length(SZ) 
        for i = 1:length(SZ)
            cnt = cnt + 1;
            drawResults.big_categ{cnt} = SZ(i);
            drawResults.result(cnt) = CorL(i);
            drawResults.colors{cnt} = num2str(type);
        end
    end 
end
g = gramm('x',drawResults.big_categ','y',drawResults.result', 'color', drawResults.colors');
g.geom_point()
g.set_point_options('base_size',5);
g.stat_glm()
g.axe_property('box', 'on');
g.set_text_options('interpreter', 'tex','base_size',9)
g.set_names('x', 'Flock Size (mm)', 'y', 'Correlation Length (mm)','color', " ",'label',{}) 
g.set_order_options('x',0);
figure;
figSize_L = 11;
figSize_W = 5;
set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
g.draw();