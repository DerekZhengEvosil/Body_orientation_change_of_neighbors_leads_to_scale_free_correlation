%% Empirical Analysis on the U-turn trajectory, corresponding to Figure1-(a-d) in the main text
clear;clc;
addpath("../Utility/")
addpath("../Data/")
addpath("../Data/Empirical Data/")
load("data_emipircal_analysis.mat")
%% Draw the Basic property of Collective U-turn
figure();
set(gcf,'position',[50 10 780 500]);
traj = axes('Posi',[0.1 0.1 0.3 0.8]);
cuv_op = axes('Posi',[0.5 0.7 0.37 0.20]);
cuv_sin = axes('Posi',[0.5 0.4 0.37 0.20]);
stem_turning_rank = axes('Posi',[0.5 0.1 0.37 0.20]);
axes(traj)
% caclulate op
op = (nanmean(U_turn_vel_x(1:end,:),2).^2 + ......
            nanmean(U_turn_vel_y(1:end,:),2).^2).^(0.5);
% turning rank
[sorted_turning_start, ascend_idx] = sort(turning_start);

start_pos = memory{traj_start}(1:2,:);
end_pos = memory{traj_end}(1:2,:);
before_uturn_pos_x = U_turn_traj(1:traj_start,1:3:end);
before_uturn_pos_y = U_turn_traj(1:traj_start,2:3:end);
after_uturn_pos_x = U_turn_traj(traj_end:end,1:3:end);
after_uturn_pos_y = U_turn_traj(traj_end:end,2:3:end);
U_turn_pos_x = U_turn_traj(traj_start:traj_end,1:3:end);
U_turn_pos_y = U_turn_traj(traj_start:traj_end,2:3:end);
select_U_turn_pos_x = U_turn_traj(hier_start:hier_end,1:3:end);
select_U_turn_pos_y = U_turn_traj(hier_start:hier_end,2:3:end);
line(U_turn_pos_x, U_turn_pos_y,'LineWidth', 3,'Color',[145, 145, 145]/255)
hold on 
for p = 1:Nc
    patch([select_U_turn_pos_x(:,p)' NaN], [select_U_turn_pos_y(:,p)' NaN], [[hier_start:hier_end]/(hier_end -hier_start)  NaN], ...
                                'EdgeColor','interp','MarkerFaceColor','flat','LineWidth',3)
    colormap(turbo)
    hold on
end
h = colorbar('FontSize',15);
t=get(h,'YTickLabel');
t=strcat(t,'s');
set(h,'YTickLabel',t);
set(get(h,'Title'),'string','[ T-\tau, T ]');
hold on 
text(start_pos(1,:),start_pos(2,:),num2str([1:Nc]'),'HorizontalAlignment', ...
    'center','FontSize',12,'FontWeight','normal','Color',[0,0,0])
hold on 
scatter(end_pos(1,:),end_pos(2,:),50,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]) 
hold on
box on
interwall_r = 250;
outwall_r = 360; 
rectangle('position',[0-outwall_r, 0-outwall_r, outwall_r*2, outwall_r*2], ...
    'curvature',[1,1], 'LineWidth',4, 'EdgeColor',[0 0 0]/255);
xlabel("x")
ylabel("y")
xlim([145.2241330966653,450.5744183248542])
ylim([-179.5661795468842,199.1628563950557])
set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
axes(cuv_op)
plot([1:size(op,1)] * cyctime,op,'LineWidth', 2,'Color','k')
hold on 
plot([hier_start:hier_end] * cyctime, op([hier_start:hier_end]),'LineWidth', 2,'Color','r')
hold on
arrow_start = [hier_start * cyctime, op(hier_start)];
arrow_end = [hier_start * cyctime, 0];
drawarrow(arrow_start,arrow_end)
hold on
arrow_start = [hier_end * cyctime, op(hier_end)];
arrow_end = [hier_end * cyctime, 0];
drawarrow(arrow_start,arrow_end)
xlabel("time(s)")
ylabel("\phi")
set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
axes(cuv_sin)
plot([1:size(theta_wi,2)] * cyctime, theta_wi','LineWidth',2)
xlabel("time(s)")
ylabel("sin(\theta_{w})")
set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
axes(stem_turning_rank)
stem(sorted_turning_start * cyctime,':diamondr','filled','Color',[241 41 79]/255,'LineWidth',2)
xticks([1:Nc])
xticklabels(cellstr(string(ascend_idx)))
xlabel("Fish ID sorted by ascending cost time")
ylabel("time(s)")
set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
%% Draw the leader-follower nwtwork 
figure
figSize_L = 12;
figSize_W = 10;
set(gcf,'Units','centimeter','Position',[5 5 figSize_L figSize_W]);
C_ij_mat = gen_period_C_ij(U_turn_vel_x, U_turn_vel_y, hier_end, hier_end - hier_start, cyctime);
C_ij_mat = tril(C_ij_mat);
s = [];
t = [];
weights = [];
for k = 1:size(C_ij_mat,1)
    i2j = C_ij_mat(:,k);
    % < 0 i follow j / > 0 j follow i
    tau_ij_lz = find(i2j < 0);
    s_tmp = k * ones(1,length(tau_ij_lz));
    t_tmp = tau_ij_lz';
    s = [s,s_tmp];
    t = [t,t_tmp];
    weights = [weights,abs(i2j(tau_ij_lz))'];
    tau_ij_mz = find(i2j > 0);
    s_tmp = tau_ij_mz';
    t_tmp = k * ones(1,length(tau_ij_mz));
    s = [s,s_tmp];
    t = [t,t_tmp];
    weights = [weights,abs(i2j(tau_ij_mz))'];
end
G = digraph(t,s,weights); % For each pairwise comparison, the directed edge points from the leader to the follower
if size(G.Nodes,1) ~= size(C_ij_mat,1)
    G = addnode(G,size(C_ij_mat,1) - size(G.Nodes,1));
end
[~, order] = sort(G.outdegree,'descend');
G = reordernodes(G,order);
LWidths = 5 * G.Edges.Weight/max(G.Edges.Weight);
p = plot(G,'Layout','layered','LineWidth',LWidths,'NodeLabel',{});
p.Marker = 'o';
p.MarkerSize = 30;
unique_outdegree = unique(G.outdegree);
use_hierColor = hierColor.GMTred2green(length(unique_outdegree):-1:1,:);
for order_id = 1:size(G.Nodes,1)
    color_idx = find(G.outdegree(order_id) == unique_outdegree);
    highlight(p, order_id, 'NodeColor', use_hierColor(color_idx,:));
end

p.EdgeColor = [135 135 135]/255;
p.ArrowSize = 10;
text(p.XData,p.YData,num2str(order),'HorizontalAlignment','center','FontSize',15,'FontWeight','bold','Color',[0,0,0],'FontName','helvetica')
L_t = G.outdegree / (Nc - 1);
[~, sort_order_id] = sort(order,'ascend');
L_t_rank_by_id = L_t(sort_order_id);   
M_s = nanmean(M_ij_tmp,1);
[rho, ~] = corr(L_t_rank_by_id,M_s','type','Spearman');
figure
figSize_L = 10;
figSize_W = 10;
set(gcf,'Units','centimeter','Position',[5 5 figSize_L figSize_W]);
[ha,pos] = tight_subplot(Nc,Nc); % ,gap,marg_h,marg_w
trick_mat = eye(10);
cnt = 0;
draw_b = 0.5;
a=1;
unit_M_ij_tmp = M_ij_tmp./max(M_ij_tmp);
for tmp1 = 1:Nc
    for tmp2 = 1:Nc
        cnt = cnt + 1;
        axes(ha(cnt))
        xlim([-1 1])
        ylim([-1 1])
        if trick_mat(tmp1, tmp2)
            axis off
            continue
        else
            ellipse_tmp = ellipsedraw(a/2, draw_b/2, 0, 0, pi/2, '-');
            set(ellipse_tmp,'Color',[1,1,1])
            patch(ellipse_tmp.XData, ellipse_tmp.YData,[11 92 175]/255 , 'FaceAlpha', unit_M_ij_tmp(tmp1, tmp2),...
                                                    'EdgeColor',[1,1,1])
            axis off  %  [244, 149, 10]/255
        end
    end
end 
%% Draw the Correlation between the BOC and leadership
figure
figSize_L = 12;
figSize_W = 4;
set(gcf,'Units','centimeter','Position',[5 5 figSize_L figSize_W]);
[ha,pos] = tight_subplot(1,Nc); % ,gap,marg_h,marg_w
size_adjust = 1;
theta=0:pi/100:2*pi;
expand_size = 1;
for pp = 1:Nc
    axes(ha(pp))
    rectangle('position',[0-(L_t_rank_by_id(pp)*size_adjust), 0-(L_t_rank_by_id(pp)*size_adjust),  ...
                (L_t_rank_by_id(pp))*2*size_adjust, (L_t_rank_by_id(pp))*2*size_adjust],'curvature', ...
                    [1,1],'FaceColor',[241 41 79]/255,'EdgeColor',[1,1,1]);
    xlim([-expand_size expand_size])
    ylim([-expand_size expand_size])

    hold on
    text(0,0,num2str(pp),'HorizontalAlignment', ...
    'center','FontSize',12,'FontWeight','normal','Color',[0,0,0])

    text(0,-3, num2str(round(L_t_rank_by_id(pp), 2)), 'HorizontalAlignment', ...
    'center','FontSize',12,'FontWeight','normal','Color',[0,0,0])
    axis equal
    axis off
end

figure
figSize_L = 12;
figSize_W = 4;
set(gcf,'Units','centimeter','Position',[5 5 figSize_L figSize_W]);
[ha,pos] = tight_subplot(1,Nc); % ,gap,marg_h,marg_w
size_adjust = 1;
theta=0:pi/100:2*pi;
expand_size = 1;
for pp = 1:Nc
    axes(ha(pp))
    rectangle('position',[0-(M_s(pp)./max(M_s)*size_adjust), 0-(M_s(pp)./max(M_s)*size_adjust),  ...
                (M_s(pp)./max(M_s))*2*size_adjust, (M_s(pp)./max(M_s))*2*size_adjust],'curvature', ...
                    [1,1],'FaceColor',[37 66 116]/255,'EdgeColor',[1,1,1]);
    xlim([-expand_size expand_size])
    ylim([-expand_size expand_size])

    hold on
    text(0,0,num2str(pp),'HorizontalAlignment', ...
    'center','FontSize',12,'FontWeight','normal','Color',[0,0,0])

    text(0,-3, num2str(round(M_s(pp)./max(M_s), 2)), 'HorizontalAlignment', ...
    'center','FontSize',12,'FontWeight','normal','Color',[0,0,0])
    axis equal
    axis off
end
figure
figSize_L = 12;
figSize_W = 6;
set(gcf,'Units','centimeter','Position',[5 5 figSize_L figSize_W]);
[L_t_ascend, ascend_idx] = sort(L_t_rank_by_id, 'ascend');
stem(1:length(L_t_ascend), L_t_ascend,'filled','Color',[241 41 79]/255)
xticks([1:Nc])
xticklabels(cellstr(string(ascend_idx)))
xlabel("Fish ID sorted by ascending ")
ylabel("$L_i(T, \tau) or G_i(T,\tau)$", 'interpreter','latex')
hold on
stem(1:length(L_t_ascend),M_s(ascend_idx)./max(M_s),"filled","Color",[11 92 175]/255)
set(gca, 'Fontname', 'helvetica', 'FontSize', 15)