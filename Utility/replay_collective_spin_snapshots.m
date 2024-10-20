function [heading, op, transfer_Ang] = replay_collective_spin_snapshots(folder_name,plot_figure,activateTime)
txtFiles = dir(folder_name);
G = struct;
G.cycTime = 0.2;
G.BL = 0.06; % meter
for i = 2:length(txtFiles)
    param = split(txtFiles(i).name, '_');
    if param{1} == "simData"
        robotId = str2double(erase(param{2}, ".txt")) + 1;
        G.actor{robotId}.memory = load([folder_name + '/' + txtFiles(i).name]);
    end
end
G.num = length(G.actor);
for i = 1:G.num
    len = size(G.actor{i}.memory,1);
    savingData(1:len,i,1:2)=G.actor{i}.memory(:,1:2)./ G.BL;
    savingData(1:len,i,3) = cos(G.actor{i}.memory(:,3) - pi/2);
    savingData(1:len,i,4) = sin(G.actor{i}.memory(:,3) - pi/2);
end
cyctime = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G.expNum = G.num;
G.total_step = size(savingData,1);
chk = 0;
all_pos = [];
all_vel = [];
op = [];
nnd = [];

activate_time = activateTime;
informed_id_file = load([folder_name + '/' + "informed_id.txt"]);
informed_id = informed_id_file(:,1);
informed_id = informed_id + 1; % python 中 index 从0开始
start_time_file = load([folder_name + '/' + "start_time.txt"]);
start_time = (start_time_file(1,:) - activate_time) * cyctime;

tracks = [];
for t = 1:size(savingData,1)
    G.robotsPosH = squeeze(savingData(t,:,1:4));
    all_pos = G.robotsPosH(1:G.expNum, [1:2])';
    all_vel = G.robotsPosH(1:G.expNum, [3:4])';
    op(t,:) = (nanmean(all_vel(1,:))^2 + nanmean(all_vel(2,:))^2)^(0.5);
    dist_xy = squareform((pdist(all_pos(:,[1:G.num])','euclidean'))); 
    dist_xy(logical(eye(G.num))) = NaN;
    nnd(t,:) = min(min(dist_xy));
    heading(:,t) = atan2(all_vel(2,:), all_vel(1,:)); 
end
changed_angle = 0;
informed_stop_time = 0;
for t = activate_time:size(savingData,1)-1
    cur_robotsPosH = squeeze(savingData(t,:,1:4));
    next_robotsPosH = squeeze(savingData(t+1,:,1:4));
    cur_all_vel = cur_robotsPosH(1:G.expNum, [3:4])';
    next_all_vel = next_robotsPosH(1:G.expNum, [3:4])';
    cur_infomred_vel = cur_all_vel(:, informed_id);
    next_informed_vel = next_all_vel(:, informed_id);
    changed_angle = changed_angle + (acos(dot(cur_infomred_vel, next_informed_vel)));
    if changed_angle > 2 * pi
        informed_stop_time = t;
        break;
    end 
end

for id = 1:G.expNum
    one_frame = [];
    cnt = 0;
    for t = activate_time:informed_stop_time
        cnt = cnt + 1;
        G.robotsPosH = squeeze(savingData(t,:,1:4));
        all_pos = G.robotsPosH(1:G.expNum, [1:2])';
        all_vel = G.robotsPosH(1:G.expNum, [3:4])';
        one_frame(cnt,:) = [id, all_pos(1,id), all_pos(2,id), t*cyctime, all_vel(1,id), all_vel(2,id)];
    end
    tracks = [tracks ;one_frame];
end

robotsPosH = squeeze(savingData(1,:,1:4));
all_pos = G.robotsPosH(1:G.expNum, [1:2])';
informed_pos = all_pos(:,informed_id);
transfer_distance = (sum((all_pos - repmat(informed_pos, 1, G.expNum)).^2, 1)).^0.5;  % unit is Body length
[start_time_new, I]=sort(start_time);
for i=1:length(I)
    Rank(I(i))=i; % = 1 is the first bird start to turn
end

T=unique(tracks(:,4));
[ascend_delay, delay_idx] = sort(start_time,'ascend');
ascend_trans_dis = transfer_distance(delay_idx);

d_i = zeros(1,G.expNum);
for k = 1:length(Rank)
    if Rank(k) == 1
        d_i(k) = 0;
    else
        ID_tmp = find(Rank <= Rank(k));
        first_rank_pos = all_pos(:, ID_tmp);
        radius_tmp = max(max(squareform((pdist(first_rank_pos(:,[1:length(ID_tmp)])','euclidean')))));
        d_i(k) = (Rank(k) .* (pi * (radius_tmp/2)^2) ./ G.expNum).^0.5;
    end
end
d_i_asc = d_i(delay_idx);
% scatter(ascend_delay, ascend_trans_dis)
% scatter(ascend_delay, d_i(delay_idx))

id = find(tracks(:,4) == T(1)); % time frame when turn start
xyz=tracks(id,2:3);
ID1 = find(Rank < 0.2 * G.expNum);
ID2 = find(Rank > 0.8 * G.expNum);
vs = mean(xyz(ID1,:),1)-mean(xyz(ID2,:),1); % opposit to propagate direction
vs = vs/sum(vs.^2)^0.5;
group_vel = [cos(nanmean(heading(:, activate_time))) sin(nanmean(heading(:, activate_time)))];
transfer_Ang = acos(sum(group_vel.*vs));


start_heading = heading(:, 1);
stop_time_heading = heading(:, informed_stop_time);
rotate_angle = zeros(1,G.expNum);
for t = activate_time:informed_stop_time + 1
    cur_robotsPosH = squeeze(savingData(t,:,1:4));
    past_robotsPosH = squeeze(savingData(t - 1,:,1:4));
    cur_all_vel = cur_robotsPosH(1:G.expNum, [3:4])';
    past_all_vel = past_robotsPosH(1:G.expNum, [3:4])';
    rotate_angle = rotate_angle + rad2deg(acos(dot(cur_all_vel, past_all_vel)));
end
err_heading = abs(stop_time_heading - start_heading);
acted_id = find((rotate_angle) > 20); % 
cmap = crameri('roma');
if plot_figure == 1
    snap_list = [1, length(T) - 20, length(T) - 10, length(T)];
    for snap_idx = snap_list
        figure
        figSize_L = 12;
        figSize_W = 10;
        set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
        id=find(tracks(:,4)==T(snap_idx)); 
        xyz=tracks(id,2:3);
        u=tracks(id,5:6);
        u=u./sqrt(sum(u.*u,2));
        phi = atan2(u(:,2), u(:,1));
        a = 4; 
        b = 1.2;
        colormap(cmap);
        start_time_normalized = (start_time - min(start_time)) / (max(start_time) - min(start_time)); % 将start_time归一化到 [0, 1] 范围内
        colormap_index = round(start_time_normalized * (size(cmap, 1) - 1)) + 1; % 计算归一化后的索引
        mapped_colors = zeros(G.expNum, 3);
        none_move_id = setdiff([1:G.expNum], acted_id);
        for c_i = 1:G.expNum
            mapped_colors(c_i, :) = interp1(1:size(cmap, 1), cmap, colormap_index(c_i));
        end
        start_time(start_time<0) = 0;
        start_time(start_time > informed_stop_time) = informed_stop_time;
        hs = scatter(xyz(:,1),xyz(:,2),300,start_time,'o','filled', 'LineWidth',1.5,'MarkerFaceAlpha', 0);
        scatter_colors = get(hs, 'CData');
        cnt = 0;
        ellipse_X = [];
        ellipse_Y = [];
        for p_idx = 1:length(acted_id)
            p = acted_id(p_idx);
            cnt = cnt + 1;
            [ellipse_X(:,cnt), ellipse_Y(:,cnt)] = ellipsepatch(a/2, b/2, xyz(p,1), xyz(p,2), phi(p));
        end
        patch(ellipse_X, ellipse_Y, scatter_colors(acted_id))
        hold on
        cnt = 0;
        ellipse_X = [];
        ellipse_Y = [];
        for p_idx = 1:length(none_move_id)
            p = none_move_id(p_idx);
            cnt = cnt + 1;
            [ellipse_X(:,cnt), ellipse_Y(:,cnt)] = ellipsepatch(a/2, b/2, xyz(p,1), xyz(p,2), phi(p));
        end
        patch(ellipse_X, ellipse_Y, [127, 127, 127]/255)
        hold on
        [info_X, info_Y] = ellipsepatch(a/2, b/2, xyz(informed_id,1), xyz(informed_id,2), phi(informed_id));
        patch(info_X, info_Y, 'r')
        hold off
        grid off;
        box on
        axis equal
        colormap(cmap);
        grid off;box on
        color_h=colorbar;
        set(get(color_h,'Title'),'string','spinning lag(step)','interpreter','latex');
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
%         title("t = " + num2str(snap_idx + activate_time))
        set(gca, 'Fontname', 'helvetica', 'FontSize', 10)
        xlabel("$x$",'interpreter','latex')
        ylabel("$y$",'interpreter','latex')
        set(gca,'linewidth',3)
    end

    figure;
    figSize_L = 24;
    figSize_W = 8;
    set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
    draw_sec = 1;
    heading = rad2deg(heading);
    group_heading = heading;
    group_heading(informed_id,:) = [];
    mean_group_heading = nanmean(group_heading,1);
    plot([1:draw_sec:informed_stop_time], group_heading(:,1:draw_sec:informed_stop_time), "LineWidth",2)
    hold on
    plot([1:draw_sec:informed_stop_time], heading(informed_id,1:draw_sec:informed_stop_time),'Color','k', "LineWidth",4)
    hold on
    ylabel("heading")
    xlabel("time(step)")
    set(gca, 'Fontname', 'helvetica', 'FontSize', 15)

    figure;
    figSize_L = 24;
    figSize_W = 8;
    set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
    draw_sec = 1;
    plot([1:draw_sec:informed_stop_time]*cyctime, op(1:draw_sec:informed_stop_time), "LineWidth",4)
    ylabel("\phi")
    xlabel("time(step)")
    set(gca, 'Fontname', 'helvetica', 'FontSize', 15)

end
end
