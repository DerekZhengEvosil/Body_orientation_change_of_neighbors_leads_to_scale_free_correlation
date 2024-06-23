clear, clc, close all; 
% The "proj.m" and "rot.m" are kindly provided in https://github.com/francelo/Visual-Servoing-IBVS-vs-PBVS
O = [0; 0; 0]; 
C = [0.0; 0.0; 0.0]; 
versor_origin = 2;
versor_camera = 2;
origin_axis = {'O';'X';'Y';'Z'};
camera_axis = {'oc','yc','xc','yc'};
plane_y = 1.5;      
plane_z = 1.5;    
f = 0.05;    
p = [1.0; 1.85; 0.8]; 
d = C - O;            
ang = [pi/2 0.0 pi/2];  
% create image plane
x = linspace(C(1)+f,C(1)+f,20);
y = linspace(C(2)-plane_y,C(2)+plane_y,20);
z = linspace(C(3)-plane_z,C(3)+plane_z,20);
[X, Y] = meshgrid(x,y);
Z = meshgrid(z);
%% projection of an object
% plot object
rotate_list = [0, 30, 60, 90];
a = 1; 
b = 0.53;
c = 0.5; 
cnt = 1;
for rotate_angle = rotate_list
    figure;
    figSize_L = 10;
    figSize_W = 10;
    set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
    rows = 1;
    cols = 2;
    width = 1 / cols * 0.8; 
    height = 1 / rows * 0.8; 
    my_color = [254,178,76
                253,187,132
                252,141,89
                239,101,72
                215,48,31]./255;
    idx = find(rotate_angle ==  rotate_list);
    center = [5, 0, 0]; 
    rotation_angles = [0, 0, rotate_angle];
    theta = linspace(0, 2*pi, 20); 
    phi = linspace(0, pi, 20); 
    [theta, phi] = meshgrid(theta, phi);
    x = a * sin(phi) .* cos(theta);
    y = b * sin(phi) .* sin(theta);
    z = c * cos(phi);
   
    rotation_matrix_x = [1, 0, 0; 0, cosd(rotation_angles(1)), -sind(rotation_angles(1)); 0, sind(rotation_angles(1)), cosd(rotation_angles(1))];
    rotation_matrix_y = [cosd(rotation_angles(2)), 0, sind(rotation_angles(2)); 0, 1, 0; -sind(rotation_angles(2)), 0, cosd(rotation_angles(2))];
    rotation_matrix_z = [cosd(rotation_angles(3)), -sind(rotation_angles(3)), 0; sind(rotation_angles(3)), cosd(rotation_angles(3)), 0; 0, 0, 1];
    
    points_matrix = [x(:), y(:), z(:)];
    
    rotated_points_matrix = (rotation_matrix_z * rotation_matrix_y * rotation_matrix_x * points_matrix')';
   
    x_rotated = reshape(rotated_points_matrix(:, 1), size(x)) + center(1);
    y_rotated = reshape(rotated_points_matrix(:, 2), size(y)) + center(2);
    z_rotated = reshape(rotated_points_matrix(:, 3), size(z)) + center(3);

    [X_surf, Y_surf, Z_surf] = deal(x_rotated, y_rotated, z_rotated);
    points = [X_surf(:), Y_surf(:), Z_surf(:)];

    proj_points = zeros(length(points), 2);
    for i = 1:length(points)
        [proj_points(i,1), proj_points(i,2)] = proj(points(i, :)', ang, d, f);
    end

    grid on;
    quiver3([O(1);O(1);O(1)],[O(2);O(2);O(2)],[O(3);O(3);O(3)],[versor_origin;0;0],[0;versor_origin;0],[0;0;versor_origin])

    axis equal
    hold on;
    quiver3([C(1);C(1);C(1)],[C(2);C(2);C(2)],[C(3);C(3);C(3)],[versor_camera;0;0],[0;versor_camera;0],[0;0;versor_camera])
    link_points = find(points(:,1) < 5);
    for i = 1:length(points)
        plot3(linspace(C(1),points(i,1)), linspace(C(2),points(i,2)), linspace(C(3),points(i,3)),'-','Color', [my_color(idx,:),0.2])
        alpha(.5)
    end
    hold on;
    s = surf(X,Y,Z,'FaceAlpha',0.3);
    s.EdgeColor = 'none';
    s = surf(X_surf,Y_surf,Z_surf,'Facealpha',1,'FaceColor',my_color(idx,:));
    s.EdgeColor = [37,37,37]./255;
    xlabel("x")
    ylabel("y")
    zlabel("h")
    ylim([-1.5, 1.5])
    view(135, 30)
    hold on
    box on
    proj_points = proj_points.* 100;
    k = convhull(proj_points(:,1), proj_points(:,2));
    fill3(zeros(1,length(k)), proj_points(k,1), proj_points(k,2), my_color(idx,:), 'FaceAlpha', 1,'EdgeColor','none'); % 填充凸包
    proj_x = proj_points(:,1);
    proj_y = proj_points(:,2);
    max_proj_x = max(proj_x);
    min_proj_x = min(proj_x);
    max_proj_y = max(proj_y);
    min_proj_y = min(proj_y);
    box_center = nanmean(proj_points);
    x_left = box_center(2) - (max_proj_y - min_proj_y)/2;
    x_right = box_center(2) + (max_proj_y - min_proj_y)/2;
    z_top = box_center(1) + (max_proj_x - min_proj_x)/2;
    z_bottom = box_center(1) - (max_proj_x - min_proj_x)/2;

    points_rect = [z_top, x_left;
              z_top, x_right;
              z_bottom, x_right;
              z_bottom, x_left;
              z_top, x_left]; 

    plot3(zeros(1, size(points_rect, 1)), points_rect(:, 1), points_rect(:, 2), 'k-', 'LineWidth', 3,'Color',[1,0,0]);
    hold on;
    set(gca, 'Fontname', 'helvetica', 'FontSize', 15)
end