clear, clc, close all; 
% The "proj.m" and "rot.m" are kindly provided in https://github.com/francelo/Visual-Servoing-IBVS-vs-PBVS
O = [0; 0; 0];        
C = [0.0; 0.0; 0.0];  
versor_origin = 0.4;
versor_camera = 0.2;
origin_axis = {'O';'X';'Y';'Z'};
camera_axis = {'oc','zc','xc','yc'};
plane_y = 1;      
plane_z = 1;       
f = 0.05;             
p = [1.0; 1.85; 0.8]; 
d = C - O;            
ang = [pi/2 0.0 pi/2];  
%% projection of an object
rotate_list = [0:0.5:90]; 
my_colormap = jet(length(rotate_list));
a = 1; 
b = 0.53; 
c = 0.5; 
cnt = 1;
for rotate_angle = rotate_list
    idx = find(rotate_angle ==  rotate_list);
    center = [5, 0, 0]; 
    rotation_angles = [0, 0, rotate_angle];
    theta = linspace(0, 2*pi, 10); 
    phi = linspace(0, pi, 10);
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
    hold on
    projection = proj_points;
    hold on
    proj_x = proj_points(:,1);
    proj_y = proj_points(:,2);
    max_proj_x = max(proj_x);
    min_proj_x = min(proj_x);
    max_proj_y = max(proj_y);
    min_proj_y = min(proj_y);
    box_area(idx) = ((max_proj_x - min_proj_x) * (max_proj_y - min_proj_y));
end
box_area = box_area * 10000;
figure
figSize_L = 16;
figSize_W = 30;
set(gcf, 'Units', 'centimeter','Position', [5 5 figSize_L figSize_W])
drawSec = 1;
plot(rotate_list([1:drawSec:end, length(box_area)]), box_area([1:drawSec:end, length(box_area)]), ...
                    '-','MarkerFaceColor',[69, 126, 180]./255,'MarkerEdgeColor', ...
                [69, 126, 180]./255,'LineWidth',3,'Color',[69, 126, 180]./255,'MarkerSize',10)
zeros_idx = find(rotate_list == 0);
thirty_idx = find(rotate_list == 30);
sixty_idx = find(rotate_list == 60);
ninty_idx = find(rotate_list == 90);
hold on
plot(rotate_list([zeros_idx, thirty_idx, sixty_idx, ninty_idx]), ...
            box_area([zeros_idx, thirty_idx, sixty_idx, ninty_idx]),'o','Color', [228,26,28]./255,'LineWidth',2,'MarkerSize',20)
hold on
plot(rotate_list([zeros_idx, thirty_idx, sixty_idx, ninty_idx]), ...
            box_area([zeros_idx, thirty_idx, sixty_idx, ninty_idx]),'o', ...
            'Color', [228,26,28]./255,'LineWidth',2,'MarkerSize',10, 'MarkerFaceColor',[228,26,28]./255)
xlabel("orientation respect to axis x")
ylabel("Area of projected image")
xlim([0,90])
grid on
box on
set(gca, 'Fontname', 'helvetica', 'FontSize', 25)