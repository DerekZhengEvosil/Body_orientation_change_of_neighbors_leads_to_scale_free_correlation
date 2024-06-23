load("./example_data.mat")
[r, C, GroupSZ] = CorrFunction(all_pos',all_vel');
x = unique(r, 'stable');
y = unique(C, 'stable');
nan_idx = find(isnan(x));
x(nan_idx) = [];
y(nan_idx) = [];
interp_method = 'spline';
interp_x = linspace(min(x), max(x), 1000);
interp_y = interp1(x, y, interp_x, interp_method);
cross_x = interp_x(find(interp_y < 0, 1));
cross_y = 0;  
CL = cross_x;
fprintf('Correlation_length = %.2f, Group Size = %.2f \n', cross_x, GroupSZ);
plot(x, y, 'o', 'DisplayName', 'data points'); 
hold on;
plot(interp_x, interp_y, 'r');  
hold on
plot(cross_x, cross_y, 'ro', 'MarkerSize', 10, 'DisplayName', 'cross_point'); 
hold on
yline(0,'--','lineWidth',1)