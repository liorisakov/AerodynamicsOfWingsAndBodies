% ---------------------------------- %
%  AERODYNAMICS OF WINGS AND BODIES  %
%             HOMEWORK 2             %
%  Lior Isakov          Eden Shazar  %
%   200928372            204378608   %
% ---------------------------------- %

%% Fresh start
clc
close all
clear variables

%% Parameters
% user input
airfoil = NACA('2412');
N = 50;                % number of vortex line vortices
U_oo = 100;            % [m/s] - uniform flow speed
rho = 1.225;           % [kg/m^3] - air density (chosen as ISA at sea level)
alpha = deg2rad(8);    % [rad] - angle of attack

% calculated parameters
x = linspace(0, 1, N+2);    % vortex line sectioning; extra two points as
                            % LE and TE points don't have vortices


%% Airfoil
% vortex line definition
vortex.t = x(2:end-1);                     % vortex points
vortex.p = ((vortex.t + x(3:end))/2).';    % control points

% camber line
dz_dx = airfoil.camber_slope(vortex.p).';    % slope
zi = airfoil.camber_line(vortex.p);          % z coordinate of control points
zj = airfoil.camber_line(vortex.t);          % z coordinate of vortex points

%% Solver
% parameters
xi = vortex.p;    % x coordinate of control points
xj = vortex.t;    % x coordinate of vortex points
delta2 = atan((zi-zj)./(xi-xj));    % angles from vortices to control points
delta3 = delta2 - atan(dz_dx);      % angles between velocities induced by
                                    % vortices and local camber line normal

% solve A*mj = b
A = 1/(2*pi) * cos(delta2).*cos(delta3)./(xi - xj);
b = U_oo * (alpha - atan(dz_dx)).';
vortex.Gamma = (A\b).';    % resuling magnitudes of vortices

%% Aerodynamic equations
% flow components
r_func = @(x, y) sqrt((x - xj).^2 + (y - zj).^2);    % distance to all vortices
v_theta_func = @(x, y) - 1./(2*pi*r_func(x, y)) .* vortex.Gamma;
u_func = @(x, y) cos(alpha)*U_oo + ...
                 sum((cos(acos((x - xj)./r_func(x, y)) + pi/2) .* ...
                     v_theta_func(x, y)));
v_func = @(x, y) sin(alpha)*U_oo + ...
                 sum((sin(acos((x - xj)./r_func(x, y)) + pi/2) .* ...
                     v_theta_func(x, y)));

% pressure coefficient
Delta_x = xj(2) - xj(1);
Delta_Cp = 2*vortex.Gamma/(U_oo*Delta_x);
Cp_func = @(u, v) 1 - (u.^2 + v.^2)./U_oo^2;

%% Plots
% plot parameters
lw = 1.2;    % line width
fs = 15;     % font size
ms = 7;      % marker size
load('colors.mat')


% figure
% hold on
% % airfoil.plot_self()
% % plot(xi, zi, '.', 'Color', colors.blue)
% % plot(xj, zj, '- .', 'Color', colors.red)
% plot(xj, Delta_Cp, '- .', 'Color', colors.blue)
% % axis image
% y_limits = ylim;
% ylim([y_limits(1), 4]);
% grid on
% hold off


% plotted results
x_plot = linspace(1e-3, 1, 100);
x_surf = linspace(-0.5, 1.5, 100);
y_surf = linspace(-1, 1, 100);
[X_grid, Y_grid] = meshgrid(x_surf, y_surf);
X_grid = X_grid.';
Y_grid = Y_grid.';
[u, v, Cp] = deal(zeros(size(X_grid)));
[u_plus, v_plus] = deal(zeros(size(x_plot)));
[u_min, v_min, Cp_min] = deal(1);
[u_max, v_max, Cp_max] = deal(0);
for i = 1:length(x_surf)
    for j = 1:length(y_surf)
        u(i,j) = u_func(x_surf(i), y_surf(j));
        v(i,j) = v_func(x_surf(i), y_surf(j));
        Cp(i,j) = Cp_func(u(i,j), v(i,j));
        [yU, yL] = airfoil.generate_geometry(x_surf(i));
        if x_surf(i) <= 0 || x_surf(i) >= 1 || ...
           y_surf(j) >= yU || y_surf(j) <= yL
            u_min = min([u(i,j), u_min]);     
            u_max = max([u(i,j), u_max]);
            v_min = min([v(i,j), v_min]);
            v_max = max([v(i,j), v_max]);
            Cp_min = min([Cp(i,j), Cp_min]);
            Cp_max = max([Cp(i,j), Cp_max]);
        end
    end
end
% for i = 1:length(x_plot)
%     u_plus(i) = u_func(x_plot(i), airfoil.yU(x_plot(i)));
%     v_plus(i) = v_func(x_plot(i), airfoil.yU(x_plot(i)));
% end
% Cp_plus = Delta_Cp(u_plus, v_plus);

% FIGURE 1: airfoil geometry and vortex line
figure
hold on
airfoil.plot_self()
plot(xi, zi, '.', 'Color', colors.yellow)
plot(xj, zj, '.', 'Color', colors.red)
% title({sprintf('NACA $00%02d$ airfoil geometry ', airfoil.xx), ...
%        'and defined vortex line'}, 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{y}{c}$', 'FontSize', fs)
% legend([plots.geometry, plots.vortex], ...
%        sprintf('NACA $00%02d$ geometry', airfoil.xx), ...
%        'Doublet line sections', ...
%        'Location', 'Northeast');
axis equal
grid on
hold off

% FIGURES 2-3: u flow component
figure
hold on
contour(X_grid, Y_grid, u/U_oo, 300)
caxis([u_min/U_oo, u_max/U_oo])
% contour(X_grid, Y_grid, u/U_oo, [1, Inf], 'ShowText', 'on', ...
%         'LineWidth', lw, 'LineColor', colors.dark_red)
% surf(X_grid, Y_grid, u/U_oo, 'LineStyle', 'none', 'FaceColor', 'interp')
airfoil.plot_self()
title('Flow component: $\frac{u}{U_\infty}$', 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{y}{c}$', 'FontSize', fs)
view(2)
colorbar
axis image
hold off

% figure
% hold on
% airfoil.plot_self()
% plot(x_plot, u_plus/U_oo, 'LineWidth', lw, 'Color', colors.red)
% title('Horizontal flow component along upper surface', 'FontSize', fs)
% xlabel('$\frac{x}{c}$', 'FontSize', fs)
% ylabel('$\frac{u}{U_{\infty}}$', 'FontSize', fs)
% axis image
% axis manual
% y_limits = ylim;
% ylim([y_limits(1), 1.4]);
% grid on
% hold off

% FIGURES 4-5: v flow component
figure
hold on
contour(X_grid, Y_grid, v/U_oo, 300)
caxis([v_min/U_oo, v_max/U_oo])
airfoil.plot_self()
title('Flow component: $\frac{v}{U_\infty}$', 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{y}{c}$', 'FontSize', fs)
view(2)
colorbar
axis image
hold off

% figure
% hold on
% airfoil.plot_self()
% plot(x_plot, v_plus/U_oo, 'LineWidth', lw, 'Color', colors.red)
% title('Vertical flow component along upper surface', 'FontSize', fs)
% xlabel('$\frac{x}{c}$', 'FontSize', fs)
% ylabel('$\frac{v}{U_{\infty}}$', 'FontSize', fs)
% axis image
% axis manual
% y_limits = ylim;
% ylim([y_limits(1), 0.8]);
% grid on
% hold off

% FIGUREs 5-6: pressure coefficient 
figure
hold on
contour(X_grid, Y_grid, Cp, 300)
caxis([Cp_min, Cp_max])
airfoil.plot_self()
title('Pressure coefficient', 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{y}{c}$', 'FontSize', fs)
view(2)
colorbar
axis image
hold off

figure
hold on
plot(xj, Delta_Cp, '.', 'LineWidth', lw, 'MarkerSize', ms)
title({'Pressure coefficient difference', 'along camber line'}, ...
      'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\Delta C_p$', 'FontSize', fs)
y_limits = ylim;
axis([0, 1, y_limits(1), 4]);
grid on
hold off

% FIGURE 7: stream lines
lines = 80;
starty1 = linspace(min(min(Y_grid)), max(max(Y_grid)), lines);
% starty2 = linspace(-0.03, 0.1, 4);
startx1 = min(min(X_grid))*ones(1, lines);
% startx2 = max(max(X_grid))*ones(1, 4);
figure
hold on
streamline(X_grid.', Y_grid.', u.', v.', startx1, starty1)
airfoil.plot_self(false)
% streamline(X_grid.', Y_grid.', -u.', -v.', startx2, starty2)
title('Stream lines', 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{y}{c}$', 'FontSize', fs)
xlim([min(min(X_grid)), max(max(X_grid))]);
ylim([min(min(Y_grid)), max(max(Y_grid))]);
axis image
grid on
hold off

function plots = plot_NACA_00xx(airfoil, vortex, plot_doublet_line)
% plot parameters
lw = 1.2;    % line width
ms = 7;      % marker size
load('colors.mat')

if nargin == 2
    plot_doublet_line = false;
end

geometry_x = linspace(0, 1);
geometry_y = airfoil.Y(geometry_x);

patch([geometry_x, fliplr(geometry_x)], ...
      [geometry_y, -fliplr(geometry_y)], ...
      colors.blue, 'FaceAlpha', 0.2)
if plot_doublet_line
    plots.vortex = plot(vortex.tj, zeros(length(vortex.tj), 1), '- .', ...
                         'LineWidth', lw, 'MarkerSize', ms, 'Color', colors.red);
end
plots.geometry = plot([geometry_x, fliplr(geometry_x)], ...
                      [geometry_y, -fliplr(geometry_y)], ...
                      'LineWidth', lw, 'Color', colors.blue);
end
