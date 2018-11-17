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
% xs_offset = 1e-3;    % forward offset of start of vortex line from origin
% xf_offset = 1e-3;    % backward offset of end of vortex line from cord end
N = 50;     % number of vortex line vortices
U_oo = 1;            % horizontal uniform flow speed
alpha = deg2rad(8);    % angle of attack

% calculated parameters
x = linspace(0, 1, N+2);    % vortex line sectioning;
                            % LE and TE points don't have vortices

% vortex.xs = xs_offset;        % start of vortex line
% vortex.xf = 1 - xf_offset;    % end of vortex line


%% Wing
% vortex line definition
% vortex.c_star = vortex.xf - vortex.xs;        % vortex line length
vortex.t = x(2:end-1);    % vortex points
vortex.p = ((vortex.t + x(3:end))/2).';    % control points

% camber line
dz_dx = airfoil.camber_slope(vortex.p).';
zi = airfoil.camber_line(vortex.p);
zj = airfoil.camber_line(vortex.t);

%% Solver
% parameters
xi = vortex.p;
xj = vortex.t;
delta2 = atan((zi-zj)/(xi-xj));
delta3 = delta2 + pi/2 - atan(dz_dx);

% solve A*mj = b
A = 1/(2*pi) * cos(delta2).*cos(delta3)./(xi - xj);
b = U_oo * (alpha - atan(dz_dx)).';
vortex.Gamma = (A\b).';

%% Aerodynamic equations
% flow components
% u_eq = @(x, y) U_oo * (1 + ...
%                        sum(vortex.mj.*((vortex.tj(2:end) - x) ./ ...
%                                         (((vortex.tj(2:end) - x).^2) + y^2) - ...
%                                         (vortex.tj(1:end-1) - x) ./ ...
%                                         (((vortex.tj(1:end-1) - x).^2) + y^2))));
% v_eq = @(x, y) - U_oo * ...
%                (sum(vortex.mj.*(y./((vortex.tj(2:end) - x).^2 + y^2) - ...
%                                  y./((vortex.tj(1:end-1) - x).^2 + y^2))));

% pressure coefficient
Delta_x = x(2) - x(1);
Delta_Cp = 2*vortex.Gamma/(U_oo*Delta_x);

%% Plots
% plot parameters
lw = 1.2;    % line width
fs = 15;     % font size
ms = 7;      % marker size
load('colors.mat')

% plotted results
x_plot = linspace(1e-3, 1, 400);
x_surf = linspace(-0.5, 1.5, 400);
y_surf = linspace(-1, 1, 400);
[X_grid, Y_grid] = meshgrid(x_surf, y_surf);
X_grid = X_grid.';
Y_grid = Y_grid.';
[u, v] = deal(zeros(size(X_grid)));
[u_plus, v_plus] = deal(zeros(size(x_plot)));
for i = 1:length(x_surf)
    for j = 1:length(y_surf)
        if x_surf(i) <= 0 || x_surf(i) >= 1 || ...
           y_surf(j) >= airfoil.(x_surf(i)) || y_surf(j) <= - airfoil.Y(x_surf(i))
            u(i,j) = u_eq(x_surf(i), y_surf(j));
            v(i,j) = v_eq(x_surf(i), y_surf(j));
        else
            [u(i,j), v(i,j)] = deal(NaN);
        end
    end
end
Cp = Delta_Cp(u, v);
for i = 1:length(x_plot)
    u_plus(i) = u_eq(x_plot(i), airfoil.Y(x_plot(i)));
    v_plus(i) = v_eq(x_plot(i), airfoil.Y(x_plot(i)));
end
Cp_plus = Delta_Cp(u_plus, v_plus);

% FIGURE 1: airfoil geometry and vortex line
figure
hold on
plots = plot_NACA_00xx(airfoil, vortex, true);
title({sprintf('NACA $00%02d$ airfoil geometry ', airfoil.xx), ...
       'and defined vortex line'}, 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{y}{c}$', 'FontSize', fs)
legend([plots.geometry, plots.vortex], ...
       sprintf('NACA $00%02d$ geometry', airfoil.xx), ...
       'Doublet line sections', ...
       'Location', 'Northeast');
axis equal
grid on
hold off

% FIGURES 2-3: u flow component
figure
hold on
surf(X_grid, Y_grid, u/U_oo, 'LineStyle', 'none', 'FaceColor', 'interp')
plot_NACA_00xx(airfoil, vortex);
title('Flow component: $\frac{u}{U_\infty}$', 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{y}{c}$', 'FontSize', fs)
view(2)
colorbar
axis image
hold off

figure
hold on
plot_NACA_00xx(airfoil, vortex);
plot(x_plot, u_plus/U_oo, 'LineWidth', lw, 'Color', colors.red)
title('Horizontal flow component along upper surface', 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{u}{U_{\infty}}$', 'FontSize', fs)
axis image
axis manual
y_limits = ylim;
ylim([y_limits(1), 1.4]);
grid on
hold off

% FIGURES 4-5: v flow component
figure
hold on
surf(X_grid, Y_grid, v/U_oo, 'LineStyle', 'none', 'FaceColor', 'interp')
plot_NACA_00xx(airfoil, vortex);
title('Flow component: $\frac{v}{U_\infty}$', 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{y}{c}$', 'FontSize', fs)
view(2)
colorbar
axis image
hold off

figure
hold on
plot_NACA_00xx(airfoil, vortex);
plot(x_plot, v_plus/U_oo, 'LineWidth', lw, 'Color', colors.red)
title('Vertical flow component along upper surface', 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{v}{U_{\infty}}$', 'FontSize', fs)
axis image
axis manual
y_limits = ylim;
ylim([y_limits(1), 0.8]);
grid on
hold off

% FIGUREs 5-6: pressure coefficient 
figure
hold on
surf(X_grid, Y_grid, Cp, 'LineStyle', 'none', 'FaceColor', 'interp')
plot_NACA_00xx(airfoil, vortex);
title('Pressure coefficient', 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{y}{c}$', 'FontSize', fs)
view(2)
colorbar
axis image
hold off

figure
hold on
plot_NACA_00xx(airfoil, vortex);
plot(x_plot, 1 - Cp_plus, 'LineWidth', lw, 'Color', colors.red)
title('Pressure coefficient along upper surface', 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{u^2+v^2}{U_{\infty}^2}$', 'FontSize', fs)
axis image
axis manual
y_limits = ylim;
ylim([y_limits(1), 1.5]);
grid on
hold off

% FIGURE 7: stream lines
figure
hold on
plot_NACA_00xx(airfoil, vortex);

starty = linspace(min(min(Y_grid)), max(max(Y_grid)), 50);
startx = min(min(X_grid))*ones(size(starty));
streamline(X_grid.', Y_grid.', u.', v.', startx, starty)

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
