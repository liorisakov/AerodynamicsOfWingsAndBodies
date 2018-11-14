clc
close all
clear variables

%% Parameters
xx = 15;          % last two NACA airfoil digits
c = 1;            % cord
t = xx/100;       % max thickness as fraction of cord

x_s = 1e-3;        % start of doublet line
x_f = c - 1e-3;    % end of doublet line
N = 50;            % number of doublet line sections

U_oo = 1;

%% Wing
% doublet line
theta_j = linspace(0, pi, N);
c_star = x_f - x_s;
t_j = c_star/2 * (1 - cos(theta_j)) + x_s;
x_i = ((t_j(1:end-1) + t_j(2:end))/2).';

% geometry
Y = @(x) 5*t * (0.2969 * sqrt(x) - ...
                0.1260 * x - ...
                0.3516 * x.^2 + ...
                0.2843 * x.^3 - ...
                0.1015 * x.^4);

%% Solver
A = 1./Y(x_i).*(atan((t_j(2:end) - x_i)./Y(x_i)) - ...
                atan((t_j(1:end-1) - x_i)./Y(x_i)));
% for i = 1:N-1
%     for j = 1:N-1
%         l = (t_j(j+1) - x_i(i))/Y(x_i(i));
%         m = (t_j(j) - x_i(i))/Y(x_i(i));
%         A2(i,j) = (atan(l) - atan(m))/Y(x_i(i));
%         fprintf('%f ', l)
%     end
%     fprintf('\n')
% end
b = ones(N-1, 1);

m_j = A\b;

u = @(x, y) U_oo*(1 + ...
                  sum(m_j.' .* ((t_j(2:end) - x)./(((t_j(2:end) - x).^2) + y^2) - ...
                                (t_j(1:end-1) - x)./(((t_j(1:end-1) - x).^2) + y^2))));
v = @(x, y) - U_oo*(sum(m_j.' .* (y./((t_j(2:end) - x).^2 + y^2) - ...
                                  y./((t_j(1:end-1) - x).^2 + y^2))));
% psi = @(x, y) U_oo*Y(x) - 1/(2*pi) * ...
%               sum(m_j.' .* (atan((t_j(2:end) - x)./Y(x)) - ...
%                             atan((t_j(1:end-1) - x)./Y(x))));
%% Plots
% plot parameters
lw = 1.2;    % line width
fs = 15;     % font size
ms = 7;      % marker size
load('colors.mat')

% temp
x = linspace(-1, 2, N);
y = linspace(-0.5, 0.5, N);
[X_grid, Y_grid] = meshgrid(x, y);
X_grid = X_grid.';
Y_grid = Y_grid.';
for i = 1:length(x)
    for j = 1:length(y)
        U(i,j) = u(x(i), y(j));
        V(i,j) = v(x(i), y(j));
    end
end

figure
hold on
geometry_x = linspace(0, c);
geometry_y = Y(geometry_x);
patch([geometry_x, fliplr(geometry_x)], ...
      [geometry_y, -fliplr(geometry_y)], ...
      colors.blue, 'FaceAlpha', 0.2)
plot(t_j, zeros(length(t_j), 1), '.', ...
     'LineWidth', lw, 'MarkerSize', ms, 'Color', colors.red)
plot(t_j, zeros(length(t_j), 1), ...
     'LineWidth', lw, 'Color', colors.red)
plot([geometry_x, fliplr(geometry_x)], ...
     [geometry_y, -fliplr(geometry_y)], ...
     'LineWidth', lw, 'Color', colors.blue)
axis equal
grid on
hold off

figure
hold on
surf(X_grid, Y_grid, U, 'LineStyle', 'none')
view(2)
colorbar
grid on
hold off

figure
hold on
surf(X_grid, Y_grid, V, 'LineStyle', 'none')
view(3)
colorbar
grid on
hold off

figure
hold on
geometry_x = linspace(0, c);
geometry_y = Y(geometry_x);
patch([geometry_x, fliplr(geometry_x)], ...
      [geometry_y, -fliplr(geometry_y)], ...
      colors.blue, 'FaceAlpha', 0.2)
plot([geometry_x, fliplr(geometry_x)], ...
     [geometry_y, -fliplr(geometry_y)], ...
     'LineWidth', lw, 'Color', colors.blue)

starty = linspace(min(min(Y_grid)), max(max(Y_grid)), 50);
startx = -1*ones(size(starty));
streamline(X_grid.', Y_grid.', U.', V.', startx, starty)
axis equal
ylim([min(min(Y_grid)), max(max(Y_grid))]);

grid on
hold off
