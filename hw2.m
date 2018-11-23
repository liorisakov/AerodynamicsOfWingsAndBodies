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
digits = 2412;         % NACA airfoil four-digit designation
N = 50;                % number of vortex line vortices
U_oo = 50;             % [m/s] - uniform flow speed
rho = 1.225;           % [kg/m^3] - air density (chosen as ISA at sea level)
alpha = deg2rad(8);    % [rad] - angle of attack

% calculated parameters
airfoil = NACA(num2str(digits));    % airfoil object
q_oo = 1/2*rho*U_oo^2;      % uniform flow dynamic pressure
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
delta2 = atan((zi - zj)./(xi - xj));    % angles from vortices to control points
delta3 = delta2 - atan(dz_dx);    % angles between velocities induced by
                                  % vortices and local camber line normal

% solve A*mj = b
A = 1/(2*pi) * cos(delta2).*cos(delta3)./(xi - xj);
b = U_oo * (alpha - atan(dz_dx)).';
vortex.Gamma = (A\b).';    % resuling magnitudes of vortices

%% Aerodynamic equations
% flow components
r_func = @(x, y) sqrt((x - xj).^2 + (y - zj).^2);    % distances to all vortices
v_theta_func = @(x, y) - 1./(2*pi*r_func(x, y)) .* vortex.Gamma;
u_func = @(x, y) cos(alpha)*U_oo + ...
                 sum((cos(acos((x - xj)./r_func(x, y)) + pi/2) .* ...
                     v_theta_func(x, y)));
v_func = @(x, y) sin(alpha)*U_oo + ...
                 sum((sin(acos((x - xj)./r_func(x, y)) + pi/2) .* ...
                     v_theta_func(x, y)));

% pressure coefficient
Delta_x = xj(2) - xj(1);
Delta_Cp = rho*U_oo*vortex.Gamma/(q_oo*Delta_x);
Cp_func = @(u, v) 1 - (u.^2 + v.^2)./U_oo^2;

% lift and moment coefficients
L_prime = rho*U_oo*sum(vortex.Gamma);
Cl = 2*sum(vortex.Gamma)/U_oo;
Cm_025c = 2*sum(vortex.Gamma.*(0.25 - xj))/U_oo;

%% Numerical results
% numerical output and Abbot & von Doenhoff comparison
load('NACA_2412_data.mat', 'abbott_data')
Cl_err = abs(Cl - abbott_data.Cl)/abbott_data.Cl;
Cm_025c_err = abs(Cm_025c - abbott_data.Cm_025c)/abbott_data.Cm_025c;
fprintf('Lift per unit wing span: %.0f [N/m]\n\n', L_prime)
fprintf('Lift coeffcient:\n')
fprintf('    Calculated: %.3f\n', Cl)
fprintf('    Abbott & von Doenhoff: %.3f\n', abbott_data.Cl)
fprintf('    Absolute error: %.1f%%\n\n', Cl_err*100)
fprintf('Moment coeffcient about quarter chord:\n')
fprintf('    Calculated: %.4f\n', Cm_025c)
fprintf('    Abbott & von Doenhoff: %.4f\n', abbott_data.Cm_025c)
fprintf('    Absolute error: %.1f%%\n\n', Cm_025c_err*100)

% Kuethe & Chow comparison
load('NACA_2412_data.mat', 'kuethe_data')
Delta_Cp_comp = interp1(xj, Delta_Cp, kuethe_data.x, 'pchip');
Delta_Cp_err = abs(Delta_Cp_comp - kuethe_data.Delta_Cp)./kuethe_data.Delta_Cp;
T = table(kuethe_data.x.', ...
          kuethe_data.Delta_Cp.', Delta_Cp_comp.', Delta_Cp_err.', ...
          'VariableNames', {'x', ...
          'Delta_Cp_Kuethe', 'Delta_Cp_calculation', 'Delta_Cp_error'});
writetable(T, 'Kuethe & Chow comparison.xlsx')

%% Plots
% plot parameters
lw = 1.2;    % line width
fs = 15;     % font size
ms = 7;      % marker size
load('colors.mat')

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

% FIGURE 1: airfoil geometry and vortex line
figure
hold on
h1 = airfoil.plot_self();
h2 = plot(xj, zj, '.', 'Color', colors.red);
h3 = plot(xi, zi, '.', 'Color', colors.yellow);
title({sprintf('NACA $%d$ airfoil geometry,', digits), ...
       'vortices and control points'}, 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{y}{c}$', 'FontSize', fs)
legend([h1, h2, h3], sprintf('NACA $%d$ geometry', digits), ...
       'Vortex points', 'Control points', 'Location', 'Northeast');
axis equal
grid on
hold off

% FIGURE 2: u flow component
figure
hold on
contour(X_grid, Y_grid, u/U_oo, 300)
caxis([u_min/U_oo, u_max/U_oo])
airfoil.plot_self();
title('Flow component: $\frac{u}{U_\infty}$', 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{y}{c}$', 'FontSize', fs)
view(2)
colorbar
axis image
hold off

% FIGURE 3: v flow component
figure
hold on
contour(X_grid, Y_grid, v/U_oo, 300)
caxis([v_min/U_oo, v_max/U_oo])
airfoil.plot_self();
title('Flow component: $\frac{v}{U_\infty}$', 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{y}{c}$', 'FontSize', fs)
view(2)
colorbar
axis image
hold off

% FIGURES 4-5: pressure coefficient 
figure
hold on
contour(X_grid, Y_grid, Cp, 300)
caxis([Cp_min, Cp_max])
airfoil.plot_self();
title('Pressure coefficient', 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{y}{c}$', 'FontSize', fs)
view(2)
colorbar
axis image
hold off

figure
hold on
plot(xj, Delta_Cp, '- .', 'LineWidth', lw, 'MarkerSize', ms)
plot(kuethe_data.x, kuethe_data.Delta_Cp, '- .', ...
     'LineWidth', lw, 'MarkerSize', ms)
title({'Pressure coefficient difference', 'along camber line'}, ...
      'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\Delta C_p$', 'FontSize', fs)
legend('Calculation', 'Kuethe \& Chow 120 panel data', ...
       'Location', 'Northeast');
y_limits = ylim;
axis([0, 1, y_limits(1), 4]);
grid on
hold off

% FIGURE 7: stream lines
lines = 80;
starty = linspace(min(min(Y_grid)), max(max(Y_grid)), lines);
startx = min(min(X_grid))*ones(1, lines);
figure
hold on
streamline(X_grid.', Y_grid.', u.', v.', startx, starty)
airfoil.plot_self(true);
airfoil.plot_camber();
plot(xj, zj, 'Color', colors.red, 'LineWidth', lw)
title('Stream lines', 'FontSize', fs)
xlabel('$\frac{x}{c}$', 'FontSize', fs)
ylabel('$\frac{y}{c}$', 'FontSize', fs)
xlim([min(min(X_grid)), max(max(X_grid))]);
ylim([min(min(Y_grid)), max(max(Y_grid))]);
axis image
grid on
hold off