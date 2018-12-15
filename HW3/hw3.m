% ---------------------------------- %
%  AERODYNAMICS OF WINGS AND BODIES  %
%             HOMEWORK 3             %
%  Lior Isakov          Eden Shazar  %
%   200928372            204378608   %
% ---------------------------------- %

%% Fresh start
clc
close all
clear variables

%% Parameters
load('boeing 737-300 geometry data.mat')
load('alpha for required CL.mat')

% user input
N = [20, 50];          % number of tiles used per side of wing for VLM
% N = 20;              % uncomment to draw representative geometry
uoo = 50;              % [m/s] - uniform flow speed
rho = 1.225;           % [kg/m^3] - air density (chosen as ISA at sea level)
alpha = alphaForRequiredCL;    % [rad] - angle of attack
% alpha = deg2rad(8);    % uncomment to input angle of attack manually

% calculated parameters
qoo = 1/2*rho*uoo^2;

%% Wing
wing = FlatWing(chordLeadingEdges, chordLengths, N);

%% Solver
% initialization
xTileCount = size(wing.tiles, 1);
yTileCount = size(wing.tiles, 2);
tileCount = xTileCount*yTileCount;
A = zeros(tileCount/2, tileCount/2);

% sum downwash at every control point on positive y half of wing
% to build matrix A
k = 1;
for i = 1:xTileCount
    for j = 1:yTileCount/2
        % current symmetric tile pair, at positive and negative y
        yPosTile = wing.tiles{i,j};
        yNegTile = wing.tiles{i,j+yTileCount/2};
        
        % contributions of both symmetric horseshoes at current control point
        yPosDownwash = yPosTile.Downwash(wing.controlPointsX, ...
                                         wing.controlPointsY, ...
                                         'solve for gamma', true).';
        yNegDownwash = yNegTile.Downwash(wing.controlPointsX, ...
                                         wing.controlPointsY, ...
                                         'solve for gamma', true).';
        
        A(:,k) = yPosDownwash(:) + yNegDownwash(:);

        k = k + 1;
    end
end

b = - uoo*sin(alpha)*ones(tileCount/2, 1);

gamma = A\b;
gamma = reshape(gamma, yTileCount/2, xTileCount).';

for i = 1:xTileCount
    for j = 1:yTileCount/2
        wing.tiles{i,j}.horseshoe.gamma = gamma(i,j);
        wing.tiles{i,j+yTileCount/2}.horseshoe.gamma = gamma(i,j);
    end
end

%% Aerodynamic equations
% calculation points along span
y = linspace(-wing.span/2, wing.span/2);

% aerodynamic coefficients
[cl, cmApex] = wing.AeroCoeffsAtY(y, uoo);
[clTimesC, cmApexTimesCSquared] = ...
    wing.AeroCoeffsAtY(y, uoo, 'cl times c', 'cm times c squared');
[cL, cMApex] = wing.AeroCoeffs(y, uoo);
cLAlpha = cL/alpha;
cMApexAlpha = cMApex/alpha;
xCpNondimensional = - cMApexAlpha/cLAlpha;
xCp = xCpNondimensional*wing.mac;

% lift and moment about apex
lPrime = qoo * clTimesC;
cmApexPrime = qoo * cmApex;
L = qoo * wing.area * cL;
MApex = qoo * wing.area * wing.mac * cMApex;

% angle of attack required for CL = 0.5
alphaForRequiredCL = 0.5/cLAlpha;

% mean aerodynamic chord quarter chord x coordinate
macQuarterChord = wing.macLeadingEdgeX + 0.25*wing.mac;

% xCp as a percentage of mean aerodynamic chord
xCpMacPercentage = (xCp - wing.macLeadingEdgeX)/wing.mac;

%% Numerical results
% Kuethe & Chow comparison at stations along y
% for over 20 tiles in the x direction only
if wing.N(1) >= 20
    load('Kuethe and Chow data.mat')
    clComparison = wing.AeroCoeffsAtY(kuethe.yStations, uoo, ...
                                      'return separate x values');
    deltaCpComparison.station2 = clComparison(:,1).' * xTileCount;
    deltaCpComparison.station4 = clComparison(:,2).' * xTileCount;
    deltaCpComparison.station6 = clComparison(:,3).' * xTileCount;
    deltaCpComparison.station8 = clComparison(:,4).' * xTileCount;

    deltaCpError.station2 = abs(deltaCpComparison.station2 - ...
                                kuethe.deltaCp.station2) ./ ...
                            kuethe.deltaCp.station2;
    deltaCpError.station4 = abs(deltaCpComparison.station4 - ...
                                kuethe.deltaCp.station4) ./ ...
                            kuethe.deltaCp.station4;
    deltaCpError.station6 = abs(deltaCpComparison.station6 - ...
                                kuethe.deltaCp.station6) ./ ...
                            kuethe.deltaCp.station6;
    deltaCpError.station8 = abs(deltaCpComparison.station8 - ...
                                kuethe.deltaCp.station8) ./ ...
                            kuethe.deltaCp.station8;

    % table output
    x = linspace(0, 1, 20);
    T = table(x.', kuethe.deltaCp.station2.',  ...
              deltaCpComparison.station2.', deltaCpError.station2.', ...
              'VariableNames', {'xOverC', 'deltaCpKuethe', ...
              'deltaCpCalculation', 'deltaCpError'});
    writetable(T, 'Kuethe & Chow comparison.xlsx', 'Sheet', 'Station 2')
    T = table(x.', kuethe.deltaCp.station4.',  ...
              deltaCpComparison.station4.', deltaCpError.station4.', ...
              'VariableNames', {'xOverC', 'deltaCpKuethe', ...
              'deltaCpCalculation', 'deltaCpError'});
    writetable(T, 'Kuethe & Chow comparison.xlsx', 'Sheet', 'Station 4')
    T = table(x.', kuethe.deltaCp.station6.',  ...
              deltaCpComparison.station6.', deltaCpError.station6.', ...
              'VariableNames', {'xOverC', 'deltaCpKuethe', ...
              'deltaCpCalculation', 'deltaCpError'});
    writetable(T, 'Kuethe & Chow comparison.xlsx', 'Sheet', 'Station 6')
    T = table(x.', kuethe.deltaCp.station8.',  ...
              deltaCpComparison.station8.', deltaCpError.station8.', ...
              'VariableNames', {'xOverC', 'deltaCpKuethe', ...
              'deltaCpCalculation', 'deltaCpError'});
    writetable(T, 'Kuethe & Chow comparison.xlsx', 'Sheet', 'Station 8')
end

% Kuethe % Chow numerical comparison
fprintf('Results for N = [%d,%d] and alpha = %.3f [deg]:\n\n', ...
        xTileCount, yTileCount/2, rad2deg(alpha))
fprintf('Wing geometric properties:\n')
fprintf('b = %.2f [m]\n', wing.span)
fprintf('S = %.2f [m^2]\n', wing.area)
fprintf('AR = %.3f\n', wing.aspectRatio)
fprintf('MAC = %.3f [m]\n\n', wing.mac)
fprintf('Aerodynamic coefficients:\n')
fprintf('cL = %.3f\n', cL)
fprintf('cM about apex = %.3f\n', cMApex)
fprintf('cLAlpha = %.3f [1/rad]\n', cLAlpha)
fprintf('cMAlpha about apex = %.3f [1/rad]\n\n', cMApexAlpha)
fprintf('Center of pressure\n')
fprintf('xCp = %.3f [m]\n', xCp)
fprintf('xCp as percentage along MAC = %.1f%%\n\n', xCpMacPercentage*100)

%% Plots
% plot parameters
lw = 1.2;    % line width
fs = 15;     % font size
ms = 7;      % marker size
load('colors.mat')

% FIGURE 1: wing and method geometry (for less than 20 total tiles only)
if prod(N) <= 20
    figure
    hold on
    grid on
    handles = wing.PlotSelf('tiles', 'horseshoes', 'control points');
    title({'Boeing $737$-$300$ wing geometry', ...
           'and VLM visualization'}, 'FontSize', fs)
    xlabel('$x$', 'FontSize', fs)
    ylabel('$y$', 'FontSize', fs)
    axis image
    xlim([-10, 15])
    legend([handles.outline, handles.tileHandles.tile, ...
            handles.tileHandles.horseshoe, ...
            handles.tileHandles.controlPoint], ...
           'Wing geometry', 'VLM tiles', ...
           'Vortex horseshoes', 'Control Points', ...
           'Location', 'Southwest');
    hold off
end

% FIGURE 2: downwash around planar wing (for less than 20 total tiles only)
if prod(N) <= 20
    [X, Y] = meshgrid(linspace(-8, 17, 500), linspace(-18, 18, 500));
    w = 0;
    for i = 1:xTileCount
        for j = 1:yTileCount
            w = w + wing.tiles{i,j}.Downwash(X, Y);
        end
    end
    
    figure
    hold on
    surf(X, Y, w, 'LineStyle', 'None');
    handles = wing.PlotSelf('no fill');
    set(handles.outline, 'Color', colors.pink);
    view(2)
    c = colorbar;
    caxis([-20, 20])
    cLabels = get(c, 'TickLabels');
    cLabels{1} = '$\le-20$';
    cLabels{end} = '$\ge20$';
    set(c, 'TickLabelInterpreter', 'latex', 'FontSize', 10, ...
        'TickLabels', cLabels);
    title({'Induced vertical velocity in the vicinity', ...
           'of the wing in $\frac{m}{s}$'}, 'FontSize', fs)
    xlabel('$x$', 'FontSize', fs)
    ylabel('$y$', 'FontSize', fs)
    legend(handles.outline, 'Wing geometry', 'Location', 'Southwest')
    axis image
    hold off
end

% FIGURE 2: cl(y)
figure
hold on
grid on
plot(y/(wing.span/2), cl, ...
     'LineWidth', lw, 'Color', colors.red)
title({'Lift coefficient of wing section', ...
       'as a function of $y$'}, 'FontSize', fs)
xlabel('$\frac{y}{b}$', 'FontSize', fs)
ylabel('$C_l(y)$', 'FontSize', fs)
axis image
yTicks = yticks;
middleY = (max(yTicks) + min(yTicks))/2;
wing.PlotSelf('rotate', 'normalize', 'origin', [0, middleY]);
hold off

% FIGURE 3: cL(y)*c(y)
figure
hold on
grid on
plot(y/(wing.span/2), clTimesC, ...
     'LineWidth', lw, 'Color', colors.red)
title({'Lift coefficient of wing section', ...
       'times chord as a function of $y$'}, 'FontSize', fs)
xlabel('$\frac{y}{b}$', 'FontSize', fs)
ylabel('$C_l(y)\cdot c(y)\ [m]$', 'FontSize', fs)
axis image
yTicks = yticks;
middleY = (max(yTicks) + min(yTicks))/2;
wing.PlotSelf('rotate', 'normalize', 'origin', [0, middleY]);
hold off

% FIGURE 4: cmApex(y)
figure
hold on
grid on
plot(y/(wing.span/2), cmApex, ...
     'LineWidth', lw, 'Color', colors.red)
title({'Moment coefficient of wing section', ...
       'about apex as a function of $y$'}, 'FontSize', fs)
xlabel('$\frac{y}{b}$', 'FontSize', fs)
ylabel('$C_{m,apex}(y)$', 'FontSize', fs)
axis image
yTicks = yticks;
middleY = (max(yTicks) + min(yTicks))/2;
wing.PlotSelf('rotate', 'normalize', 'origin', [0, middleY]);
hold off

% FIGURE 5: cmApex(y)*c(y)^2
figure
hold on
grid on
plot(y/(wing.span/2), cmApexTimesCSquared, ...
     'LineWidth', lw, 'Color', colors.red)
title({'Moment coefficient of wing section about apex', ...
       'times chord squared as a function of $y$'}, 'FontSize', fs)
xlabel('$\frac{y}{b}$', 'FontSize', fs)
ylabel('$C_{m,apex}(y)\cdot c(y)^2\ [m^2]$', 'FontSize', fs)
axis image
yTicks = yticks;
middleY = (max(yTicks) + min(yTicks))/2;
wing.PlotSelf('rotate', 'normalize', 'origin', [0, middleY]);
hold off

% Kuethe & Chow comparison (for over 20 tiles in the x direction only)
if wing.N(1) == 20    % if there are enough tiles in the x direction
    load('Kuethe and Chow data.mat')
    x = linspace(0, 1, 20);
    y = kuethe.yStations;
	
	% FIGURE 6: DeltaCp comparison with Kuethe & Chow
    figure
    sgtitle('$\Delta C_p$ comparison at four stations along span', ...
            'FontSize', fs)
    subplot(2, 2, 1)
    hold on
    grid on
    plot(x, deltaCpComparison.station2, 'LineWidth', lw)
    plot(x, kuethe.deltaCp.station2, 'LineWidth', lw)
    title(strcat('Station $2$: $\frac{y}{0.5b}=', sprintf('%.3f$', ...
          kuethe.yStations(1)/(wing.span/2))), 'FontSize', fs-4);
    xlabel('$\frac{x}{c(y)}$', 'FontSize', fs-4)
    ylabel('$\Delta C_p$', 'FontSize', fs-4)
    hold off
    subplot(2, 2, 2)
    hold on
    grid on
    plot(x, deltaCpComparison.station4, 'LineWidth', lw)
    plot(x, kuethe.deltaCp.station4, 'LineWidth', lw)
    title(strcat('Station $4$: $\frac{y}{0.5b}=', sprintf('%.3f$', ...
          kuethe.yStations(2)/(wing.span/2))), 'FontSize', fs-4);
    xlabel('$\frac{x}{c(y)}$', 'FontSize', fs-4)
    ylabel('$\Delta C_p$', 'FontSize', fs-4)
    legend('Calculation', 'Kuethe $\&$ Chow data')
    hold off
    subplot(2, 2, 3)
    hold on
    grid on
    plot(x, deltaCpComparison.station6, 'LineWidth', lw)
    plot(x, kuethe.deltaCp.station6, 'LineWidth', lw)
    title(strcat('Station $6$: $\frac{y}{0.5b}=', sprintf('%.3f$', ...
          kuethe.yStations(3)/(wing.span/2))), 'FontSize', fs-4);
    xlabel('$\frac{x}{c(y)}$', 'FontSize', fs-4)
    ylabel('$\Delta C_p$', 'FontSize', fs-4)
    hold off
    subplot(2, 2, 4)
    hold on
    grid on
    plot(x, deltaCpComparison.station8, 'LineWidth', lw)
    plot(x, kuethe.deltaCp.station8, 'LineWidth', lw)
    title(strcat('Station $8$: $\frac{y}{0.5b}=', sprintf('%.3f$', ...
          kuethe.yStations(4)/(wing.span/2))), 'FontSize', fs-4);
    xlabel('$\frac{x}{c(y)}$', 'FontSize', fs-4)
    ylabel('$\Delta C_p$', 'FontSize', fs-4)
    hold off
    
    % FIGURE 7: DeltaCp comparison errors
    figure
    sgtitle('$\Delta C_p$ comparison errors at four stations along span', ...
            'FontSize', fs)
    subplot(2, 2, 1)
    hold on
    grid on
    plot(x, deltaCpError.station2, 'LineWidth', lw)
    title(strcat('Station $2$: $\frac{y}{0.5b}=', sprintf('%.3f$', ...
          kuethe.yStations(1)/(wing.span/2))), 'FontSize', fs-4);
    xlabel('$\frac{x}{c(y)}$', 'FontSize', fs-4)
    ylabel('$\Delta C_p$ absolute errors', 'FontSize', fs-4)
    hold off
    subplot(2, 2, 2)
    hold on
    grid on
    plot(x, deltaCpError.station4, 'LineWidth', lw)
    title(strcat('Station $4$: $\frac{y}{0.5b}=', sprintf('%.3f$', ...
          kuethe.yStations(2)/(wing.span/2))), 'FontSize', fs-4);
    xlabel('$\frac{x}{c(y)}$', 'FontSize', fs-4)
    ylabel('$\Delta C_p$ absolute errors', 'FontSize', fs-4)
    hold off
    subplot(2, 2, 3)
    hold on
    grid on
    plot(x, deltaCpError.station6, 'LineWidth', lw)
    title(strcat('Station $6$: $\frac{y}{0.5b}=', sprintf('%.3f$', ...
          kuethe.yStations(3)/(wing.span/2))), 'FontSize', fs-4);
    xlabel('$\frac{x}{c(y)}$', 'FontSize', fs-4)
    ylabel('$\Delta C_p$ absolute errors', 'FontSize', fs-4)
    hold off
    subplot(2, 2, 4)
    hold on
    grid on
    plot(x, deltaCpError.station8, 'LineWidth', lw)
    title(strcat('Station $8$: $\frac{y}{0.5b}=', sprintf('%.3f$', ...
          kuethe.yStations(4)/(wing.span/2))), 'FontSize', fs-4);
    xlabel('$\frac{x}{c(y)}$', 'FontSize', fs-4)
    ylabel('$\Delta C_p$ absolute errors', 'FontSize', fs-4)
    hold off
end

% FIGURE 8: mean aerodynamic chord and xCp
figure
hold on
grid on
handles = wing.PlotSelf('MAC');
handles.xCp = plot(xCp*[1, 1], wing.span/2*[-1, 1], '--', ...
                   'LineWidth', lw, 'Color', colors.blue);
title({'$x_{C_p}$ compared to quarter length', ...
       'of mean aerodynamic chord'}, 'FontSize', fs)
xlabel('$x\ [m]$', 'FontSize', fs)
ylabel('$y\ [m]$', 'FontSize', fs)
legend([handles.mac, handles.quarterMac, handles.xCp], ...
       'Mean aerodynamic chord', ...
       'Mean aerodynamic quarter chord', '$x_{C_p}$', ...
       'Location', 'EastOutside')
axis image
hold off
