% ---------------------------------- %
%  AERODYNAMICS OF WINGS AND BODIES  %
%             HOMEWORK 6             %
%  Lior Isakov          Eden Shazar  %
%   200928372            204378608   %
% ---------------------------------- %

%% Fresh start
clc
close all
clear variables

%% Parameters
% user input
N = 300;               % number of tiles used per side of wing for VLM
uoo = 1;             % [m/s] - uniform flow speed
rho = 1.225;           % [kg/m^3] - air density (chosen as ISA at sea level)
alpha = deg2rad(6);    % [rad] - angle of attack

% givens
load('kfir_data.mat', 'chordLeadingEdges', 'chordLengths')

% calculated parameters
qoo = 1/2*rho*uoo^2;

%% Wing definition and calculations
kfir = FlatWing(chordLeadingEdges, chordLengths, N);
kfir.CalculateGammaValues(uoo, alpha);

rootChord = kfir.chordLengths(1);
tipChord = kfir.chordLengths(2);
span = kfir.span;
sweep = atan((rootChord - tipChord)/(span/2));
clear rootChord tipChord span

[cL, cMApex] = kfir.AeroCoeffs(uoo);
cLAlpha = cL/alpha;
cT = kfir.LeadingEdgeSuction(uoo, alpha);

%% Polhamus method
polhamus.kP = cLAlpha/cos(alpha)^2;
polhamus.kV = cT/(sin(alpha)^2*cos(sweep));
polhamus.cL = @(alpha) polhamus.kP*sin(alpha).*cos(alpha).^2 + ...
                       polhamus.kV*sin(alpha).^2.*cos(alpha);
polhamus.cDi = @(alpha) polhamus.cL(alpha).*tan(alpha);

%% Numeical results
load('wind_tunnel_data.mat', 'windTunnelData')

% error calculation
cLError = abs((polhamus.cL(deg2rad(windTunnelData.alpha)) - ...
              windTunnelData.cL) ./ ...
          windTunnelData.cL);
cDiError = abs((polhamus.cDi(deg2rad(windTunnelData.alpha)) - ...
               windTunnelData.cDi) ./ ...
          windTunnelData.cDi);
      
% excel sheets of data and errors
T = table(polhamus.cL(windTunnelData.alpha).', ...
          windTunnelData.cL.', ...
          cLError.'*100, ...
          polhamus.cDi(windTunnelData.alpha).', ...
          windTunnelData.cDi.', ...
          cDiError.'*100, ...
          'VariableNames', ...
          {'cL_polhamus', 'cL_wind_tunnel', 'cL_error_percent', ...
           'cDi_polhamus', 'cDi_wind_tunnel', 'cDi_error_percent'});
writetable(T, 'wind tunnel data comparison.xlsx')

%% Plots
% plot parameters
lw = 1.2;      % line width
fs = 15;       % font size
fss = 8.2;     % small font size
msx = 8;       % x marker size
msd = 11;      % dot marker size
load('colors.mat')

% FIGURE 1: geometry of wing and VLM tiles
figure
hold on
grid on
title('Kfir wing geometry and $VLM$ tiles', 'FontSize', fs)
xlabel('$x\ [m]$', 'FontSize', fs)
ylabel('$y\ [m]$', 'FontSize', fs)
wingColors = {colors.light_blue, colors.blue, colors.royal_blue colors.purple};
kfir.PlotSelf('tiles');
axis image
hold off

% FIGURE 2: cLAlpha comparison with wind tunnel data
alphaDegLinspace = linspace(-10, 60);
alphaRadLinspace = deg2rad(alphaDegLinspace);

figure
hold on
grid on
title('Polhamus method $C_L$ against wind tunnel data', 'FontSize', fs)
xlabel('$\alpha\ [deg]$', 'FontSize', fs)
ylabel('$C_L$', 'FontSize', fs)
plot(alphaDegLinspace, polhamus.cL(alphaRadLinspace), ...
     'LineWidth', lw, 'Color', colors.blue);
yLimits = ylim;
fill([15, 15, 30, 30], ...
     [yLimits(1), yLimits(2), yLimits(2), yLimits(1)], ...
     colors.yellow, 'LineStyle', 'none', 'FaceAlpha', 0.1)
plot([15, 15], yLimits, '--k')
plot([30, 30], yLimits, '--k')
plots(1) = plot(alphaDegLinspace, polhamus.cL(alphaRadLinspace), ...
                'LineWidth', lw, 'Color', colors.blue);
plots(2) = plot(windTunnelData.alpha, windTunnelData.cL, '.', ...
                'LineWidth', lw, 'MarkerSize', msd, 'Color', colors.red);
text(22.5, 0.17, {'Vortex breakdown', 'travels along chord'}, ...
     'FontSize', fss, 'HorizontalAlignment', 'center')
text(14.5, 1.3, {'Vortex breakdown', 'at trailing edge'}, ...
     'FontSize', fss, 'HorizontalAlignment', 'right')
text(30.5, 0.65, {'Vortex breakdown', 'at apex'}, ...
     'FontSize', fss, 'HorizontalAlignment', 'left')
legend(plots, ...
       'Polhamus method', 'Wing tunnel data', ...
       'Location', 'Southeast');
hold off

% FIGURE 3: cDi comparison with wind tunnel data
figure
hold on
grid on
title('Polhamus method $C_{D_i}$ against wind tunnel data', 'FontSize', fs)
xlabel('$\alpha\ [deg]$', 'FontSize', fs)
ylabel('$C_{D_i}$', 'FontSize', fs)
plot(alphaDegLinspace, polhamus.cDi(alphaRadLinspace), ...
     'LineWidth', lw);
yLimits = ylim;
fill([15, 15, 30, 30], ...
     [yLimits(1), yLimits(2), yLimits(2), yLimits(1)], ...
     colors.yellow, 'LineStyle', 'none', 'FaceAlpha', 0.1)
plot([15, 15], yLimits, '--k')
plot([30, 30], yLimits, '--k')
plots(1) = plot(alphaDegLinspace, polhamus.cDi(alphaRadLinspace), ...
                'LineWidth', lw, 'Color', colors.blue);
plots(2) = plot(windTunnelData.alpha, windTunnelData.cDi, '.', ...
                'LineWidth', lw, 'MarkerSize', msd, 'Color', colors.red);
text(22.5, 1.8, {'Vortex breakdown', 'travels along chord'}, ...
     'FontSize', fss, 'HorizontalAlignment', 'center')
text(14.5, 2.7, {'Vortex breakdown', 'at trailing edge'}, ...
     'FontSize', fss, 'HorizontalAlignment', 'right')
text(30.5, 2.7, {'Vortex breakdown', 'at apex'}, ...
     'FontSize', fss, 'HorizontalAlignment', 'left')
legend(plots, ...
       'Polhamus method', 'Wing tunnel data', ...
       'Location', 'Southeast');
hold off
