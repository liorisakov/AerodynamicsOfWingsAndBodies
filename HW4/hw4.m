% ---------------------------------- %
%  AERODYNAMICS OF WINGS AND BODIES  %
%             HOMEWORK 4             %
%  Lior Isakov          Eden Shazar  %
%   200928372            204378608   %
% ---------------------------------- %

%% Fresh start
clc
close all
clear variables

%% Parameters
% user input
% N = [20, 50];          % number of tiles used per side of wing for VLM
N = 150;               % number of tiles used per side of wing for VLM
uoo = 50;              % [m/s] - uniform flow speed
rho = 1.225;           % [kg/m^3] - air density (chosen as ISA at sea level)
alpha = deg2rad(6);    % [rad] - angle of attack
rootChord = 1;         % [m] - chord length at wing root

% givens
AR = [1, 3, 6, 10];    % [-] - wing aspect ratios
lambda = 0.2;          % [-] - trapezoidal wing taper ratio

% calculated parameters
qoo = 1/2*rho*uoo^2;
tipChord = lambda*rootChord;
span = @(AR) AR*(rootChord + tipChord)/2;
sweep = @(span) atan((rootChord - tipChord)./(span/2));

wingCount = length(AR);
chordLeadingEdges = cell(1, wingCount);
chordLengths = cell(1, wingCount);
for i = 1:wingCount
    b = span(AR(i));
    Lambda(i) = sweep(b);
    
    chordLeadingEdges{i} = [0, 0;
                            (1 - lambda)*rootChord, b/2];
    chordLengths{i} = [rootChord, tipChord];
end
clear b

%% Wing definition and calculations
% initialization
wings = cell(1, wingCount);
[cL, cLAlpha, xCpMacPercentage, xCpCStarPercentage, cT, cDi] = ...
    deal(zeros(1, wingCount));

for i = 1:wingCount
    % define wing and perform VLM calculations
    wings{i} = FlatWing(chordLeadingEdges{i}, chordLengths{i}, N);
    wings{i}.CalculateGammaValues(uoo, alpha);
    
    % lift and moment coefficients
    [cL(i), cMApex] = wings{i}.AeroCoeffs(uoo);
    cLAlpha(i) = cL(i)/alpha;
    cMApexAlpha = cMApex/alpha;
    
    % center of pressure
    xCpNondimensional = - cMApexAlpha./cLAlpha(i);
    xCp = xCpNondimensional*wings{i}.mac;
    xCpMacPercentage(i) = (xCp - wings{i}.macLeadingEdgeX)/wings{i}.mac;
    xCpCStarPercentage(i) = xCp/(rootChord - tipChord);
    
    % leading edge suction and induced drag
    cT(i) = wings{i}.LeadingEdgeSuction(uoo, alpha);
    cDi(i) = cL(i)*alpha - cT(i);
end
clear cMApex cMApexAlpha xCpNondimensional macQuarterChord

%% Numeical results
e = Oswald(AR, Lambda);    % Oswald efficiency number
cLFunc = @(cLAlpha) cLAlpha*alpha;
cDiFunc = @(AR, cL) cL.^2./(pi*AR);

% finite wing theory functions
finite.cLAlphaFunc = @(Lambda, AR, e) ...
                     2*pi*cos(Lambda)./(1 + 2*cos(Lambda)./(AR.*e));

% finite wing theory results
finite.cLAlpha = finite.cLAlphaFunc(Lambda, AR, e);
finite.cL = cLFunc(finite.cLAlpha);
finite.cDi = cDiFunc(AR, finite.cL);
finite.xCpMacPercentage = 0.25;

% errors in relation to finite wing theory
finite.cLAlphaError = abs(cLAlpha - finite.cLAlpha)./finite.cLAlpha;
finite.cDiError = abs(cDi - finite.cDi)./finite.cDi;
finite.xCpMacPercentageError = ...
    abs(xCpMacPercentage - finite.xCpMacPercentage) / ...
    finite.xCpMacPercentage;

% slender wing theory functions
slender.cLAlphaFunc = @(AR) pi*AR/2;

% slender wing theory results
slender.cLAlpha = slender.cLAlphaFunc(AR);
slender.cL = cLFunc(slender.cLAlpha);
slender.cDi = cDiFunc(AR, slender.cL);
slender.xCpCStarPercentage = 2/3;

% errors in relation to slender wing theory
slender.cLAlphaError = abs(cLAlpha - slender.cLAlpha)./slender.cLAlpha;
slender.cDiError = abs(cDi - slender.cDi)./slender.cDi;
slender.xCpCStarPercentageError = ...
    abs(xCpCStarPercentage - slender.xCpCStarPercentage) / ...
    slender.xCpCStarPercentage;

% excel sheets of theoretical results and errors
T = table(finite.cLAlpha.', ...
          finite.cDi.', ...
          finite.xCpMacPercentage*ones(4,1), ...
          finite.cLAlphaError.', ...
          finite.cDiError.', ...
          finite.xCpMacPercentageError.', ...
          'VariableNames', {'cLAlpha', 'cDi', 'xCp', ...
                            'cLAlpha_error', 'cDi_error', 'xCp_error'});
writetable(T, 'theory comparison.xlsx', 'Sheet', 'Finite wing theory')
T = table(slender.cLAlpha.', ...
          slender.cDi.', ...
          slender.xCpCStarPercentage*ones(4,1), ...
          slender.cLAlphaError.', ...
          slender.cDiError.', ...
          slender.xCpCStarPercentageError.', ...
          'VariableNames', {'cLAlpha', 'cDi', 'xCp', ...
                            'cLAlpha_error', 'cDi_error', 'xCp_error'});
writetable(T, 'theory comparison.xlsx', 'Sheet', 'Slender wing theory')
clear T

% text output
fprintf('=================================================\n')
fprintf('VLM results for N = [%d,%d] and alpha = %.1f [deg]\n', ...
        wings{1}.xTileCount, wings{1}.yTileCount/2, rad2deg(alpha))
fprintf('=================================================\n')
for i = 1:wingCount
    fprintf('Wing %d\n', i)
    fprintf('------\n')
    fprintf('AR = %d\n', AR(i))
    fprintf('cLAlpha = %.3f\n', cLAlpha(i))
    fprintf('xCp as percentage of MAC = %.3f%%\n', ...
            xCpMacPercentage(i)*100)
    fprintf('xCp as percentage of (cR - cT) = %.3f%%\n', ...
            xCpCStarPercentage(i)*100)
    fprintf('cT = %.3e\n', cT(i))
    fprintf('cDi = %.3e\n\n', cDi(i))
end

%% Plots
% plot parameters
lw = 1.2;    % line width
fs = 15;     % font size
msX = 8;      % marker size
load('colors.mat')

% FIGURE 1: geometry of wings and VLM tiles
figure
hold on
grid on
title({'Geometry of wings and $VLM$ tiles:', ...
       'various aspect ratios'}, 'FontSize', fs)
xlabel('$x$', 'FontSize', fs)
ylabel('$y$', 'FontSize', fs)
zlabel('$AR$', 'FontSize', fs)
wingColors = {colors.light_blue, colors.blue, colors.royal_blue colors.purple};
for i = 1:wingCount
    wings{i}.PlotSelf('color', wingColors{i}, 'tiles', 'z offset', i-1);
end
zticks([0, 1, 2, 3])
zticklabels(string(AR))
view(-124, 36)
axis image
hold off

% FIGURE 2: cLAlpha comparison with theory
figure
hold on
grid on
ARLinspace = linspace(min(AR), max(AR));
LambdaLinspace = sweep(span(ARLinspace));
eLinspace = Oswald(ARLinspace, LambdaLinspace);
title('$c_{L,\alpha}$ comparison with theoretical results', 'FontSize', fs)
xlabel('$AR$', 'FontSize', fs)
ylabel('$c_{L,\alpha}\ [\frac{1}{rad}]$', 'FontSize', fs)
plot(ARLinspace, ...
     finite.cLAlphaFunc(LambdaLinspace, ARLinspace, eLinspace), ...
    'LineWidth', lw, 'Color', colors.green)
plot(ARLinspace, slender.cLAlphaFunc(ARLinspace), '--', ...
     'LineWidth', lw, 'Color', colors.yellow)
plot(AR, cLAlpha, 'x', ...
     'LineWidth', lw, 'MarkerSize', msX, 'Color', colors.red)
fill([6, 6, 10, 10], [0, 6, 6, 0], colors.green, ...
     'LineStyle', 'none', 'FaceAlpha', 0.1)
fill([1, 1, 1.5, 1.5], [0, 6, 6, 0], colors.yellow, ...
     'LineStyle', 'none', 'FaceAlpha', 0.1)
legend('Finite wing theory', 'Slender wing theory', '$VLM$ results', ...
       'Finite wing theory region', 'Slender wing theory region', ...
       'Location', 'Southeast');
xticks(AR)
ylim([0, 6])
hold off

% FIGURE 3: cDi comparison with theory
figure
hold on
grid on
cDicLSquared = 1./(pi*ARLinspace);
title('$c_{D_i}$ comparison with theoretical results', 'FontSize', fs)
xlabel('$AR$', 'FontSize', fs)
ylabel('$c_{D_i}$', 'FontSize', fs)
plot(ARLinspace, cDicLSquared, ...
    'LineWidth', lw, 'Color', colors.blue)
plot(AR, cDi./cL.^2, 'x', ...
     'LineWidth', lw, 'MarkerSize', msX, 'Color', colors.red)
legend('Theoretical curve', '$VLM$ results', ...
       'Location', 'Northeast');
xticks(AR)
hold off

% FIGURE 4: xCp comparison with finite wing theory
figure
hold on
grid on
cDicLSquared = 1./(pi*ARLinspace);
title('$x_{C_p}$ comparison with finite wing theory', 'FontSize', fs)
xlabel('$AR$', 'FontSize', fs)
ylabel('$\frac{x_{C_p}-x_{MAC}}{c_{MAC}}$', 'FontSize', fs)
plot([min(AR), max(AR)], finite.xCpMacPercentage*[1, 1], ...
    'LineWidth', lw, 'Color', colors.blue)
plot(AR, xCpMacPercentage, 'x', ...
     'LineWidth', lw, 'MarkerSize', msX, 'Color', colors.red)
xticks(AR)
ylim([0.1, 0.4])
legend('Location according to finite wing theory', '$VLM$ location', ...
       'Location', 'Southeast');
hold off

% FIGURE 5: xCp comparison with slender wing theory
figure
hold on
grid on
title('$x_{C_p}$ comparison with slender wing theory', 'FontSize', fs)
xlabel('$AR$', 'FontSize', fs)
ylabel('$\frac{x_{C_p}}{c_{root}-c_{tip}}$', 'FontSize', fs)
plot([min(AR), max(AR)], slender.xCpCStarPercentage*[1, 1], ...
    'LineWidth', lw, 'Color', colors.blue)
plot(AR, xCpCStarPercentage, 'x', ...
     'LineWidth', lw, 'MarkerSize', msX, 'Color', colors.red)
xticks(AR)
ylim([0.5, 0.8])
legend('Location according to finite wing theory', '$VLM$ location', ...
       'Location', 'Southeast');
hold off
