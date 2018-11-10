clc
close all
clear variables

%% Parameters
xx = 15;    % last two NACA airfoil digits
t = xx/100;    % max thickness as fraction of cord
xs = 0;
xf = 1;

%% Wing
yt = @(x) 5*t * (0.2969 * sqrt(x) - ...
                 0.1260 * x - ...
                 0.3516 * x.^2 + ...
                 0.2843 * x.^3 - ...
                 0.1015 * x.^4);

%% Solver

%% Plots
