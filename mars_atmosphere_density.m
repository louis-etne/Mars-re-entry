clear all;
close all;
clc;

%% Load data
addpath('Data');
addpath('Functions');

load('constants.mat', 'hs', 'rho0');

%% Plot Mars atmosphere density
h = 0:120000;

figure;
ax = axes;
plot(h, MarsDensity(h));
xlabel('Altitude (m)');
ylabel('Density \rho (kg.m^{-3})');
title('Mars density versus Altitude using an exponential model');
set(ax, 'Xdir', 'reverse');
saveas(ax, 'Plots/density_vs_altitude.png');