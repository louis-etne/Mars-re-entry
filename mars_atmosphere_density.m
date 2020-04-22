clear all;
close all;
clc;

%% Load data
addpath('Data');
addpath('Functions');

%% Plot Mars atmosphere density
h = 0:120000; % m - Altitudes at which the density will be calculated

figure;
ax = axes;
plot(h, MarsDensity(h));
xlabel('Altitude (m)');
ylabel('Density \rho (kg.m^{-3})');
title('Mars density versus Altitude using an exponential model');
set(ax, 'Xdir', 'reverse');